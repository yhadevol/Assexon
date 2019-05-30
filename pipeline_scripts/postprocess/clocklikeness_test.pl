#!/usr/bin/env perl

use warnings;
use strict;
use List::Util qw/shuffle/;
use Statistics::Distributions qw< chisqrprob >;
use Parallel::ForkManager;
use Getopt::Long;   # include the module for input

# input dir of aligned loci, output dir of loci following molecular clock assumption, concatenated ML tree
my ($dir, $clocklike, $besttree, $help);
my $outfile = "likelihood_list.txt"; # List of likelihood statistics for each gene 
my $run_paup = "run_paup"; # running dir
my $score = "score"; # paup output
my $log = "log"; # log file
my $conf_level = 0.05; # conf_level of chi-square test
my $pdis_trsd = 0; # Minimum p-distance required to conduct molecular clock test
my $cpu = 1;

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'indir:s', \$dir,
                      'besttree:s', \$besttree,
					  'clocklike:s', \$clocklike,
                      'likelihood_list:s', \$outfile,
                      'pdis:s', \$pdis_trsd,
                      'conf_level:f', \$conf_level,
                      'cpu:i', \$cpu,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check essential options
my @essential = qw/indir besttree clocklike/;
check_option(\@essential, \@ARGVoptions);                    

# check paup
`paup -n` or die "\nERROR: Cannot find paup in \$PATH. It may just not named as 'paupaxxx' rather than paup\n\n"; # check paup 

# check input dir and get file name
opendir (DIR, $dir) or die "\nERROR: Cannot find input directory \"$dir\", please check --indir ($!)\n\n"; #open folder containing all fasta data file
my @infile = grep {$_ !~ /^\./} readdir DIR;
@infile = grep {$_ =~ /fa$|fas$|fasta$/i} @infile;
closedir DIR;

# check and read provided tree
open BESTTREE, "$besttree" or die "\nERROR: Cannot find concatenated ML tree \"$besttree\", please check --besttree ($!)\n\n";
chomp(my $tree = <BESTTREE>);
close BESTTREE;

# get all taxa in tree
my (@all_taxa) = $tree =~ /([^:|^\(|^\)|^\,]+)\:/g;

# replace name
my (%all_taxa, @replaced_name);
my $order = 1;
foreach my $taxa (@all_taxa) {
my $replaced_name = "taxon$order";
$tree =~ s/$taxa\:/$replaced_name\:/;
$all_taxa{$taxa} = $replaced_name;
push @replaced_name, $replaced_name;
$order++;
}

# copy provided tree, and add ":0.0;"
my ($longest_str) = $tree =~ /(\(\S+\))/;
$longest_str = sentific_to_normal($longest_str);
open BESTTREE, ">besttree";
print BESTTREE "$longest_str:0.0;\n";
close BESTTREE;

# create running and output dir
if (!(-e $clocklike)) {
mkdir $clocklike or die "\nERROR: Cannot create output directory \"$clocklike\", please check --clocklike ($!)\n\n"; # clock like gene
}
mkdir $run_paup; # running dir
mkdir $score; # paup output
mkdir $log; # log file

# prepare output file and header
my ($OUT_FILE);
open ($OUT_FILE, ">$outfile") or die "\nERROR: Cannot write list of statistics for each locus \"$outfile\", please check --likelihood_list ($!)\n\n";
print $OUT_FILE "#Null=Molecular clock constrained\n#Alt=Molecular clock not constrained\n";
print $OUT_FILE "Gene Name\tNumber of Taxa\tLength of Alignment (bp)\tLikelihood Ratio\tdf\tP-value\t-In(Alt)\t-In(Null)\n";

# split dataset into $split part
my $split = $cpu;
my $splited_array = split_array(\@infile, $split);

# initiate multi-process
my $pm = Parallel::ForkManager->new(int($cpu/2));

# find clocklike loci in parallel
DATA_LOOP: while (my $array = shift @$splited_array) {
$pm->start and next DATA_LOOP;

	LOOP: while (my $genename = shift @$array) {
		next if ($genename =~ /^\./);
		
		# get loci name
		my ($gene) = $genename =~ /(\S+)\.fa\S*$/;
		
	    # check input file is wrongly formatted or have any strange characters 
		my $genefile = "$dir/$genename";
		my ($seqhash, $input_order, $aln_len) = readfas($genefile);
		
		my @exist_taxa;
		foreach my $taxon (@$input_order) {
			if (exists $all_taxa{$taxon}) {
			push @exist_taxa, $all_taxa{$taxon};
			}
		}
				
		# skip loci with too few taxa
		my $taxa_num = @exist_taxa;
		if ($taxa_num < 3) {
		print STDERR "WARNING: At least 3 taxa are required, and these taxa must also exist in provieded tree. There are only $taxa_num in $genefile. Skip this file\n";
		next LOOP;
		}
		
		# skip not aligned loci
		if ($aln_len eq "na") {
		print STDERR "WARNING: Length of sequences in $genefile are not the same. It may be not aligned. Skip this file\n";
		next LOOP;
		}
		
		# hash of sequence
		my %sequence;
		foreach my $taxon (sort keys %$seqhash) {
		$sequence{$all_taxa{$taxon}} = $seqhash->{$taxon};
		}
		
		# taxa, alignment length will written in paup input
		my @taxa = @exist_taxa;
		my $numofchar = $aln_len;
	
		# calculate num of different nucleotide between sequence
		my $numoftaxa = scalar @taxa;
		my $diff = 0;
		my $compare = 0;
		for (my $i=0; $i<($numoftaxa-1); $i++) {
		my $taxa1seq = $sequence{$taxa[$i]};
			for (my $j=$i+1; $j<@taxa; $j++) {
			my $taxa2seq = $sequence{$taxa[$j]};
			$diff += diff_nucleo($taxa1seq, $taxa2seq);
			$compare++;
			}
		}
		
		# do not conduct molecular clock test if p-dis is too close
		my $ave_pdis = sprintf("%.3f" , $diff/$compare);
		if ($ave_pdis <= $pdis_trsd) {
		print STDERR "WARNING: Average P-distance among $genename is equal to or lower than $pdis_trsd, so it won't be included in molecular clock test\n";
		next;
		}
		
		# get the tree only containing taxa in provided fas file
		my $reftree = subtree (\@taxa); #trim the tree
		
		# generate paup file
		generate_paup_file(\%sequence, $reftree, $gene, $numofchar);

		# run paup
		`paup -f -n $run_paup/$gene.paup`;
	
		my $lra = 0;
		#parse the log likelihood value for the paup result files
		open (my $INFILE, "$score/$gene.score") or die "\nERROR: Cannot find $score/$gene.score ($!)\n\n";
		while (my $line = readline ($INFILE)){
			chomp $line;
			if ((my $ln0) = $line =~ /^1\s+(\S+)\s+[\d\.]+/){ # get the likelihood when molecular clock is not constrained
			readline ($INFILE);
			my $line4 = readline ($INFILE);
			(my $ln1) = $line4 =~ /^1\s+(\S+)\s+\S+/; # get the likelihood when molecular clock is constrained
			$lra = 2 * ($ln1 - $ln0); # likelihood ratio
			my $degrees_of_freedom = $numoftaxa-2; # degree of freedom
			my $probability = sprintf("%0.3f", chisqrprob($degrees_of_freedom, $lra)); # probability of that there's no difference between alt and null
			$lra = sprintf("%0.3f", $lra);
			$ln0 = sprintf("%0.3f", $ln0);
			$ln1 = sprintf("%0.3f", $ln1);
			print $OUT_FILE "$gene\t$numoftaxa\t$numofchar\t$lra\t$degrees_of_freedom\t$probability\t$ln0\t$ln1\n"; # output statistics
			`cp $dir/$genename $clocklike/$genename` if ($probability > $conf_level); # mv gene to output dir if there's no difference between alt and null
			}
		}
		close $INFILE;
	}
	
$pm->finish;
}
$pm->wait_all_children;	

close $OUT_FILE;

# remove dir and copied file
`rm -rf $run_paup`;
`rm -rf $score`;
`rm -rf $log`;
`rm -f besttree`;

#####################################################
# Subroutines
#####################################################

# subroutine to trim the tree and remove a list of taxa
sub subtree {
my ($taxaref) = shift @_; # taxa want to kept in tree

my @commontaxa = sort @$taxaref; #make the array for the list of taxa want to be kept
my %commontaxa = map {$_ , ""} @commontaxa;

#get the best tree
my $treefile = "besttree";
open (my $TREEFILE, "<$treefile") or die "\nERROR: Can't find $treefile\n\n";
my $besttree = readline ($TREEFILE);
chomp $besttree;

# get taxa not exist in current gene
my @taxa;
foreach my $taxa (@replaced_name) {
push @taxa, $taxa if (! exists $commontaxa{$taxa});
}
my @sortedtaxa = sort @taxa;

foreach my $taxon (@sortedtaxa) {
my ($part1, $part2) = split /$taxon\:/, $besttree; #get the string before and after the taxon, split with ":" in case of wrong split;
$part2 = ":$part2"; # add the ":"

my $parenth = 1; #initial the parentheses variable
		
if ($besttree =~ s/\($taxon\:[\d\.]+\,//) {# remove the taxon on left side and affiliated junk
my $len = length ($part2);
	my $i;
	for ($i = 0; $i < $len; $i ++) {#loop through the string, find the paired parenthesis
		my $char = substr ($part2, $i, 1);#get the current character
		if ($parenth == 0) {#check if the total number of parentheses matches
		last;
		}
		elsif ($char eq "\(") {
		$parenth ++;
		}
		elsif ($char eq "\)") {
		$parenth --;
		}
	} 
	
	my $substr = substr($part2, $i + 1);
	
	if ($substr ne "0\.0;") {
	my ($branch2) = $substr =~ /(^[\d\.]+)/;
	my ($branch1) = $part2 =~ /\:([\d\.]+)\)\:$branch2/;
	
	my $newbranch = sprintf("%.20f", $branch1 + $branch2);
	$besttree =~ s/$branch1\)\:$branch2/$newbranch/;
	}else {
	$besttree =~ s/\)\:0\.0\;$/\;/;
	}
} elsif ($besttree =~ /\,$taxon\:[\d\.]+\)\:[\d\.]+/) {# remove the taxon on right side and affiliated junk
my ($branch1, $branch2) = $besttree =~ /\:([\d\.]+)\,$taxon\:[\d\.]+\)\:([\d\.]+)/;
my $newbranch = sprintf("%.20f", $branch1 + $branch2);

	my $len = length ($part1);
	my $i;
	for ($i = 0; $i < $len; $i ++) {#loop through the string, find the paired parenthesis
		my $char = substr ($part1, $len - $i - 1, 1);#get the current character
		if ($parenth == 0) {#check if the total number of parentheses matches
		last;
		}
		elsif ($char eq "\)") {
		$parenth ++;
		}
		elsif ($char eq "\(") {
		$parenth --;
		}
	}
	 
	my $substr = substr ($part1, $len - $i + 1);#get the string at the matched left parenthesis
	my ($sub) = $substr =~ /(^[^:]+:)/;
	
	if (my ($paren, $partsub) = $sub =~ /(^\(+)(\S+)/) {#if the $sub has leading (
		my $cnt = length ($paren);
		my $temp = "\\\(";
		for ($i = 1; $i < $cnt; $i ++) {
		$temp = $temp . "\\" . "\(";
		}
		$besttree =~ s/\($temp$partsub/$sub/; #delete the (
	}
	else { 
	$besttree =~ s/\($sub/$sub/; #delete the (
	}
	
	$besttree =~ s/\:[\d\.]+\,$taxon\:[\d\.]+\)\:[\d\.]+/\:$newbranch/; #delete the taxon and replace with new branch length
}
	#if the outgroup taxon was deleted, unroot the tree
	if ($besttree !~ /\S+\:0.0\;$/) {
	my ($temp2) = $besttree =~ /(^.+);$/;
	$besttree = "\(" . "$temp2" . "\):0.0\;";
	}
}

# remove branch length at root
my ($inner_group) = $besttree =~ /\((.+)\)/;

# get the position of outermost parenthese
my (@parenth_pair_left, @parenth_pair_right);
my $parenth = 0; #initial the parentheses variable
my $len = length ($inner_group);
for (my $i = 0; $i < $len; $i ++) {#loop through the string, find the paired parenthesis
	my $char = substr ($inner_group, $i, 1);#get the current character
	
	if (($parenth == 0)&&($char eq "\(")) {
	push @parenth_pair_left, $i;
	}
	
	if ($char eq "\(") {
	$parenth ++;
	}
	elsif ($char eq "\)") {
	$parenth --;
	}
	
	if (($parenth == 0)&&($char eq "\)")) {
	push @parenth_pair_right, $i;
	}
} 

# remove content of outermost group
my @substr;
my $substract_len = 0;
for (my $i = 0; $i < @parenth_pair_left; $i++) {
my $left_pos = $parenth_pair_left[$i] - $substract_len;
my $right_pos = $parenth_pair_right[$i] - $substract_len;
my $length_longest_str = length($inner_group);

my $substr_left = substr $inner_group, 0, $left_pos;
my $substr_right = substr $inner_group, $right_pos+1, $length_longest_str-$right_pos-1;
$inner_group = $substr_left . $substr_right;

my $len = $right_pos-$left_pos+1;
$substract_len += $len;
}

# num of outermost group
my @comma = $inner_group =~ /\,/g;

if (@comma == 0){
($inner_group) = $besttree =~ /\((.+)\)/;
($besttree) = $inner_group =~ /(\(.+\))/;
$besttree .= ":0.0;";
}

# return trimmed tree
return $besttree;
}

# subroutine generate paup file
sub generate_paup_file {
my $sequence = shift; # sequence
my $reftree = shift; # trimmed tree
my $gene = shift; # gene name 
my $numofchar = shift; # alignment length

my $file = "$run_paup/$gene.paup"; # paup file
my $score = "../$score/$gene.score"; # paup output
my $logfile = "../$log/$gene.log"; # log file

# write paup file
my $FILE;
open ($FILE, ">$file") or die "\nERROR: Can't open the $file ($!)\n\n";
my $numoftaxa = scalar keys %$sequence;
# sequence block
print $FILE "#NEXUS\nBEGIN DATA;\ndimensions ntax=$numoftaxa nchar=$numofchar;\nformat missing=?\ndatatype=DNA gap= -;\n\nmatrix\n\n[$gene]\n\n";
foreach my $taxon (sort keys %$sequence) {print $FILE "$taxon\t$sequence->{$taxon}\n";}
# tree block
print $FILE ";\nend;\n\nbegin trees;\n\ttree mytree=[&U]$reftree\nend;\n\nbegin paup;\n";
# compute log likelihood
print $FILE "\tlog file =$logfile replace =yes start=yes;\n";
print $FILE "\tset root=outgroup criterion=likelihood warnroot=no taxlabels=full MONITOR = NO;\nend;\n\nbegin paup;\n";
print $FILE "\tlscores 1/ nst=6 rates=gamma rmatrix=estimate shape=estimate clock=no scorefile=$score replace=yes;\n\troottrees;\n";
print $FILE "\tlscores 1/ nst=6 rates=gamma rmatrix=previous shape=estimate clock=yes scorefile=$score append=yes;\nend;\n";
close $FILE;
}

# find the number of different nucleotide between sequences
sub diff_nucleo {
my($string1, $string2) = @_; # seq1, se12

# we assume that the strings have the same length
my($length) = length($string1);
my($position);
my($count) = 0;
	
	# compare nucleotide position by position
	for ($position=0; $position < $length ; ++$position) {
		if ((substr($string1,$position,1) =~ /[A-Z?*]{1}/) && (substr($string2,$position,1) =~ /[A-Z?*]{1}/)){#ignore gaps and missing data
			if(uc substr($string1,$position,1) ne uc substr($string2,$position,1)) {
			++$count;
			}
		}
	}

# return number of different nucleotide 
return ($count);
}

# subroutine to split array in to several sub-array
sub split_array {
my $array = shift;
my $split = shift;
	
	$split++ if ($split == 0);
	my $block_num = int(@$array/$split);
	$block_num++ if ($block_num == 0);
	
	my @splited_array;
	my $size_cnt = 0;
	my $block_cnt = 0;
	
	while (my $file = shift @$array) {
		if ($size_cnt >= $block_num) {
		$size_cnt = 0;
		$block_cnt ++;
		}
	push @{$splited_array[$block_cnt]}, $file;
	$size_cnt ++;
	}

return(\@splited_array);
}

# subroutine to reformat input fasta file
sub readfas {
my $file = shift;

open INFILE, $file or die "\nERROR: Cannot find input file \"$file\" ($!)\n\n";

my (%seq, %seqlength, @input_order, $seq, $taxon, $lasttaxon, %taxa_cnt);

while (my $line = <INFILE>) {

	# remove enter
	$line =~ s/\r//g;
	chomp $line;
	
	# if find >
	if ($line =~ /^>(\S+)/) {
		$taxon = $1;
		
		# if find previous taxon name, and it hasn't recorded before
		if ($lasttaxon) {
			if (!(exists $seq{$lasttaxon})) {
				$seq =~ s/\s//g;
				$seq = uc $seq;
				my @strange_char = $seq =~ /[^A-Z?*\.\-]/g;
				my $length = length($seq);
				if (($length > 0)&&(@strange_char == 0)) {
				$seq{$lasttaxon} = $seq;
				push @input_order, $lasttaxon;
				$seqlength{$length} = ""; 
				} else {
					if ($length == 0) {
					print STDERR "WARNING: No nucleotide found sequence of $lasttaxon in $file. This taxon will be discarded.\n";
					}
					if (@strange_char > 0) {
					my $strange_char = join " ", @strange_char;
					print STDERR "WARNING: Found strange character \"$strange_char\" in sequence of $lasttaxon in $file. This taxon will be discarded.\n";					
					}
				} 
			} else {
			$taxa_cnt{$lasttaxon}++
			}
		}
		
		$seq = "";
		$lasttaxon = $taxon;
	} else {
	$seq .= $line;
	}
}

if ($lasttaxon) {
	if (!(exists $seq{$lasttaxon})) {
		$seq =~ s/\s//g;
		$seq = uc $seq;
		my @strange_char = $seq =~ /[^A-Z?*\.\-]/g;
		my $length = length($seq);
		if (($length > 0)&&(@strange_char == 0)) {
		$seq{$lasttaxon} = $seq;
		push @input_order, $lasttaxon;
		$seqlength{$length} = ""; 
		} else {
			if ($length == 0) {
			print STDERR "WARNING: No nucleotide found sequence of $lasttaxon in $file. This taxon will be discarded.\n";
			}
			if (@strange_char > 0) {
			my $strange_char = join " ", @strange_char;
			print STDERR "WARNING: Found strange character \"$strange_char\" in sequence of $lasttaxon in $file. This taxon will be discarded.\n";					
			}
		} 
	} else {
	$taxa_cnt{$lasttaxon}++
	}
}

close INFILE;

foreach my $taxon (sort keys %taxa_cnt) {
my $taxon_num = $taxa_cnt{$taxon};
$taxon_num++;
print STDERR "WARNING: Found $taxon_num \"$taxon\" in $file. Only first sequence named in \"$taxon\" will be kept\n";
}

my @seqlength = sort keys %seqlength;

my $aln_len;
if (@seqlength == 1) {
$aln_len = $seqlength[0];
} else {
$aln_len = "na";
}

return (\%seq, \@input_order, $aln_len);
}

# Subroutine to check missing option
sub check_option {
my $essential = shift;
my $ARGVoptions = shift;

	# get all arguments
	my %args;
	foreach my $args (@$ARGVoptions) {
		if ($args =~ /^--*(\S+)/) {
		$args{$1} = "";
		}
	}
	
	# find all missing arguments
	my @missing_option;
	foreach my $ess_arg (@$essential) {
		if (!(exists $args{$ess_arg})) {
		push @missing_option, "--$ess_arg";
		}
	}
	
	# print out all missing arguments
	if (@missing_option >= 1) {
	my $missing_option = join ", ", @missing_option;
	die "\nERROR: Missing option $missing_option\n\nUse -h or --help for more information on usage\n\n";
	}
}

sub sentific_to_normal {
my $tree = shift;

my @sentific = $tree =~ /(\d+\.\d+E\-\d+)/ig;

foreach my $sentific (@sentific) {
my ($float, $index) = $sentific =~ /\d+\.(\d+)E\-(\d+)/i;
my $float_len = length($float);
$float_len += $index;
my $normal = sprintf("%.${float_len}f", $sentific);
$tree =~ s/$sentific/$normal/;
}

return ($tree);
}

sub usage {
print STDERR "
Script name: clocklikeness_test.pl

This is a script to pick out loci following molecular clock hypothesis in fasta format.

Dependencies: 
(1) paup 4.0a (build 161) or higher (rename 'paupa***' as 'paup and put it under \$PATH before use). If current version is expired, download a lastest version from 'http://phylosolutions.com/paup-test/'
(2) perl module: 
	1. Statistics::Distributions
	2. Parallel::ForkManager

Example usage:
(1) Alignments under 'nf_aligned' and concatenated ML tree 'besttree.tre' are used to calculate likelihoods whicg are recorded in 'likelihood_list.txt'. Loci following molecular clock hypothesis are written to 'clocklike_dir'. Run script in 4 process:

	perl clocklikeness_test.pl --indir nf_aligned --besttree besttree.tre --clocklike clocklike_dir --cpu 4

Input files: 
(1) nf_aligned
(2) besttree.tre

Output files:
(1) clocklike_dir
(2) likelihood_list.txt

Options:
--indir
  Directory containing aligned nucleotide sequences
--besttree
  Unrooted concatenated maximum likelihood tree 
--clocklike
  Directory containing genes which follows molecular clock hypothesis
--likelihood_list
  List of statistics for each locus, named as '$outfile' in default 
--pdis
  Minimum p-distance required to conduct molecular clock test, $pdis_trsd in default
--conf_level
  Confidence level of molecular clock test, $conf_level in default, molecular clock hypothesis will be rejected if p-value of the gene is lower than this threshold
--cpu
  Limit the number of CPUs, $cpu in default
--help , -h
 Show this help message and exit

Author: Hao Yuan                                                                     
        Shanghai Ocean University                                               
        Shanghai, China, 201306                                                               
                                                                                         
Created by: Nov 20, 2018                                                              
                                                                                         
Last modified by: 
";
exit;
}
