#!/usr/bin/env perl

use Getopt::Long; 
use warnings;
use strict;
use Parallel::ForkManager;

my ($dir, $mono_constrain, $ml_tree_dir, $constrained_tree_dir, $help);
my $flag = "monophyletic_loci"; # dir containing genes which follow monophyletic assumption
my $run_paup_dir = "run_paup"; # dir for running paup
my $run_paup_log_dir = "log"; # dir for log file
my $tree_copy_path = "copy_tree";
my $likelihood_list = "likelihood_list.txt"; # list recording likelihood under alternative and null assumption
my $topology_test = "SH"; # topology test used
my $bootstrap = 1000; # times of rell bootstrap
my $conf_level = 0.05; # p-value
my $pdis_trsd = 0; # minimum p-dis required for topology test
my $cpu = 1;

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'indir:s', \$dir,
					  'mono_constrain:s', \$mono_constrain,
					  'ml_tree_dir:s', \$ml_tree_dir,
					  'constrained_tree_dir:s', \$constrained_tree_dir,
                      'monophyly_dir:s', \$flag,
                      'likelihood_list:s', \$likelihood_list,
                      'bootstrap:s', \$bootstrap,
                      'pdis:s', \$pdis_trsd,
                      'topology_test:s', \$topology_test,
                      'conf_level:f',\$conf_level,
                      'cpu:i', \$cpu,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# if essential options are missed
my @essential = qw/indir ml_tree_dir constrained_tree_dir/;

# check missing option
check_option(\@essential, \@ARGVoptions);

# check paup
`paup -n` or die "\nERROR: Cannot find paup in \$PATH. It may just not named as 'paupaxxx' rather than paup \n\n";

# check whether name of topology test is correct
$topology_test = uc $topology_test;
if (($topology_test ne "AU") && ($topology_test ne "SH")) {
die "\nERROR: Only \"AU\" and \"SH\" test is avaliable, please specify a correct name of topology test \n\n";
}

my $tree_suffix;
if (($ml_tree_dir)&&($constrained_tree_dir)) { # if both ml and constrained tree is provided

	# check existence of ml and constrained tree dir
	if (!(-e $ml_tree_dir)) { 
	die "\nERROR: Cannot find directory containing unconstrained ML tree \"$ml_tree_dir\", please check --ml_tree_dir ($!)\n";
	}
	if (!(-e $constrained_tree_dir)) {
	die "\nERROR: Cannot find directory containing constrained ML tree \"$constrained_tree_dir\", please check --constrained_tree_dir ($!)\n";
	}
	
	# check whether number of tree is the same
	my $ml_tree_num = `ls $ml_tree_dir|wc -l` =~ /(\d+)/;
	my $constrained_tree_num = `ls $constrained_tree_dir|wc -l` =~ /(\d+)/;
	if ($ml_tree_num != $constrained_tree_num) {
	die "\nERROR: File number under $ml_tree_dir and $constrained_tree_dir is not the same\n\n";
	}
	
	# check whether tree under ml and constrained tree dir are the same
	opendir ML_TREE_DIR, $ml_tree_dir;
	my @ml_tree = grep {$_ !~ /^\./} readdir ML_TREE_DIR;
	closedir ML_TREE_DIR;
	opendir CONSTRAINED_TREE_DIR, $constrained_tree_dir;
	my @constrained_tree = grep {$_ !~ /^\./} readdir CONSTRAINED_TREE_DIR;
	closedir CONSTRAINED_TREE_DIR;
	my %constrained_tree = map {$_, ""} @constrained_tree;
	foreach my $tree (@ml_tree) {
	($tree_suffix) = $tree =~ /\S+\.(\S+)$/;
		if (! exists $constrained_tree{$tree}) {
		die "\nERROR: Cannot find $tree under $constrained_tree_dir. Tree name under $constrained_tree_dir and $ml_tree_dir must be the same \n\n";
		}
	}
} else { # die if one of ml or constrained tree is not provided
die "\nERROR: Please specify --ml_tree_dir and --constrained_tree_dir if you have prepared ML and constrained tree\n\n";
}

opendir DIR, $dir or die "\nERROR: Cannot find input directory $dir, please check --indir ($!)\n";
my @infile = grep {$_ !~ /^\./} readdir DIR;
@infile = grep {$_ =~ /fa$|fas$|fasta$/i} @infile;
closedir DIR;

# run dir, log dir and dir for non-monophyletic gene
mkdir $run_paup_dir;
mkdir $run_paup_log_dir;

if (!(-e $flag)) {
mkdir $flag or die "\nERROR: Cannot create directory containing loci following the defined monophyletic group, please check --monophyly_dir ($!)\n\n";
}

# print header
open LIKELIHOOD_LIST, ">$likelihood_list" or die "\nERROR: Cannot write list of statistics for each locus \"$likelihood_list\", please check --likelihood_list ($!)\n\n";
print LIKELIHOOD_LIST "Gene Name\tNumber of Taxa\tLength of Alignment (bp)\tP-value\t-In(ML)\t-In(Constrained)\n";

$tree_copy_path = absolute_path($tree_copy_path);
mkdir $tree_copy_path;

my $split = $cpu;
my $splited_array = split_array(\@infile, $split);

my $pm = Parallel::ForkManager->new(int($cpu/2));

DATA_LOOP: while (my $array = shift @$splited_array) {
	
	LOOP: while (my $genename = shift @$array) {
		next if ($genename =~ /^\./);
		my ($gene) = $genename =~ /(\S+)\.fa\S*$/;

		my $genefile = "$dir/$genename"; #get the name of fasta file
		my ($seqhash, $input_order, $aln_len) = readfas($genefile);
	
		my $seq_num = @$input_order;
		if ($seq_num < 4) {
		print STDERR "WARNING: At least 4 taxa are required. There are only $seq_num in $genefile. Skip this file\n";
		next LOOP;
		}
	
		if ($aln_len eq "na") {
		print STDERR "WARNING: Number of characters among sequences are not the same in $genefile. It may be not aligned\n";
		next LOOP;
		}
	
		my %sequence = %$seqhash;
		my @taxa = @$input_order;
		my $numofchar = $aln_len;

		my $numoftaxa = @taxa;
	
		# calculate the number of difference between pair of sequences
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
	
		# do not conduct monophyly test if p-distance is too low
		my $ave_pdis = sprintf("%.3f" , $diff/$compare);
		if ($ave_pdis <= $pdis_trsd) {
		print STDERR "Average P-distance among $genename is equal to or lower than $pdis_trsd, so it won't be included in monophyly test\n";
		next;
		}
	
		# path to constrained and ml tree
		my $constrained_path = "$constrained_tree_dir/$gene.$tree_suffix";
		my $ml_path = "$ml_tree_dir/$gene.$tree_suffix";
		
		# conduct topology test if we find both constrained and ml tree
		if ((-e $ml_path) && (-e $constrained_path)) {
		
		# copy tree and replace taxa name as taxon 1.2.3..
		
		# save replaced name
		my $cnt = 0;
		my %replaced_seq;
		foreach my $taxon (sort keys %sequence) {
		my $replaced_taxa = "taxon$cnt";
		$replaced_seq{$taxon}->{seq} = $sequence{$taxon};
		$replaced_seq{$taxon}->{replaced_taxa} = $replaced_taxa;
		$cnt++;
		}
		
		# replace taxa name in ML tree
		open ML_TREE, $ml_path;
		chomp(my $ml_tree = <ML_TREE>);
		foreach my $taxon (sort keys %replaced_seq) {
		my $replaced_taxa = $replaced_seq{$taxon}->{replaced_taxa};
		$ml_tree =~ s/$taxon\:/$replaced_taxa\:/;
		}
		close ML_TREE;
		
		my $copied_ML_tree = "$tree_copy_path/$gene.ML.tree";
		open REPLACED_TREE, ">$copied_ML_tree";
		print REPLACED_TREE "$ml_tree\n";
		close REPLACED_TREE;
		
		# replace taxa name in constrained tree
		open CONSTRAINED_TREE, $constrained_path;
		chomp(my $constrained_tree = <CONSTRAINED_TREE>);
		foreach my $taxon (sort keys %replaced_seq) {
		my $replaced_taxa = $replaced_seq{$taxon}->{replaced_taxa};
		$constrained_tree =~ s/$taxon\:/$replaced_taxa\:/;
		}
		close CONSTRAINED_TREE;
		
		my $copied_constrained_tree = "$tree_copy_path/$gene.constrained.tree";
		open REPLACED_TREE, ">$copied_constrained_tree";
		print REPLACED_TREE "$constrained_tree\n";
		close REPLACED_TREE;
		
		# topo test
		topology_test_by_paup($copied_ML_tree, $copied_constrained_tree, \%replaced_seq, $gene, $numofchar, $topology_test); # generate paup file
		`paup -f -n $run_paup_dir/$gene.topotest.paup`; # run paup
			open LOG, "$run_paup_log_dir/$gene.topotest.log"; # get likelihood value from log
			while (my $line = <LOG>) {
				if ($line =~ /\s+Tree\s+-ln\s+L\s+Diff\s+-ln\s+L\s+[SH|AU].*/) {
				<LOG>;
				my $best = <LOG>;
				my $subopt = <LOG>;
				my @best = $best =~ /(\S+)/g;
				my @subopt = $subopt =~ /(\S+)/g;
				my %ln;
				$ln{$best[0]} = sprintf("%0.3f", $best[1]); # tree with optimal likelihood
				$ln{$subopt[0]} = sprintf("%0.3f", $subopt[1]); # tree with sub optimal likelihood
				my ($p_value) = $subopt[3] =~ /\~*([^\*|^\~]+)\**/; # get p-value
				$p_value = sprintf("%0.3f", $p_value);
				print LIKELIHOOD_LIST "$gene\t$numoftaxa\t$numofchar\t$p_value\t$ln{1}\t$ln{2}\n"; # print value
					if ($best[0] == 1) {
					`cp $dir/$genename $flag/$genename` if ($p_value > $conf_level); # move gene which following provided monophyly group
					}
				}
			}
		} else {
		print STDERR "WARNING: Cannot found either directory containing ML trees \"$ml_path\" or constrained ML trees \"$constrained_path\". \n";
		}
		
	}
	
$pm->finish;
}
$pm->wait_all_children;

close LIKELIHOOD_LIST;

# remove run dir, log dir and  copied tree path
`rm -rf $run_paup_dir`;
`rm -rf $run_paup_log_dir`;
`rm -rf $tree_copy_path`;

#####################################################
# Subroutines
#####################################################

# count the number of different nucleotide between pair of sequences
sub diff_nucleo {
my($string1, $string2) = @_;

# we assume that the strings have the same length
my($length) = length($string1);
my($position);
my($count) = 0;

	for ($position=0; $position < $length ; ++$position) {
		if ((substr($string1,$position,1) =~ /[A-Z?*]{1}/) && (substr($string2,$position,1) =~ /[A-Z?*]{1}/)){#ignore gaps and missing data
			if(uc substr($string1,$position,1) ne uc substr($string2,$position,1)) {
			++$count;
			}
		}
	}

return ($count);
}

# subroutine to generate paup file for topology test
sub topology_test_by_paup {
my $ml_path = shift;
my $constrained_path = shift;
my $sequence = shift;
my $gene = shift;
my $numofchar = shift;
my $topology_test = shift;

my $numoftaxa = scalar keys %$sequence;

# log, ml tree, constrained tree and paup file
my $log = "../$run_paup_log_dir/$gene.topotest.log";
my $ml_tree = "$ml_path";
my $constrain_tree = "$constrained_path";
my $paup = "./$run_paup_dir/$gene.topotest.paup";

open (my $FILE, ">$paup") or die "\nERROR: Cannot write paup file $paup ($!)\n\n";

# print header
print $FILE 
"#NEXUS
BEGIN DATA;
dimensions ntax=$numoftaxa nchar=$numofchar;
format missing=?\ndatatype=DNA gap= -;

matrix

[$gene]

";

# print sequence
foreach my $taxon (sort keys %$sequence) {
print $FILE "$sequence->{$taxon}->{replaced_taxa}\t$sequence->{$taxon}->{seq}\n";
}

# end of sequence block
print $FILE 
";
end;

";

# print tree block
print $FILE 
"Begin paup;
set autoclose=yes warntree=no warnreset=no increase=no ErrorBeep=no;
log start file=$log replace;

gettrees file=$ml_tree;
gettrees file=$constrain_tree mode=7;
DerootTrees;\n";

	# topology test
	if ($topology_test eq "SH") {
	print $FILE "lscores 1-2/ shtest=yes rell=yes bootReps=$bootstrap;\n";
	} elsif ($topology_test eq "AU") {
	print $FILE "lscores 1-2/ autest=yes rell=yes bootReps=$bootstrap;\n";
	}

print $FILE "end;"; # end of topology test block

close ($FILE);
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

# subroutine to translate relative path into absolute path
sub absolute_path {
my $path = shift;

my @str = $path =~ /([^\/]+)/g;
my $suffix = $str[-1];
$path =~ s/$suffix$//;
my $prefix = $path;

if (! $prefix) {
$prefix = "./";
}
my ($absolute_path) = `cd $prefix;pwd` =~ /(\S+)/;
$absolute_path .= "/$suffix";

return $absolute_path;
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

# subroutine to print usage
sub usage {
print STDERR "
Script name: monophyly_test.pl

This is a script to pick out aligned loci in fasta format which topology are not congurence with provided monophyly group.

Dependencies: 
(1) paup 4.0a (build 161) or higher (rename 'paupa***' as 'paup' and put it under \$PATH before use)

Example usage:
(1) ML tree and ML trees constrained by monophyletic group are user-provided. Then SH test is performed based on these tree:

	perl monophyly_test.pl --indir nf_aligned --topology_test SH --ml_tree_dir ml_tree_dir --constrained_tree_dir constrained_tree_dir

Input files:
(1) nf_aligned
(2) ml_tree_dir (if --user_tree is specified)
(3) constrained_tree_dir (if --user_tree is specified)

Output files:
(1) monophyletic_loci

Options:
--indir
  Directory containing aligned nucleotide sequences
--ml_tree_dir
  Directory containing ML tree, number and name of trees must the same as files under --constrained_tree_dir
--constrained_tree_dir
  Directory containing constrained trees, which can be produced by construct_tree.pl. Number and name of trees must the same as files under --constrained_tree_dir
--monophyly_dir
  Directory containing loci which topologies are following provided monophyly group, dir named as '$flag' in default
--likelihood_list
  List of statistics for each locus, named as '$likelihood_list' in default 
--bootstrap
  Number of rell bootstrap replicates, $bootstrap in default
--pdis
  Minimum p-distance required to conduct topology test, $pdis_trsd in default
--topology_test
  Topology tests avaliable including 'SH' (Shimodaira & Hasegawa Test) and 'AU' (Approximately Unbiased Test), 'SH' in default
--conf_level
  Confidence level of topology test, $conf_level in default, loci with p-value lower than this threshold will be considered that provided monophyletic group is not formed in this gene
--cpu:
  Limit the number of CPUs, $cpu in default
--help , -h
  Show this help message and exit

Author: Hao Yuan                                                                     
        Shanghai Ocean University                                               
        Shanghai, China, 201306                                                               
                                                                                         
Created by: June 27, 2018                                                              
                                                                                         
Last modified by:
";
exit;
}