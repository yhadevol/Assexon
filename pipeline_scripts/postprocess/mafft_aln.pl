#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use Parallel::ForkManager;
use Bio::Seq;

my ($dna_unaligned_dir, $dna_aligned_dir, $non_codon_aln, $keep_aa_aln, $help);
my $aa_aligned_dir = "aa_aligned"; # dir for aligned amino acid
my $cpu = 1;

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'dna_unaligned:s', \$dna_unaligned_dir,
					  'dna_aligned:s', \$dna_aligned_dir,
					  'non_codon_aln!', \$non_codon_aln,
					  'keep_aa_aln!', \$keep_aa_aln,
					  'cpu:i', \$cpu,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/dna_unaligned/;
check_option(\@essential, \@ARGVoptions);

# check mafft
`mafft --version 2> mafft_check`;
my $mafft_check = `grep "command not found" mafft_check`;
unlink "mafft_check";
if ($mafft_check) {
die "\nERROR: Cannot find mafft in \$PATH, It may just not named as \"mafft\"\n\n";
}

if (($non_codon_aln)&&($keep_aa_aln)) {
say STDERR "WARNING: No amino acid sequences will be generated if --non_codon_aln is specified";
}

# check input nucleotide dir
if (!(-e $dna_unaligned_dir)) {
die "\nERROR: Cannot find $dna_unaligned_dir ($!)\n\n";
}

# get all fasta file
opendir (DIR, $dna_unaligned_dir) or die "\nERROR: Cannot find input directory \"$dna_unaligned_dir\", please check --dna_unaligned ($!)\n\n";
my @unaligned = grep {$_ !~ /^\./} readdir DIR;
@unaligned = grep {$_ =~ /fa$|fas$|fasta$/i} @unaligned;
closedir DIR;

# automatically named as xx_aligned if output name of aligned nucleotide in not provided
$dna_aligned_dir = $dna_unaligned_dir . "_aligned" if (!($dna_aligned_dir)); 
if (!(-e $dna_aligned_dir)) {
mkdir $dna_aligned_dir or die "\nERROR: Cannot create output directory \"$dna_aligned_dir\", please check --dna_aligned ($!)\n\n";
}

my $split = $cpu;
my $splited_array = split_array(\@unaligned, $split);

my $pm = Parallel::ForkManager->new(int($cpu/2));

# 2 modes of alignment
if (! $non_codon_aln) { # align in codon
	say STDOUT "All sequences will be aligned in codon";
	
	DATA_LOOP: while (my $unaligned = shift @$splited_array) {
	$pm->start and next DATA_LOOP;
	codon_aln($unaligned);
	$pm->finish;	
	}
	$pm->wait_all_children;
	
	# keep amino acid alignment if --keep_aa_aln is specified
	if (! $keep_aa_aln) {
	`rm -rf $aa_aligned_dir`;
	}
	
} else { # simply align
	say STDOUT "All sequences will be simply aligned (not aligned in codon)";

	DATA_LOOP: while (my $unaligned = shift @$splited_array) {
	$pm->start and next DATA_LOOP;
	normal_aln($unaligned);
	$pm->finish;	
	}
	$pm->wait_all_children;
}

#####################################################
# Subroutines
#####################################################

# subroutine to align sequence in codon
sub codon_aln {
my $unaligned = shift;

mkdir "$aa_aligned_dir"; # create aligned amino acid dir

NEXT_FILE: while (my $file = shift @$unaligned) {
	my ($gene) = $file =~ /(\S+)\.fa\S*$/;

	# write unaligned amino acid sequences first
	my $translated_seq_num = 0;

	my $dna_unaligned_file = "$dna_unaligned_dir/$file";
	my $dna_unaligned_file_copy = "$aa_aligned_dir/$gene.dna.fas";
	my $aa_unaligned_file = "$aa_aligned_dir/$gene.aa.fas";
	
	my ($seqhash, $input_order, $aln_len) = readfas($dna_unaligned_file);

	my $seq_num = @$input_order;
	if ($seq_num == 0) {
	print STDERR "WARNING: Nothing found in $dna_unaligned_file. Skip this file.\n";
	next NEXT_FILE;
	}
	
	if ($seq_num == 1) {
	print STDERR "WARNING: Only 1 sequence found in $dna_unaligned_file. This file will be copied to $dna_aligned_dir\n";
	`cp $dna_unaligned_file $dna_aligned_dir/$file`;
	next NEXT_FILE;
	}
	
	open DNA_UNALIGNED_COPY, ">$dna_unaligned_file_copy";
	open AA_UNLALIGNED, ">$aa_unaligned_file";
	NEXT_SEQ: foreach my $taxon (@$input_order) {
		my $dna_seq = $seqhash->{$taxon};
		
		if ($dna_seq =~ /-/) { # skip if sequences have been aligned
		print STDERR "WARNING: Sequences of $dna_unaligned_file has already be aligned. Gaps will be removed and sequences will be realigned\n";
		$dna_seq =~ s/-//g;
		}
		if (length($dna_seq)%3 != 0)  { # skip if length of sequences cannot exact divided by 3
		print STDERR "WARNING: Length of sequence of $taxon in $dna_unaligned_file cannot be divided by 3. This sequence will be discarded before aligning\n";
		next NEXT_SEQ;
		}
		
		my $aa_seq = translate($dna_seq, 1); # translate sequence
		
		print DNA_UNALIGNED_COPY ">$taxon\n$dna_seq\n"; # write out
		print AA_UNLALIGNED ">$taxon\n$aa_seq\n"; # write out
		
		$translated_seq_num++;
	}
	close AA_UNLALIGNED;
	close DNA_UNALIGNED_COPY;
	
	if ($translated_seq_num == 0) {
	unlink $dna_unaligned_file_copy;
	unlink $aa_unaligned_file;
	next NEXT_FILE;
	}

	# align amino acid by mafft
	my $aa_aligned_file = "$aa_aligned_dir/$gene.aa_aln.fas";
	mafft($aa_unaligned_file, $aa_aligned_file);

	# translate back to nucleotide alignment
	my $dna_aligned_file = "$dna_aligned_dir/$file";
	translate_back($aa_aligned_file, $dna_unaligned_file_copy, $dna_aligned_file);

	unlink($dna_unaligned_file_copy);
	unlink($aa_unaligned_file);
}

}

# subroutine simply align sequences
sub normal_aln {
my $unaligned = shift;

NEXT_FILE: while (my $file = shift @$unaligned) {

	my $dna_unaligned_file = "$dna_unaligned_dir/$file"; # unaligned nucleotide
	my $dna_unaligned_nogap_file = "$dna_aligned_dir/$file.nogap"; # unaligned nucleotide
	my $dna_aligned_file = "$dna_aligned_dir/$file"; # aligned nucleotide

	my ($seqhash, $input_order, $aln_len) = readfas($dna_unaligned_file);

	my $seq_num = @$input_order;
	if ($seq_num == 0) {
	print STDERR "WARNING: Nothing found in $dna_unaligned_file. Skip this file.\n";
	next NEXT_FILE;
	}
	
	if ($seq_num == 1) {
	print STDERR "WARNING: Only 1 sequence found in $dna_unaligned_file. This file will be copied to $dna_aligned_dir\n";
	`cp $dna_unaligned_file $dna_aligned_dir/$file`;
	next NEXT_FILE;
	}
	
	open DNA_UNALIGNED_NOGAP, ">$dna_unaligned_nogap_file";
	NEXT_SEQ: foreach my $taxon (@$input_order) {
		my $dna_seq = $seqhash->{$taxon};
		if ($dna_seq =~ /-/) { # skip if sequences have been aligned
		print STDERR "WARNING: Sequences of $dna_unaligned_file have already be aligned. Gaps will be removed and sequences will be realigned\n";
		$dna_seq =~ s/-//g;
		}
		print DNA_UNALIGNED_NOGAP ">$taxon\n$dna_seq\n"; # write out
	}
	close DNA_UNALIGNED_NOGAP;

	mafft($dna_unaligned_nogap_file, $dna_aligned_file); # align by mafft
	unlink $dna_unaligned_nogap_file;
}

}

sub translate{
my $seq = shift @_;
my $start_codon = shift @_;
$start_codon --;
my $seq_obj=Bio::Seq->new(-seq=>$seq,-alphabet=>'dna');
my $pro=$seq_obj->translate(-frame => $start_codon, -terminator => "*", -unknown => "X");
$pro=$pro->seq; 
return($pro);
}

sub mafft {
my $unaligned = shift;
my $aligned = shift;

# align sequence by mafft
my $aln_temp_path = $aligned . ".temp";
`mafft --auto --preservecase --quiet --bl 80 "$unaligned" > "$aln_temp_path"`;

# unwrap the output
unwrap($aln_temp_path, $aligned);
unlink $aln_temp_path;
}

# subroutine to unwrap the output of mafft
sub unwrap {
my $in = shift;
my $unwrap = shift;

open IN, $in;
open UNWRAP, ">$unwrap";

my $seq;
my $line1 = <IN>;
my ($name_lag) = $line1 =~ />(\S+)/;
while (my $line = <IN>) {
chomp $line;
	if ($line =~ />(\S+)/) { # if find the header 
	print UNWRAP ">$name_lag\n$seq\n"; # print sequence
	$name_lag = $1;
	$seq = "";
	} else { # save sequence into $seq
	$seq .= $line;
	}
}
print UNWRAP ">$name_lag\n$seq\n"; # print last sequence

close IN;
close UNWRAP;
}

# subroutine to translate aligned amino acid back to nucleotide sequences
sub translate_back {
my $inpath_aa = shift; # aligned amino acid
my $inpath_nt = shift; # unaligned nucleotides
my $outfile = shift; # outfile name

# get nucleotide sequences and save to hash
my ($FILE_NT, $FILE_AA, $OUTFILE, %dna, %aa, @inputorder);
open ($FILE_NT, "<$inpath_nt") or die "Cannot find $inpath_nt ($!)";
while (my $line = <$FILE_NT>) {
	if ($line =~ />(\S+)/) {
	chomp(my $seq = <$FILE_NT>);
    $dna{$1} = $seq;
	}
}

# get amino acid sequences and save to hash
open ($FILE_AA, "<$inpath_aa") or die "Cannot find $inpath_aa ($!)"; # open aligned aa seq file
while (my $line = <$FILE_AA>) {
	if ($line =~ />(\S+)/) {
	push @inputorder, $1;
	chomp(my $seq = <$FILE_AA>);
    $aa{$1} = $seq;
	}
}
close $FILE_NT;
close $FILE_AA;	

open ($OUTFILE, ">$outfile") or die "can't open $outfile for writing!!!"; # create outfile

	foreach my $taxon (@inputorder){
	next unless ((exists $aa{$taxon})&&(exists $dna{$taxon})); # skip nucleotide sequences not exists in amino acid alignment
	
		# get amino acid and nucleotide sequences
		my $aln_aa_seq = $aa{$taxon};
		my $nucleo_seq = $dna{$taxon};
		$nucleo_seq =~ s/-//g;		
				
		print $OUTFILE ">$taxon\n";                                    
	
		my ($chars, @chars);
		
		my $pervious_char = substr $aln_aa_seq, 0, 1;
		$chars .= $pervious_char;
		for (my $i=1; $i<length($aln_aa_seq); $i++) {
		my $cur_char = substr $aln_aa_seq, $i, 1;
			if (($pervious_char ne "-")&&($cur_char eq "-")) {
			push @chars, $chars;
			$chars = "";
			} elsif (($pervious_char eq "-")&&($cur_char ne "-")) {
			push @chars, $chars;
			$chars = "";
			}
			$chars .= $cur_char;
			$pervious_char = $cur_char;
		}
		push @chars, $chars;
		
		my $aln_nucleo_seq;
		my $pointer = 0;
		foreach my $chars (@chars) {
			if ($chars =~ /-/) {
			my $gaps = $chars x 3;
			$aln_nucleo_seq .= $gaps;
			} else {
			my $length = length($chars)*3;
			my $sub_seq = substr $nucleo_seq, $pointer, $length;
			$aln_nucleo_seq .= $sub_seq;
			
			$pointer += $length;
			}
		}
			
		print $OUTFILE "$aln_nucleo_seq\n";		
		
	}
	
close $OUTFILE;
}

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
	if ($taxon_num > 1) {
	print STDERR "WARNING: Found $taxon_num \"$taxon\" in $file. Only first sequence named in \"$taxon\" will be kept\n";
	}
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

# subroutine to print out usage
sub usage {
print STDERR "
Script name: mafft_aln.pl

This is a script to either codon or nomarlly align nucleotide sequences in fasta format.

Dependencies: 
(1) BioPerl v1.007001 or higher
(2) Mafft v7.294b or higher (rename it as 'mafft' and put it under \$PATH)
(3) Perl module:
	1. Parallel::ForkManager

Example usage:
(1) Align sequences under 'nf' in codon. Write output to 'nf_aligned'. Run script with 4 process:

	perl mafft_aln.pl --dna_unaligned nf --dna_aligned nf_aligned --cpu 4

(2) Align sequences under 'nf' normally. Write output to 'nf_aligned'. Run script with 4 process:

	perl mafft_aln.pl --dna_unaligned nf --dna_aligned nf_aligned --non_codon_aln --cpu 4

(3) Align sequences under 'nf' in codon. Write output to 'nf_aligned'. Keep directory containing aligned amino acid sequences called 'aa_aligned'. Run script with 4 process:

	perl mafft_aln.pl --keep_aa_aln --dna_unaligned nf --dna_aligned nf_aligned --cpu 4
	
Input files:
(1) nf

Output files:
(1) nf_aligned
(2) aa_aligned (if --aa_aligned is specified)

Options:
--dna_unaligned
  Directory containing unaligned nucleotide sequences
--dna_aligned
  Directory containing aligned nucleotide sequences, named as 'xx_aligned' if this option is not specified
--non_codon_aln
  Do not align nucleotide sequences in codon. This option is turned off by default
--keep_aa_aln
  Keep directory containing aligned amino acid sequences named as 'aa_aligned', which is deleted by default
--cpu
  Limit the number of CPUs, 1 in default
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