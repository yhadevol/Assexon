#!/usr/bin/env perl

use Getopt::Long; 
use warnings;
use strict;

my ($bam_dir, $help);

my $outfile = "map_statistics.txt";
# my $max_dep = 8000; # max depth allowed for a position, this option is disabled in v0.1.19
my $dep_dir = "dep_dir";

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'marked_bam:s', \$bam_dir,
					  'outfile:s', \$outfile,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check essential options
my @essential = qw/marked_bam/;
check_option(\@essential, \@ARGVoptions);

# check samtools
`samtools 2> samtools_check`;
my $samtools_check = `grep "command not found" samtools_check`;
unlink "samtools_check";
if ($samtools_check) {
die "\nERROR: Cannot find samtools in \$PATH \n\n";
}

# get file name
opendir BAM, $bam_dir or die "\nERROR: Cannot find input directory \"$bam_dir\", please check --bam ($!)\n\n";
my @bam = grep {$_ =~ /marked\.bam$/} readdir BAM;
closedir BAM;

mkdir $dep_dir;

open OUTFILE, ">$outfile" or die "\nERROR: Cannot write output to \"$outfile\", please check --outfile ($!)\n\n";

# header
print OUTFILE "Sample\tTotal reads\tMapped reads\tMapped PCR duplicates\tUniquely mapped reads\tOn-target reads(\%)\tSeq. depth\n";

# summary statistics
foreach my $bam (@bam) {
	my ($prefix) = $bam =~ /(.+)\.marked.bam$/;

	my $path = "$bam_dir/$bam";

	# calculate depth
	my $depth_file = "$dep_dir/$prefix.dep";
	`samtools depth $path > $depth_file`;
	
	# average depth
	my $pos_num = 0;
	my $total_depth = 0;
	open DEPTH, $depth_file;
	while (my $line = <DEPTH>) {
	my @line = $line =~ /\S+/g;
	$total_depth += $line[-1];
	$pos_num++;
	}
	close DEPTH;
	
	if ($pos_num == 0) {
	say STDERR "WARNING: No reads aligned in $path";
	next;
	}
	
	my $ave_dep = sprintf("%.2f", $total_depth/$pos_num);
	
	# use flagstat total, duplicated and mapped reads
	my $flagstat = `samtools flagstat $path`;
	
	# total reads
	my ($tot1, $tot2) = $flagstat =~ /(\d+)\s\+\s(\d+)\sin\stotal/;
	my $tot = $tot1+$tot2;

	# duplicated reads
	my ($dup1, $dup2) = $flagstat =~ /(\d+)\s\+\s(\d+)\sduplicates/;
	my $dup = $dup1+$dup2;
	
	# mapped reads
	my ($map1, $map2) = $flagstat =~ /(\d+)\s\+\s(\d+)\smapped/;
	my $map = $map1+$map2;
	
	# uniquely mapped 
	my $uniq_map = $map-$dup; 
	
	# on target rate
	my $on_target = sprintf("%.2f", $uniq_map/$tot*100);

	print OUTFILE "$prefix\t$tot\t$map\t$dup\t$uniq_map\t$on_target\t$ave_dep\n";
}

close OUTFILE;

`rm -rf $dep_dir`;

#####################################################
# Subroutines
#####################################################

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

sub usage {
say STDERR 
"Script name: map_statistics.pl

This is a script to summary statistics from bam file after marking duplicates (from pipeline of GATK), including:
(1) Number of total reads
(2) Number of mapped reads
(3) Number of mapped duplicated reads (JUST MAPPED PCR DUPLICATES. NOT ALL PCR DUPLICATES!)
(4) Number of uniquely mapped reads (mapped reads without duplicates)
(5) Divide number of uniquely mapped reads by number of all reads
(6) Average sequencing depth

Dependencies:
(1) samtools v0.1.19 or higher

Example usage:
(1) Summary statistics from 'xxx.marked.bam' files under 'bam_dir'. Write output to 'map_statistics.txt':
	
	perl map_statistics.pl --marked_bam bam_dir

Input:
(1) bam_dir

Output:
(1) map_statistics.txt

Option:
--marked_bam
  Input directory containing duplicates-marked bam file, bam file must named as 'xxx.smarked.bam'. xxx is name of sample
--outfile
  Tab delimited tab of summarized statistics for each sample, named in map_statistics.txt in default
--help, h
  Show this help message and exit

Author: Hao Yuan                                                                     
        Shanghai Ocean University                                               
        Shanghai, China, 201306                                                               
                                                                                         
Created by: June 27, 2018                                                              
                                                                                         
Last modified by:
";
exit;
}