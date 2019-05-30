#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long; 

# input dir of reads, output file of statistics
my ($dir, $help);
my $outfile = "reads_bases_cnt.txt";

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'indir:s', \$dir,
					  'outfile:s', \$outfile,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";
                       
# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check essential options
my @essential = qw/indir/;
check_option(\@essential, \@ARGVoptions);

# get name of input reads
opendir DIR, $dir or die "\nERROR: Cannot find input directory \"$dir\", please check --indir ($!)\n\n";
my @reads = grep {$_ =~ /fq$|fastq$/} readdir DIR;
closedir DIR;

# write output
open OUTFILE, ">$outfile" or die "\nERROR: Cannot write output to \"$outfile\", please check --outfile ($!)\n\n";

# header
print OUTFILE "Sample\tNumber of reads\tNumber of bases(bp)\n";

my %counted;
foreach my $reads (@reads) {
	next if (exists $counted{$reads}); # skip read reads
	
	my $suffix;
	if ($reads =~ /\S+\.fq/) {
	$suffix = "fq";
	} elsif ($reads =~ /\S+\.fastq/) {
	$suffix = "fastq";
	}
	
	# base and reads counter
	my $base_cnt = 0;
	my $reads_cnt = 0;
	
	my $prefix;
	if ($reads =~ /(\S+)\_R[1|2]\.$suffix/i) { # if reads are paired
		$prefix = $1;
		
		# R1 and R2
		my $path_R1 = "$dir/${prefix}_R1.$suffix";
		my $path_R2 = "$dir/${prefix}_R2.$suffix";
		
		# check existence
		if (!(-e $path_R1)) {
	    die "\nERROR: Cannot find R1 of $prefix \"$path_R1\" ($!)\n\n";
		}
		
		if (!(-e $path_R2)) {
		die "\nERROR: Cannot find R2 of $prefix \"$path_R1\" ($!)\n\n";
		}
		
		# count bases and reads of R1
		my ($R1_bases, $R1_reads) = reads_bases_cnt($path_R1);
		$base_cnt += $R1_bases;
		$reads_cnt += $R1_reads;
		
		# count bases and reads of R2
		my ($R2_bases, $R2_reads) = reads_bases_cnt($path_R2);
		$base_cnt += $R2_bases;
		$reads_cnt += $R2_reads;
		
		# save read sample
		$counted{"${prefix}_R1.$suffix"} = "";
		$counted{"${prefix}_R2.$suffix"} = "";
	} elsif ($reads =~ /(\S+)\.$suffix/i) {
		$prefix = $1;
		
		my $path = "$dir/${prefix}.$suffix";
		
		# check existence
		if (!(-e $path)) {
		die "\nERROR: Cannot find R2 of $prefix \"$path\" ($!)\n\n";
		}
		
		# count bases and reads of current sample
		my ($bases, $reads) = reads_bases_cnt($path);
		$base_cnt += $bases;
		$reads_cnt += $reads;	
		
		# save read sample
		$counted{$reads} = "";
	}
	
	# print output
	print OUTFILE "$prefix\t$reads_cnt\t$base_cnt\n";	
}
close OUTFILE;

#####################################################
# Subroutines
#####################################################

# subroutine to count reads and bases
sub reads_bases_cnt {
my $file = shift; # reads file

my $base_cnt = 0;
my $read_cnt = 0;

# count reads and bases
open READS, $file;
while (my $line = <READS>){
	if ($line =~ /^\@\S+\s*\S*/){ # beginning of the reads
	chomp(my $seq = <READS>); # seq
	<READS>; # plus
	<READS>; # quality
	$base_cnt += length($seq); # add up bases
	$read_cnt ++; # add up reads
	}
}
close READS;

# return number of bases and reads
return ($base_cnt, $read_cnt);
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
Script name: count_reads_bases.pl

This is a script to count number of reads and bases pairs for each file in fastq format.

Example usage:
(1) Count the number of bases of fastq under 'raw_reads'. Write output to 'reads_bases_cnt.txt'

	perl count_bases.pl --indir	raw_reads --outfile reads_bases_cnt.txt
	
Input files:
(1) raw_reads

Output files:
(1) reads_bases_cnt.txt

Options:
--indir
  Directory containing reads in fastq format. Either paired or single file can be counted. Suffix of file must be fq or fastq, and paired files must named as xxx_R1.fq and xxx_R2.fq.
--outfile
  Tab delimited table recording number of bases and reads in each paired or single file
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