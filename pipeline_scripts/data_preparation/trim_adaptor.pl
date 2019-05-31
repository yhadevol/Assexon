#!/usr/bin/env perl

use Getopt::Long; 
use warnings;
use strict;
use Parallel::ForkManager;

# dir of raw reads, dir of trimmed reads, R1 adaptor, R2 adaptor,path to cutadapt
my ($dir, $outdir, $oria, $oria2, $cutadapt, $help); 
$oria = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
$oria2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
my $reads_bases_cnt = "trimmed_reads_bases_count.txt"; # file summarized with number of bp and reads in raw and trimmed reads
my $process = 1; # number of process used in default

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'raw_reads:s', \$dir,
                      'adapter_R1:s', \$oria,
                      'adapter_R2:s', \$oria2,
                      'cutadapt_path:s', \$cutadapt,
                      'trimmed:s', \$outdir,
                      'cpu:s', \$process,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential;
@essential = qw/raw_reads trimmed/;
check_option(\@essential, \@ARGVoptions);

# check whether trim_galore and cutadapt properly installed
`trim_galore -v` or die "\nERROR: Cannot find trim_galore under \$PATH\n\n";
`cutadapt --version` or die "\nERROR: Cannot find cutadapt under \$PATH\n\n" if (! $cutadapt);

# initiate multi-process
my $pm = Parallel::ForkManager->new($process);

# get the name of read R1
opendir DIR, $dir or die ("\nERROR: Cannot find input directory \"$dir\", please check --raw_reads ($!)\n\n");
my @reads = grep {$_ !~ /^\./} readdir DIR;
@reads = grep {$_ =~ /\S+_R1.fq/} @reads;
closedir DIR;

# mkdir if output dir don't exist
if (!(-e $outdir)){
mkdir $outdir or die "\nERROR: Cannot create output directory \"$outdir\", please check --trimmed ($!)\n\n";
}

# trim adaptor
adapt($dir, $outdir, $cutadapt, $oria, $oria2);

# mv log to trimming report
if (!(-e "trimming_report")){
mkdir "trimming_report" or die "\nERROR: Cannot create output directory of trimming_report \"trimming_report\"\n\n";
}
`mv -f ./$outdir/*trimming_report* ./trimming_report`;

# rename outfile
my @outfile = glob("$outdir/*");
map {
my $ori = $_;
$_ =~ s/_val_\d//;
`mv $ori $_`;
} @outfile;

# summary sample statistics from trimming_report
open READS_BASES_COUNT, ">$reads_bases_cnt";

print READS_BASES_COUNT "Sample\tReads num. of raw reads\tBases num. of raw reads (bp)\tReads num. of cleaned reads\tBases num. of cleaned reads (bp)\n";

# get name of R1 reads in output dir
opendir OUTDIR, $outdir;
my @outdir_reads = grep {$_ !~ /^\./} readdir OUTDIR;
@outdir_reads = grep {$_ =~ /\S+_R1.fq/} @outdir_reads;
closedir OUTDIR;

# summary statistics for each reads file
foreach my $R1 (@outdir_reads) {
my ($prefix) = $R1 =~ /(\S+)\_R1.fq/;

# summary number and bp of raw reads from trimming report
my $trimming_report_R1 = "trimming_report/${prefix}_R1.fq_trimming_report.txt";
my $trimming_report_R2 = "trimming_report/${prefix}_R2.fq_trimming_report.txt";

my ($R1_all_bases, $R1_all_reads) = all_reads_bases($trimming_report_R1);
my ($R2_all_bases, $R2_all_reads) = all_reads_bases($trimming_report_R2);

# summary number and bp of trimmed reads from trimmed reads
my $path_R1 = "$outdir/${prefix}_R1.fq";
my $path_R2 = "$outdir/${prefix}_R2.fq";

my ($R1_trimmed_bases, $R1_trimmed_reads) = trimmed_reads_bases($path_R1);
my ($R2_trimmed_bases, $R2_trimmed_reads) = trimmed_reads_bases($path_R2);

# sum up number of read and bases
my $all_bases = $R1_all_bases+$R2_all_bases;
my $all_reads = $R1_all_reads+$R2_all_reads;

my $trimmed_bases = $R1_trimmed_bases+$R2_trimmed_bases;
my $trimmed_reads = $R1_trimmed_reads+$R2_trimmed_reads;

print READS_BASES_COUNT "$prefix\t$all_reads\t$all_bases\t$trimmed_reads\t$trimmed_bases\n";
}
close READS_BASES_COUNT;

######################################################
# Subroutines
######################################################

# subroutine to cut adaptor
sub adapt {
	my $dir = shift; # dir containing raw reads
	my $outdir = shift; # output dir
	my $cutadapt = shift; # path to cutadapt
	my $a = shift; # R1 adaptor
	my $a2 = shift; # R2 adaptor
		
		# cut adaptor samples in parallel
		DATA_LOOP: foreach my $read (@reads) {
		$pm->start and next DATA_LOOP;
		my ($speciesname) = $read =~ /(\S+)_R1.fq/;
		
		say STDOUT "Start to trim adaptors and low quality bases in reads of $speciesname";

			if ($cutadapt) { # specify --path_to_cutadapt if have $cutadapt
			`trim_galore -a $a -a2 $a2 --paired $dir/${speciesname}_R1.fq $dir/${speciesname}_R2.fq --path_to_cutadapt $cutadapt -o $outdir`;
			} else {
			`trim_galore -a $a -a2 $a2 --paired $dir/${speciesname}_R1.fq $dir/${speciesname}_R2.fq -o $outdir`;
			}
		
		say STDOUT "Reads of $speciesname has been trimmed";
		
		$pm->finish;
		}
		$pm->wait_all_children();
}

# count number of reads and bases from trimming report
sub all_reads_bases {
my $trimming_report = shift; # trimming report

my ($all_reads, $all_bases);
open TRIMMING_REPORT, $trimming_report;
while (my $all_reads_line = <TRIMMING_REPORT>) {
	if ($all_reads_line =~ /Processed\sreads:/) { # find the line with "Processed reads:"
	my @all_reads = $all_reads_line =~ /\S+/g;
	$all_reads = $all_reads[2]; # the second char is the number of reads
	
	my $all_bases_line = <TRIMMING_REPORT>; # next line havinh number of bases
	my @all_bases = $all_bases_line =~ /\S+/g;
	$all_bases = $all_bases[2]; # the second char is the number of bases
	
	last;
	}
}
close TRIMMING_REPORT;

# return number of bases and reads
return($all_bases, $all_reads);
}

# count number of reads and bases from trimmed reads
sub trimmed_reads_bases {
my $file = shift; # trimmed reads file

# initial number of read and base
my $base_cnt = 0;
my $read_cnt = 0;

# count number of reads and bases
open READS, $file;
while (my $line = <READS>){
	if ($line =~ /^\@\S+\s*\S*/){ # if find @, which is the beginning line of a read
	chomp(my $seq = <READS>); # sequence
	<READS>; # plus
	<READS>; # quality
	$base_cnt += length($seq); # add number of bases
	$read_cnt ++; # add number of reads
	}
}
close READS;

# return number of bases and reads
return ($base_cnt, $read_cnt);
}

# Subroutine to check missing option
sub check_option { 
my $essential = shift; # mandatory options
my $ARGVoptions = shift; # input options

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
print STDERR "
Script name: trim_adaptor.pl

This is a script to trim low quality bases and adaptor from 3\' end of reads

Dependencies:
(1) trim_galore v0.4.1 or higher
(2) cutadapt v1.2.1 or higher

Example usage:

	perl trim_adaptor.pl --raw_reads raw_reads --trimmed trimmed	

Input files: 
(1) raw_reads

Output files:
(1) trimmed
(2) trimmed_reads_bases_count.txt (file summarized number of reads and bases in raw and trimmed reads)

Options:
--raw_reads
  Directory containing expanded raw reads
--trimmed
  Directory containing cleaned reads
--adapter_R1
  Adaptor to be trimmed on R1 side, 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' in default
--adapter_R2
  Adaptor to be trimmed on R2 side, 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT' in default
--cutadapt_path
  Specify a path to the Cutadapt executable, else it is assumed that Cutadapt is in \$PATH
--cpu
  Limit the number of CPUs, $process in default
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