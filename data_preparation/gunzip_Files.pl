#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use Parallel::ForkManager;

# dir with gzipped raw reads, output dir
my ($dirname, $outdir, $help);
my $process = 1; # number of process in default

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'gzip:s', \$dirname,
                      'gunzipped:s', \$outdir,
                      'cpu:s', \$process,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/gzip gunzipped/;
check_option(\@essential, \@ARGVoptions);

# initiate multi-process 
my $pm = Parallel::ForkManager->new($process); 

# open dir of gzipped raw reads
opendir(my $DIR1, $dirname) or die "ERROR: Cannot find gzipped directory \"$dirname\", please check --gzip ($!)";

# create outdir if it doesn't exist
if (!(-e $outdir)) {
mkdir $outdir or die "ERROR: Cannot create ouput directory \"$outdir\", please check --gunzipped ($!)";
}

# expand each sample in parallel
DATA_LOOP: while (my $foldname = readdir($DIR1)){
next if ($foldname =~ /^\./);
$pm->start and next DATA_LOOP;
	
	# open subdir of sample
	my $subfolder = "$dirname/$foldname";
    opendir(my $DIR2, $subfolder) or die "ERROR: Cannot find $subfolder";
    
    # only get name of gunzipped reads
    my @gz = grep {$_ =~ /gz$/} readdir($DIR2);
    closedir($DIR2);
    
    # expand reads
    if (@gz % 2 == 0) { # if reads are paired
    	print "Start expanding $foldname\n";
		foreach my $gz (@gz) {
		my ($prefix) = $gz =~ /(\S+).gz/;
		`gunzip -cf $subfolder/$prefix.gz > $outdir/$prefix`;
		}
		print "$foldname have been expanded\n";
    } else { # if reads are not paired
	print STDERR "Reads under $subfolder is not paired\n";    
    }
    
$pm->finish;
}
$pm->wait_all_children();
closedir($DIR1);

######################################################
# Subroutines
######################################################

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
Script name: gunzip_Files.pl

This is a script to gunzip raw data(*.gz)

Example usage:

	perl gunzip_Files.pl --gzip gzipped_raw_reads --gunzipped unzipped_raw_data

Input files: 
(1) gzipped_raw_reads

Output files:
(1) unzipped_raw_data

Options:
--gzip
  Directory containing gzipped raw data, the structure of directory must be:
  
  gzipped_raw_reads
      | 
      └── test1
      |     ├── test1_R1.fq.gz
      |     └── test1_R2.fq.gz
      └── test2
            ├── test2_R1.fq.gz
            └── test2_R2.fq.gz
  
--gunzipped
  Directory containing expanded raw data
--cpu
  Limit the number of CPUs. $process in default
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

