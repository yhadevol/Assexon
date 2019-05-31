#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my ($infile, $outfile, $help);

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'infile:s', \$infile,
					  'outfile:s', \$outfile,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/infile outfile/;
check_option(\@essential, \@ARGVoptions);

# check input and output name
if ($infile eq $outfile) {
die "\nERROR: Name of input and output file must be different\n\n";
}

if ($infile !~ /fa$|fas$|fasta$/i) {
die "\nERROR: Input sequences must be in fasta format\n\n";
}

my $unixfile = "$infile.unix";
ConvertToUnix($infile, $unixfile); # substitute line break
unwrap($unixfile, $outfile); # unwrap sequences
unlink $unixfile;

#####################################################
# Subroutines
#####################################################

# subroutine to change line breaks into unix style
sub ConvertToUnix {
my ($infile, $unixFile) = @_;
open (IN, $infile) or die "\nERROR: Cannot find input fasta file \"$infile\", please check --infile ($!)\n\n";
open (OUT, ">$unixFile") or die "\nERROR: Cannot write file with unix style line break \"$unixFile\" ($!)\n\n";
my @buffer = <IN>;
close IN;
my $line = "";
foreach my $element (@buffer) {
	$line .= $element;
}
if ($line =~ /\r\n/) {
	$line =~ s/\r//g;
}elsif ($line =~ /\r/) {
	$line =~ s/\r/\n/g;
}
print OUT $line;	
close OUT;	
}

# subroutine to unwrap folded sequences
sub unwrap {
my $in = shift;
my $unwrap = shift;

open IN, $in or die "\nERROR: Cannot find file with unix style line break \"$in\" ($!)\n\n";
open UNWRAP, ">$unwrap" or die "\nERROR: Cannot find write output \"$outfile\", please check --outfile ($!)\n\n";

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
Script name: unixlb_unwarp.pl

This is a script to change link break into unix(LF) and unwrap folded sequences in fasta format.

Example usage:
(1) Substitute link break as unix(LF) and unwrap folded sequences in 'locus.fas'. Write reformatted file to 'locus.unwrap.fas'

	perl unixlb_unwarp.pl --infile locus.fas --outfile locus.unwrap.fas

Input files:
(1) locus.fas

Output files:
(1) locus.unwrap.fas

Options:
--infile
  Input fasta file
--outfile
  Line break substituted and unfolded outfile
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