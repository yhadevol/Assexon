#!/usr/bin/env perl

use Getopt::Long;
use warnings;
use strict;

my ($indir, $outdir, $essential_dir, $essential_taxa, $help);
my $min_seq_num = 2;

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'indir:s', \$indir,
                      'outdir:s', \$outdir,
                      'min_seq:i', \$min_seq_num,
                      'essential_dir:s', \$essential_dir,
                      'essential_taxa:s', \$essential_taxa,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/indir outdir/;
check_option(\@essential, \@ARGVoptions);

# get name of input dir
my @indir = $indir =~ /\S+/g;

foreach my $dir (@indir) { # check existence
	if (!(-e $dir)) {
	die "\nERROR: Cannot find input directory \"$dir\", please check --indir ($!)\n\n";
	}
}

my (@essential_taxa, %essential_taxa, $essential_taxa_num);
if ($essential_taxa) {
@essential_taxa = $essential_taxa =~ /\S+/g;
%essential_taxa = map {$_, ""} @essential_taxa;
$essential_taxa_num = scalar @essential_taxa;
}

my (@essential_dir, %essential_dir, $essential_dir_num);
if ($essential_dir) {
	@essential_dir = $essential_dir =~ /\S+/g;

	my @error_dir;
	my %indir = map {$_, ""} @indir;
	foreach my $essential_dir (@essential_dir) {
		if (! exists $indir{$essential_dir}) {
		push @error_dir, $essential_dir;
		}
	}
	if (@error_dir > 0) {
	my $error_dir = join " ", @error_dir;
	die "\nERROR: Unrecognized name of directory found (\"$error_dir\") after --essential_dir \n\n";
	}
	
	%essential_dir = map {$_, ""} @essential_dir;
	$essential_dir_num = scalar @essential_dir;
}

my %hash;
foreach my $dir (@indir) {
	opendir DIR, $dir;
	
	LOOP: while (my $infile = readdir DIR) {
		next if ($infile =~ /^\./);
		next if ($infile !~ /fa$|fas$|fasta$/i);
		
		my ($prefix, $suffix);
		if ($infile =~ /(\S+)\.(aln\.fa\S*$)/) {
		$prefix = $1;
		$suffix = $2;
		} elsif ($infile =~ /(\S+)\.(fa\S*$)/) {
		$prefix = $1;
		$suffix = $2;
		}
	
		# get name of each header
		my $path = "$dir/$infile";
		my ($seq_array, $input_order, $aln_len) = readfas($path);
		
		if (@$input_order == 0) {
		print STDERR "WARNING: Nothing found in $path. Skip this file. \n";
		next LOOP;
		}
		
		my ($taxon, $seq);
		for (my $i=0; $i<@$seq_array; $i++) {
			my $element = $seq_array->[$i];
			if ($i % 2 == 0) {
			$taxon = $element;
			$hash{$prefix}->{taxon}->{$taxon}++;
			push @{$hash{$prefix}->{element}->{$dir}}, ">$taxon";
			} else {
				$seq = $element;
				if ($seq =~ /-/) { # warning if found gap in sequence
				print STDERR "WARNING: Found gaps in $taxon of $path. Gaps will be removed. \n";
				$seq =~ s/-//g;
				}
				push @{$hash{$prefix}->{element}->{$dir}}, $seq;
			}
		}
				
		push @{$hash{$prefix}->{suffix}}, $suffix; # save suffix
	}
	closedir DIR;
}

# merge outfile
if (!(-e $outdir)) {
mkdir $outdir or die "\nERROR: Cannot create output directory \"$outdir\", please check --outdir ($!)\n\n";
}

foreach my $gene (sort keys %hash) {
	# set suffix of outfile
	my @suffix = @{$hash{$gene}->{suffix}};
	my %suffix = map {$_, ""} @suffix;
	my $suffix_num = scalar keys %suffix;
	
	my $suffix;
	# if ($suffix_num == 1) {
	# $suffix = $suffix[0];
	# } else {
	$suffix = "fas";
	# }
	
	# warning if the same taxa exist multiple time
	my %taxa_list = %{$hash{$gene}->{taxon}};
	foreach my $taxon (sort keys %taxa_list) {
		my $taxon_num = $taxa_list{$taxon};
		if ($taxon_num > 1) { # Warning if found same taxa
		print STDERR "WARNING: Found $taxon_num \"$taxon\" in loci \"$gene\". Please remove them later if you accidentally add it\n";
		}
	}
	
	# output
	my @output_order;
	
	my %dir = %{$hash{$gene}->{element}};
	foreach my $dir (@indir) {
		if (exists $dir{$dir}) {
		push @output_order, $dir;
		}
	}
	
	if ($essential_dir) {
	my $output_order_ref = essential_dir(\@output_order, \%dir);
	@output_order = @$output_order_ref;
	
	next if (@output_order == 0);
	}
	
	my $all_element = 0;
	my $essential_taxa_cnt = 0;
	foreach my $dir (sort keys %dir) {
		my @element = @{$dir{$dir}};
		$all_element += @element;
		
		if ($essential_taxa) {
			for (my $i=0; $i<@element; $i++) {
				if ($i % 2 == 0) {
				my ($taxon) = $element[$i] =~ />(\S+)/;
				$essential_taxa_cnt++ if (exists $essential_taxa{$taxon});
				$i++;
				}
			}
		}
	}
	next if (($essential_taxa)&&($essential_taxa_cnt < $essential_taxa_num));
	
	my $seq_num = $all_element/2;
	next if ($seq_num < $min_seq_num);

	output(\%dir, \@output_order, "$outdir/$gene.$suffix");
}

#####################################################
# Subroutines
#####################################################

# subroutine to find whether loci from essential dir exists
sub essential_dir {
my $output_order = shift;
my $dir = shift;

my $essential_dir_cnt = 0;
foreach my $essential (sort keys %$dir) {
$essential_dir_cnt++ if (exists $essential_dir{$essential});
}

my @essential_output_order;
if ($essential_dir_cnt == $essential_dir_num) {
@essential_output_order = @$output_order;
} else {
@essential_output_order = ();
}

return (\@essential_output_order);
}

# subroutine to write output
sub output {
my $element = shift;
my $input_order = shift;
my $outpath = shift;

open OUTFILE, ">$outpath";

my %element = %$element;
foreach my $dir (@$input_order) {
	my @element = @{$element{$dir}};
	foreach my $element (@element) {
	print OUTFILE "$element\n";
	}
}
close OUTFILE;
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

# subroutine to reformat input fasta file, permit multiple same taxa 
sub readfas {
my $file = shift;

open INFILE, $file or die "\nERROR: Cannot find input file \"$file\" ($!)\n\n";

my (@seq, %seqlength, @input_order, $seq, $taxon, $lasttaxon);

while (my $line = <INFILE>) {

	# remove enter
	$line =~ s/\r//g;
	chomp $line;
	
	# if find >
	if ($line =~ /^>(\S+)/) {
		$taxon = $1;
		
		# if find previous taxon name, and it hasn't recorded before
		if ($lasttaxon) {
			my $length = length($seq);
			if ($length > 0) {
			push @input_order, $lasttaxon;
			push @seq, $lasttaxon;
			push @seq, $seq;
			$seqlength{length($seq)} = "";
			} else {
			print STDERR "WARNING: No nucleotide found sequence of $lasttaxon in $file. This taxon will be discarded.\n";			
			}
		}
		
		$seq = "";
		$lasttaxon = $taxon;
	} else {
	$seq .= $line;
	}
}

if ($lasttaxon) {
	my $length = length($seq);
	if ($length > 0) {
	push @input_order, $lasttaxon;
	push @seq, $lasttaxon;
	push @seq, $seq;
	$seqlength{length($seq)} = "";
	} else {
	print STDERR "WARNING: No nucleotide found sequence of $lasttaxon in $file. This taxon will be discarded.\n";			
	}
}

close INFILE;

my @seqlength = sort keys %seqlength;

my $aln_len;
if (@seqlength == 1) {
$aln_len = $seqlength[0];
} else {
$aln_len = "na";
}

return (\@seq, \@input_order, $aln_len);
}

sub usage {
print STDERR "
Script name: merge_loci.pl

This is a script to merge unaligned sequences in fasta format under several directories from the same loci.

Example usage:
(1) Merge fasta files with the same prefix from 'dir1' 'dir2'. Only merged files with at least 2 sequences will be written to 'outdir':

	perl merge_loci.pl --indir 'dir1 dir2' --outdir outdir --min_seq 2

Input files:
(1) dir1
(2) dir2

Output files:
(1) outdir

Options:
--indir
  Space delimited list of directories containing sequences
--outdir
  Directory containing merged gene files
--essential_dir 
  Space delimited list of directories. Only output loci if it contains sequences from these directory
--essential_taxa
  Space delimited list of taxa. Only output loci have these taxa
--min_seq
  Minimum sequences required in merged file, $min_seq_num in default
--help , -h
  Show this help message and exit

Author: Hao Yuan                                                                     
        Shanghai Ocean University                                               
        Shanghai, China, 201306                                                               
                                                                                         
Created by: June 27, 2018                                                              
                                                                                         
Last modified by: 2020.7.7

The output files are only suffixed with \"fas\", instead of identifying suffix from input files

";
exit;
}