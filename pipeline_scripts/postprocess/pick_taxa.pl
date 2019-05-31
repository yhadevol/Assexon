#!/usr/bin/env perl

use Getopt::Long;
use warnings;
use strict;

my ($dir, $outdir, $essential_taxa, $deselected_specieslist, $selected_specieslist, $rm_dup_taxa, $help);

my $minnumseq = 2; #the minimum number of sequences for keeping a locus

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'indir:s', \$dir,
                      'outdir:s', \$outdir,
                      'min_seq:s', \$minnumseq,
                      'essential_taxa:s', \$essential_taxa,
                      'selected_taxa:s', \$selected_specieslist,
                      'deselected_taxa:s', \$deselected_specieslist,
                      'rm_dup_taxa!', \$rm_dup_taxa,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/indir outdir/;
check_option(\@essential, \@ARGVoptions);

# check dir
if (!(-e $dir)) {
die "\nERROR: Cannot find input directory \"$dir\", please check --indir ($!)\n\n";
}

# check co-existence of --deselected_species and --selected_species
if (($selected_specieslist)&&($deselected_specieslist)) {
die "\nERROR: Cannot specify --deselected_species and --selected_species at the same time\n\n";
}

# check co-existence of --essential_taxa with --deselected_taxa or --selected_species
if (($essential_taxa)&&(($selected_specieslist)||($deselected_specieslist))) {
die "\nERROR: Cannot specify --essential_taxa with either --deselected_taxa or --selected_species at the same time\n\n";
}
 
if (!(-e $outdir)) {                
mkdir $outdir or die "\nERROR: Cannot create output directory \"$outdir\", please check --outdir ($!)\n\n";
}

# get sample names of --selected_species
my (@selected_list, %selected_list);
if ($selected_specieslist) {
@selected_list = $selected_specieslist =~ /(\S+)/g;
%selected_list = map {$_, ""} @selected_list;
}

# get sample names of --deselected_species
my (@deselected_list, %deselected_list);
if ($deselected_specieslist) {
@deselected_list = $deselected_specieslist =~ /(\S+)/g;
%deselected_list = map {$_, ""} @deselected_list;
}

my (@essential_taxa, %essential_taxa, $essential_taxa_num);
if ($essential_taxa) {
@essential_taxa = $essential_taxa =~ /(\S+)/g;
%essential_taxa = map {$_, ""} @essential_taxa;
$essential_taxa_num = scalar @essential_taxa;
}

opendir DIR, $dir;

LOOP: while (my $file = readdir(DIR)) { #read all fasta files under the folder
	next if ($file =~ /^\./);
	next if ($file !~ /fa$|fas$|fasta$/i); #skip files beginning with.
	
	my $infilename = "$dir/$file"; #get the name of fasta file
    
    my ($seq_array, $input_order, $aln_len) = readfas($infilename);

	if (@$input_order == 0) {
	print STDERR "WARNING: Nothing found in $infilename. Skip this file\n";
	next LOOP;
	}
    
    my $countseq = 0;
	my (@seq, %taxa_list);
    if ($selected_specieslist) { # if --selected_species is specified, keep samples in the list
    	for (my $i=0; $i<@$seq_array; $i++) {
			if ($i % 2 == 0) {
				my $taxon = $seq_array->[$i];
				
				if (exists $selected_list{$taxon}) { # keep samples in the list
					$i++;
					my $seq = $seq_array->[$i];
					
					if ($seq =~ /-/) {
					print STDERR "WARNING: Found gaps in $taxon of $infilename. Gaps will be removed. \n";
					$seq =~ s/-//g;
					}
					
					push @seq, ">$taxon";
					push @seq, "$seq";
					
					$taxa_list{$taxon}++;
				}
			}
    	}
    } elsif ($deselected_specieslist) { # if --deselected_species is specified, discard samples in the list
		for (my $i=0; $i<@$seq_array; $i++) {
			if ($i % 2 == 0) {
				my $taxon = $seq_array->[$i];
				
				if (!(exists $deselected_list{$taxon})) { # keep samples in the list
					$i++;
					my $seq = $seq_array->[$i];
					
					if ($seq =~ /-/) {
					print STDERR "WARNING: Found gaps in $taxon of $infilename. Gaps will be removed. \n";
					$seq =~ s/-//g;
					}
					
					push @seq, ">$taxon";
					push @seq, "$seq";
					
					$taxa_list{$taxon}++;
				}
			}
    	}
    } else { # if --selected_species and --deselected_species are both not specified 
    	for (my $i=0; $i<@$seq_array; $i++) {
			if ($i % 2 == 0) {
				my $taxon = $seq_array->[$i];
				
				$i++;
				my $seq = $seq_array->[$i];
				
				if ($seq =~ /-/) {
				print STDERR "WARNING: Found gaps in $taxon of $infilename. Gaps will be removed. \n";
				$seq =~ s/-//g;
				}
				
				push @seq, ">$taxon";
				push @seq, "$seq";
				
				$taxa_list{$taxon}++;
			}
    	}
    }
    
    if (! $rm_dup_taxa) {
		foreach my $taxon (sort keys %taxa_list) {
			my $taxon_num = $taxa_list{$taxon};
			if ($taxon_num > 1) { # Warning if found same taxa
			print STDERR "WARNING: Found $taxon_num \"$taxon\" in selected taxa of $infilename. Please remove them using option \"--rm_dup_taxa\" if you accidentally add it\n";
			}
		}
    }
    
    # print sequences with ay least $minnumseq
    my $seq_num = @seq/2;
    if ($seq_num >= $minnumseq) {
    	my $essential_cnt = 0;
		my $outfile = "$outdir/$file";
		
		open OUTFILE, ">$outfile";
		
		if ($rm_dup_taxa) {
			my %taxa_cnt;
			for (my $i=0; $i<@seq; $i++) {
				my $line = $seq[$i];
				if ($i%2 == 0) {
				my ($taxon) = $line =~ />(\S+)/;
				$essential_cnt++ if (exists $essential_taxa{$taxon});
				$taxa_cnt{$taxon}++;
				print OUTFILE "$line\n" if ($taxa_cnt{$taxon} == 1);
				} else {
				my ($taxon) = $seq[$i-1] =~ />(\S+)/;
				print OUTFILE "$line\n" if ($taxa_cnt{$taxon} == 1);
				}
			}
		} else {
			for (my $i=0; $i<@seq; $i++) {
				my $line = $seq[$i];
				if ($i%2 == 0) {
				my ($taxon) = $line =~ />(\S+)/;
				$essential_cnt++ if (exists $essential_taxa{$taxon});
				} 
				print OUTFILE "$line\n";
			}
		}
		
		close OUTFILE;
		
		unlink $outfile if (($essential_taxa)&&($essential_cnt < $essential_taxa_num));
    }
}
closedir DIR;

#####################################################
# Subroutines
#####################################################

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
Script name: pick_taxa.pl

This is a script to pick out needed taxa or discard unneeded taxa from unaligned sequences in fasta format.

Example usage:
(1) Preserve the files in 'nf' with at least 2 sequences. Write filtered files to 'selected':

	perl pick_taxa.pl --indir nf --outdir selected --min_seq 2

(2) Only keep 'taxon1', 'taxon2' and 'taxon3' for each file in 'nf'. Write filtered files to 'selected':

	perl pick_taxa.pl --indir nf --outdir selected --selected_taxa 'taxon1 taxon2 taxon3'

(3) Discard 'taxon1', 'taxon2' and 'taxon3' for each file in 'nf'. Write filtered files to 'selected':

	perl pick_taxa.pl --indir nf --outdir selected --deselected_taxa 'taxon1 taxon2 taxon3'

(4) Discard the file do not have 'taxon1' and 'taxon2' in 'nf'. Write filtered files to 'selected':

	perl pick_taxa.pl --indir nf --outdir selected --essential_taxa 'taxon1 taxon2'

(5) Remove sequences of duplicated taxa. Write deduplicated files to 'deduplicated':

	perl pick_taxa.pl --indir nf --outdir deduplicated --rm_dup_taxa


Input files:
(1) nf

Output files:
(1) selected

Options:
--indir
  Directory containing unaligned sequences
--outdir
  Directory containing sequences of selected taxon
--min_seq
  Keep the file with at least --min_seq sequences, $minnumseq in default
--essential_taxa
  List of taxa must exist in the loci, discard the loci if there are no these taxa in it,  each taxon is delimited by space
--selected_taxa
  List of taxa want to keep, each taxon is delimited by space
--deselected_taxa
  List of taxa want to be discarded, each taxon is delimited by space
--rm_dup_taxa
  Remove sequences of duplicated taxa
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
