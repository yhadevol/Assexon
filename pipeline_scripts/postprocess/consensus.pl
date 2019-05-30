#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

# input dir of aligned loci, majority consensus reference
my ($dir, $outfile, $help);

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'indir:s', \$dir,
                      'outfile:s', \$outfile,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/indir outfile/;
check_option(\@essential, \@ARGVoptions);

# open input dir
opendir (DIR, $dir) or die "\nERROR: Cannot find input directory \"$dir\", please check --indir ($!)\n\n";

# write outfile
open (my $OUT_FILE, ">$outfile") or die "\nERROR: Cannot write consensus sequences \"$outfile\", please check --outfile\n\n";

my ($IN_FILE, $infilename);
LOOP: while (my $file = readdir DIR) {
    next if ($file =~ /^\./);
    next if ($file !~ /fa$|fas$|fasta$/i);
	
	# get the name of fasta file
    $infilename = "$dir/$file"; 
    
    # check input file is wrongly formatted or have any strange characters 
    my ($seqhash, $input_order, $aln_len) = readfas($infilename);
    
    # skip if no taxa in it
    if (@$input_order == 0) {
    print STDERR "WARNING: Nothing is found in $infilename. Skip this file\n";
	next LOOP;
    }
    
    # skip if loci are not aligned
    if ($aln_len eq "na") {
    die "\nERROR: Number of characters among sequences are not the same in $infilename. It may be not aligned\n\n";
    }
    
    # split sequences into array of single nucleotide
    my %sequence;
	foreach my $taxon (sort keys %$seqhash) {
	my @ac = split//, $seqhash->{$taxon};
    push @{$sequence{$taxon}->{ac}}, @ac;
	}
    
    # print gene name as header
    my ($genename) = $file =~ /^(\S+)\.fa\S*$/;
    print $OUT_FILE ">$genename\n"; #print the gene name

    my @keys = keys %sequence;
    my @sortedkeys = sort @keys;
    
    # get majority consensus nucleotide for each column and print it to outfile sequentially
    for (my $i = 0 ; $i < $aln_len ; ++$i) {
        my $conchar;
        my @chars; # all nucleotide in current column
        foreach my $key (@sortedkeys) {        
        push @chars, @{$sequence{$key}->{ac}}[$i];
        }
    
        $conchar = consensus(\@chars); # majority consensus nucleotide
        print $OUT_FILE "$conchar"; # print it to outfile
    }
    
    print $OUT_FILE "\n";
}

closedir (DIR);
close ($OUT_FILE);

######################################################
# Subroutines
######################################################

# Subroutine to find the majority rule consensus nucleotides
sub consensus {
my $arrayref = shift @_; # col of nucleotide

my $numofa = 0;
my $numoft = 0;
my $numofc = 0;
my $numofg = 0;
my $numofgap = 0;
foreach my $char (@$arrayref) { # count the occurence of each nucleotide and gap
	if ($char eq "A") {
	$numofa ++;
	} elsif ($char eq "T") {
	$numoft ++;
	} elsif ($char eq "C") {
	$numofc ++;
	} elsif ($char eq "G") {
	$numofg ++;
	} else {
	$numofgap ++;
	}
}

if (($numofa >= $numoft) && ($numofa >= $numofc) && ($numofa >= $numofg)){ # majority nucleotide is A 
return "A";
} elsif (($numoft >= $numofa) && ($numoft >= $numofc) && ($numoft >= $numofg)){ # majority nucleotide is T 
return "T";
} elsif (($numofc >= $numofa) && ($numofc >= $numoft) && ($numofc >= $numofg)){ # majority nucleotide is C
return "C";
} elsif (($numofg >= $numofa) && ($numofg >= $numoft) && ($numofg >= $numofc)){ # majority nucleotide is G
return "G";
} else { # majority char is gap
return "?";
}
    
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
Script name: consensus.pl

This is a script to read the fasta sequences from each file and make a majority consensus sequence.

Example usage:
(1) Read aligned sequences from 'nf_aligned'. Write consensus sequences to 'consensus.fas':

	perl consensus.pl --indir nf_aligned --outfile consensus.fas

Input files:
(1) nf_aligned

Output files:
(1) consensus.fas

Options:
--indir
  Directory containing aligned nucleotide sequences
--outfile
  File containing majority consensus sequences for each gene
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

