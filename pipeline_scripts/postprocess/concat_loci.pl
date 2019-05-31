#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;   # include the module for input

my $dir; # Dir containing aligned nucleotide sequences
my $outfile; # concatenated genes
my $first_name; # Name of first sequence
my $help; # help message
my $minnumseq = 1; #the minimum num of seq for a gene being selected

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'indir:s', \$dir,
                      'outfile:s', \$outfile,
                      'min:i', \$minnumseq,
                      'first_name:s', \$first_name,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/indir outfile/;
check_option(\@essential, \@ARGVoptions);  

# get all fasta file
opendir DIR, $dir or die "\nERROR: Cannot find input directory \"$dir\", please check --indir ($!)\n\n";
my @infile = grep {$_ !~ /^\./} readdir DIR;
@infile = grep {$_ =~ /fa$|fas$|fasta$/i} @infile;
closedir DIR;

my (%taxon_list, %skip);
LOOP: foreach my $infile (@infile) {
	# read the first sequence
	my $path = "$dir/$infile";
	
	my ($seqhash, $input_order, $aln_len) = readfas($path);
	
	if (@$input_order == 0) {
	print STDERR "WARNING: Nothing is found in $path. Skip this file\n";
	$skip{$infile} = "";
	next LOOP;
	}
	
	if ($aln_len eq "na") {
	die "\nERROR: Number of characters among sequences are not the same in $path. It may be not aligned\n\n";
	}
	
	my $first_taxon = $input_order->[0];
	
	# die if first name is gene name 
	if (! $first_name) {
		if ($infile =~ /^$first_taxon$/) {
		die "\nERROR: Name of first sequence is similar to its file ($path) which may be the name of gene. Please use --first_name or rename the first sequence, or script cannot determine name of the first concatenated sequence\n\n";
		}
	}
	
	# summary enriched gene number for each taxon
	if ($first_name) { # deem taxon of first sequence as $first_name if --first_name is specified
	$taxon_list{$first_name}++;
	} else {
	$taxon_list{$first_taxon}++;
	}
	
	# read rest of file and summary enriched gene number for each taxon
	foreach (1..(@$input_order-1)) {
	$taxon_list{$input_order->[$_]}++;
	}
}

# name of output
my $nex_outfile = $outfile . ".nex"; # nexus file
my $fas_outfile = $outfile . ".fas"; # fasta file
my $phy_outfile = $outfile . ".phy"; # phylip file
my $LOGFILE = $dir . ".numofgenescaptured.txt"; # number of capture gene

open OUTFILE, ">$nex_outfile" or die "\nERROR: Cannot write outfile in nex format \"$nex_outfile\", please check --outfile ($!)\n\n"; # write out nexus file first

my (%genenum, %fasta, %taxa_gap);
my $gap_num = 0; # gap count
my $all_num = 0; # char count
my $start = 1;
my $end = 1;
my $numselect = 0; #num of gene selected according to num of seq

foreach my $infile (@infile) {
	next if (exists $skip{$infile});
	
	#get the name of fasta file
	my $infile_path = "$dir/$infile"; 

	my %seq;
	
	# check input file is wrongly formatted or have any strange characters 
	my ($seqhash, $input_order, $aln_len) = readfas($infile_path, "nowarn");
	
	foreach my $taxon (sort keys %$seqhash) {
	$seq{$taxon} = $seqhash->{$taxon};
	$genenum{$taxon}++;
	}
	open INFILE, $infile_path or die "\nERROR: Cannot find $infile_path ($!)\n\n";
    
    # print it to nexus output if there are at least $minnumseq
    my $numofhash = scalar keys %seq;
    if ($numofhash >= $minnumseq) {
    $numselect ++;
     
    my $gaps = "-" x $aln_len; # fill gap if there's no taxon in this gene

	my ($seqid) = $infile =~ /(\S+)\.fa\S*$/; #get the sequence name
	$end = $start + $aln_len - 1; # end position
	print OUTFILE "[$seqid $aln_len bp $start-$end bp]\n"; # header 
	$start = $end + 1; # new start position
	
		# print sequence or gap and save sequence to %fasta prepare for fasta output
        foreach my $taxon (sort keys %taxon_list) {
            print OUTFILE "$taxon\t";
            if (exists $seq{$taxon}){
            my $seq = $seq{$taxon};
            
			# count gap and all char num
			my @gaps = $seq =~ /-/g;
			$gap_num += @gaps;
			$all_num += $aln_len;
			
			$taxa_gap{$taxon} += @gaps;
            
            print OUTFILE "$seq\n";
            $fasta{$taxon} .= $seq{$taxon};
            } else {
            # count gap and all char num
            $gap_num += $aln_len;
			$all_num += $aln_len;
            
            $taxa_gap{$taxon} += $aln_len;
            
            print OUTFILE "$gaps\n";    
            $fasta{$taxon} .= $gaps;                 
            } 
        }
    }
    print OUTFILE "\n";
}
close OUTFILE;

# print number of genes captured for each taxon
open LOGFILE, ">$LOGFILE" or die "\nERROR: Cannot write log file \"$LOGFILE\" ($!)\n\n";
foreach my $taxon (sort keys %taxon_list) {
print LOGFILE "$taxon\t$genenum{$taxon}\n";
}
close LOGFILE;

# write fasta output 
open FASTA_OUTFILE, ">$fas_outfile" or die "\nERROR: Cannot write outfile in fasta format \"$fas_outfile\", please check --outfile ($!)\n\n";
foreach my $taxon (sort keys %taxon_list) { #open files for output
print FASTA_OUTFILE ">$taxon\n$fasta{$taxon}\n";
}
close FASTA_OUTFILE;

# write phylip output
my $taxa_num = scalar keys %taxon_list;
my $total_len = $end;
open PHY_OUTFILE, ">$phy_outfile" or die "\nERROR: Cannot write outfile in phy format \"$phy_outfile\", please check --outfile ($!)\n\n";
print PHY_OUTFILE "$taxa_num $total_len\n";
foreach my $taxon (sort keys %taxon_list) { #open files for output
print PHY_OUTFILE "$taxon\t$fasta{$taxon}\n";
}
close PHY_OUTFILE;

# report number of enriched gene
print "\nThere are $numselect genes included in the output files!\n";

# report percentage of gap
my $gappct = sprintf("%.2f", $gap_num/$all_num*100);
print "\nThere is ${gappct}% gap in concatenated alignment\n";

# report percentage of gap in each taxon
my $taxon_aln_len = $all_num/$taxa_num;
print "\nLength of concatenated alignment is $taxon_aln_len bp\n\n";

foreach my $taxon (sort keys %taxa_gap) {
my $gappct = sprintf("%.2f", $taxa_gap{$taxon}/$taxon_aln_len*100);
print "$taxon has ${gappct}% gaps\n";
}

#####################################################
# Subroutines
#####################################################

# subroutine to reformat input fasta file
sub readfas {
my $file = shift;
my $nowarn = shift;

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
					if (! $nowarn) {
						if ($length == 0) {	
						print STDERR "WARNING: No nucleotide found sequence of $lasttaxon in $file. This taxon will be discarded.\n";
						}
						if (@strange_char > 0) {
						my $strange_char = join " ", @strange_char;
						print STDERR "WARNING: Found strange character \"$strange_char\" in sequence of $lasttaxon in $file. This taxon will be discarded.\n";					
						}
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
			if (! $nowarn) {
				if ($length == 0) {	
				print STDERR "WARNING: No nucleotide found sequence of $lasttaxon in $file. This taxon will be discarded.\n";
				}
				if (@strange_char > 0) {
				my $strange_char = join " ", @strange_char;
				print STDERR "WARNING: Found strange character \"$strange_char\" in sequence of $lasttaxon in $file. This taxon will be discarded.\n";					
				}
			}
		} 
	} else {
	$taxa_cnt{$lasttaxon}++
	}
}

close INFILE;

if (! $nowarn) {
	foreach my $taxon (sort keys %taxa_cnt) {
	my $taxon_num = $taxa_cnt{$taxon};
	$taxon_num++;
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

sub usage {
print STDERR "
Script name: concat_loci.pl

This is a script to read alignments in fasta format, then concatenate them into nexus, fasta and phylip file

Example usage:
(1) Concatenate all loci in directory 'nf_aligned'

	perl concat_loci.pl --indir nf_aligned --outfile concat

(2) Concatenate all loci in directory 'nf_aligned', and only include loci with at least 2 sequences:

	perl concat_loci.pl --indir nf_aligned --outfile concat --min 2

(3) Concatenate all loci in directory 'nf_aligned', replace the name of first sequence (it\'s probably a reference) as 'ref_species':

	perl concat_loci.pl --indir nf_aligned --outfile concat --first_name ref_species

Input files:
(1) nf_aligned

Output files:
(1) concat.nex
(2) concat.fas
(3) concat.phy

Options:
--indir
  Directory containing aligned nucleotide sequences
--outfile
  Name prefix of outfile. If --outfile is 'concat', the outfile name would be 'concat.nex', 'concat.fas' and 'concat.phy'
--min
  Only include loci with at least --min sequences, $minnumseq in default
--first_name
  Name of first sequence. Please use this option, if first sequence is a reference and name of it in each file is different 
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