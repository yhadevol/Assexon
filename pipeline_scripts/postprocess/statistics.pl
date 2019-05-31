#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use Bio::AlignIO;
use Bio::Align::DNAStatistics;
use Parallel::ForkManager;

my ($non_flank, $flank, $nf_out, $f_out, $indi_out, $help);
$nf_out = "coding_summary.txt";
$f_out = "flanking_summary.txt";
$indi_out = "sample_summary.txt";

my $flank_num_cnt_trsd = 20;
my $cpu = 1;

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'coding_aligned:s', \$non_flank,
                      'flanking_aligned:s', \$flank,
                      'coding_stat:s', \$nf_out,
                      'flanking_stat:s', \$f_out,
                      'sample_stat:s', \$indi_out,
                      'cpu:s', \$cpu,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/coding_aligned flanking_aligned/;
check_option(\@essential, \@ARGVoptions);

# get file name of coding sequences and flanking sequences
opendir NF_DIR, $non_flank or die "\nERROR: Cannot find input directory containing aligned coding sequences \"$non_flank\", please check --coding_aligned ($!)\n\n";
my @nf_infile = grep {$_ !~ /^\./} readdir NF_DIR;
@nf_infile = grep {$_ =~ /fa$|fas$|fasta$/i} @nf_infile;
closedir NF_DIR;

# make copy file dir
my $nf_dir_copy = "${non_flank}_copy";
mkdir $nf_dir_copy;

my $f_dir_copy = "${flank}_copy";
mkdir $f_dir_copy;

# multi process
my $split = $cpu;
my $splited_array = split_array(\@nf_infile, $split);

my $pm = Parallel::ForkManager->new(int($cpu/2));

my (%indi_result, %gene_nf_result, %gene_f_result);

# callback of hash containing number of reads
$pm->run_on_finish(
	sub {
		my $info = $_[-1];
		
		foreach my $taxon (sort keys %{$info->[0]}) {
		$indi_result{$taxon} = $info->[0]->{$taxon};
		}
	
		foreach my $loci (sort keys %{$info->[1]}) {
		$gene_nf_result{$loci} = $info->[1]->{$loci};
		}
	
		foreach my $loci (sort keys %{$info->[2]}) {
		$gene_f_result{$loci} = $info->[2]->{$loci};
		}
	}
);

DATA_LOOP: while (my $nf_infile_array = shift @$splited_array) {
$pm->start and next DATA_LOOP;

	# summary statistics
	my (%individual, %gene_nf, %gene_f);
	LOOP: while (my $nf_infile = shift @$nf_infile_array) {
		
		my ($gene) = $nf_infile =~ /(\S+)\.fa\S*$/; # get gene name 
		my $nf_infile_path = "$non_flank/$nf_infile"; # set file path
			
		# summary several statistics
		my $totalgap = 0;
	
		my ($seq_nf_hash, $input_order, $aln_len) = readfas($nf_infile_path);
	
		if (@$input_order == 0) {
		print STDERR "WARNING: Nothing found in $nf_infile_path. Skip this file\n";
		next LOOP;
		}
	
		if ($aln_len eq "na") {
		print STDERR "WARNING: Number of characters among sequences are not the same in $nf_infile_path. It may be not aligned. Skip this file\n";
		next LOOP;
		}
	
		my $all_nf_gc = 0;
		my $all_nf_base = 0;
	
		my $num_of_sample = @$input_order;
	
		my $copy_nf_fas = "$nf_dir_copy/$nf_infile";
		open COPY_NF_FAS, ">$copy_nf_fas";
	
		foreach my $taxon (@$input_order) {
		my $seq_nf = $seq_nf_hash->{$taxon};
	
		print COPY_NF_FAS ">$taxon\n$seq_nf\n";
	
		my @gap = $seq_nf =~ /-/g;
		$totalgap += @gap; # total number of gap in alignment
	
		$seq_nf =~ s/-//g;
		my $seq_nf_len = length($seq_nf);
	
		my @gc = $seq_nf =~ /([G|C])/g;
		$all_nf_gc += @gc;
		$all_nf_base += $seq_nf_len;

		$individual{$taxon}->{gc} += @gc;
		$individual{$taxon}->{all_base} += $seq_nf_len;
		$individual{$taxon}->{num} ++;
		}
	
		close COPY_NF_FAS;
	
		my $gap_pct = sprintf("%.2f", $totalgap/(@$input_order*$aln_len)*100); # gap percentage in alignment
		my $gc_content = sprintf("%.2f", $all_nf_gc/$all_nf_base*100); # gc content in alignment
	
		# calculate p-distance
		my $pdis = ave_pdis($copy_nf_fas);
		
		$gene_nf{$gene}->{num} = $num_of_sample;
		$gene_nf{$gene}->{nf_aln} = $aln_len;
		$gene_nf{$gene}->{nf_gc} = $gc_content;
		$gene_nf{$gene}->{nf_gap} = $gap_pct;
		$gene_nf{$gene}->{nf_pdis} = $pdis;
				        
		# statistics of flanking sequences
		my $f_infile_path = "$flank/$nf_infile";
			
		if (-e $f_infile_path) {
			my ($seq_f_hash, $input_order, $aln_len) = readfas($f_infile_path, "f");
						
			if (@$input_order == 0) {
			print STDERR "WARNING: Nothing found in $f_infile_path. Skip this file\n";
			next LOOP;
			}
		
			if ($aln_len eq "na") {
			print STDERR "WARNING: Number of characters among sequences are not the same in $f_infile_path. It may be not aligned. Skip this file\n";
			next LOOP;
			}
		
			if (@$input_order > 0) {
			my ($left_block, $right_block) = cut_block($seq_nf_hash, $seq_f_hash, $gene);
			
			my $left_block_file = "$f_dir_copy/$gene.left";
			my $right_block_file = "$f_dir_copy/$gene.right";
			
			# left block
			my $left_aln_len = 0;
			my $left_print_cnt = 0;
			open LEFT_BLOCK, ">$left_block_file";
			foreach my $taxon (sort keys %$left_block) {
				my $left_seq = $left_block->{$taxon};
			
				# length without gap
				my $temp_seq = $left_seq;
				$temp_seq =~ s/-//g;
				my $left_no_gap_len = length($temp_seq);
			
				# flank length
				$individual{$taxon}->{flen} += $left_no_gap_len;
			
				# aln length
				$left_aln_len = length($left_seq);
			
				if ($left_aln_len >= $flank_num_cnt_trsd) {
					if ($left_no_gap_len > 0) {
					print LEFT_BLOCK ">$taxon\n$left_seq\n";
					$left_print_cnt++;
					}
				} else {
				last;
				}	
			}
			close LEFT_BLOCK;
			
			# right block
			my $right_aln_len = 0;
			my $right_print_cnt = 0;
			open RIGHT_BLOCK, ">$right_block_file";
			foreach my $taxon (sort keys %$right_block) {
				my $right_seq = $right_block->{$taxon};
			
				# length without gap
				my $temp_seq = $right_seq;
				$temp_seq =~ s/-//g;
				my $right_no_gap_len = length($temp_seq);
			
				# flank length
				$individual{$taxon}->{flen} += $right_no_gap_len;
			
				# aln length
				$right_aln_len = length($right_seq);
			
				if ($right_aln_len >= $flank_num_cnt_trsd) {
					if ($right_no_gap_len > 0) {
					print RIGHT_BLOCK ">$taxon\n$right_seq\n";
					$right_print_cnt++;
					}
				} else {
				last;
				}
			}
			close RIGHT_BLOCK;
			
			my $left_pdis;
			if ($left_print_cnt > 1) {
			$left_pdis = ave_pdis($left_block_file);
			} else {
			$left_pdis = "NA";
			}
			
			my $right_pdis;
			if ($right_print_cnt > 1) {
			$right_pdis = ave_pdis($right_block_file);
			} else {
			$right_pdis = "NA";
			} 
						
			$gene_f{$gene}->{left_aln_len} = $left_aln_len;
			$gene_f{$gene}->{right_aln_len} = $right_aln_len;
			$gene_f{$gene}->{left_pdis} = $left_pdis;
			$gene_f{$gene}->{right_pdis} = $right_pdis;
			}
		}
	}
	
$pm->finish(0, [\%individual, \%gene_nf, \%gene_f]);	
}
$pm->wait_all_children;

# print NF OUT
open NON_FLANK_OUT, ">$nf_out" or die "\nERROR: Cannot write summary statistics table of coding region \"$nf_out\", please check --coding_stat ($!)\n\n";
print NON_FLANK_OUT "Locus\tNum. of samples\tAlignment length(bp)\tGC content(%)\tMissing data(%)\tPairwise distance\n";
foreach my $loci (sort keys %gene_nf_result) {
print NON_FLANK_OUT "$loci\t$gene_nf_result{$loci}->{num}\t$gene_nf_result{$loci}->{nf_aln}\t$gene_nf_result{$loci}->{nf_gap}\t$gene_nf_result{$loci}->{nf_pdis}\n";
}
close NON_FLANK_OUT;

# print F OUT
open FLANK_OUT, ">$f_out" or die "\nERROR: Cannot write summary statistics table of flanking region \"$f_out\", please check --flanking_stat ($!)\n\n";
print FLANK_OUT "Locus\tAlignment length of left flank(bp)\tAlignment length of right flank(bp)\tPairwise distance of left flank (>= $flank_num_cnt_trsd bp)\tPairwise distance of right flank (>= $flank_num_cnt_trsd bp)\n";
foreach my $loci (sort keys %gene_f_result) {
print FLANK_OUT "$loci\t$gene_f_result{$loci}->{left_aln_len}\t$gene_f_result{$loci}->{right_aln_len}\t$gene_f_result{$loci}->{left_pdis}\t$gene_f_result{$loci}->{right_pdis}\n";
}
close FLANK_OUT;

# statistics for each sample
open INDIVIDUAL, ">$indi_out" or die "\nERROR: Cannot write table of summary statistics for each sample \"$indi_out\", please check --sample_stat ($!)\n\n";
print INDIVIDUAL "Sample\tNum. of enriched loci\tGC content(%)\tAver. length of flanks(bp)\n"; #header

foreach my $taxon (sort keys %indi_result) {
	my $inaverlen;
	if (exists $indi_result{$taxon}->{flen}) {
	$inaverlen = sprintf("%.2f", $indi_result{$taxon}->{flen}/$indi_result{$taxon}->{num}); # ave whole length
	} else {
	$inaverlen = "NA";
	}
	my $inavergc = sprintf("%.2f", $indi_result{$taxon}->{gc}/$indi_result{$taxon}->{all_base});  # ave gc
	my $genenum = $indi_result{$taxon}->{num}; # enriched gene num
	print INDIVIDUAL "$taxon\t$genenum\t$inavergc\t$inaverlen\n";
}
close INDIVIDUAL;

`rm -rf $nf_dir_copy
rm -rf $f_dir_copy`;

#####################################################
# Subroutines
#####################################################

# subroutine to calculate average
sub average {
my $array = shift;

my @array = @$array;

my $total;
map {$total += $_} @array;
my $aver = $total/@array;
}

# subroutine to reformat input fasta file
sub readfas {
my $file = shift;
my $f = shift;

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

# subroutine to calculate average pdis
sub ave_pdis {
my $aln_file = shift;

my $stats = Bio::Align::DNAStatistics->new();
$stats->verbose(-1);

my $alignin = Bio::AlignIO->new(-format => 'fasta',
								-file   => $aln_file);
							
my $aln = $alignin->next_aln;
my $Uncorrected_matrix = $stats->distance(-align => $aln, 
								-method => 'Uncorrected');

my @names = @{$Uncorrected_matrix->{_names}}; 
my @distance = @{$Uncorrected_matrix->{_values}};

my $sum = 0;
my $entry_cnt = 0;
for (my $row = 0; $row < @distance; $row++) {
my @row_value = @{$distance[$row]};
my $row_name = $names[$row];

	for (my $col = $row+1; $col < @row_value; $col++) {
	my $col_name = $names[$col];
	my $entry_value = $row_value[$col];	
		if ($entry_value >= 0) {
		$sum += $entry_value;
		$entry_cnt++;
		}	
	}
}

my $ave;
if ($entry_cnt > 0) {
$ave = sprintf("%.2f", $sum/$entry_cnt);
} else {
$ave = "NA";
}

return $ave;
}

# Subroutine to cut off flanking sequences blocks
sub cut_block {
my $seq_nf = shift;
my $seq_f = shift;
my $gene = shift; 

my (%left_block, %right_block);
my $left_max_len = 0;
my $right_max_len = 0;
foreach my $taxon (sort keys %$seq_nf) {
	my $nfseq = $seq_nf->{$taxon};
	my $fseq = $seq_f->{$taxon};
	
	if (($nfseq)&&($fseq)) {
		my ($left_f_len, $right_f_len) = left_right_f_len($nfseq, $fseq);
	
		my $left_seq = substr $fseq, 0, $left_f_len;
		my $temp_seq = $left_seq;
		$temp_seq =~ s/-//g;
	
		if (! $temp_seq) {
		$left_block{$taxon} = "";
		} else {
		$left_block{$taxon} = $left_seq;
		}
	
		my $right_seq = substr $fseq, length($fseq)-$right_f_len, $right_f_len;
		$temp_seq = $right_seq;
		$temp_seq =~ s/-//g;
	
		if (! $temp_seq) {
		$right_block{$taxon} = "";
		} else {
		$right_block{$taxon} = $right_seq;
		}
	
		$left_max_len = length($left_block{$taxon}) if (length($left_block{$taxon}) > $left_max_len);
		$right_max_len = length($right_block{$taxon}) if (length($right_block{$taxon}) > $right_max_len);
	}
}

foreach my $taxon (sort keys %left_block) {
my $left_seq = $left_block{$taxon};
my $gap_num = $left_max_len - length($left_seq);
my $gap_chain = "-" x $gap_num;
$left_seq = $left_seq . $gap_chain;
$left_block{$taxon} = $left_seq;
}

foreach my $taxon (sort keys %right_block) {
my $right_seq = $right_block{$taxon};
my $gap_num = $right_max_len - length($right_seq);
my $gap_chain = "-" x $gap_num;
$right_seq = $gap_chain . $right_seq;
$right_block{$taxon} = $right_seq;
}

return (\%left_block, \%right_block);
}

# subroutine to find coding-flank boundary
sub countpos {
my $seq = shift;
my $len = shift;

my $nucleocount = 0;
my $allcount = 0;
my @seq = split //, $seq;
foreach my $nucleo (@seq) {
$nucleocount++ if ($nucleo ne "-");
$allcount++;
last if ($nucleocount == $len+1);
}

return $allcount;
}

# subroutine to find boundary of coding-flanking
sub left_right_f_len {
my $nfseq = shift;
my $fseq = shift;

my $nftemp = $nfseq;
$nftemp =~ s/-//g;
my $ftemp = $fseq;

$ftemp =~ s/-//g;
my $stpos = index $ftemp, $nftemp;
my $leftlen = $stpos;
my $rightlen = length($ftemp)-length($nftemp)-$stpos;

my $leftst;
if ($leftlen > 0) {
$leftst = countpos($fseq, $leftlen);
$leftst--;
} else {
$fseq =~ /^(\-*)[A-Z]+.+/;
$leftst = length($1);
}

my $righted;
my $rfseq = reverse($fseq);
if ($rightlen > 0) {
$righted = countpos($rfseq, $rightlen);
$righted--;
} else {
$rfseq =~ /^(\-*)[A-Z]+.+/;
$righted = length($1);
}

return($leftst, $righted);
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
Script name: statistics.pl

This script is used to summary statistics for each locus and sample from their alignments in fasta format. 

Summarized statistics for coding region of each locus including:
(1) Number of enriched samples
(2) Alignment length
(3) GC content
(4) Percentage of Missing data
(6) Average pairwise distance among sequences

Summarized statistics for flanking region of each locus including:
(1) Alignment length of left flanking region
(2) Alignment length of right flanking region
(3) Average pairwise distance in left flanking region (only calculated for flanking region >= 20 bp)
(4) Average pairwise distance in right flanking region (only calculated for flanking region >= 20 bp)

Summarized statistcis for each sample including:
(1) Number of enriched loci
(2) GC content
(3) Average length of flanking region (left + right)

Dependencies: 
(1) Perl module:
	1. Bio::AlignIO (included in Bioperl)
	2. Bio::Align::DNAStatistics (included in Bioperl)
	3. Parallel::ForkManager

Example Usage:

	perl statistics.pl --coding_aligned nf_aligned_filtered --flanking_aligned f_aligned_filtered --cpu 12

Input files:
(1) nf_aligned_filtered
(2) f_aligned_filtered

Output files:
(1) $nf_out (Summarized statistics for coding region)
(2) $f_out (Summarized statistics for flanking region)
(3) $indi_out (Summarized statistcis for each sample)

Option:
--coding_aligned
  Directory comprising aligned coding sequences
--flanking_aligned
  Directory comprising aligned coding sequences with flanks
--coding_stat
  Tab delimited table of summarized statistics for coding region of each locus, named as '$nf_out' in default
--flanking_stat
  Tab delimited table of summarized statistics for flanking region of each locus, named as '$f_out' in default
--sample_stat 
  Tab delimited table of summarized statistics for each sample, named as '$indi_out' in default
--cpu
  Limit the number of CPUs, $cpu in default
--help , -h
  Show this page and exit

Author: Hao Yuan                                                                     
        Shanghai Ocean University                                               
        Shanghai, China, 201306                                                               
                                                                                         
Created by: June 27, 2018                                                              
                                                                                         
Last modified by:
";
exit;
}
