#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
use Bio::AlignIO;
use Bio::Align::DNAStatistics;

# non-flank, flank, outdir, name of reference, whether keep reference in output
my ($nf, $f, $outdir, $query_taxa, $query_output, $help);

my $deviation_trsd = 3; # if flank sequences of a sample is shorter 3 bp than flanking gaps of reference, sequences of this sample will be discarded for the alignment of coding sequences is problematic

# parameters for removing unevenly enriched sequences at two sides of the alignment
my $few_col_len = 5; # length of consecutive columns with enough number and percentage of nucleotides
my $few_col_num = 1; # minimum number of nucleotides required in a columns
my $few_col_pct = 0.5; # minimum percentage of nucleotides required in a columns 

# remove taxa with too long unique insertion which may mis-assembly or intron insertion
my $gap_col_taxon_num_trsd = 5; # minimum sequences required to remove misassembly by detecting consecutively and uniquely gap-aligned sequences 
my $unique_col_len_trsd = 10; # if length of consecutively and uniquely gap-aligned sequences are longer than $unique_col_len_trsd, the sequences will be removed

# remove too variable flanks
my $pdis_high_trsd = 0.4; # maximum required p-dis allowed to keep a sequence
my $step_size = 10; # window step size
my $window_size = 20; # window size
my $pdis_num_trsd = 0; # if more than $pdis_num_trsd window have distant flanks, sequences will trimmed from distant window to the end 

# remove poor flanking block
my $min_block_len = 10; # minimum length of flanking block required
my $cov_trsd = 0.65; # remove low coverage sequences in the block
my $min_len = 10; # minimum number of nucleotide require in a flanking sequence
my $seq_pct_low_trsd = 0.65; # minimum percentage of non-gap sequences required to keep the block
my $seq_num_trsd = 1; # minimum number of non-gap sequence required to keep the sequence

# if pct of gaps in left and right block are lower than this value, this loci won't be aligned by mafft
my $block_gap_pct_trsd = 0.8;

my $block_dir = "block_dir"; # temporary dir

my $cpu = 1;

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'flank:s', \$f,
					  'nonflank_filtered:s', \$nf,
					  'flank_filtered:s', \$outdir,
					  'ref_taxa:s', \$query_taxa,
					  'keep_ref!', \$query_output,
					  'few_col_len:i', \$few_col_len,
					  'few_col_num:i', \$few_col_num,
					  'few_col_pct:f', \$few_col_pct,
					  'unique_col_len:i', \$unique_col_len_trsd,
					  'pdis:f', \$pdis_high_trsd,
					  'low_pdis_num:i', \$pdis_num_trsd,
                      'window_size:i', \$window_size,
                      'step_size:i', \$step_size,
					  'block_len:i', \$min_block_len,
					  'seq_len:i', \$min_len,
					  'len_pct:f', \$cov_trsd,
					  'seq_pct:f', \$seq_pct_low_trsd,
					  'seq_num:i', \$seq_num_trsd,
					  'cpu:i', \$cpu,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/flank nonflank_filtered flank_filtered ref_taxa/;
check_option(\@essential, \@ARGVoptions);

# check mafft
`mafft --version 2> mafft_check`;
my $mafft_check = `grep "command not found" mafft_check`;
unlink "mafft_check";
if ($mafft_check) {
die "\nERROR: Cannot find mafft in \$PATH, It may just not named as \"mafft\"\n\n";
}

# check whether --ref_taxa is specified
my %query_taxa;
if ($query_taxa) {
my @query_taxa = $query_taxa =~ /\S+/g;
%query_taxa = map {$_, ""} @query_taxa;
} else {
die "\nERROR: Please provide space delimited list of reference taxa by --ref_taxa \n\n";
}

# get name of loci in nf dir
my %filename;
opendir NF_DIR, $nf or die "\nERROR: Cannot find directory to aligned coding sequences \"$nf\", please check --nonflank_filtered ($!)\n\n";
while (my $nf_file = readdir NF_DIR) {
next if ($nf_file =~ /^\./);
my $loci_name = extract_loci_name($nf_file);
$filename{$loci_name}->{nf} = "$nf/$nf_file";
}
closedir NF_DIR;

# get name of loci in f dir
opendir F_DIR, $f or die "\nERROR: Cannot find directory to aligned sequences with flanking regions \"$f\", please check --flank ($!)\n\n";
while (my $f_file = readdir F_DIR) {
next if ($f_file =~ /^\./);
my $loci_name = extract_loci_name($f_file);
$filename{$loci_name}->{f} = "$f/$f_file";
}
closedir F_DIR;

# create outdir if it do not exists
if (!(-e $outdir)) {
mkdir $outdir or die "\nERROR: Cannot create output directory \"$outdir\", please check --flank_filtered ($!)\n\n";
}

# create dir of temporary files if it do not exists
if (!(-e $block_dir)) {
mkdir $block_dir or die "\nERROR: Cannot create directory \"$block_dir\" for temporary files, which has already exists\n\n";
}

# quite STDERR and STDOUT from Bio::Align::DNAStatistics
my $stats = Bio::Align::DNAStatistics->new();
$stats->verbose(-1);

# split dataset into several parts
my @filename = sort keys %filename;

my $split = $cpu;
my $splited_array = split_array(\@filename, $split);

# initiate multi-process
my $pm = Parallel::ForkManager->new(int($cpu/2));

# filter in parallel
DATA_LOOP: while (my $array = shift @$splited_array) {
$pm->start and next DATA_LOOP;

	LOOP: while (my $loci_name = shift @$array) {
		my $key_num = scalar keys %{$filename{$loci_name}};
		
		# start filter if both sequences with or without flanks exists
		if ($key_num == 2) {
			my $nf_file = $filename{$loci_name}->{nf};
			my $f_file = $filename{$loci_name}->{f};
			
			# extract information of coding sequences
			my (%seq, @input_order, %seq_len);
			open NF, $nf_file;
			while (my $line = <NF>) {
				if ($line =~ />(\S+)/) {
				my $header = $1;
				push @input_order, $header;
		
				chomp(my $nfseq = <NF>);
				$seq_len{length($nfseq)} = "";
				$seq{$header}->{nf} = $nfseq;
				}
			}
			close NF;
			
			# skip file with too few sequences
			my $len_key_num = scalar keys %seq_len;
			if ($len_key_num > 1) {
			print STDERR "WARNING: Sequences of $nf_file are not aligned. Skip this locus\n";
			next LOOP;
			}
			
			# extract information of flanking sequences
			%seq_len = ();
			open F, $f_file;
			while (my $line = <F>) {
				if ($line =~ />(\S+)/) {
				my $header = $1;		
				chomp(my $fseq = <F>);
				$seq_len{length($fseq)} = "";
				$seq{$header}->{f} = $fseq;
				}
			}
			close F;
			
			# skip file with too few sequences
			$len_key_num = scalar keys %seq_len;
			if ($len_key_num > 1) {
			print STDERR "WARNING: Sequences of $f_file are not aligned. Skip this locus\n";
			next LOOP;
			}
			
			# check number of sequences co-exist in nf and f files
			my $all_taxa_num = scalar keys %seq;
				
			my $co_exist_cnt = 0;
			foreach my $taxa (sort keys %seq) {
				my $key_num = scalar keys %{$seq{$taxa}};
				if ($key_num < 2) {
				delete $seq{$taxa} 
				} else {
				$co_exist_cnt++;
				}
			}
		
			next LOOP if ($co_exist_cnt < ($seq_num_trsd+1));
			
			# cut off left and right flanking block
			my $nf_block_unchanged = 0;
			$nf_block_unchanged = 1 if ($co_exist_cnt == $all_taxa_num);
				
			my ($left_len, $right_len, $left_block, $right_block, $nf_block, $query_cnt, $deviation_change) = cutblock(\%seq, \%query_taxa);

			if ($left_len eq "no_query") {
			print STDERR "WARNING: Cannot find reference sequences in $loci_name. Skip this file\n";
			next LOOP;
			}
			
			# trim unevenly enriched sequences from left
			if (($left_len >= $few_col_len)&&($left_len > 0)) {
			($left_block, $left_len) = few_col_base($left_block, $left_len, $query_cnt);
			}
			
			# trim unevenly enriched sequences from right
			if (($right_len >= $few_col_len)&&($right_len > 0)) {
			$right_block = reverse_block($right_block);
			($right_block, $right_len) = few_col_base($right_block, $right_len, $query_cnt);
			$right_block = reverse_block($right_block);
			}
			
			# filter poorly aligned sequences of left flanks
			if (($left_len >= $min_block_len)&&($left_len > 0)) {
			$left_block = reverse_block($left_block);
			($left_block, $left_len) = filter($left_block, $left_len, \%query_taxa, $query_cnt, $loci_name);
			$left_block = reverse_block($left_block);
			}
			
			# filter poorly aligned sequences of right flanks
			if (($right_len >= $min_block_len)&&($right_len > 0)) {
			($right_block, $right_len) = filter($right_block, $right_len, \%query_taxa, $query_cnt, $loci_name);
			}
			
			# filter too short block or block with too few sequences			
			my $left_block_gap_pct = 0;
			if ($left_len > 0) {
			$left_block = poor_block_filter($left_block, \%query_taxa, $query_cnt);
			}
			
			my $right_block_gap_pct = 0;
			if ($right_len > 0) {	
			$right_block = poor_block_filter($right_block, \%query_taxa, $query_cnt);
			}

			my %left_block = %$left_block;
			my %right_block = %$right_block;
			my %nf_block = %$nf_block;
			
			# concatenate flank and coding sequences	
			my %concat_block;
			foreach my $taxa (@input_order) {
				if ((exists $left_block{$taxa})&&(exists $right_block{$taxa})&&(exists $nf_block{$taxa})) {
				$concat_block{$taxa} = $left_block{$taxa} . $nf_block{$taxa} . $right_block{$taxa};
				}
			}		
			
			# keep query or not
			if (! $query_output){
				foreach my $query_taxa (sort keys %query_taxa) {
				delete $concat_block{$query_taxa} if (exists $concat_block{$query_taxa});
				}
			}
			
			# check whether too many gaps in alignment
			my $gapful;
			if (! $deviation_change) {
			my $left_gapful = gapless($left_block);
			my $right_gapful = gapless($right_block);
			my $nf_gapful = gapless($nf_block);
			$gapful = $left_gapful + $right_gapful + $nf_gapful;
			}
			
			# output if there is not many gaps
			if ((! $deviation_change)&&(! $gapful)) {
				my ($reduced_concat_block, $reduced_concat_block_len) = trim_gap_col(\%concat_block);
				my %reduced_concat_block = %$reduced_concat_block;
		
				open OUTFILE, ">$outdir/$loci_name.fas";
				foreach my $taxa (@input_order) {
					if (exists $reduced_concat_block{$taxa}) {
					print OUTFILE ">$taxa\n$reduced_concat_block{$taxa}\n";
					}
				}
				close OUTFILE;	
			} else { # realign if there are too many gaps
				my $seq_cnt = 0;
				my $unrealigned = "$outdir/$loci_name.fas.unrealigned";
				
				open UNREALIGNED, ">$unrealigned";
				foreach my $taxa (@input_order) {
					if (exists $concat_block{$taxa}) {
					my $seq = $concat_block{$taxa};
					$seq =~ s/-//g;
					print UNREALIGNED ">$taxa\n$seq\n";
					$seq_cnt++
					}
				}
				close UNREALIGNED;

				if ($seq_cnt > 1) {
				mafft($unrealigned, "$outdir/$loci_name.fas");
				unlink ($unrealigned);
				} elsif ($seq_cnt == 1) {
				`mv $unrealigned $outdir/$loci_name.fas`;
				}
			}	
		} else { # warn if file of the coding or coding sequences with flanks is missing 
			if (!(exists $filename{$loci_name}->{f})) {
			print STDERR "WARNING: Cannot find aligned sequences with flanking region of $loci_name under $f. This locus will be skipped\n";
			} elsif (!(exists $filename{$loci_name}->{nf})) {
			print STDERR "WARNING: Cannot find aligned coding sequences of $loci_name under $nf. It may be filtered by filter.pl. This locus will be skipped\n";
			}
		}
	}

$pm->finish;	
}
$pm->wait_all_children;

# remove temp dir
`find $block_dir -delete`;

#####################################################
# Subroutines
#####################################################

# subroutine to get loci name
sub extract_loci_name {
my $name = shift;

my $loci_name;
if ($name =~ /(\S+)\.al[i]*[g]*n[e]*[d]*\.fas\S*/) {
$loci_name = $1;
} elsif ($name =~ /(\S+)\.fas\S*/) {
$loci_name = $1;
}

return $loci_name;
}

# Subroutine to cut off flanking sequences blocks
sub cutblock {
my $seq = shift;
my $query_taxa = shift;

my %seq = %$seq;
my %query_taxa = %$query_taxa;

my @query_taxa = grep {exists $query_taxa{$_}} sort keys %seq;
my $query_cnt = @query_taxa;

if ($query_cnt > 0) {
	my $one_of_query = $query_taxa[0];

	my (%left_block, %right_block, %nf_block);
	my $longest_left_block = 0;
	my $longest_right_block = 0;

	my $query_nfseq = $seq{$one_of_query}->{nf};
	my $reverse_query_nfseq = reverse($query_nfseq);
	my ($query_left_gap) = $query_nfseq =~ /^(\-*)\S+/;
	my ($query_right_gap) = $reverse_query_nfseq =~ /^(\-*)\S+/;
	
	my $query_left_gap_len = length($query_left_gap);
	my $query_right_gap_len = length($query_right_gap);

	my $query_fseq = $seq{$one_of_query}->{f};
	my ($query_left_f_len, $query_right_f_len) = left_right_f_len($query_nfseq, $query_fseq);
	
	$left_block{$one_of_query} = substr $query_fseq, 0, $query_left_f_len;
	$right_block{$one_of_query} = substr $query_fseq, length($query_fseq)-$query_right_f_len, $query_right_f_len;
	$nf_block{$one_of_query} = substr $query_fseq, $query_left_f_len, length($query_fseq)-$query_left_f_len-$query_right_f_len;
	
	$longest_left_block = $query_left_f_len if ($query_left_f_len > $longest_left_block); 
	$longest_right_block = $query_right_f_len if ($query_right_f_len > $longest_right_block); 
		
	my %left_right_gap_nf;
	my $deviation_change = 0;
	foreach my $taxa (sort keys %seq) {
		next if ($taxa eq $one_of_query);
		my $nfseq = $seq{$taxa}->{nf};
		my $fseq = $seq{$taxa}->{f};
		my $reverse_nfseq = reverse($nfseq);
		my ($left_gap) = $nfseq =~ /^(\-*)\S+/;
		my ($right_gap) = $reverse_nfseq =~ /^(\-*)\S+/;
		
		my $left_nf_len = length($left_gap);	
		my $right_nf_len = length($right_gap);
		my ($left_f_len, $right_f_len) = left_right_f_len($nfseq, $fseq);
				
		my $nf_left_deviation = $left_nf_len - $query_left_gap_len;
		my $nf_right_deviation = $right_nf_len - $query_right_gap_len;
		
		my $f_left_deviation = $left_f_len - $query_left_f_len;
		my $f_right_deviation = $right_f_len - $query_right_f_len;
		
		next if (!(($f_left_deviation >= -$deviation_trsd)&&($f_right_deviation >= -$deviation_trsd)));
		
		if ($f_left_deviation > 0) {
		my $left_seq_btw_query_sub = substr $fseq, $query_left_f_len, $f_left_deviation;
		$left_seq_btw_query_sub =~ s/-//g;
		$deviation_change++ if (length($left_seq_btw_query_sub) > 0);
		}
		
		if ($f_right_deviation > 0) {
		my $right_seq_btw_query_sub = substr $fseq, length($fseq)-$right_f_len, $f_right_deviation;
		$right_seq_btw_query_sub =~ s/-//g;
		$deviation_change++ if (length($right_seq_btw_query_sub) > 0);
		}
			
		$left_block{$taxa} = substr $fseq, 0, $left_f_len;
		$right_block{$taxa} = substr $fseq, length($fseq)-$right_f_len, $right_f_len;
		$nf_block{$taxa} = substr $fseq, $left_f_len, length($fseq)-$left_f_len-$right_f_len;
		
		$longest_left_block = $left_f_len if ($left_f_len > $longest_left_block); 
		$longest_right_block = $right_f_len if ($right_f_len > $longest_right_block); 
	}
	
	if ($deviation_change > 0) {
		foreach my $taxa (sort keys %left_block) {
		my $seq = $left_block{$taxa};
		my $deviation = $longest_left_block-length($seq);
		my $gaps = "-" x $deviation;
		$left_block{$taxa} = $seq . $gaps;
		}

		foreach my $taxa (sort keys %right_block) {
		my $seq = $right_block{$taxa};
		my $deviation = $longest_right_block-length($seq);
		my $gaps = "-" x $deviation;
		$right_block{$taxa} = $gaps . $seq;
		}
	} else {
		foreach my $taxa (sort keys %seq) {
		next if ($taxa eq $one_of_query);
		my $fseq = $seq{$taxa}->{f};
		$left_block{$taxa} = substr $fseq, 0, $query_left_f_len;
		$right_block{$taxa} = substr $fseq, length($fseq)-$query_right_f_len, $query_right_f_len;
		$nf_block{$taxa} = substr $fseq, $query_left_f_len, length($fseq)-$query_left_f_len-$query_right_f_len;
		}
		
		$longest_left_block = $query_left_f_len;
		$longest_right_block = $query_right_f_len;
	}
	
	return ($longest_left_block, $longest_right_block, \%left_block, \%right_block, \%nf_block, $query_cnt, $deviation_change);
} else {
return ("no_query");
}

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

# Subroutine to count number of flanks with nucleotides
sub no_gap_seq_cnt {
my $block = shift;
my %block = %$block;

my %no_gap_seq;
foreach my $taxa (sort keys %block) {
my $seq = $block{$taxa};
$seq =~ s/-//g;
$no_gap_seq{$taxa} = "" if (length($seq) > 0);
}

return (\%no_gap_seq);
}

# Subroutine to remove unevenly enriched flanks
sub few_col_base {
my $block = shift;
my $block_len = shift;
my $query_cnt = shift;

my %block = %$block;

my $few_col_continue = 0;
my ($cut_col, $cut_pass, $new_block_len);
for (my $col = 0; $col<$block_len; $col++) {
	my $base_str; 
	foreach my $sample (sort keys %block) {
	my $seq = $block{$sample};
	my $base = substr $seq, $col, 1;
	$base_str .= $base;
	}
	
	my @base = $base_str =~ /[A-Z]/g;
	my $base_num = @base;
	my $col_num = length($base_str);
	$col_num -= $query_cnt;
	
	my $base_col_pct = $base_num/$col_num;
	if (($base_col_pct >= $few_col_pct)&&($base_num >= $few_col_num)) {
	$few_col_continue++;
	} else {
	$few_col_continue = 0;
	}
	
	if ($few_col_continue >= $few_col_len) {
	$cut_col = $col-$few_col_continue+1;
	$cut_pass = 1;
	last;
	}
}
	
if ($cut_pass) {
	foreach my $sample (sort keys %block) {
	my $sub_seq = substr $block{$sample}, $cut_col;
	$block{$sample} = $sub_seq;
	$new_block_len = length($sub_seq);
	} 
} else {
	foreach my $sample (sort keys %block) {
	$block{$sample} = "";
	$new_block_len = 0;
	} 
} 

return(\%block, $new_block_len);	
}

# Subroutine to reverse sequences of flank block
sub reverse_block {
my $block = shift;

	my %block = %$block;
	foreach my $sample (sort keys %block) {
	$block{$sample} = reverse $block{$sample};
	}
	
return (\%block);
}

# Subroutine to filter poorly aligned sequences
sub filter {
my $block = shift;
my $block_len = shift;
my $query_taxa = shift;
my $query_cnt = shift;
my $loci_name = shift;

my %block = %$block;
my %query_taxa = %$query_taxa;

# remove sequence with long and unique insertion
my $filtered_block = cov_col(\%block, $block_len, $query_taxa, $query_cnt);
($filtered_block, $block_len) = trim_gap_col($filtered_block);

# remove too variable flanks to other flanks
if (($block_len >= $min_block_len)&&($block_len > 0)) {
$filtered_block = var_seq($filtered_block, $block_len, $loci_name, $query_taxa);

# filter sequence with high p-distance to consensus
$filtered_block = pdis_window($filtered_block, $block_len, $query_taxa);
}

return ($filtered_block, $block_len);
}

# subroutine to remove sequence with uniq insertion
sub cov_col {
my $seq = shift;
my $block_len = shift;
my $query_taxa = shift;
my $query_cnt = shift;

my $gaps = "-" x $block_len;

my %seq = %$seq;
my %query_taxa = %$query_taxa;

my $gap_col = 0;
my @taxon = keys %seq;
my $taxon_num = scalar @taxon;
$taxon_num -= $query_cnt;
my $seqlen = length($seq{$taxon[0]});

if ($taxon_num >= $gap_col_taxon_num_trsd) { # remove sequence with uniq align if taxa number is greater than $gap_col_taxon_num_trsd
	# split sequence and push into %seq_matrix
	my %seq_matrix;
	foreach my $taxon (sort keys %seq) {
	my @seq = split //, $seq{$taxon};
	push @{$seq_matrix{$taxon}}, @seq;
	}
	
	my %unique_taxon;
	
	# foreach columns
	for (my $col=0; $col<$seqlen; $col++) {
		# count the number of nucleotide in each columns
		my @nucleo_taxon;
		my $nucleo_cnt = 0;
		foreach my $taxon (sort keys %seq_matrix) {
			if ($seq_matrix{$taxon}->[$col] ne "-") {
			push @nucleo_taxon, $taxon;
			$nucleo_cnt++;
			}
		}
		
		# save the column number if there's only one nucleotide found
		if ($nucleo_cnt == 1) {
			foreach my $nucleo_taxon (@nucleo_taxon) {
			push @{$unique_taxon{$nucleo_taxon}}, $col;
			}
		}	
	}
	
	# find consecutively and uniquely gap-aligned block
	my %unique_taxon_block;
	foreach my $taxon (sort keys %unique_taxon) {
		my @col = @{$unique_taxon{$taxon}};
		my $col_lag = shift @col;
		my @continue;
		push @continue, $col_lag;
		for (my $col_order=0; $col_order<@col; $col_order++) {
		my $cur_col = $col[$col_order];
			if ($cur_col-$col_lag > 1) { # refresh @continue if block is not consecutive
			push @{$unique_taxon_block{$taxon}}, scalar @continue;
			@continue = ();
			}
			push @continue, $cur_col;
			$col_lag = $cur_col;
		}
		push @{$unique_taxon_block{$taxon}}, scalar @continue;
	}
	
	foreach my $taxon (sort keys %unique_taxon_block) {
	next if (exists $query_taxa{$taxon}); # skip reference
		my @unique_taxon_len = sort {$a<=>$b} @{$unique_taxon_block{$taxon}};
		if ($unique_taxon_len[-1] > $unique_col_len_trsd) { # delete file if consecutive block is longer than $unique_col_len_trsd
		$seq{$taxon} = $gaps;
		}
	}
}

return(\%seq);
}

# subroutine to remove sequence with too high p-distance to consensus in window
sub pdis_window {
my $seq = shift;
my $block_len = shift;
my $query_taxa = shift;

my %seq = %$seq;
my %query_taxa = %$query_taxa;

my $gaps = "-" x $block_len;

my $no_gap_seq = no_gap_seq_cnt($seq);
my %no_gap_seq = %$no_gap_seq;
my $no_gap_seq_num = scalar keys %no_gap_seq;

my %skip_taxa;
foreach my $query (sort keys %query_taxa) {
$skip_taxa{$query} = "";
}
foreach my $sub (sort keys %seq) {
$skip_taxa{$sub} = "" if (!(exists $no_gap_seq{$sub}));
}

my $consensus_seq = consensus(\%seq);
	
# count the window which got sequences with too high p-distance to reference
my $outlier_taxa = outlier_window($seq, $consensus_seq, \%skip_taxa);

foreach my $taxon (sort keys %$outlier_taxa) {
	next if (exists $skip_taxa{$taxon}); # skip reference
	my @outlier_st_point = @{$outlier_taxa->{$taxon}};
	
	if (@outlier_st_point > $pdis_num_trsd) {
	@outlier_st_point = sort {$a<=>$b} @outlier_st_point;
	my $first_outlier_st_point = $outlier_st_point[0];
	my $whole_seq = $seq{$taxon};
	my $pre_seq = substr $whole_seq, 0, $first_outlier_st_point;
	my $append_seq = "-" x ($block_len-$first_outlier_st_point);
	$seq{$taxon} = $pre_seq . $append_seq;
	}
}

return (\%seq);
}

# subroutine to generate consensus sequence
sub consensus {
my $seq = shift;
my %seq = %$seq;

# split sequence and push into @split_aln
my @split_aln;
foreach my $taxa (sort keys %seq) {
my $seq = $seq{$taxa};
my @split_seq = split //, $seq;
push @split_aln, \@split_seq;
}

# get alignment length
my $seq_len = scalar @{$split_aln[0]};

# generate majority consensus nucleotide for each columns
my @consensus_seq;
for (my $pos=0 ;$pos<$seq_len ;$pos++) {
	my %nucleo;
	foreach my $split_seq (@split_aln) {
	$nucleo{$split_seq->[$pos]}++;
	}
	delete $nucleo{"-"};
	my @sort_nucleo = sort {$nucleo{$a}<=>$nucleo{$b}} keys %nucleo;
	if (@sort_nucleo > 0) {
	push @consensus_seq, $sort_nucleo[-1];
	} else {
	push @consensus_seq, "-";
	}
}

# join all consensus nucleotide
my $consensus_seq = join "", @consensus_seq;

return($consensus_seq);
}

# subroutine to find window got sequence with too high pdistance to consensus
sub outlier_window {
my $seq = shift;
my $consensus_seq = shift;
my $skip_taxa = shift;

my %skip_taxa = %$skip_taxa;

my %seq = %$seq;
my $aln_len = length($consensus_seq);

my %outlier_cnt;
for (my $st_point=0; 0<1 ;$st_point+=$step_size) {
my $rest_len = $aln_len-$st_point;
last if (($rest_len < $window_size)&&($st_point != 0));
	
	# compute the window size, it may alter when move the end of the alignment
	my $cur_window_size;
	if ($rest_len >= ($window_size+$step_size)) { # if rest alignment length is enough for another move
	$cur_window_size = $window_size; # $window_size is window_size 
	} else { # if rest alignment length is not enough for another move
	$cur_window_size = $rest_len; # rest length is window size
	}
	
	# cut sequence according to window
	my (%window_seq);
	foreach my $taxa (sort keys %seq) {
	next if (exists $skip_taxa{$taxa}); # skip reference
	my $seq = $seq{$taxa};
	my $window_seq = substr $seq, $st_point, $cur_window_size;
	$window_seq{$taxa} = $window_seq;
	}
	
	# cut consensus sequence according to window
	my $consensus_window_seq = substr $consensus_seq, $st_point, $cur_window_size;
	
	foreach my $taxa (sort keys %window_seq) {
		my $pdis = pdis($consensus_window_seq, $window_seq{$taxa});
		if ($pdis ne "NaN") { 
		push @{$outlier_cnt{$taxa}}, $st_point if ($pdis > $pdis_high_trsd); # cnt ++ if find a sequence got too high p-distance with reference
		}
	}
}

return (\%outlier_cnt);
}

# subroutine to calculate p-distance
sub pdis {
my $seq1 = shift;
my $seq2 = shift;

my $count = 0;
my $lengthcount = 0;

# calculate pdistance
for (my $position=0; $position < length($seq1); $position++) {
	if ((substr($seq1,$position,1) =~ /[A-Z]{1}/) && (substr($seq2,$position,1) =~ /[A-Z]{1}/) ){#ignore gaps and missing data
	$lengthcount++;
		if(uc substr($seq1,$position,1) ne uc substr($seq2,$position,1)) {
		$count++;
		}
	}
}


if ($lengthcount >= ($window_size/2)) {
return sprintf("%.3f", $count/$lengthcount);
} else { # if alignment with nucleotide is too short, do not return
return "NaN";
}
}

# Subroutine to filter short block or block with too few sequences
sub poor_block_filter {
my $block = shift;
my $query_taxa = shift;
my $query_cnt = shift;

my %query_taxa = %$query_taxa;

	my ($reduced_block, $block_len) = trim_gap_col($block);
	my %block = %$reduced_block;
	
	my $gaps = "-" x $block_len;
	
	if (($block_len >= $min_block_len)&&($block_len > 0)) {
		foreach my $sample (sort keys %block) {
		next if (exists $query_taxa{$sample});
			my @nucleo = $block{$sample} =~ /[A-Z]/g;
			my $nucleo_num = @nucleo;
			my $nucleo_cov = $nucleo_num/$block_len;
						
			if (!(($nucleo_num >= $min_len)&&($nucleo_cov >= $cov_trsd))) {
			$block{$sample} = $gaps;
			}
		}
		
		my $key_num = scalar keys %block;

		my $no_gap_seq = no_gap_seq_cnt(\%block);
		my $sample_seq_cnt = scalar keys %$no_gap_seq;
		
		$key_num -= $query_cnt;
		my $seq_pct = $sample_seq_cnt/$key_num;

		if (!(($seq_pct >= $seq_pct_low_trsd)&&($sample_seq_cnt >= $seq_num_trsd))) {
			foreach my $sample (sort keys %block) {
			$block{$sample} = "";
			}
		}
	} else {
		foreach my $sample (sort keys %block) {
		$block{$sample} = "";
		}
	}
	
return (\%block);
}

# Subroutine to check whether flanks filtered from blocks
sub check_block_change {
my $seq_num_bf_filter = shift;
my $seq_num_af_filter = shift;
my $cut_block_change = shift;
	
	if ($cut_block_change) {
		if ($seq_num_af_filter == 0) {
		return 1;
		} else {
		return 0;
		}
	} else {
		if ($seq_num_bf_filter == $seq_num_af_filter) {
		return 1;
		} else {
			if ($seq_num_af_filter == 0) {
			return 1;
			} else {
			return 0;
			}
		}
	}
}

# Subroutine to trim sequences in columns
sub trim_gap_col {
my $block = shift;

my %seq = %$block;
my @taxa_order = sort keys %seq;

	my (@split_aln, %block);
	foreach my $taxa (@taxa_order) {
	$block{$taxa} = "";

	my $seq = $seq{$taxa};
	my @split_seq = split //, $seq;
		for (my $pos=0; $pos<@split_seq; $pos++) {
		push @{$split_aln[$pos]}, $split_seq[$pos];
		}
	}

	my $block_len = scalar @split_aln;

	my %gap_col;
	for (my $pos=0; $pos<$block_len ;$pos++) {
	my @nucleo_col = @{$split_aln[$pos]};
	my $col_bases = join "", @nucleo_col;
	my @all_gaps = $col_bases =~ /-/g;

		if (@all_gaps != length($col_bases)) {
			for (my $i=0; $i<@taxa_order; $i++) {
			$block{$taxa_order[$i]} .= $nucleo_col[$i];
			}
		}
	}
	
	$block_len = length($block{$taxa_order[0]});

return (\%block, $block_len);
}

# subroutine to remove too variable sequences to others
sub var_seq {
my $block = shift;
my $block_len = shift;
my $loci_name = shift;
my $query_taxa = shift;

my %seq = %$block;
my %query_taxa = %$query_taxa;

my $gaps = "-" x $block_len;
 
my %qua_seq;
foreach my $taxa (sort keys %seq) {
next if (exists $query_taxa{$taxa});
	my $seq = $seq{$taxa};
	my @nucleo = $seq =~ /[A-Z]/g;
	if (@nucleo >= $min_len) {
	$qua_seq{$taxa} = $seq;
	}
}

my $qua_seq_num = scalar keys %qua_seq;

if ($qua_seq_num > 1) {
	my $qua_seq_path = "$block_dir/$loci_name.fas";
	open QUA_SEQ, ">$qua_seq_path"; 
	foreach my $taxa (sort keys %qua_seq) {
	print QUA_SEQ ">$taxa\n$qua_seq{$taxa}\n";
	}
	close QUA_SEQ;
	
	my $alignin = Bio::AlignIO->new(-format => 'fasta',
                                    -file   => "$qua_seq_path");
	my $aln = $alignin->next_aln;
	my $Uncorrected_matrix = $stats->distance(-align => $aln, 
									-method => 'Uncorrected');
	
	my @names_order = @{$Uncorrected_matrix->{_names}};
	my @value = @{$Uncorrected_matrix->{_values}};
	my $seq_num = scalar @value;
	
	my (%high_pdis, @high_pdis);
	for (my $i=0; $i<(@value-1); $i++) {
		for (my $j=$i+1; $j<@value; $j++) {
			my $pdis = sprintf("%.2f", $value[$i]->[$j]);
			if ($pdis > $pdis_high_trsd) {
			push @high_pdis, "$i $j";
			$high_pdis{$i} = "";
			$high_pdis{$j} = "";
			}
		}
	}
			
	my @bad_seq;
	foreach my $high_pdis_pair (@high_pdis) {
	my ($col1, $col2) = $high_pdis_pair =~ /(\S+)\s(\S+)/;
		if ((exists $high_pdis{$col1})&&(exists $high_pdis{$col2})) {
			my $pdis_all_col1 = 0;
			for (my $i=0; $i<@value; $i++) {
			next if ($col1 == $i);
			my $pdis = $value[$i]->[$col1];
			$pdis_all_col1 += $pdis;
			}
		
			my $pdis_all_col2 = 0;
			for (my $i=0; $i<@value; $i++) {
			next if ($col2 == $i);
			my $pdis = $value[$i]->[$col2];
			$pdis_all_col2 += $pdis;
			}
			if ($pdis_all_col1 > $pdis_all_col2) {
			push @bad_seq, $col1;
			delete $high_pdis{$col1};
			} else {
			push @bad_seq, $col2;
			delete $high_pdis{$col2};
			}
		}
	}
	
	foreach my $bad_seq_order (@bad_seq) {
	my $taxa = $names_order[$bad_seq_order];
	$seq{$taxa} = $gaps;
	}
	
	unlink $qua_seq_path;
}

return (\%seq);
}

# subroutine to align sequence with mafft
sub mafft {
my $unaligned = shift;
my $aligned = shift;

# align sequence by mafft
my $aln_temp_path = $aligned . ".temp";
`mafft --auto --preservecase --quiet --bl 80 "$unaligned" > "$aln_temp_path"`;

# unwrap the output
unwrap($aln_temp_path, $aligned);
unlink $aln_temp_path;
}

# subroutine to unwrap the output of mafft
sub unwrap {
my $in = shift;
my $unwrap = shift;

open IN, $in;
open UNWRAP, ">$unwrap";

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

# subroutine to check whether too many gaps in alignment
sub gapless {
my $block = shift;
my %block = %$block;

my $gapful = 0;
foreach my $taxa (sort keys %block) {
	my $seq = $block{$taxa};
	my (@nucleo) = $seq =~ /[A-Z]/g;
	my $nulceo_num = scalar @nucleo;
	if ($nulceo_num > 0) {
	$gapful++ if ($nulceo_num/length($seq) < $block_gap_pct_trsd);
	}
}

if ($gapful) {
return(1);
} else {
return(0);
}
}

# subroutine to split input data into several sets
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
Script name: flank_filter.pl

This is a script to remove sequences from alignments with flanking regions in fasta format including:
(1) Unevenly enriched sequences at two sides of the alignments
(2) Unique insertion in one taxa (> 10bp), which is potential misassembly
(3) Too variable sequences (> 0.4 pdis with other sequences)
(4) Sliding window sifts through flanks. If there are too many windows distant from their consensus sequences, then sequences will be trimmed from windows.
(5) Short sequences and poorly enriched flanking blocks

Dependencies: 
(1) Mafft v7.294b or higher (rename it as 'mafft' and put it under \$PATH)
(2) Perl module:
	1. Parallel::ForkManager
	2. Bio::AlignIO
    3. Bio::Align::DNAStatistics

Example usage:
(1) Filter sequences with flanking regions in 'f_aligned', and only files and sequences co-exist in 'f_aligned' and 'nf_filtered' will be written to 'f_filtered'. Sequence of 'taxon1' is reference. Run script with 4 process: 

	perl flank_filter.pl --flank f_aligned --nonflank_filtered nf_filtered --flank_filtered f_filtered --ref_taxa taxon1 --cpu 4

(2) Not only filter sequences, but also keep sequences of reference (taxon) in output

	perl flank_filter.pl --flank f_aligned --nonflank_filtered nf_filtered --flank_filtered f_filtered --ref_taxa taxon1 --keep_ref --cpu 4

Input files:
(1) f_aligned nf_filtered

Output files:
(1) f_filtered

Options:

Input and output:
--flank:
  Name of input directory containing aligned sequences with flanking regions
--nonflank_filtered:
  Name of input directory containing well aligned coding sequences. Poorly aligned coding sequences can be filtered by filter.pl
--flank_filtered:
  Name of output directory containing well aligned sequences with flanking regions

Taxa of reference sequences:
--ref_taxa:
  A space delimit list of reference taxa. This option is mandatory, because reference sequences will be used to find boundary of coding and noncoding sequences
--keep_ref:
  Reference sequences will be kept in output, if this option is specified. Reference sequences are discarded by default

Remove unevenly enriched sequences at two sides of the alignments:
--few_col_len:
  Minimum length of consecutive columns with enough number and percentage of nucleotides, $few_col_len in default
--few_col_num:
  Minimum number of nucleotides required in a columns, $few_col_num in default
--few_col_pct:
  Minimum percentage of nucleotides required in a columns, reference is not included, $few_col_pct in default

Filter potential misassemblies or long insertion:
--unique_col_len:
  Minimum length of unique insertion existing in only one taxa which may misassemblies or long insertion, $unique_col_len_trsd in default

Filter too variable flanks:
--pdis:
  Maximum pairwise distance allowed between flanks or consensused sequences and flanks, $pdis_high_trsd in default
--window_size:
  Size of window to sift through the alignment, $window_size bp in default
--step_size
  Step size of window movement, $step_size bp in default
--low_pdis_num
  If too many distant windows were found, sequences will be trimmed from these windows to the end of sequences. Maximum number of window is $pdis_num_trsd in default

Filter short sequences and poorly enriched flanking block:
--len_pct:
  Minimum number of nucleotides/alignment length allowed to keep a sequence, $cov_trsd in default
--block_len:
  Minimum length of flanking block, $min_block_len in default
--seq_len
  Minimum number of nucleotides in flanking sequences, $min_len in default
--seq_pct
  Minimum percentage of taxa in flanking blocks, reference is not included, $seq_pct_low_trsd in default
--seq_num
  Minimum number of sequences in flanking blocks, $seq_num_trsd in default

Other parameters:
--cpu:
  Limit the number of CPUs, 1 in default
--help or -h:
  Show this help message and exit

Author: Hao Yuan                                                                     
        Shanghai Ocean University                                               
        Shanghai, China, 201306                                                               
                                                                                         
Created by: Nov 20, 2018                                                              
                                                                                         
Last modified by: 
";
exit;
}