#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use Bio::Seq;
use Parallel::ForkManager;

# input dir, output dir, taxa of reference, whether coding or not, pdis, whether trim poor column, trim how may column in a unit, output query or not
my ($dir, $outdir, $query_taxa, $non_coding_aln, $pdis_trsd, $trim_gap, $trim_unit, $rm_query, $help);

# remove sequences with too long indels
my $indel_len_trsd = 10; # maximum length of internal gap allowed

# remove misassembly by detecting part of sequence where are impossibly distant from other sequences
my $query_consensus_trsd = 5; # minimum query number required for make a query consensus instead of randomly choose one
my $pdis_num_trsd = 0; # maximum number of low pdis window found to remove a sequence
my $window_size = 50; # window size
my $step_size = 25; # window step size

# sequences number and length of alignment
my $seq_num_trsd = 2; # minimum sequence required in an alignment
my $block_len_trsd = 100; # minimum length of alignment block
my $cov_trsd = 0.8; # minimum coverage required for a sequence

# remove sequences in column
my $col_gap_trsd = 0.5; # maximum col of gap allowed

# if pct of gaps in the block are lower than this value, this loci won't be realigned by mafft
my $block_gap_pct_trsd = 0.8;

my $cpu = 1;

my $realign_dir = "realign"; # number of realign dir
my $block_dir = "block_dir";

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'indir:s', \$dir,
                      'filtered:s', \$outdir,
                      'ref_taxa:s', \$query_taxa,
                      'seq_num:i', \$seq_num_trsd,
                      'aln_len:i', \$block_len_trsd,
                      'len_pct:f', \$cov_trsd,
                      'pdis:f', \$pdis_trsd,
                      'low_pdis_num:i', \$pdis_num_trsd,
                      'window_size:i', \$window_size,
                      'step_size:i', \$step_size,
                      'indel_len_trsd:i', \$indel_len_trsd,
                      'trim_gap!', \$trim_gap,
                      'trim_unit:i', \$trim_unit,
                      'colgap_pct:f', \$col_gap_trsd,
                      'non_coding_aln!', \$non_coding_aln,
                      'remove_ref!', \$rm_query,
                      'cpu:i', \$cpu,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/indir filtered ref_taxa/;
check_option(\@essential, \@ARGVoptions);

# check mafft
`mafft --version 2> mafft_check`;
my $mafft_check = `grep "command not found" mafft_check`;
unlink "mafft_check";
if ($mafft_check) {
die "\nERROR: Cannot find mafft in \$PATH, It may just not named as \"mafft\"\n\n";
}

# check the minimum seq num threshold
die "\nERROR: At least 2 sequences are required to keep an alignment\n\n" if ($seq_num_trsd < 2);

# co-open or close of --trim_gap and --colgap_pct
if ($trim_gap) {
print STDERR "CAUTION: --trim_gap is active right now. --colgap_pct (0.5) are also turned on\n";
} else {
print STDERR "CAUTION: --trim_gap is disabled right now. --colgap_pct are turned off\n";
}

# set trim unit and pdis
if (! $non_coding_aln) { # if coding
$trim_unit = 3 if (!($trim_unit));
$pdis_trsd = 0.3 if (!($pdis_trsd));
} else { # non-coding
$trim_unit = 1 if (!($trim_unit));
$pdis_trsd = 0.4 if (!($pdis_trsd));
}

# whether query taxa are provided
if (! $query_taxa) {
die "\nERROR: Name of references must be provided\n\n";
}

# get query taxa;
my @query_taxa = $query_taxa =~ /\S+/g;
my %query_taxa = map {$_, ""} @query_taxa;

# get all file name
opendir DIR, $dir or die "\nERROR: Cannot find input directory \"$dir\", please check --indir ($!)\n\n";
my @infile = grep {$_ !~ /^\./} readdir DIR;
@infile = grep {$_ =~ /fa$|fas$|fasta$/i} @infile;
closedir DIR;

# if dir for realignment haven't created
if (!(-e $realign_dir)) {
mkdir $realign_dir or die "\nERROR: Cannot create directory \"$realign_dir\" for realignment ($!)\n\n";
}

# if dir for outdir haven't created
if (!(-e $outdir)) {
mkdir $outdir or die "\nERROR: Cannot create output directory \"$outdir\", please check --filtered ($!)\n\n";
}

# split files into several parts
my $split = $cpu;
my $splited_array = split_array(\@infile, $split);

# initiate multi-process
my $pm = Parallel::ForkManager->new(int($cpu/2));

# filter in parallel
DATA_LOOP: while (my $array = shift @$splited_array) {
$pm->start and next DATA_LOOP;
	
	LOOP: while (my $infile = shift @$array) {
		my (%seq, $filtered_seq);
				
		my $path = "$dir/$infile";
		my ($gene) = $infile =~ /(\S+)\.fa\S*$/;
		
		# check input file is wrongly formatted or have any strange characters 
		my ($seqhash, $input_order, $aln_len, $origin_taxon_num) = readfas($path);
				
		my @inputorder = @$input_order;
		
		# skip file if alignment with too few sequences were filtered
		my $seq_num = @inputorder;
		if ($seq_num < $seq_num_trsd) {
		print STDERR "WARNING: At least $seq_num_trsd sequences are required. There are only $seq_num in $path. Skip this file\n";
		next LOOP;
		}
		
		# skip file if sequences are not aligned
		if ($aln_len eq "na") {
		print STDERR "WARNING: Number of characters among sequences are not the same in $path. It may be not aligned. Skip this file\n";
		next LOOP;
		}
		
		# skip file if length of alignment is not multiple of 3
		if (! $non_coding_aln) { 
			if ($aln_len % 3 != 0) {
			print STDERR "WARNING: Length of alignment in $path cannot exactly divided by 3. Included sequences may not full-coding. Skip this file\n";
			next LOOP;
			}
		}
		
		foreach my $taxon (sort keys %$seqhash) {
			my $nucleo_seq = $seqhash->{$taxon};
			
			# skip file if length of sequences is not multiple of 3
			if (! $non_coding_aln) {
				my $temp_seq = $nucleo_seq;
				$temp_seq =~ s/-//g;
				if (length($temp_seq)%3 != 0) {
				print STDERR "WARNING: Length of sequence of $taxon in $path cannot exactly divided by 3. It may not full-coding. Skip this file\n";
				next LOOP;
				}
				
				# generate AA seq
				my $aa_seq = translate($nucleo_seq, 1);
				$seq{$taxon}->{aa} = $aa_seq;
			}
			$seq{$taxon}->{nucleo} = $nucleo_seq;
		}
		
		# check whether length can be exact divided by $trim_unit
		if ($trim_gap) {
			if ($aln_len % $trim_unit != 0) {
			print STDERR "WARNING: Length of alignment in $path can't be exactly divided by $trim_unit which is specified by --trim_unit. Skip this file\n";
			next LOOP;
			}
		}
	
		# skip alignments that too short
		next LOOP if ($aln_len < $block_len_trsd);
	
		# check whether query seq is exist
		my $query_cnt = 0;
		foreach my $taxon (sort keys %seq) {
		$query_cnt++ if (exists $query_taxa{$taxon});
		}
	
		if ($query_cnt == 0) { # die if no query taxa found
		print STDERR "WARNING: Cannot find provided taxa of reference in $path. Skip this file\n";
		next LOOP;
		}		

		# remove sequences with too long internal gap
		$filtered_seq = long_indel(\%seq);
		
		# next if filtered taxa num is too low
		my $internal_gap_filtered_taxon_num = scalar keys %$filtered_seq;
		next LOOP if ($internal_gap_filtered_taxon_num < $seq_num_trsd);
		
		$filtered_seq = gapremove($filtered_seq);
		next LOOP if (! $filtered_seq);
		
		# filter sequence with high p-distance
		$filtered_seq = pdis_window($filtered_seq);
		
		# next if filtered taxa num is too low
		my $pdis_filtered_taxon_num = scalar keys %$filtered_seq;
		next LOOP if ($pdis_filtered_taxon_num < $seq_num_trsd);
		
		$filtered_seq = gapremove($filtered_seq);
		next LOOP if (! $filtered_seq);
		
		# filter sequence with low coverage
		$filtered_seq = poor_align_filter($filtered_seq);
		
		# keep query or not
		if ($rm_query) {
			foreach my $taxa (sort keys %query_taxa) {
			delete $filtered_seq->{$taxa};
			}
		}
			
		# next if filtered taxa num is too low
		my $poor_align_filtered_taxon_num = scalar keys %$filtered_seq;
		next LOOP if ($poor_align_filtered_taxon_num < $seq_num_trsd);
		
		$filtered_seq = gapremove($filtered_seq);
		next LOOP if (! $filtered_seq);
		
		my $gapful = gapless($filtered_seq);
		
		# realign the sequence if one or more sequences was deleted
		if ($gapful) {
		$filtered_seq = realign($gene, $filtered_seq);
		}
		
		# trim gap if --trim_gap is specified
		if ($trim_gap) {
		$filtered_seq = gapremove($filtered_seq, 1);
		}
	
		# output if alignments get through the filter
		if ($filtered_seq) {
			open OUTFILE, ">$outdir/$infile";
			foreach my $taxon (@inputorder) {
				if (exists $filtered_seq->{$taxon}) {
				print OUTFILE ">$taxon\n$filtered_seq->{$taxon}->{nucleo}\n";
				}
			}
			close OUTFILE;
		}
	}
	
$pm->finish;	
}
$pm->wait_all_children;

# remove dir containing realigned sequences
`rm -rf $realign_dir`;

#####################################################
# Subroutines
#####################################################

# remove sequences with too long indels
sub long_indel {
my $seq = shift; # hash of sequence

my %seq = %$seq;

# get one of the query sequences
my @cur_query_taxa = grep {exists $seq{$_}} sort keys %query_taxa;
my $cur_query_taxa = $cur_query_taxa[0];
my $query_seq = $seq{$cur_query_taxa}->{nucleo};

foreach my $taxon (sort keys %seq) {
next if (exists $query_taxa{$taxon});
	
	# get the length of indel except beginning and ending deletion 
	my $nucleo_seq = $seq{$taxon}->{nucleo};
	my $indel = indel($query_seq, $nucleo_seq);
	
	my @indel_len;
	foreach my $indel (@$indel) {
	my ($indel_len) = $indel =~ /\S+\s(\S+)\s\S+/;
	push @indel_len, $indel_len;
	}
	
	# get the longest indel
	@indel_len = sort {$a<=>$b} @indel_len;
	
	# remove taxa if it have long indel
	if (@indel_len > 0) {
	delete $seq{$taxon} if ($indel_len[-1] > $indel_len_trsd);
	}
}

# return filter hash of sequence
return (\%seq);
}

# subroutine to find indel in sequences
sub indel {
my $seq1 = shift; # seq1 
my $seq2 = shift; # seq2

# split nucleotide into array
my @seq1 = split //, $seq1;
my @seq2 = split //, $seq2;

# get columns with at least one base
my (@gapless_seq1, @gapless_seq2);
for (my $pos=0; $pos<length($seq1); $pos++) {
	my $base1 = $seq1[$pos];
	my $base2 = $seq2[$pos];
	
	next if (($base1 eq "-")&&($base1 eq $base2)); # skip columns with 2 gap
	
	push @gapless_seq1, $base1;
	push @gapless_seq2, $base2;
}

# array of indel, length of insertion, length of deletion, length of column with 2 bases, 
my (@indel, $insert_len, $del_len, $same_len);

for (my $pos=0; $pos<@gapless_seq1; $pos++) {
	my $base1 = $gapless_seq1[$pos];
	my $base2 = $gapless_seq2[$pos];
	
	my $last_pos = $pos-1; # position of last base in indel
	
	if (($base1 =~ /[A-Z\*\?]/)&&($base2 =~ /[A-Z\*\?]/)) { # if seq1 and seq2 are both bases
		if ($del_len) { # if deletion in seq2 goes end 
		push @indel, "del $del_len " . $last_pos; # push info into @indel
		$del_len = 0; # empty scalar
		} elsif ($insert_len) { # if insertion in seq1 goes end  
		push @indel, "insert $insert_len " . $last_pos;  # push info into @indel
		$insert_len = 0; # empty scalar
		}
		$same_len++; # increment of length count of columns with 2 bases
	} elsif (($base1 eq "-")&&($base2 ne "-")) { # if there is inserted base in seq2
		if ($same_len) { # if columns with 2 bases goes end
		$same_len = 0; # empty scalar
		} elsif ($del_len) { # if deletion in seq2 goes end 
		push @indel, "del $del_len " . $last_pos; # push info into @indel
		$del_len = 0; # empty scalar
		}
		$insert_len++; # increment of length count of insertion columns
	} elsif (($base1 ne "-")&&($base2 eq "-"))  { # if there is deletion in seq2
		if ($same_len) { # if columns with 2 bases goes end
		$same_len = 0; # empty scalar
		} elsif ($insert_len) { # if insertion in seq1 goes end 
		push @indel, "insert $insert_len " . $last_pos; # push info into @indel
		$insert_len = 0; # empty scalar
		}
		$del_len++; # increment of length count of deletion columns
	}
}

# same for last indel
my $last_pos = @gapless_seq1-1;
if ($del_len) {
push @indel, "del $del_len " . $last_pos;
} elsif ($insert_len) {
push @indel, "insert $insert_len " . $last_pos;
}

# remove first and last deletion which caused by incomplete enrichment
if (@indel > 0) {
	my ($first_indel, $first_indel_len, $first_indel_last_pos) = $indel[0] =~ /(\S+)\s(\S+)\s(\S+)/;
	my ($last_indel, $last_indel_last_pos) = $indel[-1] =~ /(\S+)\s\S+\s(\S+)/;

	# do not consider first and last deletion
	if ($first_indel eq "del") {
	my $first_indel_start_pos = $first_indel_last_pos+1-$first_indel_len;
	shift @indel if ($first_indel_start_pos == 0);
	}

	if ($last_indel eq "del") {
	pop @indel if ($last_indel_last_pos == $last_pos);
	}
}

# return internal indel
return (\@indel);
}

# subroutine to remove sequence with too high p-distance in window
sub pdis_window {
my $seq = shift; # hash of all seq

my %all_seq;
foreach my $taxon (sort keys %$seq) {
$all_seq{$taxon} = $seq->{$taxon}->{nucleo};
}

# get consensus seq
# hash of reference seq, consensus seq, length of alignment
my (%query_seq, $consensus_seq, $aln_len);
foreach my $taxon (sort keys %query_taxa) {
	if (exists $all_seq{$taxon}) {
	$query_seq{$taxon} = $all_seq{$taxon};
	}
}

my @query_seq = sort keys %query_seq;
if (@query_seq >= $query_consensus_trsd) { # if taxa num reach $query_consensus_trsd, make query consensus
$consensus_seq = consensus(\%query_seq);
$aln_len = length($consensus_seq);
} else { # use first one if sequence is not enough
$consensus_seq = $query_seq{$query_seq[0]};
$aln_len = length($consensus_seq);
}

# count the window which got sequences with too high p-distance to reference
my $outlier_cnt = outlier_window(\%all_seq, $aln_len, $consensus_seq);

foreach my $taxon (sort keys %$outlier_cnt) {
next if (exists $query_taxa{$taxon}); # skip reference
	my $outlier_num = $outlier_cnt->{$taxon};
	if ($outlier_num > $pdis_num_trsd) { # delete taxon if found at least $pdis_num_trsd outlier window
	delete $seq->{$taxon};
	}
}

# return hash of filtered seqs
return ($seq);
}

# subroutine to generate consensus sequence
sub consensus {
my $seq = shift; # hash of seq

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

# return consensus seq
return($consensus_seq);
}

# subroutine to find window got sequence with too high pdistance to reference
sub outlier_window {
my $seq = shift; # hash of seq
my $aln_len = shift; # alignment length
my $consensus_seq = shift; # consensus reference

my %seq = %$seq;

# find seq distant from reference in each window
my %outlier_cnt;
for (my $st_point=0; 0<1 ;$st_point+=$step_size) {
my $rest_len = $aln_len-$st_point; # length of rest column
last if (($rest_len < $window_size)&&($st_point != 0)); # jump out loop if screen to last window
	
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
	my $seq = $seq{$taxa};
	my $window_seq = substr $seq, $st_point, $cur_window_size;
	$window_seq{$taxa} = $window_seq;
	}
	
	# cut consensus sequence according to window
	my $consensus_window_seq = substr $consensus_seq, $st_point, $cur_window_size;
	
	foreach my $taxa (sort keys %window_seq) {
	next if (exists $query_taxa{$taxa}); # skip reference
		my $pdis = pdis($consensus_window_seq, $window_seq{$taxa});
		if ($pdis ne "NaN") { 
		$outlier_cnt{$taxa}++ if ($pdis > $pdis_trsd); # cnt ++ if find a sequence got too high p-distance with reference
		}
	}
}

# hash of distant windows in each taxa
return (\%outlier_cnt);
}

# subroutine to calculate p-distance
sub pdis {
my $seq1 = shift; # seq1
my $seq2 = shift; # seq2

my $count = 0; # number of different bases columns
my $lengthcount = 0; # number of columns with 2 bases

	# calculate pdistance
	for (my $position=0; $position < length($seq1); $position++) {
		if ((substr($seq1,$position,1) =~ /[A-Z]{1}/) && (substr($seq2,$position,1) =~ /[A-Z]{1}/) ){#ignore gaps and missing data
		$lengthcount++;
			if(uc substr($seq1,$position,1) ne uc substr($seq2,$position,1)) {
			$count++;
			}
		}
	}
	
	# return p-distance
	if ($lengthcount >= ($window_size/2)) { # if number of columns with 2 bases reach half of the window
	return sprintf("%.3f", $count/$lengthcount);
	} else { # if alignment with nucleotide is too short, do not return
	return "NaN";
	}
}

# subroutine to realign sequences
sub realign {
my $gene = shift; # loci name
my $seq = shift; # hash of seqs

# choose aln by codon or simple align
my $realign_out;
if (! $non_coding_aln) { # if input is coding sequence
$realign_out = codon_aln($gene, $seq);
} else { # if input is not coding sequence
$realign_out = simple_aln($gene, $seq);
}

# save sequence into %realignment for further analysis
my %realignment;
open REALIGN_OUT, $realign_out;
while (my $line = <REALIGN_OUT>) {
	if ($line =~ /^>(\S+)/) {
	my $taxon = $1;
	chomp(my $seq = <REALIGN_OUT>);
	$realignment{$taxon}->{nucleo} = $seq;
	}
}
close REALIGN_OUT;

# hash of realigned sequences
return(\%realignment);
}

# subroutine to simply align sequence
sub simple_aln {
my $gene = shift; # loci name
my $seq = shift; # hash of sequence

# remove gap
my $realign = "$realign_dir/$gene.realign.fas";
open REALIGN, ">$realign";
foreach my $taxon (sort keys %$seq) {
my $nucleo_seq = $seq->{$taxon}->{nucleo};
$nucleo_seq =~ s/-//g;
print REALIGN ">$taxon\n$nucleo_seq\n";
}
close REALIGN;

# align by mafft
my $realign_out = "$realign_dir/$gene.realigned.fas";
mafft($realign, $realign_out);

unlink $realign;

# name of realigned files
return($realign_out);
}

# subroutine to aln sequence in codon
sub codon_aln {
my $gene = shift; # loci name
my $seq = shift; # hash of seqs

# remove gap
my $realign = "$realign_dir/$gene.realign.fas";
my $realign_aa = "$realign_dir/$gene.realign.aa.fas";
open REALIGN, ">$realign";
open REALIGN_AA, ">$realign_aa";
foreach my $taxon (sort keys %$seq) {
my $nucleo_seq = $seq->{$taxon}->{nucleo};
$nucleo_seq =~ s/-//g;
my $aa_seq = $seq->{$taxon}->{aa};
$aa_seq =~ s/-//g;
print REALIGN ">$taxon\n$nucleo_seq\n";
print REALIGN_AA ">$taxon\n$aa_seq\n";
}
close REALIGN;
close REALIGN_AA;

# align by mafft
my $realign_aa_aln = "$realign_dir/$gene.realign.aa_aln.fas";
my $realign_out = "$realign_dir/$gene.realigned.fas";
mafft($realign_aa, $realign_aa_aln);

# translate back to nucleotide
translate_back($realign_aa_aln, $realign, $realign_out);

unlink $realign;
unlink $realign_aa_aln;
unlink $realign_aa;

# name of realigned files
return($realign_out);
}

# subroutine to remove columns of with too much gap
sub gapremove {
my $seqhash = shift; # hash of aligned seqs
my $gap_trim = shift; # flag to determine whether trim columns with too few bases or just trim columns with only gaps

my %seq = %$seqhash; # hash of aligned seqs
my @taxon = keys %seq; # array of taxa name
my $seqnum = scalar @taxon; # number of sequences
my $seqlen = length($seq{$taxon[0]}->{nucleo}); # alignment length

# split sequence and push into %gaptrim
my %gaptrim;
foreach my $sub (sort keys %seq) {
my @seq = split //, $seq{$sub}->{nucleo};
push @{$gaptrim{$sub}}, @seq;
}

# read through each col
my @trim;
for (my $i=0; $i<$seqlen; $i=$i+$trim_unit) {
	my @pos;
	map {push @pos, $i+$_} (0..($trim_unit-1)); # trim in $trim_unit
	
	# if find col with too much gap, push all @pos into @trim
	my @no_trim_gap;
	foreach my $col (@pos) {
		my $gapcnt = 0;
		foreach my $sub (sort keys %gaptrim) {
		$gapcnt++ if ($gaptrim{$sub}->[$col] eq "-");
		}
		my $gap_pct = $gapcnt/$seqnum;
		
		if ($gap_trim) {
			if ($gap_pct > $col_gap_trsd) {
			push @trim, @pos;
			last;
			}
		} else {
			if ($gap_pct == 1) {
			push @no_trim_gap, $col;
			}
		}
	}
	
	if (! $gap_trim) {
	push @trim, @pos if (@no_trim_gap == $trim_unit);
	}
}

# calculate the left alignment length
my $left_len = $seqlen - @trim;
if ($left_len >= $block_len_trsd) { # do not output if alignment length is too short
	# remove col
	foreach my $sub (sort keys %gaptrim) {
		foreach my $pos (@trim) {
		$gaptrim{$sub}->[$pos] = "";
		}
		my $trimmed_seq = join "", @{$gaptrim{$sub}};
		delete $gaptrim{$sub};
		$gaptrim{$sub}->{nucleo} = $trimmed_seq;
		$gaptrim{$sub}->{aa} = translate($trimmed_seq, 1) if (! $non_coding_aln);
	}
	# return hash of trimmed sequences
	return (\%gaptrim);	
} else {
return (0);	# return nothing if entire alignment is discarded
}
	
}

# subroutine to translate dna into aa
sub translate{
my $seq = shift @_; #import the sequences
my $start_codon = shift @_; #import the start codon
$start_codon --; #minus 1 for 1st 2nd 3rd codon in Bio::seq is 0 1 2
my $seq_obj=Bio::Seq->new(-seq=>$seq,-alphabet=>'dna'); #import the seq and claim that it's dna
my $pro=$seq_obj->translate(-frame => $start_codon, -terminator => "*", -unknown => "X"); #import the start codon, the sign for terminator and unknown codon, then translate
$pro=$pro->seq; 
return($pro); #return the aa
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

# subroutine to translate aligned amino acid back to nucleotide sequences
sub translate_back {
my $inpath_aa = shift; # aligned amino acid
my $inpath_nt = shift; # unaligned nucleotides
my $outfile = shift; # outfile name

# get nucleotide sequences and save to hash
my ($FILE_NT, $FILE_AA, $OUTFILE, %dna, %aa, @inputorder);
open ($FILE_NT, "<$inpath_nt") or die "Cannot find $inpath_nt ($!)";
while (my $line = <$FILE_NT>) {
	if ($line =~ />(\S+)/) {
	chomp(my $seq = <$FILE_NT>);
    $dna{$1} = $seq;
	}
}

# get amino acid sequences and save to hash
open ($FILE_AA, "<$inpath_aa") or die "Cannot find $inpath_aa ($!)"; # open aligned aa seq file
while (my $line = <$FILE_AA>) {
	if ($line =~ />(\S+)/) {
	push @inputorder, $1;
	chomp(my $seq = <$FILE_AA>);
    $aa{$1} = $seq;
	}
}
close $FILE_NT;
close $FILE_AA;	

open ($OUTFILE, ">$outfile") or die "can't open $outfile for writing!!!"; # create outfile
	
# translate back aa seqs of each taxa
foreach my $taxon (@inputorder){
next unless ((exists $aa{$taxon})&&(exists $dna{$taxon})); # skip nucleotide sequences not exists in amino acid alignment

	# get amino acid and nucleotide sequences
	my $aln_aa_seq = $aa{$taxon};
	my $nucleo_seq = $dna{$taxon};
	$nucleo_seq =~ s/-//g;	
			
	print $OUTFILE ">$taxon\n";                                    

	my ($chars, @chars);
	
	# split sequences into serial of gaps or bases, push these string into array
	my $pervious_char = substr $aln_aa_seq, 0, 1;
	$chars .= $pervious_char;
	for (my $i=1; $i<length($aln_aa_seq); $i++) {
	my $cur_char = substr $aln_aa_seq, $i, 1;
		if (($pervious_char ne "-")&&($cur_char eq "-")) {
		push @chars, $chars;
		$chars = "";
		} elsif (($pervious_char eq "-")&&($cur_char ne "-")) {
		push @chars, $chars;
		$chars = "";
		}
		$chars .= $cur_char;
		$pervious_char = $cur_char;
	}
	push @chars, $chars;
	
	# add gaps and string from unaligned dna
	my $aln_nucleo_seq;
	my $pointer = 0;
	foreach my $chars (@chars) {
		if ($chars =~ /-/) { # if string is gaps
		my $gaps = $chars x 3;
		$aln_nucleo_seq .= $gaps; # add 3 times of gaps into aligned dna string
		} else { # if string is bases
		my $length = length($chars)*3; 
		my $sub_seq = substr $nucleo_seq, $pointer, $length; # substr charaters with 3 times of length of AA string from unaligned dna
		$aln_nucleo_seq .= $sub_seq; # add to string into aligned dna 
		
		$pointer += $length;
		}
	}
	
	# print aligned dna
	print $OUTFILE "$aln_nucleo_seq\n";
}
	
close $OUTFILE;
}

# subroutine to split array in to several sub-array
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

# subroutine to delete columns with low coverage
sub poor_align_filter {
my $seq = shift; # hash of aligned seqs

my %seq = %$seq;

foreach my $taxon (sort keys %seq) {
next if (exists $query_taxa{$taxon}); # skip reference
	my $seq = $seq{$taxon}->{nucleo};
	my $coverage = coverage($seq);
	if ($coverage < $cov_trsd) { # delete taxon if coverage is too low
	delete $seq{$taxon};
	}
}

# return filtered seqs
return(\%seq);
}

# subroutine to calculate coverage
sub coverage {
my $seq = shift; # sequences

# compute coverage
my @gaps = $seq =~ /-/g; # number of gaps
my $seq_len = length($seq); # length of aligned seq
my $coverage = sprintf("%.3f", ($seq_len-@gaps)/$seq_len); # pct of bases in aligned seq

# return coverage
return $coverage;
}

# check whether filtered alignment have many gaps, alignment with too many gaps need to be realigned
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

if ($gapful) { # if there are many gaps in alignment
return(1);
} else { # if there are not many gaps in alignment
return(0);
}
}

# subroutine to reformat input fasta file
sub readfas {
my $file = shift;

open INFILE, $file or die "\nERROR: Cannot find input file \"$file\" ($!)\n\n";

my (%seq, %seqlength, @input_order, $seq, $taxon, $lasttaxon, %taxa_cnt);
my $origin_seq_num = 0;

while (my $line = <INFILE>) {

	# remove enter
	$line =~ s/\r//g;
	chomp $line;
	
	# if find >
	if ($line =~ /^>(\S+)/) {
		$origin_seq_num++;
	
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
	if ($taxon_num > 1) {
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

return (\%seq, \@input_order, $aln_len, $origin_seq_num);
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
Script name: filter.pl

This is a script to remove sequences or columns from alignments in fasta format including:
(1) sequences with too long insertion or deletion to reference 
(2) Sliding window sift through alignments. If too many windows are distant from reference, sequences will be removed
(3) short and poorly enriched sequences and alignments
(4) trim columns with few nucleotides

Dependencies: 
(1) BioPerl v1.007001 or higher
(2) Mafft v7.294b or higher (rename it as 'mafft' and put it under \$PATH)
(3) Perl module:
	1. Parallel::ForkManager

Example usage:
(1) Filter aligned full-coding sequences in 'nf_aligned', and write output to 'filtered'. Only remove poorly aligned sequences. Sequence of 'taxon1' and 'taxon2' are reference. Run script with 4 process: 

	perl filter.pl --indir nf_aligned --filtered filtered --ref_taxa taxon1 --cpu 4

Input files:
(1) nf_aligned

Output files:
(1) filtered

Options:

Input and output:
--indir:
  Name of input directory containing unfiltered alignments
--filtered:
  Name of output directory containing filtered alignments
  
Taxa of reference sequences:
--ref_taxa:
  A space delimit list of reference taxa. This option is mandatory
--remove_ref:
  Remove references sequences from alignment if this option is specified 

Coding sequences or not:
--non_coding_aln
  If non-coding sequence is input, Please specify this option 

Filter sequences distant from reference:
--pdis
  Maximum pairwise distance allowed between reference and other sequences, 0.3 for coding sequences (when --non_coding_aln is not specified), 0.4 for non-coding sequences (when --non_coding_aln specified) in default
--window_size
  Size of window to sift through the alignment, $window_size bp in default
--step_size
  Step size of window movement, $step_size bp in default
--low_pdis_num
  If more than --low_pdis_num distant windows to reference are detected, sequences will be removed, $pdis_num_trsd in default

Filter sequences with too long indel:
--indel_len_trsd:
  Maximum length of consecutive gaps between nucleotides allowed, $indel_len_trsd in default

Filter short sequences and alignments:
--seq_num:
  Maximum number of sequences allowed to keep a alignment, $seq_num_trsd in default
--aln_len:
  Minumum alignment length allowed to keep a alignment, $block_len_trsd in default
--len_pct:
  Minimum number of nucleotides/alignment length allowed to keep a sequence, $cov_trsd in default
                      
Trim sequences in columns:
--trim_gap:
  Trim columns with too many gaps, this option is turned on with --colgap_pct
--trim_unit:
  Number of column trimmed at once. In default, if coding sequences are input (when --non_coding_aln is not specified), --trim_unit is 3. Non-coding sequences are input (when --non_coding_aln is specified), --trim_unit is 1
--colgap_pct:
  Maximum percentage of gaps allowed in a column, $col_gap_trsd in default

Other parameters:
--cpu:
  Limit the number of CPUs, 1 in default
--help , -h:
  Show this help message and exit

Author: Hao Yuan                                                                     
        Shanghai Ocean University                                               
        Shanghai, China, 201306                                                               
                                                                                         
Created by: Nov 20, 2018                                                              
                                                                                         
Last modified by: 
";
exit;
}