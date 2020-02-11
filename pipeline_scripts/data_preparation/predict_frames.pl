#!/usr/bin/env perl

use warnings;
use strict;
use Bio::Seq;
use Getopt::Long;

my @ARGVoptions = @ARGV;

my $stop_trsd = 1; # maximum stop codon
my $id_trsd = 80; # minimum identity between protein and baits
my $cov_trsd = 90; # minimum percentage of baits covered on reference protein
my $frame_info_table = "frame_result.txt"; # frame correct result
my $uncorrected_loci_list = "uncorrected_loci.txt"; # uncorrected loci

my $temp_dir = "temp"; # temp dir to do pairwise alignment

# baits, onehitCDSmarker, reference protein, dna output, aa output, flag for onehitCDSmarker header, help
my ($baits, $cds, $proteome, $dna_out, $aa_out, $cds_header, $help);
my $opt = GetOptions( 'baits:s', \$baits,
					  'cds:s', \$cds,
					  'ref_prot:s', \$proteome,
					  'dna_out:s', \$dna_out,
            'aa_out:s', \$aa_out,
					  'stop_codon:i', \$stop_trsd,
					  'id:i', \$id_trsd,
					  'cov:i', \$cov_trsd,
					  'cds_header!', \$cds_header,
            'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/baits cds ref_prot/;
check_option(\@essential, \@ARGVoptions);

# check existence of baits file
if (!(-e $baits)) {
die "\nERROR: Cannot find file of baits sequences \"$baits\", please check --baits ($!)\n\n";
}

# check existence of cds file
if (!(-e $cds)) {
die "\nERROR: Cannot find file of frame predicted sequences \"$cds\", please check --cds ($!)\n\n";
}

# check existence of reference protein file
if (!(-e $proteome)) {
die "\nERROR: Cannot find file of reference protein \"$proteome\", please check --ref_prot ($!)\n\n";
}

# check usearch
`exonerate` or die "\nERROR: Cannot find exonerate in \$PATH. It may just not named as 'exonerate' \n\n";

# baits file format
my $prefix;
if ($baits =~ /fa$|fas$|fasta$|fna$/i) {
($prefix) = $baits =~ /(\S+)\.[^\.]+$/;
} else {
die "\nERROR: Baits sequences must be in fasta format ($!)\n\n";
}

# reference protein format
if ($proteome !~ /fa$|fas$|fasta$|fna$/i) {
die "\nERROR: Reference protein must be in fasta format ($!)\n\n";
}

# specify output name if --dna_out and --aa_out are omitted
if ((! $dna_out) || (! $aa_out)) {
$dna_out = "$prefix.dna.fas";
$aa_out = "$prefix.aa.fas";
}

# check format of fasta file, and get name and sequence of each loci ($seqhash)
my ($seqhash, $input_order, $aln_len) = readfas($baits);

# generate hash of baits sequences and loci name:
# loci position->{header} = original loci name
# 			  |->{seq} = sequence
my $baits_seq = baits_hash($seqhash);

# get hash of loci and corresponding gene and hash of gene present in onehitCDS file
# loci position = {gene name}
# gene name = ""
my ($loci_gene, $all_gene) = cds($cds);

# find frame by pairwise alignment and output
my ($frame_info, $uncorrected_loci) = findframes($baits_seq, $proteome, $loci_gene, $all_gene);

# print frame result for each loci
open FRAME_INFO, ">$frame_info_table";
print FRAME_INFO "Loci\tFrame result\n";
foreach my $loci (sort keys %$frame_info) {
print FRAME_INFO "$loci\t$frame_info->{$loci}\n";
}
close FRAME_INFO;

# print uncorrected loci
open UNCORRECTED_LOCI, ">$uncorrected_loci_list";
foreach my $loci (@$uncorrected_loci) {
print UNCORRECTED_LOCI "$loci\n";
}
close UNCORRECTED_LOCI;

# remove temp dir
`rm -rf $temp_dir`;

######################################################
# Subroutines
######################################################

# subroutine to generate hash of baits sequences and loci name:
sub baits_hash {
my $seqhash = shift;

my %baits_seq;
foreach my $header (sort keys %$seqhash) {
	my @header_fields = $header =~ /[^\.]+/g;
	my $loci_position = $header_fields[-3] . "." . $header_fields[-2] . "." . $header_fields[-1];

	my $seq = $seqhash->{$header};
	$seq = uc($seq);

	if (!(exists $baits_seq{$loci_position})) {
	$baits_seq{$loci_position}->{header} = $header;
	$baits_seq{$loci_position}->{seq} = $seq;
	} else {
	print STDERR "\nWARNING: $loci_position exists multiple times\n\n";
	}
}

return \%baits_seq;
}

# subroutine to get hash of loci and its corresponding gene and hash of gene present in onehitCDS file
sub cds {
my $cds = shift;

my (%loci_gene, %gene);
open CDS, $cds or die "\nERROR: Cannot open $cds for reading !\n\n";
<CDS> if ($cds_header);
while (my $line = <CDS>) {
my @columns = $line =~ /\S+/g;
my $column1 = $columns[0];
my @fields = $column1 =~ /[^\|]+/g;
my $loci_position_ori = $fields[0];
my $gene = $fields[1];

my @loci_position_ori = $loci_position_ori =~ /[^\:]+/g;
my $loci_position = $loci_position_ori[-3] . "." . $loci_position_ori[-2]. "." . $loci_position_ori[-1];

$loci_gene{$loci_position} = $gene;
$gene{$gene} = "";
}
close CDS;

return (\%loci_gene, \%gene);
}

# subroutine to find coding frame, output then return frame result
sub findframes {
my $baits_seq = shift;
my $proteome = shift;
my $loci_gene = shift;
my $all_gene = shift;

my $exist;
my %proteome;

open PROTEOME, $proteome;
while (my $line = <PROTEOME>) {
	if ($line =~ />(.+)/) {
		my $header = $1;
		my @header_field = $header =~ /\S+/g;
		my $prot_id = $header_field[0];
		my @gene_field = grep {$_ =~ /^gene\:(\S+)/} @header_field;
		my ($gene_id) = $gene_field[0] =~ /^gene\:(\S+)/;

		if (exists $all_gene->{$gene_id}) {
		$proteome{$gene_id}->{$prot_id} = "";
		$exist = [$gene_id, $prot_id];
		} else {
		$exist = "";
		}
	} else {
	chomp $line;
	$proteome{$exist->[0]}->{$exist->[1]} .= $line if ($exist);
	}
}
close PROTEOME;

mkdir $temp_dir;

open DNA, ">$dna_out";
open AA, ">$aa_out";

my (@uncorrected_loci, %frame_info);
foreach my $loci_position (sort keys %$baits_seq) {
	my $uncorrected_seq = $baits_seq->{$loci_position}->{seq};
	my $ori_name = $baits_seq->{$loci_position}->{header};
	my $ori_loci_position = $ori_name;

	my $print = 0;
	$frame_info{$ori_name} = "NA";

	# no ref correct
	my ($corrected_frame, $corrected_prot, $corrected_seq) = no_ref_correct($uncorrected_seq);

	if (($corrected_frame)&&($corrected_prot)) { # if only one frame is fine
		my $prot = translate($uncorrected_seq, $corrected_frame);
		my ($corrected_prot, $stopnum) = rm_trail_stop($prot);
		my $corrected_seq = substr $uncorrected_seq, $corrected_frame-1, length($corrected_prot)*3;

		print DNA ">$ori_name\n$corrected_seq\n";
		print AA ">$ori_name\n$corrected_prot\n";

		$frame_info{$ori_name} = $corrected_frame;
		$print = 1;
	} elsif (($corrected_frame)&&(! $corrected_prot)) { # ref correct if multiple frame left in consideration
		$loci_position =~ s/-/_/ if (! (exists $loci_gene->{$loci_position}));

		if (exists $loci_gene->{$loci_position}) {
			my $gene = $loci_gene->{$loci_position};
			if (exists $proteome{$gene}) {
				my %protseq = %{$proteome{$gene}};

				my $corrected_frame = pair_aln($uncorrected_seq, \%protseq, $ori_name);

				if ($corrected_frame) {
					my $prot = translate($uncorrected_seq, $corrected_frame);
					my ($corrected_prot, $stopnum) = rm_trail_stop($prot);

					if ($stopnum <= $stop_trsd) {
					my $corrected_seq = substr $uncorrected_seq, $corrected_frame-1, length($corrected_prot)*3;
					print DNA ">$ori_name\n$corrected_seq\n";
					print AA ">$ori_name\n$corrected_prot\n";
					$frame_info{$ori_name} = $corrected_frame;
					$print = 1;
					}
				}
			} else {
# 			print STDERR "\nWARNING: Reference protein of loci $ori_name ($gene) was absent in $proteome\n";
			}
		} else {
# 		print STDERR "\nWARNING: Cannot find reference gene of loci $ori_name in index file\n";
		}

		# select frame with fewest stop codon if it cannot be corrected by ref
		if ($print == 0) {
		my $prot = translate($uncorrected_seq, $corrected_frame);
		my ($corrected_prot, $stopnum) = rm_trail_stop($prot);

		my $corrected_seq = substr $uncorrected_seq, $corrected_frame-1, length($corrected_prot)*3;
		print DNA ">$ori_name\n$corrected_seq\n";
		print AA ">$ori_name\n$corrected_prot\n";
		$frame_info{$ori_name} = $corrected_frame;
		$print = 1;
		}
	}

	push @uncorrected_loci, $ori_loci_position if (! $print);
}

close DNA;
close AA;

return (\%frame_info, \@uncorrected_loci);
}

# subroutine to translate dna to aa
sub translate {
my $seq = shift @_; # DNA sequence
my $start_codon = shift @_; # start codon (1, 2 or 3)
$start_codon --; # translate function in bioperl only accept 0,1 or 2

# translate
my $seq_obj=Bio::Seq->new(-seq=>$seq,-alphabet=>'dna'); # create Bio::Seq object
my $pro=$seq_obj->translate(-frame => $start_codon, -terminator => "*", -unknown => "X"); # translate
$pro=$pro->seq;
return($pro);
}

# subroutine to remove trailing stop codons and get number of stop codons in aa
sub rm_trail_stop {
my $prot = shift;

my ($nost) = $prot =~ /(\S+\w+)\**$/;
my @stcodon = $nost =~ /\*/g;
my $num = @stcodon;

return ($nost, $num);
}

# subroutine to align baits with reference protein to find the frame
sub pair_aln {
my $uncorrected_seq = shift;
my $protseq = shift;
my $loci_name = shift;

open DNA_TEMP, ">$temp_dir/$loci_name.dna.fas";
print DNA_TEMP ">$loci_name\n$uncorrected_seq\n";
close DNA_TEMP;

open AA_TEMP, ">$temp_dir/$loci_name.aa.fas";
foreach my $prot_id (sort keys %$protseq) {
print AA_TEMP ">$prot_id\n$protseq->{$prot_id}\n";
}
close AA_TEMP;

my $exonerate_out = `exonerate -m protein2dna $temp_dir/$loci_name.aa.fas $temp_dir/$loci_name.dna.fas --score 0 --showvulgar FALSE --showalignment FALSE --subopt FALSE --ryo "sugar: %S %ps %tal %tl\n"`;
my @sugar = $exonerate_out =~ /\n(sugar.+)/g;

my %qualified_frame;
foreach my $sugar (@sugar) {
	my @fields = $sugar =~ /\S+/g;

	my $frame = ($fields[-7]%3)+1;
	my $strand = $fields[-5];
	my $aln_score = $fields[-4];
	my $id = $fields[-3];
	my $cov = sprintf("%.2f", $fields[-2]/$fields[-1]*100);

	if (($id >= $id_trsd)&&($cov >= $cov_trsd)&&($strand eq "+")) {
	$qualified_frame{$aln_score} = $frame;
	}
}

my $corrected_frame = 0;
if (%qualified_frame) {
my @sort_score = sort {$a<=>$b} keys %qualified_frame;
$corrected_frame = $qualified_frame{$sort_score[-1]};
}

return $corrected_frame;
}

# subroutine to find frame with least stop codon
sub no_ref_correct {
my $uncorrected_seq = shift;

my %corrected_frame;
foreach my $frame (1..3) {
	my $prot = translate($uncorrected_seq, $frame);
	my ($corrected_prot, $stopnum) = rm_trail_stop($prot);

	if ($stopnum <= $stop_trsd) {
	$corrected_frame{$frame}->{stop_codon} = $stopnum;
	$corrected_frame{$frame}->{corrected_prot} = $corrected_prot;
	} else {
	next;
	}
}

my @corrected_frame = sort {$corrected_frame{$a}->{stop_codon}<=>$corrected_frame{$b}->{stop_codon}} keys %corrected_frame;

my ($corrected_frame, $corrected_prot, $corrected_seq);
if (@corrected_frame == 1) { # one frame left in consideration
$corrected_frame = $corrected_frame[0];
$corrected_prot = $corrected_frame{$corrected_frame}->{corrected_prot};
$corrected_seq = substr $uncorrected_seq, $corrected_frame-1, length($corrected_prot)*3;
} elsif (@corrected_frame > 1) { # multiple frame left in consideration, keep the frame with fewest stop codon
$corrected_frame = $corrected_frame[0];
$corrected_prot = "";
$corrected_seq = "";
} elsif (@corrected_frame == 0) { # no frame is good
$corrected_frame = "";
$corrected_prot = "";
$corrected_seq = "";
}

return ($corrected_frame, $corrected_prot, $corrected_seq);
}

# subroutine to reformat input fasta file
sub readfas {
my $file = shift; # file name of input fasta

# open fasta seq
open INFILE, $file or die "\nERROR: Cannot find input file \"$file\" ($!)\n\n";

# hash of sequences and correspond taxa, hash of sequence length, input order of sequence, sequence, taxon, last taxon, number of taxa occur multiple time
my (%seq, %seqlength, @input_order, $seq, $taxon, $lasttaxon, %taxa_cnt);

# extract sequence from input file and check whether format is correct or strange characters in it
while (my $line = <INFILE>) {

	# remove enter in windows
	$line =~ s/\r//g;
	chomp $line;

	if ($line =~ /^>(\S+)/) { 	# if find > at the beginning of the line, after ">" is a taxon name
		$taxon = $1;

		# if find previous taxon name, and it hasn't recorded before
		if ($lasttaxon) {
			if (!(exists $seq{$lasttaxon})) {
				$seq =~ s/\s//g; # remove gap
				$seq = uc $seq; # uppercase
				my @strange_char = $seq =~ /[^A-Z?*\.\-]/g; # get all non-nucleotide char
				my $length = length($seq); # seq length
				if (($length > 0)&&(@strange_char == 0)) { # if DNA sequence looks fine
				$seq{$lasttaxon} = $seq; # save seq into hash
				push @input_order, $lasttaxon; # save input order
				$seqlength{$length} = ""; # save sequence length
				} else { # warn if sequence looks strange
					if ($length == 0) { # warn if no sequence under this taxon
					print STDERR "WARNING: No nucleotide found sequence of $lasttaxon in $file. This taxon will be discarded.\n";
					}
					if (@strange_char > 0) { # warn if strange char found in sequence
					my $strange_char = join " ", @strange_char;
					print STDERR "WARNING: Found strange character \"$strange_char\" in sequence of $lasttaxon in $file. This taxon will be discarded.\n";
					}
				}
			} else { # count number of taxon occur multiple times
			$taxa_cnt{$lasttaxon}++
			}
		}

		$seq = "";
		$lasttaxon = $taxon;
	} else { # else is sequences
	$seq .= $line;
	}
}

if ($lasttaxon) { # also need to take care sequence of the last taxon
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

# warn if multiple taxa have the same name
foreach my $taxon (sort keys %taxa_cnt) {
my $taxon_num = $taxa_cnt{$taxon};
$taxon_num++;
print STDERR "WARNING: Found $taxon_num \"$taxon\" in $file. Only first sequence named in \"$taxon\" will be kept\n";
}

# save length of sequences if all sequences are in the same length which may probably an alignment
my @seqlength = sort keys %seqlength;

my $aln_len;
if (@seqlength == 1) {
$aln_len = $seqlength[0];
} else {
$aln_len = "na";
}

# return hash with sequence and corresponding taxa, input order of taxa, alignment length
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
Script name: predict_frames.pl

This is a script to predict frame and generate coding and amino acid sequences of baits sequences using three files related to a Gene Capture dataset; (1) a bait sequence fasta, (2) a fasta file containing amino acid sequences with ENSMBL geneIDs for headers, (3) onehitCDSmarkers file generated from Evolmarker that is used to locate reference protein of each locus

Overview:
Each bait sequence is translated in three frames, if more than one frame left in consideration. Pairwise alignment between the translated coding sequences and the known protein sequence found from onehitCDSmarkers file is used to find the frame.

Dependencies:
(1) Bio::Seq (included in Bioperl)
(2) Exonerate v2.2.0 or higher

Example usage:

	perl predict_frames.pl --baits species.fas --cds species.onehitCDSmarkers.column1.txt --ref_prot species.pep.fas


Input files:
(1) species.fas
(2) species.onehitCDSmarkers.column1.txt
(3) species.pep.fas

Output files:
(1) species.dna.fas
(2) species.aa.fas
(3) frame_result.txt (frame prediction result for each locus)
(4) uncorrected_loci.txt (list of uncorrected loci)

Options:
--baits
  Fasta file containing bait names and sequences
--cds
  OnehitCDSmarker generated from Evolmarker. Only info in first columns will be used, so just input file with only first column is fine. Please specify --cds_header if OnehitCDSmarker has header in first line
--ref_prot
  Reference protein sequences mined from ENSMBL
--dna_out
  Coding sequences output in fasta format, name as 'xxx.dna.fas' in default, xxx is the name prefix of input file
--aa_out
  AA sequences output in fasta format, name as 'xxx.aa.fas' in default, xxx is the name prefix of input file
--stop_codon
  Maximum stop codons allowed in amino acid sequences, $stop_trsd in default
--id
  Minimum identity required between baits sequences and reference protein, $id_trsd in default
--cov
  Minimum percentage of baits sequences aligned with reference protein required, $cov_trsd in default
--cds_header
  Whether is a header at the first line of OnehitCDSmarker file. Specify this option if header exists at the first line of OnehitCDSmarker file, off in default
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
