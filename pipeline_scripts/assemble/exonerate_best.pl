#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;
use Parallel::ForkManager; 
no warnings 'experimental::smartmatch';

# AA sequences of reference, input dir of contig and graph files, sample list, output dir of filtered reads
my ($queryp, $subject, $samplelist, $outdir, $help);
my $error_rate= 0.05; # maximum error rate allowed to consider two sequences aligned 
my $overlaplowtrsd = 9; # minimum overlap between contig
my $overlapuptrsd = 25; # maximum overlap between contig
my $intron_trsd = 50; # maximum size of intron (unaligned right terminal of left contig+unaligned left terminal of right contig) allowed
my $simi_trsd = 60; # minimum similarity required between query and contig
my $st_trsd = 0; # maximum stop codon allowed in contig
my $dropoff = 2; # number of unaligned nucleotide allowed at the most left and right of the non-flanking region
my $connect_trsd = 5; # remove node in overlap graph if it connect more than 5 contigs with "na" (overlap detected in exonerate_best)
my $process = 1; # number of cpu used

my $opt = GetOptions( 'queryp:s', \$queryp,
                      'assembled:s', \$subject,
                      'samplelist:s', \$samplelist,
                      'filtered:s', \$outdir,
                      'error_rate:s', \$error_rate,
                      'min_overlap:s', \$overlaplowtrsd,
                      'max_overlap:s', \$overlapuptrsd,
                      'similarity:s', \$simi_trsd,
                      'stop_codons:s', \$st_trsd,
                      'cpu:s', \$process,
                      'help|h!', \$help) or usage(); 

# display help message if one of the following option is not specified
if (!($opt && $queryp && $subject && $samplelist && $outdir) || $help) {
usage();
} 

# check the existence of AA sequences of reference
if (!(-e $queryp)) {
die ("ERROR: Cannot open $queryp ($!)");
}

# check the existence of dir of contigs and graph files
if (!(-e $subject)) {
die ("ERROR: Cannot open $subject ($!)"); 
}

# get the name of samples from sample list
my %sample;
open SAMPLELIST, "$samplelist" or die "ERROR: Cannot find $samplelist ($!)";
while (my $samplename = <SAMPLELIST>) {
chomp $samplename;
$sample{$samplename} = "";
}
close SAMPLELIST;

# create dir containing AA sequences of query and construct a hash{gene} = seq
my %queryp;
if (-e "queryp") { # save hash only if queryp exists
	open QUERY, $queryp;
	while (my $line = <QUERY>) {
		if ($line =~ />(\S+)/) {
		chomp(my $seq = <QUERY>);
		$queryp{$1} = $seq;
		}
	}
	close QUERY;
} else { # save hash and create dir only if queryp not exists
	mkdir "queryp";
	open QUERY, $queryp;
	while (my $line = <QUERY>) {
		if ($line =~ />(\S+)/) {
		chomp(my $seq = <QUERY>);
		open QUERYPOUT, ">queryp/$1.fas";
		$queryp{$1} = $seq;
		print QUERYPOUT ">$1\n$seq\n";
		close QUERYPOUT;
		}
	}
	close QUERY;
}

# set the dir name of contig and asqg file
my $contigdir = "$subject/contig";
my $asqgdir = "$subject/asqg";

my $subject_list = get_filelist($contigdir); # get the name of all contig dir

# set names and create dir for exonerate output, filter contigs and overlap information
my $exonerateout = "$outdir/exonerate";
my $sortout = "$outdir/filtered";
my $overlapdir = "$outdir/overlap";

# create dir
mkdir $outdir; # output dir
mkdir $exonerateout; # exonerate result
mkdir $sortout; # filtered contigs
mkdir $overlapdir; # overlapping information

# split reads into $split chunk
my $split = $process;

function("filtering contigs");

# filter genes for each sample
foreach my $subject_sample (@$subject_list) {
next if (!(exists $sample{$subject_sample})); # only filter contigs in samplelist
print "Start filtering contigs in $subject_sample\n";

# make outdir foreach sample
mkdir "$exonerateout/$subject_sample";
mkdir "$sortout/$subject_sample";
mkdir "$overlapdir/$subject_sample";

# get gene list of sample
my $gene_list = get_filelist("$contigdir/$subject_sample");
map {$_ =~ s/-contigs\.fa//g} @$gene_list;
	
	my $splited_array = split_array($gene_list, $split); # split contig into $split part
	
	my $pm = Parallel::ForkManager->new(int($process/2));
	
	# filter contig in each split in parallel
	DATA_LOOP: while (my $array = shift @$splited_array) {
	$pm->start and next DATA_LOOP;
		
		# filter for each gene in each chunk
		while(my $gene = shift @$array) {
		# set name of contig, asqg, reference and outfile
		my $target  = $gene . "-contigs.fa";
		my $asqg = $gene . "-graph.asqg";
		my $queryp = $gene . ".fas";
		my $out = $gene . ".txt";
			
			# only analysis file with contigs in it
			if ((stat("$contigdir/$subject_sample/$target"))[7] > 0) {
			# run exonerate to align contigs against query
			`exonerate -m protein2dna queryp/$queryp $contigdir/$subject_sample/$target --score 0 --showvulgar FALSE --showalignment FALSE --subopt FALSE --ryo "sugar: %S %ps\npt: {%Pqs_%Pts_%Ps\t}\n" > $exonerateout/$subject_sample/$out`;
				
				# save contigs and its name in hash
				my %oricontig;
				open ORICONTIG, "$contigdir/$subject_sample/$target";
				while (my $line = <ORICONTIG>) {
					if ($line =~ />(\S+)/) {
					chomp (my $seq = <ORICONTIG>);
					$oricontig{$1} = $seq;
					}
				}
				close ORICONTIG;
		   	
		   		# read exonerate output
				open EXONERATE, "$exonerateout/$subject_sample/$out";
				my (%sub, $sugarlinecount);
				<EXONERATE>;
				while (my $line = <EXONERATE>) {
					if ($line =~ /sugar\:/) { # find the line got alignment information
						$sugarlinecount++;
						# get the start point on query, end point on query, contig name, start point on contig, end point on contig, align in same direction as query or reversed, alignment score and similarity with query
						my ($st_q, $ed_q, $sub, $st_s, $ed_s, $strand, $score, $similarity) = $line =~ /\S+\s\S+\s(\d+)\s(\d+)\s\S+\s(\S+)\s(\d+)\s(\d+)\s([+|-])\s(\d+)\s(\S+)/;
							
							# keep best hit if contig got multiple hit
							if (exists $sub{$sub}) {
							next if ($score < $sub{$sub}->{score});
							}
						
						# remove contigs with too low similarity with query
						next if ($similarity < $simi_trsd);
						
						# modify start and end point on contig and query into same format (start < end, 0-start)
						($ed_s, $st_s) = ($st_s , $ed_s) if ($st_s > $ed_s);
						$ed_q--;
						$ed_s--;
						
						# get the length of whole seq
						my $seq = $oricontig{$sub};
						my $wlen_s = length($seq);
						
						# modify all sequence to the same orientation
						my ($stq, $edq, $stsm, $edsm, $seqm) = do 
						{
							if ($strand eq "+") {
							$st_q, $ed_q, $st_s, $ed_s, $seq;
							} else {
							$seq =~ tr/ATCG/TAGC/;
							$seq = reverse $seq;
							$st_q, $ed_q, $wlen_s-$ed_s-1, $wlen_s-$st_s-1, $seq;
							}    
						};
						
						# get the non-flanking seq (aligned seq), remove it if it got lots of stop codon
						my $nfseq = substr $seqm, $stsm, $edsm-$stsm+1;
						my $nfseqp = translate($nfseq, 1);
						my @stcodon = $nfseqp =~ /\*/g;
					
						next if (@stcodon > $st_trsd);
					
						# record score start and end point on query and contigs and sequences into hash
						$sub{$sub}->{score} = $score;
						$sub{$sub}->{stq} = $stq;
						$sub{$sub}->{edq} = $edq;
						$sub{$sub}->{stsm} = $stsm;
						$sub{$sub}->{edsm} = $edsm;
						$sub{$sub}->{seq} = $seqm;
					}
				}
				close EXONERATE;
			
				next if (!($sugarlinecount)); # skip file if query do not have any hit
			
				# write output there are contigs hit with query
				if (scalar keys %sub > 0) {			
				open OUTFILE, ">$sortout/$subject_sample/$gene.fas";
					
					# write start and end point on query and contigs, score and sequences of each gene
					if (scalar keys %sub == 1) { # if there's only one contig in file
					my @sub = keys %sub;
					print OUTFILE ">$sub[0]\|$sub{$sub[0]}->{stq}\|$sub{$sub[0]}->{edq}\|$sub{$sub[0]}->{stsm}\|$sub{$sub[0]}->{edsm}\|$sub{$sub[0]}->{score}\n$sub{$sub[0]}->{seq}\n";
					} else {
						# get the overlap information from asqg file
						open ASQG, "$asqgdir/$subject_sample/$asqg";
						<ASQG>;
						my %sgaoverlap;
						while (my $line = <ASQG>) {
							if ($line =~ /ED/) {
								my @edgeinfo = $line =~ /ED\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s\S+\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s\S+/;
								my $sub1 = $edgeinfo[0]; # get overlapped contig 1
								my $sub2 = $edgeinfo[1]; # get overlapped contig 2
								next if (!(exists $sub{$sub1})); 
								next if (!(exists $sub{$sub2}));
								$sgaoverlap{"$sub1\|$sub2"} = $edgeinfo[3]-$edgeinfo[2]+1; # get the length of its overlap
							}
						}
						close ASQG;
						
						# get all overlapped contigs in filter contigs recorded in %sub, and return a hash containing overlapping information in $hash{contig1|contig2} = ""
						my $overlap = find_overlap(\%sub, \%sgaoverlap, $gene);
							
						# write out each overlapping contigs		
						my %delete;
						if (scalar keys %$overlap > 0) {
							open OVERLAP, ">$overlapdir/$subject_sample/$gene.txt";
							foreach (sort keys %$overlap) {
							print OVERLAP "$_\t$overlap->{$_}\n"; # write out overlap information
							my ($overlap1, $overlap2) = $_ =~ /(\S+)\|(\S+)/; # get the name of overlapping contigs
							# write out overlapping contigs and do not write them repeatly if it is flagged in %delete
							print OUTFILE ">$overlap1\|$sub{$overlap1}->{stq}\|$sub{$overlap1}->{edq}\|$sub{$overlap1}->{stsm}\|$sub{$overlap1}->{edsm}\|$sub{$overlap1}->{score}\n$sub{$overlap1}->{seq}\n" if (!(exists $delete{$overlap1})); 
							print OUTFILE ">$overlap2\|$sub{$overlap2}->{stq}\|$sub{$overlap2}->{edq}\|$sub{$overlap2}->{stsm}\|$sub{$overlap2}->{edsm}\|$sub{$overlap2}->{score}\n$sub{$overlap2}->{seq}\n" if (!(exists $delete{$overlap2}));
							$delete{$overlap1} = 1; # flag contig $overlap1 since it has been written
							$delete{$overlap2} = 1; # flag contig $overlap2 since it has been written
							}				    
							close OVERLAP;
						}
				
						# write out contigs with highest score
						my @sortsub = sort {$sub{$a}->{score}<=>$sub{$b}->{score}} keys %sub;
						if (!(exists $delete{$sortsub[-1]})){
						print OUTFILE ">$sortsub[-1]\|$sub{$sortsub[-1]}->{stq}\|$sub{$sortsub[-1]}->{edq}\|$sub{$sortsub[-1]}->{stsm}\|$sub{$sortsub[-1]}->{edsm}\|$sub{$sortsub[-1]}->{score}\n$sub{$sortsub[-1]}->{seq}\n";
						}
					}
				close OUTFILE;
				}  			
			}
		}
	$pm->finish;
	}
	$pm->wait_all_children;
	
print "Contigs of $subject_sample has been filtered and written to $outdir/filtered/$subject_sample/ \n";
print "Exonerate results of $subject_sample has been filtered and written to $outdir/exonerate/$subject_sample/ \n";
print "Names of overlapped contigs of $subject_sample has been filtered and written to $outdir/overlap/$subject_sample/ \n";
}

print "All contigs have been filtered and written to $outdir/filtered/ \n";
print "All exonerate results have been written to $outdir/exonerate/ \n";
print "All names of overlapped contigs have been filtered and written to $outdir/overlap/ \n";

# `rm -rf queryp`;

##################################################
# subroutine
##################################################

# subroutine to get filelist under $dir
sub get_filelist {
my $dir = shift; # input dir

opendir DIR, $dir;
my @filelist = grep {$_ !~ /^\./} readdir DIR;
closedir DIR;

# return list of file
return \@filelist;
}

# subroutine to find overlap between contigs based on position on query
sub find_overlap {
my $subject = shift; # hash of contigs of a sample and their alignment position
my $sgaoverlap = shift; # overlap information extracted from asqg file
my $gene = shift; # loci name

# To conveniently order contigs, we need to modify the start and end point of contigs with flanks.
# so modified start point may negative and end point is longer than end of reference. We also need 
# to avoid to take randomly aligned sequences as contigs with flanks, so we only deem contigs which
# align to the ends of the reference and have unaligned flanks as contigs having flanking sequences.
# However, contigs with several beginning or ending nucleotide cannot align to reference which could
# also have long flanks. Thus, we relaxed $dropoff bp of start and end point of alignment to reference
my $qstpos = $dropoff;
my $qedpos = length($queryp{$gene})*3-1-$dropoff;
    
    my %overlap;
    my @subject = sort keys %$subject;
    while (my $qsub = shift @subject) {
    # get the position information of seq1
	my $seq1qst = $subject->{$qsub}->{stq}; # start point on amino acid query
	my $seq1qed = ($subject->{$qsub}->{edq}+1)*3-1; # end point on amino acid query
	my $seq1sst = $subject->{$qsub}->{stsm}; # start point on contig
	my $seq1sed = $subject->{$qsub}->{edsm}; # end point on contig
	my $seq1lf = $seq1sst; # length of left unaligned sequence
	my $seq1rf = length($subject->{$qsub}->{seq})-$seq1sed-1;  # length of right unaligned sequence
	
	# save length of left and right unaligned sequences in hash
	my %lrf;
	$lrf{$qsub}->{lf} = $seq1lf;
	$lrf{$qsub}->{rf} = $seq1rf;
		
		# save the start and end point of seq1 into hash
		my %sort;
	    $sort{"queryst"} = $seq1qst;
		$sort{"queryed"} = $seq1qed;
		# modify the start and end point of seq1
	    if (($seq1qst <= $qstpos) && ($seq1lf > 0)) {
	    $sort{"queryst"} = $seq1qst-$seq1lf;
	    }
	    if (($seq1qed >= $qedpos) && ($seq1rf > 0)) {
		$sort{"queryed"} = $seq1qed+$seq1rf;
		}
	    
	    foreach my $sub (@subject) {
	    # get the position information of seq2, same as seq1
		my $seq2qst = $subject->{$sub}->{stq};
		my $seq2qed = ($subject->{$sub}->{edq}+1)*3-1;
		my $seq2sst = $subject->{$sub}->{stsm};
		my $seq2sed = $subject->{$sub}->{edsm};
		my $seq2lf = $seq2sst;
		my $seq2rf = length($subject->{$sub}->{seq})-$seq2sed-1;
		
		$lrf{$sub}->{lf} = $seq2lf;
		$lrf{$sub}->{rf} = $seq2rf;
			
			# save the start and end point of seq1 into hash
			$sort{"subst"} = $seq2qst;
			$sort{"subed"} = $seq2qed;
			# modify the start and end point of seq2 which part of them is in flanking region
			if (($seq2qst <= $qstpos) && ($seq2lf > 0)) {
			$sort{"subst"} = $seq2qst-$seq2lf;
			} 
			if (($seq2qed >= $qedpos) && ($seq2rf > 0)) {
			$sort{"subed"} = $seq2qed+$seq2rf;
			}
		    
		    # sort position and find which is left or right sequence
		    my ($cur, $nex);
			my @sort = sort {$sort{$a}<=>$sort{$b}} keys %sort; # sort the positon from small to big
			my ($left) = $sort[0] =~ /(\S+)st/; # get the name of left sequence
			my ($right) = $sort[-1] =~ /(\S+)ed/; # get the name of right sequence
			# if one of the sequence is not contained in other sequence
			if (($left ne $right)&&($sort{$sort[0]} != $sort{$sort[1]})&&($sort{$sort[2]} != $sort{$sort[3]})){
				# determine whether "query" or "sub" as left or right sequence
				($cur, $nex) = do {
					if ($left eq "query") {
					($qsub, $sub);
					} else {
					($sub, $qsub);
					}
				};
			} else {
			next;
			}

			# save the overlapping pair
			
			# if this pair has been detect by sga
			if (exists $sgaoverlap->{"$qsub\|$sub"}){
			$overlap{"$cur\|$nex"} = $sgaoverlap->{"$qsub\|$sub"};
			} elsif (exists $sgaoverlap->{"$sub\|$qsub"}) {     
			$overlap{"$cur\|$nex"} = $sgaoverlap->{"$sub\|$qsub"};
			} else { # if this pair has not been detect by sga
				# if intron sequence is not too long (length of right unaligned part of left sequence + length of left unaligned part of right sequence)
			    if (($lrf{$cur}->{rf}+$lrf{$nex}->{lf}) <= $intron_trsd) {
					my $overlap = ($subject->{$cur}->{edq}-$subject->{$nex}->{stq}+1)*3; # get the length of overlapping sequence
					my $overlapmax = $overlap+$lrf{$cur}->{rf}+$lrf{$nex}->{lf};
					
					if (($overlapmax <= $overlapuptrsd)&&($overlapmax >= $overlaplowtrsd)) { # if length of overlap is shorter than $overlapuptrsd (which should have been detected in sga), and longer than 0 
				   	
				   	my $curolp = substr $subject->{$cur}->{seq}, -$overlapmax; # get the overlap part in left seq
				    my $nexolp = substr $subject->{$nex}->{seq}, 0, $overlapmax; # get the overlap part in right seq
				   				   
				    	# kmer length
						my $kmer_len = $overlaplowtrsd;
						
						# split seq into kmer
						my $curkmer = kmer($curolp, $overlaplowtrsd);
						my $nexkmer = kmer($nexolp, $overlaplowtrsd);
										
				    	# find number of identical kmer in overlap of left and right sequence
				        my $count = 0; 
						foreach (sort keys %$curkmer) {
						$count++ if (exists $nexkmer->{$_});
						}
					 	
					 	# number of error nucleotide allowed
					 	my $error_num = int($overlapmax*$error_rate);
					 	
						# number of kmer in aligned region, which should be almost identical
					 	my $kmer_num = $overlapmax-$kmer_len+1;
					 	$kmer_num = 1 if ($kmer_num < 0);
					 	
					 	# number of least identical kmer between 2 overlap
					 	my $least_correct = $kmer_num-$kmer_len*$error_num;
						
						$least_correct = 1 if ($least_correct < 0);
										
						$overlap{"$cur\|$nex"} = "na" if ($count >= $least_correct); # keep the qualified overlap 
					}
			    }
			}
		}
	}
	
	# delete the node which got connection more than $connect_trsd
	foreach my $sub (sort keys %$subject) {
		my $connect = 0;
		my @deleteoverlap;
		foreach my $overlap (sort keys %overlap) {
			if (($overlap =~ /$sub/)&&($overlap{$overlap} eq "na")) {
			push @deleteoverlap, $overlap;
			$connect++;
			}
		}
		map {delete $overlap{$_}} @deleteoverlap if ($connect > $connect_trsd);
	}

# return overlapping information
return (\%overlap);
}

# subroutine to translate sequences
sub translate {
my $seq = shift @_;
my $start_codon = shift @_;

$start_codon --;
my $seq_obj=Bio::Seq->new(-seq=>$seq,-alphabet=>'dna');
my $pro=$seq_obj->translate(-frame => $start_codon, -terminator => "*", -unknown => "X");
$pro=$pro->seq; 
return($pro);
}

# subroutine to generate kmer
sub kmer {
my $seq = shift; # sequence
my $kmer = shift; # k-mer length

# save k-mer and their occurence
my %khash;
for (my $i=0; $i<(length($seq)-$kmer+1) ;$i++) {
my $str = substr $seq, $i, $kmer;
	if (!(exists $khash{$str})) {
	$khash{$str} = 1;
	} else {
	$khash{$str} ++;
	}
}

# numbering the same k-mer
my %hash;
foreach (sort keys %khash) {
	for (my $i=1; $i<=$khash{$_}; $i++) {
	$hash{"$_.$i"} = 1;
	}
}

# return hash of k-mer
return \%hash;
}

# subroutine to split array in to several sub-array
sub split_array {
my $array = shift; # array of loci name
my $split = shift; # parts to be splited
	
$split++ if ($split == 0); # at least split into 1 part
my $block_num = int(@$array/$split); # number of loci in each split
$block_num++ if ($block_num == 0); # each split have at least 1 part

my @splited_array;
my $size_cnt = 0; # number of loci
my $block_cnt = 0; # number of block

# split loci into $split block
while (my $file = shift @$array) {
	if ($size_cnt >= $block_num) { # if current block is full open another block
	$size_cnt = 0;
	$block_cnt ++;
	}
push @{$splited_array[$block_cnt]}, $file;
$size_cnt ++;
}

# return splited loci set
return(\@splited_array);
}

# subroutine to report the start of a step of program
sub function {
my $function = shift;

if ($process > 1){
print "Start $function in parallel\n";
} else {
print "Start $function in serial\n";
}
}

sub usage {
print STDERR "
Script name: exonerate_best.pl

This is a script to filter unqualified contigs and find contigs which might be furtherly assembled  

Dependencies:
(1) exonerate
(2) Perl Module:
	1. Bio::Seq (Included in Bioperl)
	2. Parallel::ForkManager
	
Example usage:

	perl exonerate_best.pl --queryp query.aa.fas --assembled assembled --filtered exonerateout --samplelist samplelist.txt

Input files:
(1) query.aa.fas
(2) assembled
(3) samplelist.txt

Output files:
(1) exonerateout

Options:
--queryp
  Amino acid sequences of target loci in fasta format
--assembled
  Directory containing assembled reads and corresponding graph files
--samplelist
  A list of sample name. Only samples in the list will be assembled
--filtered
  Directory containing filtered and potentially overlapped contigs
--error_rate
  Maximum error rate allowed to consider two sequences aligned, $error_rate in default
--min_overlap
  Minimum length of overlapped aligned region to target loci between 2 contigs, $overlaplowtrsd in default
--max_overlap
  Maximum length of overlapped aligned region to target loci between 2 contigs, $overlapuptrsd in default
--similarity
  Minimum similarity required between contigs and target loci, $simi_trsd in default
--stop_codons
  Maximum number of stop codon allowed in a contig, $st_trsd in default
--cpu
  Limit the number of CPUs, $process in default. 
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