#!/usr/bin/env perl

use Getopt::Long;
use Parallel::ForkManager;
# use strict;
use warnings;

my $blastunsortdir = "blastout_unsorted"; # dir containing unsorted blast result for each splited deduplicated reads
my $mergeblastoutdir = "blastout"; # dir containing merged and sorted blast result for each sample

# AA sequences of reference, DNA sequences of reference, input dir of deduplicated reads, sample list, outdir of parsed reads, word size, do not mask reads, keep merged deduplicated reads
my ($queryp, $queryn, $subdir, $samplelist, $outdir, $wordsize, $unmask, $keepmerge, $help);
my $nucleo_masker = "dust"; # nucleotide masker
my $aa_masker = "seg"; # amino acid masker
my $E = 0.0001; # default e-value
my $maxaccepts = 0; # terminate ublast search after find $maxaccepts hit, accept all if $maxaccepts == 0 and $maxrejects == 0
my $maxrejects = 0; # terminate ublast search after reject $maxrejects input reads, accept all if $maxaccepts == 0 and $maxrejects == 0
my $dbstep = 1; # Specifies that every Nth database word should be indexed
my $masktrsd = 1; # do not include reads pair in output if one of reads is $masktrsd% masked
my $threads = 1; # number of cpu used

my $opt = GetOptions( 'queryp:s', \$queryp,
		              'queryn:s', \$queryn,
                      'rmdup:s', \$subdir,
                      'samplelist:s', \$samplelist,
                      'parsed:s', \$outdir,
                      'unmask!', \$unmask,
                      'E_value_ubandpx:s', \$E,
                      'mask_percentage:s', \$masktrsd,
#                       'keepmerge:s', \$keepmerge,
                      'cpu:s', \$threads,
                      'help|h!', \$help) or usage(); 

# display help message if one of the following option is not specified
if (!($opt && ($queryp||$queryn) && $subdir && $samplelist && $outdir) || $help) {
usage();
} 

# take one of DNA/AA sequences of reference as database

# check whether AA sequences of reference exists, if --queryp is specified
my $query;
if ($queryp) {
	if (!(-e $queryp)) {
	die ("ERROR: Cannot find AA sequences of reference $queryp ($!)");
	} else {
	$query = $queryp; 
	}
}
# check whether DNA sequences of reference exists, if --queryn is specified
if ($queryn) {
	if (!(-e $queryn)) {
	die ("ERROR: Cannot find DNA sequences of reference $queryn ($!)");
	} else {
	$query = $queryn;
	}
}

# check whether dir od deduplicated reads exists
if (!(-e $subdir)) {
die ("ERROR: Cannot open dir of deduplicated reads $subdir ($!)"); 
}

# get samples which need to be parsed 
my %sample;
open SAMPLELIST, "$samplelist" or die "ERROR: Cannot find $samplelist ($!)";
while (my $samplename = <SAMPLELIST>) {
chomp $samplename;
$sample{$samplename} = "";
}
close SAMPLELIST;

# determine whether reference is AA or DNA
open QUERY, $query;
<QUERY>;
chomp(my $seq = <QUERY>);
my @nucleo = $seq =~ /[A|T|C|G]/g;
close QUERY;

# determine db type, which means we take reference as database
my $dbtype;
if (@nucleo/length($seq) >= 0.5) { # nucleotide
$dbtype = "nucleo";
} else { # protein
$dbtype = "prot";
}

# determine masker and word size
my $qmasker;
if ($dbtype eq "nucleo") {
$qmasker = $nucleo_masker;
$wordsize = 8;
} else {
$qmasker = $aa_masker;
$wordsize = 5;
}

# read each deduplicated reads, and only include sample in sample list
my $subdirfq = $subdir;

# array of all deduplicated reads need to be parsed, $sub{sample name} = array of corresponding reads
my (@fq, %sub);
opendir SUBDIRFQ, $subdirfq;
while (my $fq = readdir SUBDIRFQ) {
next if ($fq =~ /^\./);
next if ($fq =~ /masked/);
my ($sub) = $fq =~ /(\S+)\.\d+\.rmdup\.fq/;
	if (exists $sample{$sub}) {
	push @fq, $fq;
	push @{$sub{$sub}}, $fq;
	}
}
closedir SUBDIRFQ;

maskdedup() if (!($unmask)); # mask deduplicated reads if --unmask is not specified
mkquerydb(); # make database of deduplicated reads 
local_usearch(); # find reads hit with query

# make dir for merged deduplicate reads
my $mergermdupdir = $subdir . "_merge";
mkdir $mergermdupdir;

# make dir for ublast output
mkdir $mergeblastoutdir;

# merge splited deduplicated reads and its ublast output for each sample
foreach my $sub (sort keys %sub) {
	# array of splited deduplicated reads, array of splited blastout file
	my (@fqfile, @blastoutfile);
	foreach my $subfile (@{$sub{$sub}}) {
	my $fqfile = "$subdirfq/$subfile";
	my $fqfilemask = "$subdirfq/$subfile.masked";
	my $blastoutfile = "$blastunsortdir/$subfile.blast.unsort.txt";
		if (!($unmask)) {
		push @fqfile, $fqfilemask;
		} else {
		push @fqfile, $fqfile;
		}
	push @blastoutfile, $blastoutfile;
	}

	# merge splited deduplicated reads
	my $catfqfile = join " ", @fqfile;
	if (!($unmask)) {
	`cat $catfqfile > $mergermdupdir/$sub.rmdup.fq.masked`;
	} else {
	`cat $catfqfile > $mergermdupdir/$sub.rmdup.fq`;
	}

	# merge unsorted ublast output, keep reads with best hit to gene, sort it
	my $blastout_file = join " ", @blastoutfile;
	
	my $blastout_unsort = "$mergeblastoutdir/$sub.blast.unsort.txt";
	`cat $blastout_file > $blastout_unsort`;
	
	my $blastout_uniq = uniq_seq($blastout_unsort, $mergeblastoutdir, $sub);
	`sort -o $mergeblastoutdir/$sub.blast.txt $blastout_uniq`;
	
	unlink $blastout_unsort;
	unlink $blastout_uniq;
}

sub uniq_seq {
my $blastout = shift;
my $mergeblastoutdir = shift;
my $sub = shift;

open BLASTOUT, $blastout; # open blast result

# hash of reads id hit with the same loci, loci name of previous line
my %hitid;

while (my $line = <BLASTOUT>) {
# get loci and reads name
my ($geneid, $hitid, $bits) = $line =~ /(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)/;

	if (exists $hitid{$hitid}) {
	my ($pre_bits) = $hitid{$hitid} =~ /\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/;
	
	$hitid{$hitid} = $line if ($bits > $pre_bits);
	} else {
	$hitid{$hitid} = $line;
	}
}
close BLASTOUT;

my $blastout_uniq = "$mergeblastoutdir/$sub.blast.uniq.txt";
open UNIQ_OUT, ">$blastout_uniq";
foreach my $hitid (keys %hitid) {
print UNIQ_OUT "$hitid{$hitid}";
}
close UNIQ_OUT;

return $blastout_uniq;
}

# parse reads based on ublast output
if (!($unmask)) {
writemaskseq();
} else {
writeseq();
}

# remove intermediate files, unsorted ublast output and database of query
`find $subdirfq -mindepth 1|egrep -v "fq\$"|xargs rm -f`;
`rm -rf $blastunsortdir`;
unlink "$query.udb";

# remove merged deduplicated reads
# if (!($keepmerge)) {
# `rm -rf $mergermdupdir`;
# }

##################################################
# subroutine
##################################################

# subroutine to mask deduplicated reads
sub maskdedup {
function("masking deduplicated reads");

foreach my $fq (@fq) {
print "Start masking $fq\n";
my $in = "$subdirfq/$fq";
my $maskfq = "$subdirfq/$fq.masked";
`usearch -fastx_mask $in -fastqout $maskfq -qmask $nucleo_masker -threads $threads &> /dev/null`;
print "$fq has been masked\n";
}

print "All deduplicated reads has been masked\n";
}

# subroutine to make database of query
sub mkquerydb { 
    function("constructing database of query sequences");
	
	my $in = "$query";
	my $udb = "$query.udb";
	`usearch -makeudb_ublast $in -output $udb -wordlength $wordsize -dbmask $qmasker -dbstep $dbstep &> /dev/null`;
	
	print "Database of query sequences has been constructed\n";
}

# subroutine to ublast deduplicated reads against reference database
sub local_usearch {
    function("ublast search");
        
    mkdir $blastunsortdir;
    foreach my $fq (@fq) {
    print "Start to ublast $fq against $query.udb\n";
	my $querydb = "$query.udb";
	my $input = do {
		if ($unmask) {
		"$subdirfq/$fq";
		} else {
	    "$subdirfq/$fq.masked";
		}
	};
	my $blastout_unsort = "$blastunsortdir/$fq.blast.unsort.txt";
	`usearch -ublast $input -qmask user -db $querydb -evalue $E -maxaccepts $maxaccepts -maxrejects $maxrejects -userfields target+query+qlo+qhi+bits -userout $blastout_unsort -threads $threads &> /dev/null`;
	}
	
	print "Ublast search has been finished\n";
}

# subroutine to write out unmask reads
sub writeseq {
function("writing parsed reads to genes");

# initiate multi-process
my $pm = Parallel::ForkManager->new(int($threads));

# write parsed reads in parallel
mkdir $outdir;
DATA_LOOP: foreach my $sub (sort keys %sub) {
$pm->start and next DATA_LOOP;
print "Start writing parsed reads of $sub\n";

	my $fq = "$sub.rmdup.fq"; # set the name of deduplicated fastq
	my $fqpath = "$mergermdupdir/$fq"; # set the path of deduplicated fastq
	my $blastout = "$mergeblastoutdir/$sub.blast.txt"; # set the name of merged and sorted ublast output
	my $parsedreadsdir = "$outdir/$sub"; # set the name of output dir for $sub sample

	mkdir $parsedreadsdir;

	my $indexfq = indexfq($fqpath); # make index of each read (position of read in file) in $fq
	
	# write pair of reads to gene
	open BLASTOUT, $blastout; # open blast result
	open FASTQ, $fqpath; # open reads
	
	# hash of reads id hit with the same loci, loci name of previous line
	my (%hitid, $geneidlag);
	while (my $line = <BLASTOUT>) {
	# get loci and reads name
	my ($geneid, $hitid) = $line =~ /(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+/;
		
		# if loci name of current line is different from previous line
		if (($geneidlag) && ($geneidlag ne $geneid)) {
		printfq ("$parsedreadsdir/$geneidlag", \%hitid, $indexfq); 	# write reads of loci from previous line
		%hitid = (); # empty hash
		} 
		$hitid{$hitid} = ""; # save reads id
		$geneidlag = $geneid; # update loci name of previous line
	}
	printfq ("$parsedreadsdir/$geneidlag", \%hitid, $indexfq); # print reads of the last loci
	close BLASTOUT;
	close FASTQ;

print "Parsed reads of $sub have been written to $parsedreadsdir/\n";
$pm->finish();
}
$pm->wait_all_children();

print "All parsed reads has been written to $outdir/\n";
}

# subroutine to write out unmask reads
sub writemaskseq {
function("writing parsed reads to genes");

# initiate multi-process
my $pm = Parallel::ForkManager->new(int($threads));

# write parsed reads in parallel
mkdir $outdir;
DATA_LOOP: foreach my $sub (sort keys %sub) {
$pm->start and next DATA_LOOP;
print "Start writing parsed reads of $sub\n";

	my $fq = "$sub.rmdup.fq.masked"; # set the name of masked reads
	my $fqpath = "$mergermdupdir/$fq"; # set the path to masked reads
	my $blastout = "$mergeblastoutdir/$sub.blast.txt"; # set the path to merged and sorted ublast output
	my $parsedreadsdir = "$outdir/$sub"; # set the name of output dir for $sub sample 

	mkdir $parsedreadsdir;

	my $indexfq = indexfq($fqpath); # make index of each read (position of read in file) in $fq
	
	# write pair of reads to targeting loci
	open BLASTOUT, $blastout; # open blast result
	open FASTQ, $fqpath; # open reads
	
	# hash of reads id hit with the same loci, hash of reads id hit which have few lowercase reads (low complexity region), loci name of previous line
	my (%hitid, $hithash, $geneidlag);
	while (my $line = <BLASTOUT>) {
	my ($geneid, $hitid, $subst, $subed) = $line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
		if (($geneidlag) && ($geneidlag ne $geneid)) {
		$hithash = mask_filter(\%hitid, $indexfq); # filter a pair of reads if one of tham is masked more than $masktrsd %
		printfq ("$parsedreadsdir/$geneidlag", $hithash, $indexfq) if ((scalar keys %$hithash) > 0);
		%hitid = ();
		} 
		my ($readid, $orient) = $hitid =~ /(\S+)\/(\d)/;
		$hitid{$readid}->{$orient}->{subst} = $subst; # record the start of the hit of $readid
		$hitid{$readid}->{$orient}->{subed} = $subed; # record the end of the hit of $readid
		$geneidlag = $geneid; # update loci name of previous line
	}
	# do the same for last loci
	$hithash = mask_filter(\%hitid, $indexfq);
	printfq ("$parsedreadsdir/$geneidlag", $hithash, $indexfq) if ((scalar keys %$hithash) > 0);
	close BLASTOUT;
	close FASTQ;

print "Parsed reads of $sub have been written to $parsedreadsdir/\n";
$pm->finish();
}
$pm->wait_all_children();

print "All parsed reads has been written to $outdir/\n";

return \%sub;
}

# subroutine to index fastq file, return a hash with the position of each read in file
sub indexfq {
my $infile = shift; # reads

# get reads postion in $infile
open INFILE, $infile;
my ($id, %index, $len_noid);
my $file_pointer = 0;
while (my $line = <INFILE>) {
	if ($line =~ /(\@\S+\/\d)/) { # if we meet beginning position of reads
		# record position of pointer of last read
		if ($id) {
		$index{$id} = $file_pointer;
		$file_pointer = $file_pointer + length($id) + $len_noid + 1;
		} 
		$id = $1; # update reads id
	} else { # record length of other text
	$line .= <INFILE>;
	$line .= <INFILE>;
	$len_noid = length($line);
	}
}
$index{$id} = $file_pointer; # take care of last read
close (INFILE);

# hash of reads postion in $infile
return \%index;
}

# subroutine to write out reads
sub printfq {
my $path = shift; # output path
my $hitid = shift; # loci name
my $index = shift; # hash of postion of reads in fq file

# write out pair of reads
my ($R1, $R2);
open $R1, ">", "${path}_R1.fq";
open $R2, ">", "${path}_R2.fq";
foreach my $hitid (sort keys %$hitid) {
	foreach my $R1or2 (1..2) {
		seek FASTQ, $index->{"\@" . $hitid . "\/$R1or2"}, 0; # go to the position of reads
		my $header = <FASTQ>;  # first line is header, skip
		my $seq = <FASTQ>;	# second line is the whole sequence   
		$seq = uc($seq);
		my $plus = <FASTQ>;	# third line is "+"
		my $qua = <FASTQ>;	# fourth line is the quality score 
	
		if ($R1or2 == 1) {
		print $R1 $header . $seq . $plus . $qua; 
		} elsif ($R1or2 == 2) {
		print $R2 $header . $seq . $plus . $qua; 
		}
	}
}
close $R1;
close $R2;
}

# subroutine to filter reads if hit region is masked more than $masktrsd %
sub mask_filter {
my $hitid = shift; # hash of reads targeting the same loci
my $indexfq = shift; # reads postion in fq file
    
    my %maskfilter;
    foreach my $id (sort keys %$hitid) {
        my $count = 0;
		map {
		seek FASTQ, $indexfq->{"\@$id\/$_"}, 0; # go to the position of reads
		<FASTQ>;  # first line is header, skip
		chomp(my $maskseq = <FASTQ>);	# second line is the whole sequence 
		# get part of reads aligned to reference
		my $hitregion = do {
			if (exists $hitid->{$id}->{$_}) {
			my $subst = $hitid->{$id}->{$_}->{subst};
			my $subed = $hitid->{$id}->{$_}->{subed};
			substr($maskseq, $subst-1, ($subed-$subst+1)); 
			} else {
			$maskseq;
			}
		};
		# see how much low-complexity region in the reads
		my @lc = $hitregion =~ /[a-z]/g;
		$count++ if (@lc/length($hitregion) <= $masktrsd);
		} (1..2);
		# only get reads if not much low-complexity region in paired reads 
		$maskfilter{$id} = "" if ($count == 2);
	}

# return filter reads
return(\%maskfilter);
}

# subroutine to report the start of a step of program
sub function {
	my $function = shift;

	if ($threads > 1){
	print "Start $function in parallel\n";
	} else {
	print "Start $function in serial\n";
	}
}

sub usage {
print STDERR "
Script name: ubxandp.pl

This is a script to parse reads (DNA) to loci (AA or DNA) by ublast under usearch

Dependencies:
(1) USEARCH v10.0.240 or higher
(2) Perl Module:
	Parallel::ForkManager

Example usage:

	perl ubxandp.pl --queryp query.aa.fas --rmdup rmdup --parsed bandpout --samplelist samplelist.txt

Input files: 
(1) query.aa.fas
(2) rmdup
(3) samplelist.txt

Output files:
(1) bandpout

Options:
--queryp
  Amino acid sequences of target loci in fasta format, either DNA or AA sequences can be reference
--queryn
  Full coding nuleotide sequences of target loci in fasta format, either DNA or AA sequences can be reference
--rmdup
  Directory containing reads without PCR duplicates
--samplelist
  A list of sample name. Only samples in the list will be assembled
--parsed
  Directory containing parsed reads
--unmask
  Doesn\'t mask low complexity region in reads, disabled in default
--E_value_ubandpx
  Maximum e-value required for a hit, $E in default
--mask_percentage
  Percentage of masked nucleotide accecpted in a reads, $masktrsd in default
--cpu
  Limit the number of CPUs. $threads in default
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