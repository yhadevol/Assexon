#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager; 

# input dir of parsed reads, sample list, outdir of assembled contigs and string graph files
my ($dir, $samplelist, $outdir, $help);
my $minolp = 25; # minimum overlap require between contigs
my $errorate = 0.05; # maximum error rate allowed to consider two sequences aligned 
my $terminal_len = 110; # terminal branch in string graph will be trimmed if it is shorter than $terminal_len
my $filtersize = 250000; # size of file need to remove reads with low depth by kmer check
my $bigfile = 1000000; # size of big file need to run in multi-threads
my $small_timeout = 300; # process will be killed if a small gene is assembled more than $small_timeout
my $big_timeout = 1800; # process will be killed if a big gene is assembled more than $big_timeout
my $gene_per_loop = 1000; # run 1000 gene per loop
my $process = 1; # number of cpu used

my $rundir = "temp"; # running dir

my $opt = GetOptions( 'parsed:s', \$dir,
					  'samplelist:s', \$samplelist,
                      'assembled:s', \$outdir,
                      'min_overlap:s', \$minolp,
                      'error_rate:s', \$errorate,
                      'cpu:s', \$process,
                      'help|h!', \$help) or usage(); 

# display help message if one of the following option is not specified
if (!($opt && $dir && $samplelist && $outdir) || $help) {
usage();
} 

# get the sample name in sample list
my %sample;
open SAMPLELIST, "$samplelist" or die "ERROR: Cannot find $samplelist ($!)";
while (my $samplename = <SAMPLELIST>) {
chomp $samplename;
$sample{$samplename} = "";
}
close SAMPLELIST;

# set name and create running dir and outdir of contig and asqg file
mkdir $outdir;
my $contigdir = "$outdir/contig";
my $asqgdir = "$outdir/asqg";
mkdir $contigdir;
mkdir $asqgdir;
mkdir $rundir;

function("assembling reads by sga");

# assemble reads sample by sample
opendir BANDP, $dir or die "ERROR: Cannot open $dir ($!)";
while (my $bandp = readdir BANDP) {
	next if ($bandp =~ /^\./);
	next if (!(exists $sample{$bandp}));

	print "Start assembling reads of $bandp\n";

	# set and create the dir for input, contig and asqg file, also the running dir ($temp)
	my $sample_indir = "$dir/$bandp";
	my $sample_outdir = "$contigdir/$bandp";
	my $sample_olpdir = "$asqgdir/$bandp";
	my $temp = "$rundir/$bandp";
	mkdir $sample_outdir;
	mkdir $sample_olpdir;
	mkdir $temp;
	
	# number of split array is equal to $process
	my $split = $process;

	opendir BANDPDIR, $sample_indir;
	# array of small reads file name, array of big reads file name
	my (@smallR1, @bigR1);
	my $read_cnt = 0; # reads file counter
	while (my $R1 = readdir BANDPDIR) {
		if ($R1 =~ /\S+_R1\.fq/) {
			# category reads as big or small
			my $file_size = (stat("$dir/$bandp/$R1"))[7];
			if ($file_size <= $bigfile) { # if flie size smaller than $bigfile, then categorize it as small
			push @smallR1, $R1;
			$read_cnt++;
			} else {  # elsif categorize it as big
			push @bigR1, $R1;
			}
		}	
		
		# if $gene_per_loop of gene has been read, then start assemble
		if ($read_cnt == $gene_per_loop) {
		# split gene into $split part
		my $split_array = split_array(\@smallR1, $split);
		
		# assemble on batch
		sga_batch($split_array, $sample_indir, $sample_outdir, $sample_olpdir, $temp);
		
		# remove temporary files
		`find ./$temp -mindepth 1 -delete`;
		
		# initialize counter and array
		$read_cnt = 0;
		@smallR1 = ();
		}	
	}
	closedir BANDPDIR;
	
	# do not forget the loci in last array
	if ($read_cnt > 0) {
	my $split_array = split_array(\@smallR1, $split);
	sga_batch($split_array, $sample_indir, $sample_outdir, $sample_olpdir, $temp);
	}
	`find ./$temp -mindepth 1 -delete`;
	
	# assemble big files in serial
	chdir($temp);
	while (my $R1 = shift @bigR1) {
		my ($gene) = $R1 =~ /(\S+)_R1\.fq/;
	
		my $killcnt = 0;
		while () {
		my $killed = sga_assemble($gene, $sample_indir, $sample_outdir, $sample_olpdir, $process, $big_timeout); # assemble
		last if (! $killed); # jump out loop if assembly is sucessful
		$killcnt++;
			# skip this file if fail to assemble more than 3 times
			if ($killcnt > 3) {
			print STDERR "$gene is failed to be assembled by sga for 3 times!\n";
			last;
			}
		}
	}
	chdir("../../");
	
	# mv asqg files from outdir to olpdir
	`mv $sample_outdir/*.asqg $sample_olpdir/`;
	
	print "Reads of $bandp have been assembled and written to $contigdir/$bandp/ \n";
	print "String graphs of $bandp have been written to $asqgdir/$bandp/ \n";
}                      
closedir BANDP;

print "All reads have been assembled and written to $contigdir/ \n";
print "All string graphs have been assembled and written to $asqgdir/ \n";

`find ./$rundir -delete`; # remove running dir

##################################################
# subroutine
##################################################

# subroutine to split array into $split
sub split_array {
my $smallR1 = shift; # array of small reads
my $split = shift; # split in to $split parts
	
	$split++ if ($split == 0); # at least split into 1 part
	my $block_num = int(@$smallR1/$split); # each reads set have $block_num of file
	$block_num++ if ($block_num == 0); # each reads set have at least 1 file
	
	my @small_sample;
	my $size_cnt = 0; # counter of number of reads
	my $block_cnt = 0; # counter of number of reads set
	
	while (my $file = shift @$smallR1) {
		if ($size_cnt >= $block_num) { # if there are enough reads in a reads set open a new set
		$size_cnt = 0;
		$block_cnt ++;
		}
	push @{$small_sample[$block_cnt]}, $file;
	$size_cnt ++;
	}

return(\@small_sample); # return splited reads set
}

# subroutine to run sga in batch
sub sga_batch {
my $split_array = shift; # split reads set
my $sample_indir = shift; # input dir of reads
my $sample_outdir = shift; # output dir of contig
my $sample_olpdir = shift; # output dir of graph file
my $temp = shift; # run dir

chdir($temp);

# initiate multi-process
my $pm = Parallel::ForkManager->new(int($process/2));

# assemble splited reads set in parallel
DATA_LOOP: while (my $sample_per_proc = shift @$split_array) {
$pm->start and next DATA_LOOP;
	# reads from the same set
	while (my $R1 = shift @$sample_per_proc) {
		my ($gene) = $R1 =~ /(\S+)_R1\.fq/;
	
		# assemble reads, kill process if running time exceed $small_timeout
		my $killcnt = 0;
		while () {
		my $killed = sga_assemble($gene, $sample_indir, $sample_outdir, $sample_olpdir, 1, $small_timeout); # assemble
		last if (! $killed); # jump out of loop if assembly is sucessful 
		$killcnt++;
			# skip $gene if it fail to be assembled more than 3 time
			if ($killcnt >= 3) {
			print "$gene is failed to be assembled by sga for 3 times!\n";
			last;
			}
		}
	}
$pm->finish;	
}
$pm->wait_all_children;

chdir("../../");
}

# subroutine to assemble reads by sga
sub sga_assemble{
my $gene = shift; # loci name
my $sample_indir = shift; # input dir of reads
my $sample_outdir = shift; # output dir of reads
my $sample_olpdir = shift; # input dir of graph file
my $thread = shift; # number of thread
my $timeout = shift; # running time limit

# reads size
my $size = (stat("$sample_indir/${gene}_R1.fq"))[7];

# command of preprocess, construct FM index, reads correct, iconstruct FM index of corrected reads, filter, find overlap ans assemble
my ($preprocess, $index, $correct, $reindex, $filter, $overlap, $assemble);
if ($size < $filtersize) { # if reads size is smaller than $filtersize, do not use kmer check
$preprocess = "sga preprocess -o $gene.pro.fq -p 1 $sample_indir/${gene}_R1.fq $sample_indir/${gene}_R2.fq &> /dev/null";
$index = "sga index -t $thread -a ropebwt $gene.pro.fq &> /dev/null";
$correct = "sga correct -t $thread -O 2 -x 1 -o $gene.co.fa $gene.pro.fq &> /dev/null";
$reindex = "sga index -t $thread -a ropebwt $gene.co.fa &> /dev/null";
$filter = "sga filter -t $thread --homopolymer-check --low-complexity-check --no-kmer-check -o $gene.ft.fa $gene.co.fa &> /dev/null";
$overlap = "sga overlap -t $thread -m $minolp -e $errorate $gene.ft.fa &> /dev/null";
$assemble = "sga assemble -m $minolp -b 0 -x 1 -l $terminal_len -o $sample_outdir/$gene $gene.ft.asqg.gz &> /dev/null";
} else { # if reads size is exceed $filtersize, use kmer check to remove reads in low depth
$preprocess = "sga preprocess -o $gene.pro.fq -p 1 $sample_indir/${gene}_R1.fq $sample_indir/${gene}_R2.fq &> /dev/null";
$index = "sga index -t $thread -a ropebwt $gene.pro.fq &> /dev/null";
$correct = "sga correct -t $thread -O 2 -x 1 -o $gene.co.fa $gene.pro.fq &> /dev/null";
$reindex = "sga index -t $thread -a ropebwt $gene.co.fa &> /dev/null";
$filter = "sga filter -t $thread --homopolymer-check --low-complexity-check -o $gene.ft.fa $gene.co.fa &> /dev/null";
$overlap = "sga overlap -t $thread -m $minolp -e $errorate $gene.ft.fa &> /dev/null";
$assemble = "sga assemble -m $minolp -b 0 -x 1 -l $terminal_len -o $sample_outdir/$gene $gene.ft.asqg.gz &> /dev/null";
}

# assemble reads, and catch the ALRM signal if run time is too long
my $pid;
eval {
local $SIG{ALRM} = sub {die "time out"};
alarm $timeout;
$pid = open(ASSEMBLEPIPE, "$preprocess |") or die "Cannot fork when preprocess $gene ($!)";while (<ASSEMBLEPIPE>){};close ASSEMBLEPIPE;
$pid = open(ASSEMBLEPIPE, "$index |") or die "Cannot fork when index $gene ($!)";while (<ASSEMBLEPIPE>){};close ASSEMBLEPIPE;
$pid = open(ASSEMBLEPIPE, "$correct |") or die "Cannot fork when correct $gene ($!)";while (<ASSEMBLEPIPE>){};close ASSEMBLEPIPE;
$pid = open(ASSEMBLEPIPE, "$reindex |") or die "Cannot fork when $reindex $gene ($!)";while (<ASSEMBLEPIPE>){};close ASSEMBLEPIPE;
$pid = open(ASSEMBLEPIPE, "$filter |") or die "Cannot fork when filter $gene ($!)";while (<ASSEMBLEPIPE>){};close ASSEMBLEPIPE;
$pid = open(ASSEMBLEPIPE, "$overlap |") or die "Cannot fork when overlap $gene ($!)";while (<ASSEMBLEPIPE>){};close ASSEMBLEPIPE;
$pid = open(ASSEMBLEPIPE, "$assemble |") or die "Cannot fork when assemble $gene ($!)";while (<ASSEMBLEPIPE>){};close ASSEMBLEPIPE;
alarm (0);
};

if (! $@) { # if "die" is not caught, decompress the asqg file only
`gunzip -f $sample_outdir/${gene}-graph.asqg.gz`;	
} else { # if "die" is caught
	if ($@ =~ /time\sout/) { # kill process if time out signal is caught
		if (kill 0, $pid) { # if $pid exist
		kill 'TERM', $pid; # kill safely by TERM signal
		waitpid($pid, 0); # end current process 
			if (kill 0, $pid) { # if $pid not respond to TERM signal
			print "sga run of $sample_indir/$gene did not respond in ${timeout}s, and was abnormally killed\n";
			kill 'KILL', $pid; # kill $pid by KILL signal
			waitpid($pid, 0); # end current process 
			} else {
			print "sga run of $sample_indir/$gene did not respond in ${timeout}s, and was safely killed\n";
			}
		}
		close ASSEMBLEPIPE;
		return $gene;
	} else { # only close pipe if "Cannot not fork" is caught
	close ASSEMBLEPIPE;
	}
}
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

# subroutine to print usage
sub usage {
print "
Script name: sga_assemble.pl

This is a script to assemble parsed reads by SGA

Dependencies:
(1) SGA
(2) Perl Module:
	Parallel::ForkManager

Example usage:
	
	perl sga_assemble.pl --parsed bandpout --assembled sgaout --samplelist samplelist.txt

Input files:
(1) bandpout
(2) samplelist.txt

Output files:
(1) sgaout

Options:
--parsed
  Directory containing parsed reads 
--samplelist
  A list of sample name. Only samples in the list will be assembled
--assembled
  Directory containing assembled reads and string graphs
--min_overlap
  Minimum overlap required between 2 reads, $minolp in default
--error_rate
  Maximum error rate allowed to consider two sequences aligned, $errorate in default
--cpu
  Limit the number of CPUs, $process in default.
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