#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;

# ambiguous nucleotide hash
my %IUB = (
        A => [qw(A)],
        C => [qw(C)],
        G => [qw(G)],
        T => [qw(T)],
        U => [qw(U)],
        M => [qw(A C)],
        R => [qw(A G)],
        S => [qw(C G)],
        W => [qw(A T)],
        Y => [qw(C T)],
        K => [qw(G T)],
        V => [qw(A C G)],
        H => [qw(A C T)],
        D => [qw(A G T)],
        B => [qw(C G T)],
        N => [qw(A C G T)],
        X => [qw(A C G T)],
);

# input dir of trimmed reads, output dir of deduplicated reads, samplelist
my ($dir, $outdir, $samplelist, $help);
my $splitsize = 200; # files are splited into 200 MB before removing duplication
my $outputsize = 400; # size of output is around 400 MB
my $readcountfile = "rmdup_reads_bases_count.txt"; # file summarized number of original and deduplicated reads
my $merged_reads = "merged_reads"; # dir containing concatenated reads
my $process = 1; # number of cpu used

my $opt = GetOptions( 'trimmed:s', \$dir,
					  'samplelist:s', \$samplelist,
                      'rmdup:s', \$outdir,
                      'readcount:s', \$readcountfile,
                      'cpu:i', \$process,
                      'help|h!', \$help) or usage(); 

# display help message if one of the following option is not specified
if (!($opt && $dir && $samplelist && $outdir) || $help) {
usage();
} 
  
# open dir and get trimmed reads of sample which only exists in sample list                    
opendir DIR, $dir or die "ERROR: Cannot open $dir ($!)";
my @sample;
open SAMPLELIST, "$samplelist" or die "ERROR: Cannot find $samplelist ($!)";
while (my $samplename = <SAMPLELIST>) {
chomp $samplename;
push @sample, $samplename;
}
close SAMPLELIST;
closedir DIR;

# turn unit of size into Byte
$splitsize = getsize($splitsize);
$outputsize = getsize($outputsize);

my %readcount;
merge_reads(); # concatenate R1 and R2 reads together
dedup(); # remove concatenated reads
output(); # write output

# write file recoding original and deduplicated reads count
my $READCOUNT;
open $READCOUNT, ">", $readcountfile or die ("ERROR: Cannot write $readcountfile ($!)");
print $READCOUNT "Sample\tReads num. of cleaned reads\tBases num. of cleaned reads (bp)\tReads num. of deduplicated reads\tBases num. of deduplicated reads (bp)\n";
map {print $READCOUNT "$_\t$readcount{$_}\n"} sort keys %readcount;
close $READCOUNT;

# mv dir containing concatenated reads
`rm -rf $merged_reads`;

##################################################
# subroutine
##################################################

# subroutine to concatenated reads in R1 and R2
sub merge_reads {
function("concatenating reads");

mkdir "$merged_reads";

my $pm = Parallel::ForkManager->new(int($process));

# callback of hash containing number of reads
$pm->run_on_finish(
	sub {
	my $info = $_[-1];
	$readcount{$info->[0]} = "$info->[1]\t$info->[2]";
	}
);

# concatenate R1 and R2 for each sample
DATA_LOOP: foreach my $sample (@sample) {
$pm->start and next DATA_LOOP;
print "Start concatenating reads of $sample\n";

# get name of trimmed reads and concatenated reads
my $R1 = "$dir/$sample" . "_R1.fq";
my $R2 = "$dir/$sample" . "_R2.fq";
my $merged = "./$merged_reads/$sample.merge.fq";

# get the size of trimmed reads
my $filesizer1 = sprintf("%.2f", (stat($R1))[7]);
my $filesizer2 = sprintf("%.2f", (stat($R2))[7]);
my $frsize = $filesizer1+$filesizer2;
my $part = part($frsize); # determine splited file into $part part
	
	# concatenate R1 and R2
	my ($readcount, $basecount, $MERGE);
	open PAIR1, $R1 or die ("ERROR: Cannot open $R1 ($!)");
	open PAIR2, $R2 or die ("ERROR: Cannot open $R2 ($!)");
	open $MERGE, ">", $merged;
	while (my $line1 = <PAIR1>) {
		if ($line1 =~ /(^\@\S+)\s\S+/) { # if we get a read header which is start with @
		
		# read id, seq and quality sequence of R1
		my $id = $1;
		
		chomp(my $seq1 = <PAIR1>);
		$seq1 = ambiguous_sub($seq1); # remove ambiguous bases
		$readcount++;
		$basecount += length($seq1);
		
		<PAIR1>;
		chomp(my $qua1 = <PAIR1>);
		die ("ERROR: The length of reads and quality of $id are not the same in ${sample}_R1.fq!\n") if (length($seq1) != length($qua1));

		# read id, seq and quality sequence of R2
		<PAIR2>;
		
		chomp(my $seq2 = <PAIR2>);
		$seq2 = ambiguous_sub($seq2);
		$readcount++;
		$basecount += length($seq1);
		
		<PAIR2>;
		chomp(my $qua2 = <PAIR2>);
		die ("ERROR: The length of reads and quality of $id are not the same in ${sample}_R2.fq!\n") if (length($seq2) != length($qua2));
		
			my $seq2st = length($seq1);	# record the length of seq1 as the start point of seq2 in concatenated reads
			if ($part > 1) {
			print $MERGE "$seq1$seq2 $id\/$seq2st $qua1$qua2\n"; # concatenate reads in "concatenated_sequence id/start_point_of_seq2 concatenated_quality_sequence" format if reads need to be splited
			} else {
			print $MERGE "$id\/$seq2st\n$seq1$seq2\n+\n$qua1$qua2\n"; # concatenated reads normally in fastq format if reads do not need to be splited
			}
		}
	}
	close PAIR1;
	close PAIR2;
	close $MERGE;
	
	my $mergesize = sprintf("%.2f", (stat($merged))[7]); # get the size of concatenated reads
	my $maxline = int($splitsize/($mergesize/($readcount/2))); # get number of line in each splited file
	
	if ($part > 1) { # if concatenated reads need to be splited
		
		# first time split
		my $part = part($mergesize); # determine split reads in to $part part
		partition($merged, $sample, $part, $maxline); # split it
		my @partitioned = glob("./$merged_reads/$sample.*.fq"); # get the names of splited file
		my @bigfile = grep {(((stat($_))[7])/$splitsize)>1.5} @partitioned; # file whether some file need to be splited again(> 1.5 times of $splitsize)
		
		# split files until all of them is smaller than 1.5 times of $splitsize
		while () {
			foreach my $bigfile (@bigfile) {
			my $filesize = sprintf("%.2f", (stat($bigfile))[7]);
			my $part = part($filesize);
			next if ($part == 1);
			partition($bigfile, $sample, 4, $maxline);
			}
		my @partitioned = glob("./$merged_reads/$sample.*.fq");
		@bigfile = grep {(((stat($_))[7])/$splitsize)>1.5} @partitioned;
		last if (@bigfile == 0);
		}
		
		# reformat splited file as fastq style
		my $REFORMATOUT;
		my @reformat = glob("./$merged_reads/$sample.*.fq");
		foreach my $reformat (@reformat) {
		open REFORMAT, $reformat;
		open $REFORMATOUT, ">", "$reformat.reformat";
			while (<REFORMAT>) {
			my @reads = $_ =~ /(\S+)/g;
			print $REFORMATOUT "$reads[1]\n$reads[0]\n+\n$reads[2]\n";
			}
		close $REFORMATOUT; 
		close REFORMAT;
		unlink $reformat;
		`mv $reformat.reformat $reformat`;
		}
	
	print "Reads of $sample has been concatenated, sorted and divided into " . @reformat . " part\n";
	} else {
	print "Reads of $sample has been concatenated\n";
	} 
		
$pm->finish(0, [$sample, $readcount, $basecount]); # return number of reads
}
$pm->wait_all_children();

print "All reads have been concatenated\n";
}

# subroutine to remove PCR duplicates
sub dedup {
function("removing PCR duplicates");

# get name of all splited files under merged_reads
opendir MERGED_READS, "$merged_reads";
my @fq = grep {$_ !~ /^\./} readdir MERGED_READS;
closedir MERGED_READS;

# remove duplicated reads for each fq by fastx_uniques under usearch
my $pm = Parallel::ForkManager->new(int($process/2));
chdir "$merged_reads";

foreach my $sample (@sample) {
	print "Start removing PCR duplicates in $sample\n";
	
	my @partitioned_fq = glob("./$sample.*.fq");

	DATA_LOOP: foreach my $fq (@partitioned_fq) {
	$pm->start and next DATA_LOOP;
	
		uniq_fastx($fq); 
	
		# remove original splited reads
		my ($samplename) = $fq =~ /(\S+\.\S+)\.fq/;
		my @sample = glob("./$samplename.*");
		foreach (@sample) {
		unlink $_ if ($_ !~ /\.rmdup\.fq$/);
		}

	$pm->finish();
	}
	$pm->wait_all_children();

	print "$sample has been deduplicated\n";
}

chdir "..";

print "All reads have been deduplicated\n";
}

# subroutine to write output
sub output {
function("writing out deduplicated reads");    

# make output dir
mkdir $outdir;

# make output dir of fq
my $fqout = "$outdir";
mkdir $fqout;

# make temporary output dir for all splited fq
mkdir "fqtemp";

my $pm = Parallel::ForkManager->new(int($process));

# callback to get deduplicated reads counts
$pm->run_on_finish(
	sub {
	my $info = $_[-1];
	$readcount{$info->[0]} .= "\t$info->[1]\t$info->[2]";
	}
);

# write output foreach sample
DATA_LOOP: foreach my $sample (@sample) {
$pm->start and next DATA_LOOP;
print "Start writing out deduplicated reads of $sample\n";

my ($readcount, $basecount);

# get all splited reads belonging to $sample
my @samplepartition = glob("./$merged_reads/$sample.*.rmdup.fq");

# get name of splited reads
map {$_ =~ s/\.\/$merged_reads\///} @samplepartition;
map {$_ =~ s/\.fq$//} @samplepartition;

	# write output of each splited reads to temporary dir
	my $RECOVERFQ;
	foreach my $samplepartition (@samplepartition) {
	my $rmdupfq = "$samplepartition.fq";
	open MERGE, "./$merged_reads/$rmdupfq";
	open $RECOVERFQ, ">", "fqtemp/$rmdupfq"; 
		while (my $mline = <MERGE>) {
			if ($mline =~ /(\@\S+)\/(\d+)/) {
			my $mid = $1;
			my $mst = $2;
			
			chomp(my $mseq = <MERGE>);
			my $seq_r1 = substr $mseq, 0, $mst;
			my $seq_r2 = substr $mseq, $mst;
			
			$readcount += 2;
			$basecount += length($mseq); 
			
			<MERGE>;
			chomp(my $mqua = <MERGE>);
			my $qua_r1 = substr $mqua, 0, $mst;
			my $qua_r2 = substr $mqua, $mst;
			print $RECOVERFQ "$mid\/1\n$seq_r1\n\+\n$qua_r1\n";
			print $RECOVERFQ "$mid\/2\n$seq_r2\n\+\n$qua_r2\n";
			}
		}
	close $RECOVERFQ;
	close MERGE;
	}
	
	# rearrange the splited reads and merge them if it is too small
	my @category; # set array which put splited reads
	my $categorynum = 0; # set order of array element
	my $sizebox = 0; # set size box which can only contain reads smaller than $outputsize
	my @sortsamplepartition = sort {(stat("fqtemp/$a.fq"))[7] <=> (stat("fqtemp/$b.fq"))[7]} @samplepartition;
	
	foreach my $samplepartition (@sortsamplepartition) {
	my $fq = "fqtemp/$samplepartition.fq";
	my $size = (stat($fq))[7];
		$sizebox += $size;
		if ($sizebox <= $outputsize) { # if size of reads in $sizebox haven't reach $outputsize, push its name into array
		push @{$category[$categorynum]}, $samplepartition;
		} else { # if size of reads in $sizebox have reached $outputsize, push it into next element 
			my $size = (stat($fq))[7];
			$sizebox += $size;
			if ($sizebox > $outputsize) {
			$categorynum++ if (@{$category[$categorynum]} > 0); # update order of element if this element is not empty
			$sizebox = $size;
			}
			push @{$category[$categorynum]}, $samplepartition; 
		}	
	}
	
	# cat output if number of splited sample is more than 1, else just mv it into output dir
	for (my $i=0; $i<(@category); $i++) {
	my @array = @{$category[$i]};
		if (@array > 1) {
		my @joinfq = map {"./fqtemp/" . $_ . ".fq"} @array;
		my $joinfq = join " ", @joinfq;
		`cat $joinfq > $fqout/$sample.$i.rmdup.fq`;
		} elsif (@array == 1) {
		`mv ./fqtemp/$array[0].fq $fqout/$sample.$i.rmdup.fq`;
		}
	}

print "Deduplicated reads of $sample has been written to $fqout/$sample.*.rmdup.fq\n";
$pm->finish(0, [$sample, $readcount, $basecount]);
}
$pm->wait_all_children();

print "All Deduplicated reads has been written to $fqout/ \n";

`rm -rf fqtemp`;
}

# subroutine to determine how many parts split samples into 
sub part {
my $size = shift; # reads size

my $ratio = int($size/$splitsize); # get int part of reads should be splited
my $rest = ($size-($ratio*$splitsize))/$splitsize; # get ratio of rest size of reads to split size
# determine part of reads should be splited
my $part = do {
	if ($rest > 0.5) { # divide reads into one more split if rest size is too big
	$ratio++;
	$ratio;
	} else { # divide reads into one more split if rest size is too big
	$ratio;
	}
};
$part++ if ($part == 0); # at least have 1 part

return $part;
}

# subroutine to substitute ambiguous bases
sub ambiguous_sub {
my $seq = shift; # DNA sequence

# substitute ambiguous reads to a nucleotide in random
my @ambiguous = $seq =~ /([^A|T|C|G])/g;
	if (@ambiguous != 0){
		foreach my $base (@ambiguous) {
		my @bases = @{$IUB{$base}};
		my $sub = $bases[int(rand(scalar(@bases)))];
		$seq =~ s/$base/$sub/g;
		}
	}
return $seq;
}

# subroutine to sort big files
sub bigsort {
my $filepath = shift;
my $maxline = shift;
	
`split -l $maxline $filepath ${filepath}_`; # split reads into several parts with maximum of $maxline line
my @split = glob("./${filepath}_*");

# sort each splited parts
foreach my $split (@split) {
`sort -S 512M -o $split.sort $split`; 
}

# merge sorted parts
`cat ${filepath}_*.sort > $filepath.sorted`;

# remove splited parts
`rm -f ${filepath}_*`;
unlink $filepath;
}

# subroutine to split big files into several splits according to their prefix, longer prefix means more splits
sub partition {
my $filepath = shift; # reads file
my $sample = shift; # sample name
my $part = shift; # parts to be splited
my $maxline = shift; # maximum line allowed in a splited file for sorting

# determine the length of prefix by $part
$filepath =~ /\.\/$merged_reads\/\S+\.(\S+)\.fq/;
my $prefixlen = do {
	if ($1 eq "merge") { # if files haven't been splited yet, we need to define the prefix length of it
	my $len = int(log($part)/log(4)); # change number of split parts into prefix length
	$len++ if ($len == 0);
	$len;
	} else { # add prefix length to split reads into more parts
	my $len = length($1)+1;
	$len;
	}
};

# sort file
bigsort($filepath, $maxline);

# write reads same prefix to a file
my $prefix = "none";
my $PREFIX;
open SORT, "$filepath.sorted";
while (my $line = <SORT>) {
my $preline = substr $line, 0, $prefixlen; # get prefix
	if ($preline ne $prefix) { # if prefix of current line is not the same as reads file that is being written into
		close $PREFIX if ($prefix ne "none"); # close filehandle of reads file that is being written into
		$prefix = $preline; # update the prefix
		open $PREFIX, ">>", "$merged_reads/$sample.$preline.fq"; # open filehandle of file with the updated prefix
	}
	print $PREFIX $line; # write reads to the filehandle opening reads file
}
close SORT;
close $PREFIX;
unlink "$filepath.sorted"; # remove sorted reads
}

# subroutine to turn input size into unit of Byte
sub getsize {
my $size = shift; # size

my ($num, $unit) = $size =~ /(\d+\.*\d*)(\w+)*/; # get size and unit
my $sizeb = do {
	if ($unit) { # if size unit is found
		do {
			if ($unit =~ /K|KB/i) { # KB
			$num * 2**10;
			} elsif ($unit =~ /M|MB/i) { # MB
			$num * 2**20;
			} elsif ($unit =~ /G|GB/i) { # GB
			$num * 2**30;
			}
		};
	} else { # if size unit is not found, assume is MB
	$num * 2**20;
	}
};

# return size in byte
return $sizeb;
}

# subroutine to use fastx_uniques under usearch to remove PCR duplicates
sub uniq_fastx {
my $merge = shift; # concat reads

my $rmdup = $merge;
$rmdup =~ s/\.fq$/\.rmdup\.fq/; # substitute suffix of reads
`usearch -fastx_uniques $merge -fastqout $rmdup -strand both &> /dev/null`; # use fastx_uniques to remove duplicates
}

# subroutine to report the start of a step of the program
sub function {
my $function = shift; # name of step

if ($process > 1){ # if there are more than one process used
print "Start $function in parallel\n";
} else { # if there is only one process used
print "Start $function in serial\n";
}
}

sub usage {
print STDERR "
Script name: rmdup.pl

This is a script to remove PCR duplicates by fastx_uniques under usearch

Dependencies:
(1) USEARCH v10.0.240 or higher
(2) Perl Module:
	Parallel::ForkManager

Example usage:

	perl rmdup.pl --trimmed trimmed --rmdup rmdup --samplelist samplelist.txt

Input files:
(1) trimmed
(2) samplelist.txt

Output files:
(1) rmdup

Options:
--trimmed
  Directory containing reads without adaptor and low quality reads
--samplelist
  A list of sample name. Only samples in the list will be assembled
--rmdup
  Directory containing reads without PCR duplication
--readcount
  Tab delimited file summarizing number of reads before and after removing PCR duplicates, named in '$readcountfile' in default
--cpu
  Limit the number of CPUs, 1 in default.
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