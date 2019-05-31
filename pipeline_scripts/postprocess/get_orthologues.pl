#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;

# sequences of reference, db of reference species, db of subject, list of subject name, mask db of reference or not, mask db of subject or not, outdir, delete udb of subject
my ($query, $querydb, $subdb_list, $subname_list, $mask_querydb, $mask_subdb, $outdir, $delete_subdb, $help);
my $db_size = "1G"; # if input nucleotide size of db is larger than $db_size, it will be split
my $db_step = 2; # Every Nth ublast database word should be indexed
my $wordsize = 8; # seed word size 
my $masker = "dust"; # nucleotide masker
my $E = 1e-6; # e value
my $maxaccepts = 0; # terminate ublast search after find $maxaccepts hit, accept all if $maxaccepts == 0 and $maxrejects == 0
my $maxrejects = 0; # terminate ublast search after reject $maxrejects input reads, accept all if $maxaccepts == 0 and $maxrejects == 0

my $second_id_dis_trsd = 40; # difference of pdistance between 1st hit and 2nd hit
my $query_cov_trsd = 100; # minimum percentage of query sequences aligned to database required to deem it as exist in database
my $query_id_trsd = 98; # minimum similarity of between query sequences and database required to deem it as exist in database

my $paralog_id_dis_trsd = 40;

my $st_trsd = 0; # number of stop codon allowed in output sequence
my $min_seq_len = 100; # minimum seq length required

my $cpu = 1; # number of cpu used

my $qblasts = "qblasts"; # dir for blast out of query against provided database
my $qblastsout = "qblastsout"; # dir for sequence in provided database which hit with query
my $reblast_input = "reblast_input.fas"; # merged file of query + sequence in provided database
my $reblast = "reblast"; # dir for sequences in provided database which got reciprocal blast hit

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'query:s', \$query,
                      'querydb:s', \$querydb,
                      'subdb:s', \$subdb_list,
                      'subname:s', \$subname_list,
                      'mask_querydb!', \$mask_querydb,
                      'mask_subdb!', \$mask_subdb,
                      'outdir:s', \$outdir,
                      'E:s', \$E,
                      'min_len:i', \$min_seq_len,
					  'stop_codon:s', \$st_trsd,
					  'delete_subdb!', \$delete_subdb,
					  'cpu:s', \$cpu,
					  'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/query querydb subdb subname outdir/;
check_option(\@essential, \@ARGVoptions);

# check usearch
`usearch` or die "\nERROR: Cannot find usearch in \$PATH. It may just not named as 'usearch' \n\n";

# check whether output dir name is identical with intermediate dir
if (($outdir eq $qblasts)||($outdir eq $qblastsout)||($outdir eq $reblast_input)||($outdir eq $reblast)) {
die "\nERROR: Do not name --outdir as $qblasts, $qblastsout, $reblast_input or $reblast\n\n";
}

# turn the size into unit of byte
$db_size = getsize($db_size);

# outdir
my $nf = "$outdir/nf"; # coding
my $p = "$outdir/p"; # aa

# check whether reference is dna sequences, and save each sequences into hash
my %qseq;
open QUERY, $query or die "\nERROR: Cannot find file of target loci \"$query\", please check --query ($!)\n\n";
my $first_line = <QUERY>;
chomp(my $first_seq = <QUERY>);
my ($header) = $first_line =~ />(\S+)/;
$qseq{$header} = $first_seq;

my @nucleo = $first_seq =~ /[A|T|C|G]/g;
if (@nucleo/length($first_seq) < 0.5) {
die "\nERROR: --query must contain nucleotide sequences only\n\n";
}

while (my $line = <QUERY>) {
	if ($line =~ />(\S+)/) {
	chomp(my $seq = <QUERY>);
	$qseq{$1} = $seq;
	}
}
close QUERY;

# name of all input reference db
my @querydb = $querydb =~ /(\S+)/g;

# check the existence and format of input reference db
my $fascnt = 0; # fasta count
my $udbcnt = 0; # udb count
my $elsecnt = 0; # else file count
foreach my $querydb (@querydb) {
	if (!(-e $querydb)) {
	die "\nERROR: Cannot find database of species of target loci \"$querydb\", please check --querydb ($!)\n\n";
	}
	if ($querydb =~ /\S+\.fa$|\S+\.fas$|\S+\.fasta$/) {
	$fascnt++;
	} elsif ($querydb =~ /\S+\.udb$/) {
	$udbcnt++;
	} else {
	$elsecnt++;
	}
}

# die if found existence of unrecognized format
if ($elsecnt > 0) { 
die "\nERROR: Databases of reference genome must be in either fasta or udb format\n\n";
}
# die if found mixed format
if (($fascnt > 0)&&($udbcnt > 0)) { 
die "\nERROR: Databases of reference genome must be in either fasta or udb format. Mixed format is not supported\n\n";
}
# must provide fast format file, if db need to be masked
if ($mask_querydb) {
	if ($udbcnt > 0) {
	die "\nERROR: Please give reference genomes in fasta format if it need to be masked\n\n";
	}
}

# get all name of subject db
my @subdb_all; # array of all subject db

my @subdb_list = $subdb_list =~ /([^\|]+)/g; # get provided db(s) of each subject 
# check db of provided subject
foreach my $subdb_class (@subdb_list) {
	my @subdb = $subdb_class =~ /\S+/g; # get db(s) of one subject 
	foreach my $subdb (@subdb) {
		if (!(-e $subdb)) { # check existence
		die "\nERROR: Cannot find database of subjects genome \"$subdb\", please check --subdb ($!)\n\n";
		}
		
		# check input format
		my $subname;
		my ($fascnt, $udbcnt, $elsecnt) = (0, 0, 0);
		if ($subdb =~ /\S+\.fa$|\S+\.fas$|\S+\.fasta$/) {
		$fascnt++;
		} elsif ($subdb =~ /(\S+)\.udb$/) {
		$udbcnt++;
		} else { # else format
		$elsecnt++;
		}
		
		# die if found existence of unrecognized format
		if ($elsecnt > 0) {
		die "\nERROR: Databases of subject genome must be either in fasta or udb format\n\n";
		}
		# die if found mixed format
		if (($fascnt > 0)&&($udbcnt > 0)) { 
		die "\nERROR: Databases of subject genomes in one class of subject must be either in fasta or udb format. Mixed format is not supported\n\n";
		}
	}
	push @subdb_all, \@subdb;
}

# if subject name list is provided
my %subdb_all; # hash of input db of all subjects
if ($subname_list) {
	my @subname_list = $subname_list =~ /(\S+)/g;
	if (@subname_list != @subdb_all) { # die if name number is not match with db class number
	die "\nERROR: Number of subject databases and name in name list is different\n\n";
	} else { # parse name to corresponds db
		foreach (0..(@subdb_all-1)) {
		$subdb_all{$subname_list[$_]} = $subdb_all[$_];
		}
	}
} else {
die "\nERROR: Name list of subjects in extracted orthologs must be provided, please check --subname\n\n";
}

#################################################
function("construct database");
# prepare query udb
my @delete_files; # array of files to be deleted

my ($delete_splited_query, $queryudb_all); # array of splited reference db to be deleted, array of reference db
if ($fascnt > 0) { # construct udb if fasta is provided
	($queryudb_all, $delete_splited_query) = split_db(\@querydb); # split fasta file of reference db if too large
	push @delete_files, @$delete_splited_query; # array of splited fasta of reference db to be removed later
	
	if ($mask_querydb) { # mask
	$queryudb_all = mask_db($queryudb_all);
	push @delete_files, @$queryudb_all; # array of masked fasta of reference db to be removed later
	}
	
	$queryudb_all = make_udb($queryudb_all); # make db
} else {
$queryudb_all = \@querydb; 
}

# prepare subject udb
my ($subudb_all, @delete_subdb); # hash of all subject db, array of all sub udb to be deleted
foreach my $subdb_class (sort keys %subdb_all) {
	my ($subudb_class, $delete_splited_subdb); # array of subject udb of $subdb_class, array of splited subject fasta to be deleted later in $subdb_class
		
	if ($subdb_all{$subdb_class}->[0] =~ /\S+\.fa$|\S+\.fas$|\S+\.fasta$/) { # make udb if db is in fasta
		($subudb_class, $delete_splited_subdb) = split_db($subdb_all{$subdb_class}); # split big db
		push @delete_files, @$delete_splited_subdb; # array of splited fasta to be removed later in $subdb_class
	
		if ($mask_subdb) { # mask
		$subudb_class = mask_db($subudb_class);
		push @delete_files, @$subudb_class; # array of masked fasta to be removed later in $subdb_class
		}
	
		$subudb_class = make_udb($subudb_class); # make db
		push @delete_subdb, @$subudb_class; # array of subject udb to be removed later in $subdb_class
	} else {
	$subudb_class = $subdb_all{$subdb_class};
	}

	$subudb_all->{$subdb_class} = $subudb_class;
}

print "Databases have been constructed\n";

#################################################
function("blast reference against subject databases");

# mask reference sequences
my $masked = "$query.masked";
maskoutput($query, $masked, 1);
$query = $masked;

# reference sequences blast against subject db
my $index = qblasts($query, $subudb_all);

# get sequence hit with reference sequences
qblastsout($subudb_all, $index, \%qseq);
print "Blast(reference to subjects) search has been finished\n";
#################################################
function("blast subject against reference databases");
# mask merged reference+hit sequences files
my $reblast_input_masked = "$reblast_input.masked";
maskoutput($reblast_input, $reblast_input_masked);
$reblast_input = $reblast_input_masked;

# blast against reference db
sblastq($queryudb_all);
print "Blast(subjects to reference) search has been finished\n";
################################################################
function("find orthologues");

# find sequences with reciprocal best hit
ortholog();
print "Orthologues have been written to $outdir\n";
# ################################################################
# remove intermediate dir
`rm -rf $qblasts`;
`rm -rf $qblastsout`;
`rm -rf $reblast`;

unlink $reblast_input; # delete merged reference+hit sequence file
unlink $query; # delete masked reference

# delete splited fas and udb of reference and subject db
foreach $_ (@delete_files) {
unlink $_;
}

# delete subject udb if --delete_subdb is specified
if ($delete_subdb) {
	foreach (@delete_subdb) {
	unlink $_;
	}
}

################################################################
# Subroutines
################################################################

# subroutine to translate size into byte unit
sub getsize {
my $size = shift; # file size

my ($num, $unit) = $size =~ /(\d+\.*\d*)(\w+)*/; # get number and unit
my $sizeb = do { # convert into byte
	if ($unit) { # convert according to unit
		do {
			if ($unit =~ /K|KB/i) { # kb
			$num * 2**10;
			} elsif ($unit =~ /M|MB/i) { # mb
			$num * 2**20;
			} elsif ($unit =~ /G|GB/i) { # gb
			$num * 2**30;
			}
		};
	} else { # if no unit found, the unit is MB
	$num * 2**20;
	}
};

# return size in byte
return $sizeb;
}

# subroutine to split database 
sub split_db {
my $dblist = shift; # list of fasta genome

my (@splited_file, @delete_splited); # array of splited files, array of files to be deleted
foreach my $fas (@$dblist) {
my $file_size = (stat($fas))[7]; # get size

	if ($file_size > $db_size) { # if it is too large
	my $split_num = int($file_size/$db_size)+1; # split into $split_num parts
	my $splited = split_fas($fas, $file_size, $split_num); # split
	push @splited_file, @$splited;
	push @delete_splited, @$splited;
	} else { # if it is not large
	push @splited_file, $fas;
	}
}

# return array of splited files and files to be deleted
return (\@splited_file, \@delete_splited);
}

# subroutine to split fasta file
sub split_fas {
my $fas = shift; # fasta file
my $file_size = shift; # file size
my $split_num = shift; # part to split

my ($fas_name) = $fas =~ /([^\/]+$)/; # extract name of file from its path

my $seq_chunk = int($file_size/$split_num); # size of each split

# split file
my @splited_fas;
my $word_cnt = 0; # size counter
my $splited_num = 0; # split counter
open SPLITED, ">$fas_name.$splited_num"; # write first split
push @splited_fas, "$fas_name.$splited_num"; 
open FAS, $fas;
while (my $line = <FAS>) {
	if ($line =~ /^>.+/) {
		if ($word_cnt > $seq_chunk) { # write a new split if its size exceed $seq_chunk
		close SPLITED;
		$word_cnt = 0;
		$splited_num++;
		open SPLITED, ">$fas_name.$splited_num";
		push @splited_fas, "$fas_name.$splited_num";
		}
	} 
	print SPLITED $line; # write split
	$word_cnt += length($line);
}
close FAS;
close SPLITED;

if ($splited_num+1 < $split_num) { # if number of splited file cannot reach expected value
print STDERR "WARNING: Cannot split $fas into $split_num parts, because there's not enough sequences in $fas\n";
}

# return name of splited fasta
return(\@splited_fas);
}

# subroutine to mask database
sub mask_db {
my $dbarray = shift; # array of db

# mask db
my @masked_db;
foreach my $fas (@$dbarray) {
my ($fas_name) = $fas =~ /([^\/]+$)/;
my $masked = "$fas_name.masked";
maskoutput($fas, $masked, 1);
push @masked_db, $masked;
}

# return masked db
return \@masked_db;
}

# subroutine to mask input fasta file
sub maskoutput {
my ($unmask, $masked, $unlink) = @_; # unmasked file, masked file, remove original file or not
`usearch -fastx_mask $unmask -fastaout $masked -qmask $masker`; # mask
if (!($unlink)) { # delete original file if $unlink is specified
unlink $unmask;
}
}

# subroutine to construct udb 
sub make_udb {
my $dbarray = shift; # array of db

# make db
my @udb;
foreach my $fas (@$dbarray) {
my ($db_name) = $fas =~ /([^\/]+$)/;
my $udb = "$db_name.udb";
`usearch -makeudb_ublast $fas -output $udb -wordlength $wordsize -dbmask user -dbstep $db_step`;
push @udb, $udb;
}

# return array of udb
return \@udb;
}

# subroutine to blast query against subject db
sub qblasts {
my $query = shift; # reference sequences
my $subdb = shift; # hash of subject genome

mkdir "$qblasts"; # make dir for blast output

# blast against each sub
my %index;
foreach my $dbspecies (sort keys %$subdb) {
	# blast against each split db list
	my @unsort_blastout;
	foreach my $db (@{$subdb->{$dbspecies}}) {
	my $in = $query;
	my ($db_name) = $db =~ /([^\/]+$)/;
	my $unsorted_blastout = "$qblasts/$db_name.unsorted.blast.txt";
	push @unsort_blastout, $unsorted_blastout;
	`usearch -ublast $in -qmask user -db $db -strand both -evalue $E -maxaccepts $maxaccepts -maxrejects $maxrejects -userfields query+target+bits+trow+qstrand+id -userout $unsorted_blastout -threads $cpu`;
	}
	
	# concatenate result
	my $splited_blastout = join " ", @unsort_blastout;
	my $joint_unsorted_blastout = "$qblasts/$dbspecies.qblasts.txt";
		
	`cat $splited_blastout > $joint_unsorted_blastout`;
	
	# get position of each gene infile
	$index{$dbspecies} = gene_index($joint_unsorted_blastout);			
	foreach (@unsort_blastout) {
	unlink $_;
	}
}

# return position of each gene in blast result
return(\%index);
}

# generate position of gene name in file
sub gene_index {
my $unsorted_blastout = shift; # blast result

my $pointer = 0; # file pointer
my %index;
open UNSORT_BLASTOUT, $unsorted_blastout;
while (my $line = <UNSORT_BLASTOUT>) {
my ($gene) = $line =~ /(\S+)\s.+\s\S+\s\S+\s\S+\s\S+/;
push @{$index{$gene}}, $pointer;
$pointer += length($line);
}
close UNSORT_BLASTOUT;

# return position of each gene infile
return(\%index);
}

# subroutine to get sequence from blast output of query against db
sub qblastsout {
no strict; # We need to apply value to filehandle. If strict is used, it would report error

my $subdb = shift; # hash of subject db
my $index = shift; # return position of each gene in blast result
my $qseq = shift; # hash of query

# write query and hit sequence into this file
open QBLASTSOUT_ALL, ">$reblast_input";

mkdir $qblastsout; # create dir of sequences hit with reference

# open filehandle of each blast output
foreach my $sub (sort keys %$subdb) {
my $FH = "${sub}_FH";
my $sub_qblasts_result = "$qblasts/$sub.qblasts.txt";
open $FH, $sub_qblasts_result;
}

# read blast result each subject loci by loci
foreach my $gene (sort keys %$qseq) {
	my $query = $qseq->{$gene};
		
	# print loci name and reference sequences first
	open QBLASTSOUT, ">$qblastsout/$gene.fas"; # write each gene file
	print QBLASTSOUT ">$gene\n$query\n";
	
	print QBLASTSOUT_ALL ">$gene\n$query\n"; # write to merged file
	
	foreach my $sub (sort keys %$subdb) {
	my $FH = "${sub}_FH";
		if (exists $index->{$sub}->{$gene}) { # if sequences of $sub hit with loci $gene 
		my @pos = @{$index->{$sub}->{$gene}}; # get position in blast result
			my %seq;
			my $cnt = 0;
			
			# get information of each blast hit of $gene in reference to subject $sub
			foreach my $pos (@pos) {
			seek $FH, $pos, 0; # skip to that pos
			my $line = <$FH>;
			my ($target, $bits, $sequence, $strand, $id) = $line =~ /\S+\s(.+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)/; # read target name, score, sequence, strand and identity
			$sequence =~ s/-//g;
				if ($strand eq "-") { # reverse complement sequences if strand is -
				$sequence =~ tr/ATCG/TAGC/;
				$sequence = reverse $sequence;
				}
			$seq{$cnt}->{id} = $id;
			$seq{$cnt}->{bits} = $bits;
			$seq{$cnt}->{sequence} = $sequence;
			$cnt++;
			}
			
			# write sequences
			my @bits = sort {$seq{$a}->{bits}<=>$seq{$b}->{bits}} keys %seq; # sort hit
			if (@bits == 1) { # if sequence of reference have only single hit
			print QBLASTSOUT ">$gene\|$sub\n$seq{$bits[-1]}->{sequence}\n"; # write out into separate files
			print QBLASTSOUT_ALL ">$gene\|$sub\n$seq{$bits[-1]}->{sequence}\n"; # writes sequences into single file, this is input of re-blast in next step 
			} elsif (@bits > 1) { # if sequence of reference have multiple hits
				if ($seq{$bits[-1]}->{id}-$seq{$bits[-2]}->{id} >= $second_id_dis_trsd) { # if first hit is significantly better than next hit
				print QBLASTSOUT ">$gene\|$sub\n$seq{$bits[-1]}->{sequence}\n"; # write out into separate files
				print QBLASTSOUT_ALL ">$gene\|$sub\n$seq{$bits[-1]}->{sequence}\n"; # writes sequences into single file
				}
			}
		}
	}
	close QBLASTSOUT;
}

foreach my $sub (sort keys %$subdb) { # close filehandle
my $FH = "${sub}_FH";
close $FH;
}
close QBLASTSOUT_ALL;

$index = ""; # empty index
}

# subroutine to blast merged file against query
sub sblastq {
my $querydb = shift; # array of query db

mkdir $reblast; # create dir of re-blast

# blast against each query db to get unsorted blast output
my @unsorted_blastout;
foreach my $db (@$querydb) {
my $in = $reblast_input;
my ($db_name) = $db =~ /([^\/]+$)/;
my $unsorted_blastout = "$reblast/$db_name.unsorted.reblast.txt";
`usearch -ublast $in -qmask user -db $db -strand both -evalue $E -maxaccepts $maxaccepts -maxrejects $maxrejects -userfields query+target+bits+tlo+thi+qlo+qhi+ql+id -userout $unsorted_blastout -threads $cpu`;
push @unsorted_blastout, $unsorted_blastout;
}

# join and sort output
my $split_unsorted_blastout = join " ", @unsorted_blastout;
my $joint_unsorted_blastout = "$reblast/querydb.unsorted.reblast.txt";
my $joint_sorted_blastout = "$reblast/querydb.sorted.reblast.txt";

`cat $split_unsorted_blastout > $joint_unsorted_blastout`;
`sort -o $joint_sorted_blastout $joint_unsorted_blastout`;

# delete unsort output
foreach (@unsorted_blastout) {
unlink $_;
}
unlink $joint_unsorted_blastout;
}

# subroutine to get sequence with reciprocal best hit
sub ortholog {

# create output dir if it do not exist
if (!(-e $outdir)) {
mkdir $outdir or die "\nERROR: Cannot create output directory \"$outdir\", please check --outdir ($!)\n\n";
}

# create output dir of coding sequences and reference
mkdir $nf;
mkdir $p;

# get orthologs
my $sorted_reblast = "$reblast/querydb.sorted.reblast.txt"; # blast result of hit sequences from sub db->db of reference
my ($genelag, %hash);
open REBLAST, $sorted_reblast;
while (my $line = <REBLAST>) {
# get query, target, bit score, hit position in genome, position of part of sequences hit with genome, length of reference and identity 
my ($field1, $target, $bits, $pos1, $pos2, $qpos1, $qpos2, $qlen, $id) = $line =~ /^(\S+)\s(.+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)$/; 
($pos1, $pos2) = ($pos2, $pos1) if ($pos1 > $pos2);
my ($gene, $sub) = $field1 =~ /([^\|]+)\|*(.*)/; # get loci name and subject name
	
	# write output if all hits of previous gene has been read
	if (($genelag)&&($genelag ne $gene)) {
	writeout(\%hash, $genelag) if (exists $hash{$genelag});
	%hash = ();
	}
	
	if (!($sub)) { # gene name without subject name is reference
		my $query_cov = ($qpos2-$qpos1+1)/$qlen*100; # percentage of reference hit with genome
		# if percentage of reference hit with genome and identity to genome is enough, then this is a qualified reference
		if (($query_cov >= $query_cov_trsd)&&($id >= $query_id_trsd)) { 
		push @{$hash{$gene}->{bits}}, $bits; # save bit score
		push @{$hash{$gene}->{pos1}}, $pos1; # save position hit to genome
		push @{$hash{$gene}->{pos2}}, $pos2;
		push @{$hash{$gene}->{target}}, $target; # save hit chromosome
		}
	} else { # gene name with subject name is subject
		push @{$hash{$sub}->{bits}}, $bits;
		push @{$hash{$sub}->{pos1}}, $pos1;
		push @{$hash{$sub}->{pos2}}, $pos2;
		push @{$hash{$sub}->{target}}, $target;
		push @{$hash{$sub}->{id}}, $id; # save identity
	}
	$genelag = $gene; # update name of previous loci
}
writeout(\%hash, $genelag) if (exists $hash{$genelag}); # write last output
close REBLAST;
}

# subroutine to write output
sub writeout {
my ($hash, $genelag) = @_; # hash of hit, loci name

# get the target and hit position
my @query_target = @{$hash->{$genelag}->{target}};
my @query_pos1 = @{$hash->{$genelag}->{pos1}};
my @query_pos2 = @{$hash->{$genelag}->{pos2}};

my %ortho;
foreach my $sub (sort keys %$hash) {
next if ($sub eq $genelag);

	my @sub_pos1 = @{$hash->{$sub}->{pos1}}; # get hit position of contigs
	my @sub_pos2 = @{$hash->{$sub}->{pos2}};
	my @sub_target = @{$hash->{$sub}->{target}}; # get hit chromosome of contigs
	my @sub_id = @{$hash->{$sub}->{id}}; # identity
	my @sub_bits = @{$hash->{$sub}->{bits}}; # bit score
	
	my (@ortho_sub, @para_sub);
	for (my $i=0; $i<@query_target; $i++) {
		my $query_target = $query_target[$i];
		my $query_pos1 = $query_pos1[$i];
		my $query_pos2 = $query_pos2[$i];
		
		my %bits_id;
		for (my $j=0; $j<@sub_pos1; $j++) {
			my $sub_target = $sub_target[$j];
			my $sub_pos1 = $sub_pos1[$j];
			my $sub_pos2 = $sub_pos2[$j];
			my $sub_id = $sub_id[$j];
			my $sub_bits = $sub_bits[$j];
			
			# a hit is orthologous if contig hit the same target as query, and hit region is overlap with query
			if (($sub_target eq $query_target)&&($sub_pos1 <= $query_pos2)&&($sub_pos2 >= $query_pos1)) {
			push @ortho_sub, $j;
			} else { # else is a paralog
			push @para_sub, $j;
			}
		}
		
		# sort sample name according to bits score
		@ortho_sub = sort {$hash->{$sub}->{bits}->[$a]<=>$hash->{$sub}->{bits}->[$b]} @ortho_sub;
		@para_sub = sort {$hash->{$sub}->{bits}->[$a]<=>$hash->{$sub}->{bits}->[$b]} @para_sub;
		
		if ((@ortho_sub > 0) && (@para_sub == 0)) { # if there's only ortho hit, then it is an ortho
		$ortho{$sub} = "";
		} elsif ((@ortho_sub > 0) && (@para_sub > 0)) { # if there's both ortho and para hit 
			my $ortho_best_id = $hash->{$sub}->{id}->[$ortho_sub[-1]];
			my $para_best_id = $hash->{$sub}->{id}->[$para_sub[-1]];

			# if hit of ortholog is significant better than paralog 
			if ($ortho_best_id - $para_best_id >= $paralog_id_dis_trsd) {
			$ortho{$sub} = "";
			}
		}
	}
}

my $nf_path = "$qblastsout/$genelag.fas"; # read sequence hit with query
my $nf_out_path = "$nf/$genelag.fas"; # coding sequence out path
my $p_out_path = "$p/$genelag.fas"; # amino acid sequence out path
open NF, $nf_path;
open NF_OUT, ">$nf_out_path";
open P_OUT, ">$p_out_path";
<NF>;
chomp(my $query_seq = <NF>);

# write output
my $seqcnt = 0;
while (my $line = <NF>) {
	if ($line =~ />[^\|]+\|*(.*)/) {
	my $sub = $1;
	chomp(my $seq = <NF>);
		if (exists $ortho{$sub}) { # if subject got reciprocal blast hit
			my $first_codon = frame($query_seq, $seq); # find coding frame
			if ($first_codon > 0) {
				my $pseq = translate($seq, $first_codon); # translate
				my ($pseq_nost) = $pseq =~ /(\S+\w+)\**/;
				my @st_codon = $pseq_nost =~ /\*/g;
				next if (@st_codon > $st_trsd); # if number of stop codon is low enough
				my $coding_seq = substr $seq, $first_codon-1, length($pseq_nost)*3; # get coding sequence
				
				if (length($coding_seq) >= $min_seq_len) {
				print NF_OUT ">$sub\n$coding_seq\n"; # write coding sequence output
				print P_OUT ">$sub\n$pseq_nost\n"; # write amino acid sequence output
				$seqcnt++;
				}
			}
		}
	}
}
close NF;
close NF_OUT;
close P_OUT;

if ($seqcnt == 0) { # delete if no reciprocal blast hit was found
unlink $nf_out_path;
unlink $p_out_path;
}
} 

# subroutine to predict coding frame
sub frame {
my $query = shift; # frame corrected query
my $contig = shift; # uncorrected contigs

# translate query
my $quep = translate($query, 1);

# translate contigs from 3 frame
my $seqp1 = translate($contig, 1);
my $seqp2 = translate($contig, 2);
my $seqp3 = translate($contig, 3);

# extract kmer from query
my $quekmer = kmer($quep, 3);

# extract kmer from 3 frame contigs
my $seqkmer1 = kmer($seqp1, 3);
my $seqkmer2 = kmer($seqp2, 3);
my $seqkmer3 = kmer($seqp3, 3);
    
    # find the frame that got highest number of identical kmer
    my $st = 1;
	my $max = 0;
	my $maxcount = 0;
	foreach my $seqkmer ($seqkmer1, $seqkmer2, $seqkmer3) {
		# count the number of kmer
		my $count = 0; 
		foreach (sort keys %$seqkmer) {
		$count++ if (exists $quekmer->{$_});
		}
		
		# get frame with highest number of identical kmer 
		if ($count > $maxcount) {
		$max = $st;
		$maxcount = $count;
		}
	$st++;
	}

return $max;
}

# subroutine to extract the kmer from sequence
sub kmer {
my $seq = shift; # sequence
my $kmer = shift; # define kmer
	
	# extract kmer
    my %khash;
	for (my $i=0; $i<(length($seq)-$kmer+1) ;$i++) {
	my $str = substr $seq, $i, $kmer;
		if (!(exists $khash{$str})) {
		$khash{$str} = 1;
		} else {
        $khash{$str} ++;
		}
	}
	
	# tag identical kmer with different number 
	my %hash;
	foreach (sort keys %khash) {
		for (my $i=1; $i<=$khash{$_}; $i++) {
	    $hash{"$_.$i"} = 1;
		}
	}
	
# return hash of kmer
return \%hash;
}

# subroutine to translate dna to aa
sub translate{
my $seq = shift @_;
my $start_codon = shift @_;
$start_codon --;
my $seq_obj=Bio::Seq->new(-seq=>$seq,-alphabet=>'dna');
my $pro=$seq_obj->translate(-frame => $start_codon, -terminator => "*", -unknown => "X");
$pro=$pro->seq; 
return($pro);
}


# subroutine to report the start of a program chuck
sub function {
	my $function = shift;

	if ($cpu > 1){
	print "Start $function in parallel\n";
	} else {
	print "Start $function in serial\n";
	}
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
Script name: get_orthologues.pl

This is a script to find ortholog genes to target loci in multiple given genomes by reciprocal blast.

Dependencies:
(1) USEARCH v10.0.240 or higher
(2) BioPerl v1.007001 or higher
    
Example usage:
(1) Input single subject database

	perl get_orthologues.pl --query query.dna.fas --querydb query.genome.fas --subdb sub.genome.fas --subname sub --outdir orthologs --cpu 4

Input files:
1. query.dna.fas
2. query.genome.fas
3. sub.genome.fas 

Output files:
4. orthologs

(2) Input multiple subject databases

	perl get_orthologues.pl --query query.dna.fas --querydb query.genome.fas --subdb \"sub.genome.fas sub.genome1.fas|sub1.genome.fas\" --subname \"sub sub1\" --outdir orthologs --cpu 4

Input files:
1. query.dna.fas
2. query.genome.fas
3. genome of sub (sub.genome.fas, sub.genome1.fas), genome of sub1 (sub1.genome.fas)

Output files:
4. orthologs

Options:
--query
  File contains coding DNA sequences only ('--nucleo_out' of auto_assemble_pipeline/data_preparation/query_translate.pl)
--querydb
  Space delimited list of one or more DNA databases of reference species in either in FASTA or UDB format
--subdb
  List of DNA databases of subjects in FASTA or UDB format ONLY, but database format of the same subject need to consistent. Input database list of different subjects are delimited by '|', and database belonging to the same subject are delimited by space. e.g. \"sp1.genome.1.fas sp1.genome.2.fas|sp2.genome.1.udb sp2.genome.2.udb\"
--subname
  Space delimited list of subject name in output, which is one-to-one match to the list of subject databases. 
--outdir
  Name of output directory, which has 2 subfolders including 'nf' for coding sequences and 'p' for AA sequences. DO NOT NAME OUTDIR AS $qblasts, $qblastsout, $reblast_input or $reblast;
--mask_querydb
  Mask databases of reference by DUST, this option is turned off by default
--mask_subdb
  Mask databases of subject by DUST, this option is turned off by default
--E
  Maximum e-value required for a hit, $E in default
--min_len
  Minimum length of sequences required, $min_seq_len in default
--stop_codon
  Maximum number of stop codons allowed in coding sequences, $st_trsd in default
--delete_subdb
  Delete the GENERATED udb databases for subjects, this option is turned off by default 
--cpu
  Limit the number of CPUs, $cpu in default
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