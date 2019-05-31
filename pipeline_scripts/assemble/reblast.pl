#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

# input dir of merged assemblies, database type, mask db or not, input database, output dir of orthologs, name of references in output
# database masker, word size, step between seed, check whether reference sequences are existed in database, sample list
my ($merge, $dbtype, $mudb, $udb, $reblastout, $refname, $dbmasker, $wordsize, $dbstep, $check_query, $query, $samplelist, $help);
my $remove = "failed_seq"; # dir for seq failed reblast test
my $paralog_dir = "$remove/paralog";
# my $multicopy_dir = "$remove/multicopy";
my $nohit_dir = "$remove/nohit";
my $failed_hit = "failed_hit.txt";

my $allnfunmask = "allnf.unmask.fas"; # unmask fas file containing all contigs and query
my $allnfmask  = "allnf.mask.fas"; # masked fas file containing all contigs and query

my $E = 1e-6; # e-value for blast
my $maxaccepts = 0; # terminate ublast search after find $maxaccepts hit, accept all if $maxaccepts == 0 and $maxrejects == 0
my $maxrejects = 0; # terminate ublast search after reject $maxrejects input reads, accept all if $maxaccepts == 0 and $maxrejects == 0
my $db_size = "1G"; # if size of input db is bigger than $db_size, then it will be splited
my $threads = 1; # number of cpu used
my $query_cov_trsd = 100; # minimum percentage of query sequences aligned to database required to deem it as exist in database
my $query_id_trsd = 98; # minimum similarity of between query sequences and database required to deem it as exist in database
my $second_bits_pct = 20; # if bits of best paralog is al least $second_bits_pct pct of larger than bits of best ortholog, then its paralog
my $paralog_id_dis_trsd = 40; # difference of pdistance between 1st hit and 2nd hit (probably paralog)

my $opt = GetOptions( 'merged:s', \$merge,
					  'samplelist:s', \$samplelist,
                      'dbtype:s', \$dbtype,
                      'db:s',\$udb,
                      'ref_name:s', \$refname,
                      'reblastout:s', \$reblastout,
                      'mask_db!', \$mudb,
                      'failed_test:s', \$remove,
                      'E_value_reblast:s', \$E,
                      'db_step:i', \$dbstep,
                      'queryn:s', \$query,
                      'failed_hit:s', \$failed_hit,
                      'check_query!', \$check_query,
                      'cpu:s', \$threads,
                      'help|h!', \$help) or usage();

# display help message if following options is not specified
if ($check_query) { # if check query
	if ( !($opt && $dbtype && $udb && $query) || $help) { 
	usage();
	} 
} else {
	if ( !($opt && $merge && $samplelist && $dbtype && $udb && $reblastout && $refname) || $help) { 
	usage();
	} 
}

# if --check_query is not specifed
if (!($check_query)) {
	# check whether option of dir containing contigs --merge is specified
	if (!($merge)) {
	die ("ERROR: Please specify --merge for input($!)");
	}
	
	# check the existence of dir containing contigs
	if (!(-e $merge)) {
	die ("ERROR: Cannot find $merge ($!)");
	} else { # check the existence of subfolder
	my $nf = "$merge/nf";
	my $f = "$merge/f";
	my $p = "$merge/p";
		if (!(-e $nf)) {
		die ("ERROR: Cannot find $nf ($!)");
		}
		if (!(-e $f)) {
		die ("ERROR: Cannot find $f ($!)");
		}
		if (!(-e $p)) {
		die ("ERROR: Cannot find $p ($!)");
		}
	}
	
	# check whether option for output dir --reblastout is specified
	if (!($reblastout)) {
	die ("ERROR: Please specify --reblastout for output");
	}
}

# set different search parameter for different input db
if ($dbtype eq "nucleo") { # if input db is nucleotide
$wordsize = 8 if (!($wordsize)); # word size 
$dbmasker = "dust" if (!($dbmasker)); # db masker
$dbstep = 2 if (!($dbstep)); # dbstep
} elsif ($dbtype eq "prot") { # if input db is amino acid
$wordsize = 5 if (!($wordsize)); # word size 
$dbmasker = "seg" if (!($dbmasker)); # db masker
$dbstep = 1 if (!($dbstep)); # dbstep
} else { # die if --dbtype is not specified
die ("ERROR: --dbtype must be specified as 'nucleo' or 'prot'\n");
};

my @udb = $udb =~ /(\S+)/g; # get names of all provided db

# check the existence and format of input db
my $fascnt = 0; # count number of fas file
my $udbcnt = 0; # count number of udb file
my $elsecnt = 0; # count number of other file
foreach my $db (@udb) {
	if (!(-e $db)) {
	die "ERROR: Cannot find $db ($!)";
	}
	if ($db =~ /\S+\.fa$|\S+\.fas$|\S+\.fasta$/) {
	$fascnt++;
	} elsif ($db =~ /\S+\.udb$/) {
	$udbcnt++;
	} else {
	$elsecnt++;
	}
}

if ($elsecnt > 0) { # die if found existence of unrecognized format
die "Provided databases must be in either fasta or udb format";
} 
if (($fascnt > 0)&&($udbcnt > 0)) { # die if found existence of both fasta and udb format
die "Provided databases must be in either fasta or udb format";
}

# must provide fast format file, if db need to be masked
if ($mudb) {
	if ($udbcnt > 0) {
	die "ERROR: Please give genomes in fasta format if it need to be masked";
	}
}

# max thread can be set is 10 in usearch
if ($threads > 10) { 
$threads = 10;
}

# translate db size into byte
$db_size = getsize($db_size);

# prepare the db
my $delete_split;
if ($fascnt > 0) { # if fas db is provided
($udb, $delete_split) = split_db(@udb); # split db if it's bigger than $db_size 
$udb = mask_db($udb) if ($mudb); # mask db if --mask_db is specified
$udb = make_udb($udb);  # make udb database
} else { # if udb is provided
my %udb;
push @{$udb{"udb"}}, @udb;
$udb = \%udb; 
}

# check whether all genes in reference is existed in provided genome
if ($check_query) {
	if (!(($check_query)&&($query))) {
	die ("ERROR: Must specify --check_query and --queryn at the same time to check whether query sequences exist in given database \n");
	}
	if (($check_query)&&($reblastout)) {
	die ("ERROR: Don't specify --check_query and --reblastout at the same time\n");
	}
	
	# check whether --queryn is nucleotide sequence
	open QUERY, $query;
	<QUERY>;
	chomp(my $seq = <QUERY>);
	my @nucleo = $seq =~ /[A|T|C|G]/g;
	close QUERY;

	if (@nucleo/length($seq) < 0.5) {
	die "ERROR: --queryn must be nucleotide sequences";
	} 
	
	check_query($query, $udb);
	
	foreach (@$delete_split) {
	unlink $_;
	}
	
	exit;
}

# must provide name of the reference sequences
if (! $refname) {
die "ERROR: Species name of the reference sequences must be provided";
}

# get the name of sample from sample list
my %sample;
open SAMPLELIST, "$samplelist" or die "ERROR: Cannot find $samplelist ($!)";
while (my $samplename = <SAMPLELIST>) {
chomp $samplename;
$sample{$samplename} = "";
}
close SAMPLELIST;

# set the input dir name
my $mergenf = "$merge/nf"; # non-flanking
my $mergef = "$merge/f"; # flanking
my $mergep = "$merge/p"; # amino acid

# get names of all gene
opendir MERGE, $mergenf;
my @merge = grep {$_ !~ /^\./} readdir MERGE;
closedir MERGE;

# merge query and contigs together for ublast
open MERGEALL, ">$allnfunmask";
foreach my $mergefile (@merge) {
	open MERGEFILE, "$mergenf/$mergefile";
	chomp(my $query = <MERGEFILE>);
	chomp(my $qseq = <MERGEFILE>);
	print MERGEALL "$query\n$qseq\n";
	my ($gene) = $query =~ />(\S+)/;
		while (my $line = <MERGEFILE>) {
			if ($line =~ />(\S+)/) {
			my $header = $1;
			my $seq = <MERGEFILE>;
			print MERGEALL ">$gene\|$header\n$seq" if (exists $sample{$header});
			}
		}
	close MERGEFILE;
}
close MERGEALL;

# mask merged file with all query and contigs
maskoutput($allnfunmask, $allnfmask, "dust");

# start ublast
function("ublasting contigs against database");
# set the name of input file and udb
my $in = $allnfmask;
my %udb = %$udb;

# ublast merged file with all query and contigs against all db
my @unsort_blastout;
foreach my $dblist (sort keys %udb) {
	foreach my $subdb (@{$udb{$dblist}}) {
	my ($subdb_name) = $subdb =~ /([^\/]+$)/;
	my $blastresultunsort = "$subdb_name.blastout.unsort.txt";
	push @unsort_blastout, $blastresultunsort;
	`usearch -ublast $in -qmask user -db $subdb -strand both -evalue $E -maxaccepts $maxaccepts -maxrejects $maxrejects -userfields query+target+bits+tlo+thi+qlo+qhi+ql+id -userout $blastresultunsort -threads $threads &> /dev/null`;
	}
}

# join and sort blast output
my $split_unsort_blastout = join " ", @unsort_blastout;
my $joint_unsort_blastout = "$allnfmask.unsort.blastout.txt";
my $blastresult = "$allnfmask.blastout.txt";
`cat $split_unsort_blastout > $joint_unsort_blastout`;
`sort -o $blastresult $joint_unsort_blastout`;

# remove unsorted blast output
unlink $joint_unsort_blastout;
foreach (@unsort_blastout) {
unlink $_;
}
print "Ublast search has been finished\n";

# set name and create outdir 
mkdir $reblastout;
my $reblastoutnf = "$reblastout/nf"; # non-flanking
my $reblastoutf = "$reblastout/f"; # flanking
my $reblastoutp = "$reblastout/p"; # amino acid
mkdir $reblastoutnf;
mkdir $reblastoutf;
mkdir $reblastoutp;

mkdir $remove; # failed reblast test
mkdir $paralog_dir;
# mkdir $multicopy_dir;
mkdir $nohit_dir;

my (%paralog_cnt, %nohit_cnt);
# my (%paralog_cnt, %multicopy_cnt, %nohit_cnt);

# start write output
function("writing out ortholog contigs");
open BLASTRESULT, $blastresult;
my ($genelag, %hash);
while (my $line = <BLASTRESULT>) {
	# get the gene or gene|sample, target, bit score, start pos on db, end pos on db, start of alignment on query, end of alignment on query, query length, similarity to database
	my ($field1, $target, $bits, $pos1, $pos2, $qlo, $qhi, $ql, $id) = $line =~ /^(\S+)\s(.+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)$/;
	($pos1, $pos2) = ($pos2, $pos1) if ($pos1 > $pos2);
	my ($gene, $sub) = $field1 =~ /([^\|]+)\|*(.*)/; # get the gene and contig name
	
	# write output if hits of previous gene has been read
	if (($genelag)&&($genelag ne $gene)) {
		writeout(\%hash, $genelag) if (exists $hash{$genelag});
		%hash = ();
	}
	
	# save hit of $gene into hash
	if (!($sub)) { # if there's no sample name in $field1, then it is a query
		my $query_cov = ($qhi-$qlo+1)/$ql*100;
		if (($query_cov >= $query_cov_trsd)&&($id >= $query_id_trsd)) { # if query is exist in db(base on coverage and similarity)
		push @{$hash{$gene}->{bits}}, $bits;  # push all hit into hash
		push @{$hash{$gene}->{pos1}}, $pos1;
		push @{$hash{$gene}->{pos2}}, $pos2;
		push @{$hash{$gene}->{target}}, $target;
		}
	} else { # if there's sample name in $field1, then it is a contig
		push @{$hash{$sub}->{bits}}, $bits;
		push @{$hash{$sub}->{pos1}}, $pos1;
		push @{$hash{$sub}->{pos2}}, $pos2;
		push @{$hash{$sub}->{target}}, $target;
		push @{$hash{$sub}->{id}}, $id;
	}
	$genelag = $gene;
}
writeout(\%hash, $genelag) if (exists $hash{$genelag}); # write the last gene
close BLASTRESULT;

open FAILED_SUMMARIZE, ">$failed_hit";

print FAILED_SUMMARIZE "Potential paralogs\n";
foreach my $sub (sort keys %paralog_cnt) {
print FAILED_SUMMARIZE "$sub\t$paralog_cnt{$sub}\n";
}
print FAILED_SUMMARIZE "\n";

# print FAILED_SUMMARIZE "Muticopy\n";
# foreach my $sub (sort keys %multicopy_cnt) {
# print FAILED_SUMMARIZE "$sub\t$multicopy_cnt{$sub}\n";
# }
# print FAILED_SUMMARIZE "\n";

print FAILED_SUMMARIZE "No blast hit\n";
foreach my $sub (sort keys %nohit_cnt) {
print FAILED_SUMMARIZE "$sub\t$nohit_cnt{$sub}\n";
}
print FAILED_SUMMARIZE "\n";

close FAILED_SUMMARIZE;

print "All ortholog coding sequences have been written to $reblastout/nf/ \n";
print "All ortholog flanking sequences have been written to $reblastout/f/ \n";
print "All ortholog amino acid sequences have been written to $reblastout/p/ \n";

unlink $allnfmask; # remove merged file containing query and contigd

foreach (@$delete_split) { # remove split database
unlink $_;
}

##################################################
# subroutine
##################################################

# subroutine to translate size into byte unit
sub getsize {
my $size = shift; # size of input file

my ($num, $unit) = $size =~ /(\d+\.*\d*)(\w+)*/;
my $sizeb = do {
	if ($unit) {
		do {
			if ($unit =~ /K|KB/i) {
			$num * 2**10;
			} elsif ($unit =~ /M|MB/i) {
			$num * 2**20;
			} elsif ($unit =~ /G|GB/i) {
			$num * 2**30;
			}
		};
	} else {
	$num * 2**20;
	}
};

# size in byte
return $sizeb;
}

# subroutine to split database 
sub split_db {
my @dblist = @_; # database list

# split each fasta file and return the list of split file
my (%splited_file, @delete_splited);
foreach my $fas (@dblist) {
my $file_size = (stat($fas))[7];

	if ($file_size > $db_size) {
	my $split_num = int($file_size/$db_size)+1;
	my $splited = split_fas($fas, $file_size, $split_num);
	push @{$splited_file{$fas}}, @$splited;
	push @delete_splited, @$splited;
	} else {
	push @{$splited_file{$fas}}, $fas;
	}
}

# hash of splited file of each input db, array of splited name need to be deleted later
return (\%splited_file, \@delete_splited);
}

# subroutine to split fasta file
sub split_fas {
my $fas = shift; # fasta file 
my $file_size = shift; # file size of fasta file
my $split_num = shift; # parts to be splited

# size of each split
my $seq_chunk = int($file_size/$split_num);

# write the split output
my @splited_fas;
my $word_cnt = 0;
my $splited_num = 0;

open SPLITED, ">$fas.$splited_num";
push @splited_fas, "$fas.$splited_num";
open FAS, $fas;
while (my $line = <FAS>) {
	if ($line =~ /^>.+/) {
		if ($word_cnt > $seq_chunk) { # open a new split if size of current split exceed $seq_chunk
		close SPLITED;
		$word_cnt = 0;
		$splited_num++;
		open SPLITED, ">$fas.$splited_num";
		push @splited_fas, "$fas.$splited_num";
		}
	} 
	print SPLITED $line;
	$word_cnt += length($line);
}
close FAS;
close SPLITED;

# warn if cannot split to expected parts
if ($splited_num+1 < $split_num) {
print STDERR "Cannot split $fas into $split_num parts, because there's not enough sequences in $fas\n";
}

# return array of splited fast file name
return(\@splited_fas);
}

# subroutine to mask database
sub mask_db {
my $dbhash = shift; # hash of db

# mask each database
my %masked_db;
foreach my $faslist (sort keys %$dbhash) {
	foreach my $fas (@{$dbhash->{$faslist}}) {
	my $masked = "$fas.masked";
	maskoutput($fas, $masked, $dbmasker, 1);
	push @{$masked_db{$faslist}}, $masked;
	}
}

# return hash of masked database
return \%masked_db;
}

# subroutine to mask input fasta file
sub maskoutput {
# name of unmasked file, masked file, masker, whether need to delete origin file
my ($unmask, $masked, $masker, $unlink) = @_; 

# mask
`usearch -fastx_mask $unmask -fastaout $masked -qmask $masker &> /dev/null`; 

# delete origin file if $unlink is specified
if (!($unlink)) {
unlink $unmask;
}
}

# subroutine to construct udb 
sub make_udb {
my $dbhash = shift; # hash of db

# make database
my %udb;
foreach my $faslist (sort keys %$dbhash) {
	foreach my $fas (@{$dbhash->{$faslist}}) {
	my $udb = "$fas.udb";
	`usearch -makeudb_ublast $fas -output $udb -wordlength $wordsize -dbmask user -dbstep $dbstep &> /dev/null`;
	push @{$udb{$faslist}}, $udb;
	}
}

# return hash of constructed database name
return \%udb;
}

# subroutine to write output base on blast hit
sub writeout {
my ($hash, $genelag) = @_; # hash of sample sequences targeting the same loci, loci name
	
# get the target and hit position
my @query_target = @{$hash->{$genelag}->{target}};
my @query_pos1 = @{$hash->{$genelag}->{pos1}};
my @query_pos2 = @{$hash->{$genelag}->{pos2}};

# get orthologs of each loci
# hash of orthologs and paralogs
my (%ortho, %paralog);
foreach my $sub (sort keys %$hash) {
next if ($sub eq $genelag);
	
	my @sub_pos1 = @{$hash->{$sub}->{pos1}}; # get hit position of contigs
	my @sub_pos2 = @{$hash->{$sub}->{pos2}}; 
	my @sub_target = @{$hash->{$sub}->{target}}; # get hit target chromosome of contigs
	my @sub_id = @{$hash->{$sub}->{id}}; # get similarity to database
	my @sub_bits = @{$hash->{$sub}->{bits}}; # get bit score
	
	# array of orthologous sample name, array of paralogous sample name
	my (@ortho_sub, @para_sub);
	for (my $i=0; $i<@query_target; $i++) {
		# targeted chromosome and position in reference
		my $query_target = $query_target[$i];
		my $query_pos1 = $query_pos1[$i];
		my $query_pos2 = $query_pos2[$i];
		
		# a hit is qualified if contig hit the same targeted chromosome as reference, and hit region is overlap with reference
		for (my $j=0; $j<@sub_pos1; $j++) {
			my $sub_target = $sub_target[$j];
			my $sub_pos1 = $sub_pos1[$j];
			my $sub_pos2 = $sub_pos2[$j];
			my $sub_id = $sub_id[$j];
			my $sub_bits = $sub_bits[$j];
			
			# deem it as ortholog if hit is overlap with reference
			if (($sub_target eq $query_target)&&($sub_pos1 <= $query_pos2)&&($sub_pos2 >= $query_pos1)) {
			push @ortho_sub, $j;
			} else {
			push @para_sub, $j;
			}
		}
	}
	
	# sort sample name according to bits score
	@ortho_sub = sort {$hash->{$sub}->{bits}->[$a]<=>$hash->{$sub}->{bits}->[$b]} @ortho_sub;
	@para_sub = sort {$hash->{$sub}->{bits}->[$a]<=>$hash->{$sub}->{bits}->[$b]} @para_sub;
	
	if ((@ortho_sub > 0) && (@para_sub == 0)) { # if there's only ortho hit, then it is an ortho
	$ortho{$sub} = "";
	} elsif ((@ortho_sub == 0) && (@para_sub > 0)) { # if there's only para hit, then it is an para
	$paralog{$sub} = "";
	} elsif ((@ortho_sub > 0) && (@para_sub > 0)) { # if there's both ortho and para hit 
		my $ortho_best_id = $hash->{$sub}->{id}->[$ortho_sub[-1]];
		my $para_best_id = $hash->{$sub}->{id}->[$para_sub[-1]];
		
		# if id interval between orthologous and paralogous hit is large enough then it is a ortholog
		if ($ortho_best_id - $para_best_id >= $paralog_id_dis_trsd) {
		$ortho{$sub} = "";
		} else {
		$paralog{$sub} = "";
		}
	}
}

# input sequence
my $nf = "$mergenf/$genelag.fas";
my $f = "$mergef/$genelag.fas";
my $p = "$mergep/$genelag.fas";

# read all sequence into hash

# read non-flanking sequence
my %seqhash;
open NF, $nf;
chomp(my $qnfref = <NF>);
$qnfref = ">$refname" if ($refname); # change the name of query if it need to be rename
chomp(my $qnfseq = <NF>);
while (my $line = <NF>) {
my ($sub) = $line =~ />(\S+)/;
chomp(my $seq = <NF>);
$seqhash{$sub}->{nf} = $seq;
}
close NF;

# read flanking sequence
open F, $f;
chomp(my $qfref = <F>);
$qfref = ">$refname" if ($refname);
chomp(my $qfseq = <F>);
while (my $line = <F>) {
my ($sub) = $line =~ />(\S+)/;
chomp(my $seq = <F>);
$seqhash{$sub}->{f} = $seq;
}
close F;

# read amino acid sequence
open P, $p;
chomp(my $qpref = <P>);
$qpref = ">$refname" if ($refname);
chomp(my $qpseq = <P>);
while (my $line = <P>) {
my ($sub) = $line =~ />(\S+)/;
chomp(my $seq = <P>);
$seqhash{$sub}->{p} = $seq;
}
close P;

# write out ortholog and paralog sequence

# flag for open and close of the file handle
my $orthoopen = 0; # ortholog file 
my $paraopen = 0; # paralog file
my $nohitopen = 0; # no hit file
foreach my $sub (sort keys %seqhash) {
next if (!(exists $sample{$sub}));

	if (exists $ortho{$sub}) { # write to ortholog file
		if ($orthoopen == 0) {
		open ORTHONF, ">$reblastoutnf/$genelag.fas";
		open ORTHOF, ">$reblastoutf/$genelag.fas";
		open ORTHOP, ">$reblastoutp/$genelag.fas";
		print ORTHONF "$qnfref\n$qnfseq\n";
		print ORTHOF "$qfref\n$qfseq\n";
		print ORTHOP "$qpref\n$qpseq\n";
		$orthoopen = 1;
		}
		print ORTHONF ">$sub\n$seqhash{$sub}->{nf}\n";
		print ORTHOF ">$sub\n$seqhash{$sub}->{f}\n";
		print ORTHOP ">$sub\n$seqhash{$sub}->{p}\n";
	} elsif (exists $paralog{$sub}) { # write to paralog file
		if ($paraopen == 0) {
		open PARANF, ">$paralog_dir/$genelag.fas";
		$paraopen = 1;
		}
		print PARANF ">$sub\n$seqhash{$sub}->{nf}\n";
		$paralog_cnt{$sub}++;
	} else {
		if ($nohitopen == 0) { # write to nohit file
		open NOHITNF, ">$nohit_dir/$genelag.fas";
		$nohitopen = 1;
		}
		print NOHITNF ">$sub\n$seqhash{$sub}->{nf}\n";
		$nohit_cnt{$sub}++;
	}
}

# close file handle if it has been opened
if ($orthoopen == 1) {
close ORTHONF;
close ORTHOF;
close ORTHOP;
}
if ($paraopen == 1) {
close PARANF;
}
if ($nohitopen == 1) {
close NOHITNF;
}
} 

# subroutine to check whether sequences in query exist in provided database
sub check_query {
my $query = shift; # DNA sequence of reference
my $udb = shift; # database

maskoutput($query, "$query.masked",  "dust", "1"); # mask input query

my $in = "$query.masked"; # masked query
my %udb = %$udb; # list of db

# ublast query against each db
my @test_blastout;
foreach my $dblist (sort keys %udb) {
	foreach my $querydb (@{$udb{$dblist}}) {
	my $queryblastresult = "$querydb.test.masked.blastout.txt";
	push @test_blastout, $queryblastresult;
	`usearch -ublast $in -qmask user -db $querydb -strand both -evalue $E -maxaccepts $maxaccepts -maxrejects $maxrejects -userfields query+id+qlo+qhi+ql -userout $queryblastresult -threads $threads &> /dev/null`;
	}
}

# join all blast output together
my $split_test_blastout = join " ", @test_blastout;
my $joint_test_blastout = "$query.test.masked.blastout.txt";
`cat $split_test_blastout > $joint_test_blastout`;

# construct a hash containing all gene name
my %query;
open QUERY, "$query";
while (my $line = <QUERY>) {
	if ($line =~ />(\S+)/) {
	$query{$1} = "";
	}
}
close QUERY;

# delete gene name from hash if it got a qualified hit
my %multi_hit;
open QUERYBLASTOUT, $joint_test_blastout;
while (my $line = <QUERYBLASTOUT>) {
my ($query, $id, $qlo, $qhi, $ql) = $line =~ /^(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)$/; # get name of gene, similarity,  start of alignment on query, end of alignment on query, query length 
my $query_cov = ($qhi-$qlo+1)/$ql*100; # percentage of alignment of query
	if (($id >= $query_id_trsd)&&($query_cov >= $query_cov_trsd)) { # if this hit is qualified (enough coverage and identity)
	delete $query{$query} if (exists $query{$query}); # delete it from hash
	}
}
close QUERYBLASTOUT;

my $misscnt = scalar keys %query; # find the number of gene remain in hash, which is not hit with db
if ($misscnt > 0) { # print the name of gene if it is not exist in provided db
say STDOUT "#### $misscnt genes below are not found in provided database ####";
map {say STDOUT "$_"} (sort keys %query);
} else {
say STDOUT "#### All genes are found in provided database ####";
}

unlink "$query.masked";
unlink $joint_test_blastout;
foreach (@test_blastout) {
unlink $_;
}
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
Script name: reblast.pl

This is a script to remove potential paralogs by ublast under usearch

Dependencies:
(1) USEARCH v10.0.240 or higher

Example usage:
(1) Remove potential paralogs in mergeout

	perl reblast.pl --merged mergeout --db db/path --reblastout reblastout --dbtype nucleo --samplelist samplelist.txt

(2) Check whether sequences in query.fas can be found in db/path
	
	perl reblast.pl --queryn query.fas --db db/path --check_query --dbtype nucleo

Input files:
(1) mergeout
(2) db/path 
(3) samplelist.txt
(4) query.fas (if --check_query is specifed)

Output files:
(1) reblastout
(2) failed_seq.txt (tab delimited text which discribes number of multicopies, paralogs and sequences without blast hit for each sample)

Options:
--merged
  Directory containing further assembled nucleotide coding sequences(nf), whole sequences with flanking region(f), and amino acid sequences(p)
--samplelist
  A list of sample name. Only samples in the list will be assembled
--dbtype
  Type of database either \'nucleo\' or \'prot\'
--db
  Space delimited list of one or more databases in udb or fasta format. Database will be splited if its size is larger than $db_size.
--reblastout
  Directory containing ortholog nucleotide coding sequences(nf), flanking sequences(f), and amino acid sequences(p)
--ref_name
  Name of reference sequences in final assemblies
--failed_hit
  Tab delimited text file recording number of sequences failed reblast test
--mask_db
  Mask input database, disabled in default
--failed_test
  Directory containing paralogs and sequences without blast hit, named in 'failed_seq' default
--E_value_reblast
  Maximum e-value required for a hit, $E in default
--dbstep
  Every Nth ublast database word should be indexed, bigger means smaller database, 2 in default for nucleotide database, 1 in default for peptide database
--queryn
  Full coding nuleotide sequences of target loci in fasta format. It is always used with --check_query to check whether these sequences can be found in given database
--check_query
  Check query sequences existing in given database, and return list of missing query, then exit
--cpu
  Limit the number of CPUs, $threads in default
--help, -h
  Show this help message and exit

Author: Hao Yuan                                                                     
        Shanghai Ocean University                                               
        Shanghai, China, 201306                                                               
                                                                                         
Created by: Nov 20, 2018                                                              
                                                                                         
Last modified by: 
";
exit;
}
