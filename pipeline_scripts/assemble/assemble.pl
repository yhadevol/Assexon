#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

# list of sample name
my $samplelist = "samplelist.txt"; 

# statistics
my $readcount = "rmdup_reads_bases_count.txt"; # number of reads and bases in origin and deduplicated fq
my $failed_hit = "failed_hit.txt"; # number of sequences failed reblast in each sample
my $enriched_cnt = "enriched_gene.txt"; # list of number of enriched gene for sample

# intermediate output dir and file
my $run_dir = "run_dir"; # dir where run the scripts
my $rmdup = "$run_dir/rmdup"; # dir where remove PCR duplicates (rmdup.pl)
my $rmdup_super_reads_dir = "$run_dir/merged_reads"; # dir containing concatenated reads (rmdup.pl)
my $parsed = "$run_dir/parsed"; # dir where parse reads to genes (ubandpx.pl)
my $ubxandp_blastout_dir = "$run_dir/blastout"; # dir where put blast results of reads parsing (ubandpx.pl) 
my $blastunsortdir = "$run_dir/blastout_unsorted"; # dir containing unsorted blast result for each splited deduplicated reads
my $merged_rmdup = "${rmdup}_merge"; # for merged deduplicated reads
my $assembled = "$run_dir/assembled"; # dir where assemble reads into contigs (sga_assemble.pl)
my $assemble_temp_dir = "$run_dir/temp"; # run dir of sga (sga_assemble.pl)
my $filtered = "$run_dir/filtered"; # dir where filter unqualified contigs and find contigs which might be furtherly assembled (exonerate_best.pl)
my $querypdir = "$run_dir/queryp"; # dir including amino acid sequences of query (exonerate_best.pl)
my $merge = "$run_dir/merged"; # dir where assemble contigs further and get best contigs of each gene (merge.pl)
my $merge_getbesttemp = "$run_dir/getbesttemp"; # run dir for subroutine getbest in merge.pl
my $reblastout = "$run_dir/reblastout"; # dir where select orthologs (reblast.pl)
my $failed_reblast = "$run_dir/failed_seq"; # dir where put sequence failed reblast (reblast.pl)

# default parameters
my $E_value_ubandpx = 1e-4; # E value of ublast in ubandpx.pl
my $E_value_reblast = 1e-6; # E value of ublast in reblast.pl
my $similarity = 60; # minimum similarity between query and contigs
my $len_percentage = 80; # minimum percentage of reference covered in reference-contig aligment
my $stop_codons = 0; # maximum number of stop codon allowed in contigs
my $nonflank_len = 0; # minimum length required of non-flanking region

# path to scripts
my $script_path = ""; 

my @ARGVoption = @ARGV;

# apply default parameters to corresponding options
my %opthash = ('rmdup' => $rmdup,
				'parsed' => $parsed,
				'assembled' => $assembled,
				'filtered' => $filtered,
				'merged' => $merge,
				'reblastout' => $reblastout,
				'E_value_ubandpx' => $E_value_ubandpx,
				'similarity' => $similarity,
				'len_percentage' => $len_percentage,
				'stop_codons' => $stop_codons,
				'nonflank_len' => $nonflank_len,
				'E_value_reblast' => $E_value_reblast,
				'readcount' => $readcount,
				'failed_hit' => $failed_hit,
				'enriched_cnt' => $enriched_cnt,
				'script_path' => $script_path);

# all options
GetOptions(\%opthash, 
            # input dir or file
			'samplelist=s', # sample list
			'trimmed=s', # trimmed 
			'rmdup=s', # output dir of rmdup.pl
			'parsed=s', # output dir of bxandp.pl
			'assembled=s', # output dir of sga_assemble.pl
			'filtered=s', # output dir of exonerate_best.pl
			'merged=s', # output dir of merge.pl
			'reblastout=s', # output dir of reblast.pl
			'outdir=s', # output dir of nf, f, and p assemblies
			'queryp=s', # DNA coding sequences of reference
			'queryn=s', # AA sequences of reference
			'db:s', # fas or udb database of genome
			'script_path:s', # path to scripts wrapped by assemble.pl
			# parameters
			'E_value_ubandpx=f', # E value of ublast in ubandpx.pl
			'E_value_reblast=f', # E value of ublast in reblast.pl
			'similarity=i', # minimum similarity between query and contigs
			'len_percentage=i', # minimum percentage of reference covered in reference-contig aligment
			'stop_codons=i', # maximum number of stop codon allowed in contigs
			'nonflank_len=i', # minimum length required of non-flanking region 
			'dbtype=s', # database type nucleo or prot
			'mask_db!', # mask database or not
			'ref_name:s', # name of reference in final assemblies
			'cpu=i', # number of process used
			# run partial pipeline
			'stop_after_rmdup!', # stop pipeline after rmdup.pl finished
			'stop_after_ubxandp!', # stop pipeline after ubxandp.pl finished
			'stop_after_sga_assemble!', # stop pipeline after sga_assemble.pl finished
			'stop_after_exonerate_best!', # stop pipeline after exonerate_best.pl finished
			'stop_after_merge!', # stop pipeline after merge.pl finished
			'restart_from_ubxandp!', # start pipiline from ubxandp.pl
			'restart_from_sga_assemble!', # start pipiline from sga_assemble.pl
			'restart_from_exonerate_best!', # start pipiline from exonerate_best.pl
			'restart_from_merge!', # start pipiline from merge.pl
			'restart_from_reblast!', # start pipiline from reblast.pl
			# check completeness of query
			'check_query!',
			# check dependencies of assemble.pl
			'check_depends!',
			# clean intermediate output
			'clean!',
			'help|h!') or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# call help message if -h, --help or nothing is specified
if ((@ARGVoption == 0)||($opthash{"help"})) {
usage();
}

# clean intermediate output
clean() if ($opthash{"clean"});

# get absolute path of run dir
$run_dir = absolute_path($run_dir);

# check dependencies
check_depends() if ($opthash{"check_depends"}); 

# get number of avaliable cpu automatically 
if (!($opthash{"cpu"})) { 
$opthash{"cpu"} = get_cpu_num();
}

# report pipeline will be run in serial or parallel
if ($opthash{"cpu"} > 1) { # more than 1 process
say STDOUT "Wrapper will be run in $opthash{cpu} threads\n";
} elsif ($opthash{"cpu"} == 1) { # 1 process
say STDOUT "Wrapper will be run in serial\n";
} else { # illegal number of process
die "\nERROR: Number of cpu must equal to or greater than 1\n\n";
}

# ensure the path to interpreter and scripts
my $report_script_path; # scalar for script path showed in STDOUT
$script_path = $opthash{"script_path"}; # get path to wrapped scripts
if (! $script_path) { # if --script_path is not specified
$report_script_path = "\$PATH"; # report wrapped scripts are located at $PATH
} else { # if --script_path is specified
$report_script_path = $script_path; # report wrapped scripts are located at specified path

# get absolute path to script and add a "/" after it
$script_path = absolute_path($script_path); 
$script_path .= "/";
}

# get the interpreter
my @args = `ps -o args= -p $$` =~ /(\S+)/g; # get full command line in current process $$

# get arguments before current script name $0, which are interperter
my @interpreter;
foreach my $args (@args) { # get arguments before script name
	if ($args ne $0) {
	push @interpreter, $args;
	} else {
	last;
	}
}
my $report_interpreter = join " ", @interpreter; # join these arguments together

# report interpreter and path to wrapped script in current run
say STDOUT "You will use '$report_interpreter' to run the scripts under $report_script_path\n";

my $interpreter;
if ($script_path) { # if --script_path is specified, the interpreter is $report_interpreter
$interpreter = $report_interpreter;
} else { # elsif we use the interpreter specified at the beginning of the script
$interpreter = "";
}

# check whether sequences of query exists in database
check_query() if ($opthash{"check_query"});

# option in each step
my @rmdup = qw/trimmed rmdup samplelist readcount cpu/;
my @ubxandp = qw/queryp rmdup parsed E_value_ubandpx samplelist cpu/;
my @sga_assemble = qw/parsed assembled samplelist cpu/;
my @exonerate = qw/queryp assembled filtered similarity stop_codons samplelist cpu/;
my @merge = qw/queryn queryp filtered merged similarity len_percentage nonflank_len stop_codons samplelist cpu/;
my @reblast = qw/merged dbtype mask_db db reblastout ref_name E_value_reblast failed_hit samplelist cpu/;

# name of all steps
my @allprocedure = qw/rmdup ubxandp sga_assemble exonerate_best merge reblast/;
# save into hash $allprocedure{order of step} = name of step
my %allprocedure;
for (0..(@allprocedure-1)) {$allprocedure{$allprocedure[$_]} = $_;}

# generate required steps in current run
my ($startstep, $stopstep) = procedure();

# mandatory input of each step(except input from last step)
my @rmdup_depend = qw/trimmed/;
my @ubxandp_depend = qw/queryp/;
my @sga_assemble_depend = qw//;
my @exonerate_best_depend = qw/queryp/;
my @merge_depend = qw/queryp queryn/;
my @reblast_depend = qw/db dbtype ref_name/;

# hash of mandatory input
my %depend = ("rmdup" => \@rmdup_depend,
				"ubxandp" => \@ubxandp_depend,
				"sga_assemble" => \@sga_assemble_depend,
				"exonerate_best" => \@exonerate_best_depend,
				"merge" => \@merge_depend,
				"reblast" => \@reblast_depend);

# generate sample list
my $samplelist_hash = samplelist($startstep);

# check existence of all requirements
check_file_exists($startstep, $stopstep, $samplelist_hash);

# change relative path of input to absolute path
$opthash{"trimmed"} = absolute_path($opthash{"trimmed"}) if ($opthash{"trimmed"});
$opthash{"queryp"} = absolute_path($opthash{"queryp"}) if ($opthash{"queryp"});
$opthash{"queryn"} = absolute_path($opthash{"queryn"}) if ($opthash{"queryn"});
$opthash{"outdir"} = absolute_path($opthash{"outdir"}) if ($opthash{"outdir"});
$opthash{"samplelist"} = absolute_path($opthash{"samplelist"});

# change relative path of databases (probably have multiple input) to absolute path
if ($opthash{"db"}) {
	my @db = $opthash{"db"} =~ /\S+/g; # get path of each database
	my @absolute_db_path;
	foreach my $db (@db) { # get absolute path of each database
	push @absolute_db_path, absolute_path($db);
	}
	$opthash{"db"} = join " ", @absolute_db_path; # join all path into a string delimited by space
	$opthash{"db"} = "\"" . $opthash{"db"} . "\""; # add " at two sides of string
}

# change relative path of output statistics files to absolute path
$opthash{"readcount"} = absolute_path($opthash{"readcount"});
$opthash{"failed_hit"} = absolute_path($opthash{"failed_hit"});
$opthash{"enriched_cnt"} = absolute_path($opthash{"enriched_cnt"});

# running dir
if (!(-e $run_dir)) { # create running dir if it do not exist
mkdir $run_dir;
} else { # if running dir exists clean up data left before previous run
	my ($previous_num) = `find $run_dir|grep "_previous"|wc -l` =~ /(\d+)/; # number of previous data
	if ($previous_num > 0) { # delete if there are left data
	say STDOUT "Delete $run_dir/*_previous which are intermidate files left before previous time \n";
	`find $run_dir|grep "_previous"|xargs -I {} rm -rf {}`;
	}	
}

# since running dir has been created, we can change relative path of intermediate output to absolute path
$opthash{"rmdup"} = absolute_path($opthash{"rmdup"});
$opthash{"parsed"} = absolute_path($opthash{"parsed"});
$opthash{"assembled"} = absolute_path($opthash{"assembled"});
$opthash{"filtered"} = absolute_path($opthash{"filtered"});
$opthash{"merged"} = absolute_path($opthash{"merged"});
$opthash{"reblastout"} = absolute_path($opthash{"reblastout"});
$rmdup_super_reads_dir = absolute_path($rmdup_super_reads_dir);
$ubxandp_blastout_dir = absolute_path($ubxandp_blastout_dir);
$blastunsortdir = absolute_path($blastunsortdir);
$merged_rmdup = absolute_path($merged_rmdup);
$assemble_temp_dir = absolute_path($assemble_temp_dir); 
$querypdir = absolute_path($querypdir);
$merge_getbesttemp = absolute_path($merge_getbesttemp);
$failed_reblast = absolute_path($failed_reblast);

# start running
my $outdir = $opthash{"outdir"}; # specify scalar of output dir which is frequently used in following code

# change dir to running dir
chdir $run_dir;

# run accordingly from specified starting step to ending step
for ($allprocedure{$startstep}..$allprocedure{$stopstep}) {
my $step = $allprocedure[$_]; # step name
my $starttime = time(); # start time of the step
	
	# if current step is rmdup.pl
	if ($step eq "rmdup") {
		my $rmdup_dir = $opthash{"rmdup"}; # get dir
		if (-e $rmdup_dir) { # back up data if exists output from last run
		`mv $rmdup_dir ${rmdup_dir}_previous`; 
		}
		if (-e $rmdup_super_reads_dir) {
		`mv $rmdup_super_reads_dir ${rmdup_super_reads_dir}_previous`;
		}
		
		# run script of rmdup.pl
		run_script("rmdup.pl", \@rmdup);
	} elsif ($step eq "ubxandp") { # if current step is run ubxandp.pl
		my $parsed_dir = $opthash{"parsed"};
		if (-e $parsed_dir) { 
		`mv $parsed_dir ${parsed_dir}_previous 
		mv $ubxandp_blastout_dir ${ubxandp_blastout_dir}_previous`; # back up previous data
		}
		if (-e $blastunsortdir) {
		`mv $blastunsortdir ${blastunsortdir}_previous`;
		}
		if (-e $merged_rmdup) {
		`mv $merged_rmdup ${merged_rmdup}_previous`;
		}
		
		run_script("ubxandp.pl", \@ubxandp);
	} elsif ($step eq "sga_assemble") { # if current step is run sga_assemble.pl
		my $assembled_dir = $opthash{"assembled"};
		if (-e $assembled_dir) {
		`mv $assembled_dir ${assembled_dir}_previous`; # back up previous data
		}
		if (-e $assemble_temp_dir) {
		`mv $assemble_temp_dir ${assemble_temp_dir}_previous`;
		}
		
		run_script("sga_assemble.pl", \@sga_assemble);
	}  elsif ($step eq "exonerate_best") { # if current step is run exonerate_best.pl
		my $filtered_dir = $opthash{"filtered"};
		if (-e $filtered_dir) {
		`mv $filtered_dir ${filtered_dir}_previous`; # back up previous data
		}	
		
		makequerypdir($opthash{"queryp"}); # create dir of AA sequences of reference, each file contain one sequence, this is used in reference-contig alignment
		
		run_script("exonerate_best.pl", \@exonerate);
	} elsif ($step eq "merge") { # if current step is run merge.pl
		my $merged_dir = $opthash{"merged"};
		if (-e $merged_dir) {
		`mv $merged_dir ${merged_dir}_previous`; # back up previous data
		}
		if (-e $merge_getbesttemp) {
		`mv $merge_getbesttemp ${merge_getbesttemp}_previous`;
		}
		
		makequerypdir($opthash{"queryp"}); # create dir of AA sequences of reference if it not exists
		
		run_script("merge.pl", \@merge);
	} elsif ($step eq "reblast") { # if current step is run reblast.pl
		my $reblastout_dir = $opthash{"reblastout"};
		if (-e $reblastout_dir) {
		`mv $reblastout_dir ${reblastout_dir}_previous
		mv $failed_reblast ${failed_reblast}_previous`;
		}
		
		run_script("reblast.pl", \@reblast);
	
	# move output of reblast to output dir
	`mv -f $reblastout_dir $outdir`;
	
	# report output of final assemblies
	say STDOUT "Assembled contigs have been written to $outdir\n";
	}
	
	# get duration of current step
	my $endtime = time();
	my $duration = $endtime-$starttime;
	
	# change duration time into format of "xx h xx m xx s"
	my $runtime = runtime($duration);
	
	# report duration time
	say STDOUT "$step ran for $runtime\n";
}

# get out of running dir
chdir "..";

# summary number of enriched genes for each sample if run to last step
if ($stopstep eq "reblast") {
	
	# count number of enriched gene in each taxon
	my %enriched_cnt;
	opendir REBLASTOUTDIR, "$outdir/nf";
	while (my $gene = readdir REBLASTOUTDIR) {
	next if ($gene =~ /^\./);
		
		# read assmebled loci
		my $genepath = "$outdir/nf/$gene";
		open GENEFILE, "$genepath";
		
		# skip reference
		<GENEFILE>;
		<GENEFILE>;
		
		# get name of enriched taxa
		while (my $line = <GENEFILE>) {
			if ($line =~ />(\S+)/) {
			$enriched_cnt{$1}++;
			}
		}
		close GENEFILE;
	}
	closedir REBLASTOUTDIR;
	
	# write out
	open ENRICHED_CNT, ">$enriched_cnt";
	
	# get number of all loci
	my ($allgenenum) = `grep ">" $opthash{queryp}|wc -l` =~ /(\d+)/;
	print ENRICHED_CNT "total\t$allgenenum\n\n";
	
	# header
	print ENRICHED_CNT "Sample\tNum. of enriched genes\tPercentage of enriched genes(%)\n";
	
	# print taxa name, enriched loci and enriched pct of loci
	foreach my $sub (sort keys %enriched_cnt) {
	my $taxa_enriched = $enriched_cnt{$sub};
	my $enrich_pct = sprintf("%.2f", $taxa_enriched/$allgenenum*100);
	print ENRICHED_CNT "$sub\t$taxa_enriched\t$enrich_pct\n";
	}
	close ENRICHED_CNT;
}

##################################################
# subroutine
##################################################

# subroutine to print usage
sub usage {
print STDERR "
assemble.pl v1.00

This is a wrapper around several scripts to extract orthologs assemblies from reads of exon capture.

Dependencies:
Softwares (Plesase put them under \$PATH):
(1) Perl v5.18 or higher
(2) USEARCH v10.0.240 or higher
(3) SGA v0.10.15 or higher
(4) Exonerate v2.2.0 or higher
(5) Perl module: 
    1. Bio::Seq (included in Bioperl)
    2. Parallel::Forkmanager
    3. Sys::Info
(6) wrapped scripts (Put them under \$PATH. If they are located at somewhere else, please specify --script_path ):
    1. rmdup.pl
    2. ubxandp.pl
    3. sga_assemble.pl
    4. exonerate_best.pl 
    5. merge.pl
    6. reblast.pl

Before runnning this wrapper, all reads need to be expanded and cleaned by provided scripts (gunzip_Files.pl, trim_adaptor.pl) or do it by yourself. Cleaned reads need to be named and placed in following structure:

	trimmed_dir
	     ├── sample1_R1.fq
	     └── sample1_R2.fq
	     ├── sample2_R1.fq
	     └── sample2_R2.fq

Example usage: 
(1) Assemble from raw reads to orthologous assemblies

	If wrapped scripts are located under \$PATH:

	perl assemble.pl \\
	--trimmed trimmed \\
	--queryp species.aa.fas \\
	--queryn species.dna.fas \\
	--db species.fas \\
	--dbtype nucleo \\
	--ref_name species \\
	--outdir assemble_results

or 
	If wrapped scripts are located under specfied path ./xxx/xxx:
	
    perl assemble.pl \\
	--trimmed trimmed \\
	--queryp species.aa.fas \\
	--queryn species.dna.fas \\
	--db species.fas \\
	--dbtype nucleo \\
	--ref_name species \\
	--outdir assemble_results \\
	--script_path ./xxx/xxx
	
Input files: 
1) trimmed
2) species.aa.fas
3) species.dna.fas
4) species.fas

Output files:
1) assemble_results

Output statistics files:
1) rmdup_reads_bases_count.txt (file summarizing number of reads and bases of trimmed and deduplicated reads) 
2) failed_hit.txt (file summarizing number of sequences failed in reblast test in each sample)
3) enriched_gene.txt (file summarizing number and percentage of enriched loci for each sample)

(2) Check whether all sequences in DNA sequences of reference can be found in provided database:

	perl assemble.pl \\
	--check_query \\
	--queryn species.dna.fas \\
	--db species.fas \\
	--dbtype nucleo

Input files: 
1) species.dna.fas
2) species.fas

Output files:
1) Report in STDOUT

(3) Check dependencies of assemble.pl:

	perl assemble.pl \\
	--check_depends

Input files: 
1) Nothing

Output files:
1) Report in STDOUT

(4) Clean all intermediate file:
	
	perl assemble.pl \\
	--clean

Input files: 
1) Nothing

Output files:
1) Nothing
	
(5) Run partial pipeline

For example, if we run from exonerate_best.pl, only input of exonerate_best.pl and following scripts including merge.pl and reblast.pl are required. Thus, input files are:

	perl assemble.pl \\
	--queryp species.aa.fas \\ (exonerate_best.pl, merge.pl)
	--queryn species.dna.fas \\ (merge.pl)
	--db species.fas \\ (reblast.pl)
	--dbtype nucleo \\ (reblast.pl)
	--ref_name species \\ (reblast.pl)
	--outdir assemble_results \\ (reblast.pl)
	--samplelist samplelist.txt \\
	--restart_from_exonerate_best 

1) species.aa.fas
2) species.dna.fas
3) species.fas
4) ./run_dir/assembled (intermediate output of sga_assemble.pl which is last step of exonerate_best.pl)

Output files:
1) assemble_results

Output statistics files:
1) failed_hit.txt
2) enriched_gene.txt

More detailed usage please refer to 'tutorial.txt'

Frequently used options:
--trimmed
  Directory containing reads without adaptor and low quality bases
--queryp
  Amino acid sequences of reference in fasta format 
--queryn
  Nucleotide sequences of reference in fasta format
--db
  Path to database of the species of reference in either udb or fasta format. Multiple inputs are supported. Both DNA and AA databases can be input. If fasta database are input, a corresponding udb database will be automatically generated
--dbtype
  Database type either \'nucleo\' or \'prot\'
--ref_name
  Name of reference sequences in the output of last step (reblast.pl)
--outdir
  Directory containing final assemblies including coding sequence with (f) or without flanks (nf), amino acid sequences (p)
--script_path
  Path to the wrapped scripts, makes sure all scripts are placed under this path
--samplelist
  List of sample names. Only samples in the list will be assembled. A list comprising all samples will be automatically generated if wrapper is run from the first step (rmdup.pl) and this option is not specified
--cpu
  Limit the number of CPUs. Default is to use all cores available
--clean
  Remove run_dir and all intermediate output under it
--help , -h
  Show this help message and exit

Check dependencies and input:
--check_depends
  Check all dependencies for this wrapper including perl interpreter, involved perl module, wrapped scripts and softwares, then exit
--check_query
  Check DNA sequences of reference (--queryn) existing in given database, and return list of missing loci, then exit

Parameters tuning:
--E_value_ubandpx
  Maximum e-value required for a BLAST hit, $E_value_ubandpx in default
--E_value_reblast
  Maximum e-value required for a BLAST hit, $E_value_reblast in default
--similarity
  Minimum similarity required between contigs and query, $similarity in default
--len_percentage
  Percentage of reference covered in local alignment of contig and reference, $len_percentage in default
--stop_codons
  Maximum stop codon allowed in a contig, $stop_codons in default
--nonflank_len
  Minimum length of non-flanking sequences required, $nonflank_len in default
--mask_db
  Mask input database, disabled in default

Run partial pipeline:
--stop_after_rmdup
  Stop the wrapper after running rmdup.pl
--stop_after_ubxandp
  Stop the wrapper after running ubxandp.pl
--stop_after_sga_assemble
  Stop the wrapper after running sga_assemble.pl
--stop_after_exonerate_best
  Stop the wrapper after running exonerate_best.pl
--stop_after_merge
  Stop the wrapper after running merge.pl
--restart_from_ubxandp
  Start running the wrapper from ubxandp.pl
--restart_from_sga_assemble
  Start running the wrapper from sga_assemble.pl
--restart_from_exonerate_best
  Start running the wrapper from exonerate_best.pl
--restart_from_merge
  Start running the wrapper from merge.pl
--restart_from_reblast
  Start running the wrapper from reblast.pl
  
Author: Hao Yuan                                                                     
        Shanghai Ocean University                                               
        Shanghai, China, 201306                                                               
                                                                                         
Created by: Nov 20, 2018                                                              
                                                                                         
Last modified by: 
";
exit;
}

# subroutine to identify start and end step from --restart_from_* and --stop_after_*
sub procedure {
	# get specified --restart_from_*
	my @restart_procedure = qw/restart_from_ubxandp restart_from_sga_assemble restart_from_exonerate_best restart_from_merge restart_from_reblast/;
	my @restart;
	for (0..(@restart_procedure-1)) {
		if (exists $opthash{$restart_procedure[$_]}) {
		push @restart, $restart_procedure[$_];
		push @restart, $_;
		}
	}
	die "\nERROR: Only one of the --restart_from_* option can be specified\n\n" if (@restart > 2); 
	
	# get specified --stop_after_*
	my @stop_procedure = qw/stop_after_rmdup stop_after_ubxandp stop_after_sga_assemble stop_after_exonerate_best stop_after_merge/;
	my @stop;
	for (0..(@stop_procedure-1)) {
		if (exists $opthash{$stop_procedure[$_]}) {
		push @stop, $stop_procedure[$_];
		push @stop, $_;
		}
	}
	die "\nERROR: Only one of the --stop_after_* option can be specified\n\n" if (@stop > 2); 
	
	# find the start and stop step
	my ($startstep, $stopstep);
	if ((@restart == 0) && (@stop == 0)) { # if no --restart_from_* or --stop_after_* is specified
	# so start from rmdup.pl and end at reblast.pl
	$startstep = "rmdup";
	$stopstep = "reblast";
	} elsif ((@restart > 0) && (@stop > 0)){ # if one --restart_from_* and one --stop_after_* are specified
		if ($stop[1]-$restart[1] >= 1) { # end step should at least the same as start step
		($startstep) = $restart[0] =~ /restart_from_(\S+)/;
		($stopstep) = $stop[0] =~ /stop_after_(\S+)/;
		} else { # warn if end step is before start step
		die "\nERROR: Please specify a correct restart and stop option\n\n";
		}
	} elsif ((@restart > 0) && (@stop == 0)) { # if one --restart_from_* and no --stop_after_* is specified
	# start from specified step and end at reblast
	($startstep) = $restart[0] =~ /restart_from_(\S+)/;
	$stopstep = "reblast";
	} elsif ((@restart == 0) && (@stop > 0)) { # if no --restart_from_* and one --stop_after_* is specified
	# start from rmdup and end at specified step
	$startstep = "rmdup";
	($stopstep) = $stop[0] =~ /stop_after_(\S+)/;
	}

# return name of start and end step
return ($startstep, $stopstep); 
}

# subroutine to generate sample list
sub samplelist {
my $startstep = shift; # start step

my %samplelist;
if (! $opthash{"samplelist"}) { # if --samplelist is not specified
	if ($startstep eq "rmdup") { # automatically generate sample list if wrapper is start from first step
		$opthash{"samplelist"} = $samplelist; # specify value of --samplelist as default value
		
		# get sample name from trimmed reads
		if ($opthash{"trimmed"}) {
			my $trimmeddir = $opthash{"trimmed"};
			opendir TRIMMED, $trimmeddir or die "\nERROR: Cannot find directory containing reads without adaptor and low quality bases \"$trimmeddir\", please check --trimmed ($!)\n\n";
			open SAMPLELIST, ">$samplelist";
			while (my $trimmed = readdir TRIMMED) {
			next if ($trimmed =~ /^\./);
				if ($trimmed =~ /(\S+)\_R1\.fq$/) { # get sample name from R1 reads
				my $samplename = $1;
				print SAMPLELIST "$samplename\n"; # print out sample name
				$samplelist{$samplename}++; # save sample name and record number of reads file in a sample
				} elsif ($trimmed =~ /(\S+)\_R2\.fq$/) { # get sample name from R1 reads
				my $samplename = $1;
				$samplelist{$samplename}++;
				} else {
				die "\nERROR: suffix of $trimmed is neither *_R1.fq nor *_R2.fq\n\n";
				}
			}
			closedir TRIMMED;
			close SAMPLELIST;
		} else {
		die "\nERROR: Please specify --trimmed if wrapper is start from first step (rmdup.pl)\n\n";
		}
		
		# warn if reads files are not paired
		my @single_sample = grep {$samplelist{$_} == 1} (sort keys %samplelist);
		if (@single_sample > 0) {
		die "\nERROR: @single_sample only got R1 or R2\n\n";
		}
		
		say STDOUT "All samples will be assembled\n";
	} else { # list of sample must be provided if wrapper is not start from first step
	die "\nERROR: Please provided list of sample if wrapper is not start from first step (rmdup.pl)\n\n";
	}
} else { # if --samplelist is specified
	$samplelist = $opthash{"samplelist"};
	my $line_cnt = 0;
	open SAMPLELIST, "$samplelist" or die "\nERROR: Cannot find sample list \"$samplelist\", please check --samplelist ($!)\n\n";
	while (my $samplename = <SAMPLELIST>) {
	$line_cnt++;
	chomp $samplename;
	my @element = $samplename =~ /\S+/g; # get non-space words in a line
		if (@element == 1) { # if there is only one word per line
		$samplelist{$samplename} = "";
		} elsif (@element > 1) { # warn if there are multiple words
		die "\nERROR: Multiple sample names or sample name with space in line $line_cnt of \"$samplelist\"\n\n";
		}
	}
	close SAMPLELIST;
	my @samplelist = sort keys %samplelist;
	say STDOUT "Samples including @samplelist will be assembled\n";
}

# return hash of sample name
return \%samplelist;
}

# subroutine to translate relative path into absolute path
sub absolute_path {
my $path = shift; # relative path

# split path to prefix and suffix. For example xxx/xxx/aaa, xxx/xxx is prefix and aaa is suffix
my @str = $path =~ /([^\/]+)/g; # get all non-"/" word
# get suffix
my $suffix = $str[-1];
# get prefix
$path =~ s/$suffix$//;
my $prefix = $path;

# if there is no prefix, provided path is at current directory
if (! $prefix) {
$prefix = "./";
}

# get absolute path
my ($absolute_path) = `cd $prefix;pwd` =~ /(\S+)/;
$absolute_path .= "/$suffix";

# return absolute path
return $absolute_path;
}

# subroutine to run script
sub run_script {
my $script = shift; # script name
my $option = shift; # option

# generate command line
my $command = "$interpreter ${script_path}$script"; # basic command without option
my $run = generate_command($option, $command); # add options to basic command

# start running
say STDOUT "########### Start Running $script ###########";
say STDOUT "Command: $run\n";
	my $error = system("$run"); # run with system command which can return error message
	if ($error) { # if error occur
	die "\nERROR: $run unexpectedly stopped \n\n";
	}
say STDOUT "########### $script Finished ###########";
}

# subroutine to generate command line with options
sub generate_command {
my $alloption = shift; # all options
my $command = shift; # basic command
	
	# add input option one by one
	foreach my $option (@$alloption) {
		if (exists $opthash{$option}) {
		$command .= " --$option $opthash{$option}";
		}
	}

# return entire command
return $command;
}

# subroutine to check existence and correctness of input file 
sub check_file_exists {
my $startstep = shift; # start step
my $stopstep = shift; # end step
my $samplelist_hash = shift; # hash of all input sample 

	# save mandatory options of each step into @alldepends and save order of step into %step
	my (@alldepends, %steps);
	for ($allprocedure{$startstep}..$allprocedure{$stopstep}) {
	my $step = $allprocedure[$_]; # get order of step
	push @alldepends, @{$depend{$step}}; # save mandatory options of each step 
	$steps{$step} = ""; # save order of step
	}
	
	# push --outdir option if wrapper run to last step
	if ($stopstep eq "reblast") {
	push @alldepends, "outdir";
	}

	# save required option into %alldepends
	my %alldepends = map {$_, ""} (@alldepends);
	
	# save input option into %input_option
	my %input_option;
	foreach my $option (@ARGVoption) {
		if ($option =~ /^--*(\S+)/) {
		$input_option{$1} = "";	
		}
	}
	
	# check whether option missed
	my @missing_option;
	foreach my $option (sort keys %alldepends) {
		if (!(exists $input_option{$option})) {
		push @missing_option, "--$option";
		}
	}
	
	# warn if options are missed
	if (@missing_option >= 1) {
	my $missing_option = join ", ", @missing_option;
	die "\nERROR: Missing option $missing_option\n\nUse -h or --help for more information on usage\n\n";
	}
	
	# prepare the sample list
	my @samplelist = sort keys %$samplelist_hash;
	
	# check dependencies
	foreach my $depend (sort keys %alldepends) {
		if ($depend eq "trimmed") { # check dir containing trimmed reads (--trimmed)
		my $trimmed_dir = $opthash{$depend};
			if (-e $trimmed_dir) { # if $trimmed_dir exists check whether samples in it are paired
				foreach my $sample (@samplelist) {
				my $R1 = "$trimmed_dir/$sample\_R1.fq";
				my $R2 = "$trimmed_dir/$sample\_R2.fq";
					if ((!(-e $R1))&&(!(-e $R2))) {
					die "\nERROR: Cannot find R1 of $sample \"$R1\" and \"$R2\", please check files under $trimmed_dir ($!)\n\n";
					} elsif ((!(-e $R1))&&(-e $R2)) {
					die "\nERROR: Cannot find R1 of $sample \"$R1\", reads file is not paired, please check file under $trimmed_dir ($!)\n\n";
					} elsif ((-e $R1)&&(!(-e $R2))) {
					die "\nERROR: Cannot find R2 of $sample \"$R2\", reads file is not paired, please check file under $trimmed_dir ($!)\n\n";
					}
				}
			} else {
			die "\nERROR: Cannot find directory containing reads without adaptor and low quality bases \"$trimmed_dir\", please check --trimmed ($!)\n\n";
			}
		} elsif ($depend eq "queryp") { # check amino acid and nucleotide sequences of query
			my $queryp = $opthash{$depend};
			
			# check whether it is a AA sequence if char "ATCG" is below 50%
			open QUERYP, $queryp or die "\nERROR: Cannot find amino acid sequences of target loci \"$queryp\", please check --queryp ($!)\n\n";
			<QUERYP>;
			chomp(my $seqp = <QUERYP>);
			my @nucleo = $seqp =~ /[A|T|C|G]/g;
			close QUERYP;

			if (@nucleo/length($seqp) >= 0.5) {
			die "\nERROR: Please provide amino acid sequences of target loci after --queryp\n\n";
			}
		} elsif ($depend eq "queryn") { # check nucleotide sequences of query
			my $queryn = $opthash{$depend};

			# check whether it is a DNA sequence if char "ATCG" is above 50%		
			open QUERYN, $queryn or die "\nERROR: Cannot find nucleotide sequences of target loci \"$queryn\", please check --queryn ($!)\n\n";
			<QUERYN>;
			chomp(my $seqn = <QUERYN>);
			my @nucleo = $seqn =~ /[A|T|C|G]/g;
			close QUERYN;

			if (@nucleo/length($seqn) < 0.5) {
			die "\nERROR: Please provide nucleotide sequences of target loci after --queryn\n\n";
			}
		} elsif ($depend eq "db") { # check dabatase
			my $udb = $opthash{$depend};
			my @udb = $udb =~ /(\S+)/g;
			
			# check dabatase format
			my $fascnt = 0;
			my $udbcnt = 0;
			my $elsecnt = 0;
			foreach my $db (@udb) {
				if (!(-e $db)) {
				die "\nERROR: Cannot find database of the species of target loci \"$db\", please check --db ($!)\n\n";
				}
				if ($db =~ /\S+\.fa$|\S+\.fas$|\S+\.fasta$/) {
				$fascnt++;
				} elsif ($db =~ /\S+\.udb$/) {
				$udbcnt++;
				} else {
				$elsecnt++;
				}
			}
			
			if ($elsecnt > 0) {
			die "\nERROR: Provided databases must be in either fasta or udb format\n\n";
			}
			if (($fascnt > 0)&&($udbcnt > 0)) {
			die "\nERROR: Provided databases must be in either fasta or udb format\n\n";
			}
		} elsif ($depend eq "dbtype") { # check dabatase type nucleo or prot
			my $dbtype = $opthash{$depend};
			if ($dbtype) {
				if (($dbtype ne "prot")&&($dbtype ne "nucleo")) {
				die "\nERROR: --dbtype must be specified as 'nucleo' or 'prot'\n\n";
				}
			} else {
			die "\nERROR: --dbtype must be specified as 'nucleo' or 'prot'\n\n";
			}
		} elsif ($depend eq "rmdup") { # check dir containing reads without PCR duplicates
			my $rmdup_dir = $opthash{$depend};
			
			# check existence of run_dir/rmdup
			if (-e $rmdup_dir) { 
				my %sample;
				opendir RMDUP, "$rmdup_dir";
				while (my $rmdup = readdir RMDUP) {
					if ($rmdup =~ /(\S+)\.\d+\.rmdup\.fq/) {
					$sample{$1} = "";
					}
				}
				closedir RMDUP;
				
				# check existence of all sample
				foreach my $sample (@samplelist) {
					if (!(exists $sample{$sample})) {
					die "\nERROR: Cannot find deduplicated reads of $sample under \"$rmdup_dir\" \n\n";
					}
				}
			} else {
			die "\nERROR: Cannot find the directory containing deduplicated reads \"$rmdup_dir\", please check --rmdup ($!)\n\n";
			}
		} elsif ($depend eq "parsed") { # check dir parsed reads
			my $parsed_dir = $opthash{$depend};
			
			# check existence of run_dir/parsed
			if (-e $parsed_dir) {
				foreach my $sample (@samplelist) {
					# check existence of each sample
					if (!(-e "$parsed_dir/$sample")) {
					die "\nERROR: Cannot find directory containing parsed reads of $sample \"$parsed_dir/$sample\" ($!)\n\n";
					}
				}
			} else {
			die "\nERROR: Cannot find directory containing parsed reads \"$parsed_dir\", please check --parsed ($!)\n\n";
			}	
		} elsif ($depend eq "assembled") { # check dir containing contigs and asqg files
			my $assembled_dir = $opthash{$depend};
			
			# check existence of run_dir/assembled
			if (-e $assembled_dir) {
				my $contig_dir = "$assembled_dir/contig"; # contig dir
				my $asqg_dir = "$assembled_dir/asqg"; # graph dir
				
				# check existence of run_dir/assembled/contig and run_dir/assembled/asqg
				if ((-e $contig_dir)&&(-e $asqg_dir)) {
					foreach my $sample (@samplelist) {
						if (!(-e "$assembled_dir/contig/$sample")) {
						die "\nERROR: Cannot find directory containing contigs of $sample \"$assembled_dir/contig/$sample\" ($!)\n\n";
						}
						if (!(-e "$assembled_dir/asqg/$sample")) {
						die "\nERROR: Cannot find directory containing string graphs of $sample \"$assembled_dir/asqg/$sample\" ($!)\n\n";
						}
					}
				} else {
					if (!(-e $contig_dir)) {
					die "\nERROR: Cannot find directory containing contigs \"$contig_dir\" ($!)\n\n";
					}
					if (!(-e $asqg_dir)) {
					die "\nERROR: Cannot find directory containing string graphs \"$asqg_dir\" ($!)\n\n";
					}
				}
			} else {
			die "\nERROR: Cannot find directory containing contigs and string graphs \"$assembled_dir\", please check --assembled ($!)\n\n";
			}	
		} elsif ($depend eq "filtered") { # check dir containing filtered contigs including 3 subfoder "filtered", "exonerate", "overlap"
			my $filtered_dir = $opthash{$depend};
			
			# check existence of run_dir/filtered
			if (-e $filtered_dir) {
				my $contig_filtered_dir = "$filtered_dir/filtered";
				my $exonerate_dir = "$filtered_dir/exonerate";
				my $overlap_dir = "$filtered_dir/overlap";
				
				# check existence of run_dir/filtered/filtered, run_dir/filtered/exonerate and run_dir/filtered/overlap
				if ((-e $contig_filtered_dir)&&(-e $exonerate_dir)&&(-e $overlap_dir)) {
					foreach my $sample (@samplelist) {
						if (!(-e "$filtered_dir/filtered/$sample")) {
						die "\nERROR: Cannot find directory containing filtered contigs of $sample \"$filtered_dir/filtered/$sample\" ($!)\n\n";
						}
						if (!(-e "$filtered_dir/exonerate/$sample")) {
						die "\nERROR: Cannot find directory containing exonerate results of $sample \"$filtered_dir/exonerate/$sample\" ($!)\n\n";
						}
						if (!(-e "$filtered_dir/overlap/$sample")) {
						die "\nERROR: Cannot find directory containing overlapping information of $sample \"$filtered_dir/overlap/$sample\" ($!)\n\n";
						}
					}
				} else {
					if (!(-e $contig_filtered_dir)) {
					die "\nERROR: Cannot find directory containing filtered contigs \"$contig_filtered_dir\" ($!)\n\n";
					}
					if (!(-e $exonerate_dir)) {
					die "\nERROR: Cannot find directory exonerate results \"$exonerate_dir\" ($!)\n\n";
					}
					if (!(-e $overlap_dir)) {
					die "\nERROR: Cannot find directory containing overlapping information \"$overlap_dir\" ($!)\n\n";
					}
				}
			} else {
			die "\nERROR: Cannot find directory containing filtered contigs, exonerate results and overlapping information \"$filtered_dir\", please check --filtered ($!)\n\n";
			}	
		} elsif ($depend eq "merged") { # check dir containing filtered contigs including 3 subfoder "nf", "f", "p"
			my $merged_dir = $opthash{$depend};
			
			# check existence of run_dir/merged
			if (-e $merged_dir) {
				# check existence of run_dir/merged/nf, run_dir/merged/f, run_dir/merged/p
				if (!(-e "$merged_dir/nf")) {
				die "\nERROR: Cannot find directory containing full-coding nucleotide sequences \"$merged_dir/nf\" ($!)\n\n";
				}
				if (!(-e "$merged_dir/f")) {
				die "\nERROR: Cannot find directory containing whole sequences with flanking region \"$merged_dir/f\" ($!)\n\n";
				}
				if (!(-e "$merged_dir/p")) {
				die "\nERROR: Cannot find directory containing amino acid sequences \"$merged_dir/p\" ($!)\n\n";
				}
			} else {
			die "\nERROR: Cannot find directory containing further assembled sequences \"$merged_dir\", please check --merged ($!)\n\n";
			}
		} elsif ($depend eq "ref_name") {
			# check whether --ref_name is specified
			my $ref_name = $opthash{$depend};
			if (! $ref_name) {
			die "\nERROR: Please specify the species name of target loci (--ref_name)\n\n";
			}
		} elsif ($depend eq "outdir") { # check whether --outdir is specified or there's a previous existed folder
			# check whether --outdir is specified
			my $outdir = $opthash{$depend};
			
			# if it is not specified
			if (!($outdir)) {
			die "\nERROR: Please specify the name of output directory (--outdir)\n\n";
			}
			# if it have existed
			if (-e $outdir) {
			my ($outdir_path) = `cd $outdir; pwd` =~ /(\S+)/;
			die "\nERROR: Output directory \"$outdir_path\" has exists. Please specify another path\n\n";
			}
			# check specified path is correct
			mkdir $outdir or die "\nERROR: Cannot create output directory \"$outdir\". Please check --outdir\n\n";
			`rm -rf $outdir`;
		}
	}
	
	# if both exists then check whether genes in them are the same
	if ((exists $alldepends{"queryp"})&&(exists $alldepends{"queryn"})) { 
		my %query;
	
		my $queryp = $opthash{"queryp"};
		my $queryn = $opthash{"queryn"};
		
		# check occurence of loci in queryp and queryn
		open QUERYP, $queryp or die "\nERROR: Cannot find amino acid sequences of target loci \"$queryp\", please check --queryp ($!)\n\n";
		while (my $line = <QUERYP>) {
			if ($line =~ /^>(\S+)/) {
			my $header = $1;
			$query{$header}++;
			}
		}
		close QUERYP;
		open QUERYN, $queryn or die "\nERROR: Cannot find nucleotide sequences of target loci \"$queryn\", please check --queryn ($!)\n\n";
		while (my $line = <QUERYN>) {
			if ($line =~ /^>(\S+)/) {
			my $header = $1;
			$query{$header}++;
			}
		}
		close QUERYN;
		
		# warn if occurence is smaller than 2, which means some loci is not co-exist in both files
		my @single_gene = grep {$query{$_} < 2} sort keys %query;
		if (@single_gene > 0) {
		die "\nERROR: Number or name of the genes in \"$queryp\" and \"$queryn\" are not the same\n\n";
		}
	}
}

# subroutine to make dir containing AA sequences of reference
sub makequerypdir {
my $queryp = shift; # AA sequences of reference

# check whether previous exist dir got the same number of sequences as provided reference
my $flag = 0; # if $flag is bigger than 0, which means dir of AA sequences of reference has been exist 
if (-e $querypdir) {
my ($file_num) = `ls $querypdir|wc -l` =~ /(\d+)/; # count file number
my ($gene_num) = `grep ">" $queryp|wc -l` =~ /(\d+)/; # count loci number
	if ($file_num != $gene_num) { # delete origin dir if file and loci number do not consist
	`rm -rf $querypdir`;
	} else {
	$flag++;
	}
} 

# make new dir containing amino acid sequences of quer
if ($flag == 0) {
	mkdir $querypdir;
	open QUERY, "$queryp";
	while (my $line = <QUERY>) {
		if ($line =~ />(\S+)/) {
		chomp(my $seq = <QUERY>);
		open QUERYPOUT, ">$querypdir/$1.fas";
		print QUERYPOUT ">$1\n$seq\n";
		close QUERYPOUT;
		}
	}
	close QUERY;
}
}

# subroutine to check whether query exists in provided database
sub check_query {
my $queryn = $opthash{"queryn"}; # DNA sequences of reference
my $db = $opthash{"db"}; # database
my $mask_udb = $opthash{"mask_db"}; # mask db or not
my $dbtype = $opthash{"dbtype"}; # db type nucleo or prot
my $E_value_reblast = $opthash{"E_value_reblast"}; # e-value in reblast.pl
my $cpu = $opthash{"cpu"}; # number of process used
	
	# check existence of DNA sequences of reference, database and whether db type is specified
	if (!(-e $queryn)) {
	die "\nERROR: Cannot find nucleotide sequences of reference \"$queryn\", please check --queryn ($!)\n\n";
	}
	if (!(-e $db)) {
	die "\nERROR: Cannot find database of the species of reference \"$db\", please check --db ($!)\n\n";
	}
	if (!($dbtype)) {
	die "\nERROR: Please specify --dbtype as \"prot\" or \"nucleo\"\n\n";
	}
	
	# check query by running reblast.pl with --check_query option
	my $basic_command = "$interpreter ${script_path}reblast.pl --check_query --queryn $queryn --db \"$db\" --dbtype $dbtype --E_value_reblast $E_value_reblast --cpu $cpu";	
	
	# mask db if --mask_db is specified
	my $command;
	if ($mask_udb) {
	$command = $basic_command . " --mask_db";
	} else {
	$command = $basic_command;
	}
	system("$command");
	
exit;
}

# subroutine to check dependencies of wrapper
sub check_depends {	

# check interpreter of current process
my @args = `ps -o args= -p $$` =~ /(\S+)/g; # get interpreter of current process by ps command
# get argument before $0, which is script name
my @interpreter;
foreach my $args (@args) {
	if ($args ne $0) {
	push @interpreter, $args;
	} else {
	last;
	}
}
my $current_interpreter = join " ", @interpreter;

say STDOUT "Currently used interpreter is \"$current_interpreter\"\n";

# check version of interpreter
my $version_info = `$current_interpreter -v` or die "\nERROR: Please specify a correct perl interpreter\n\n";
my ($first_line) = $version_info =~ /(.+)/;
my ($version) = $first_line =~ /v(\d+\.\d+)\.\d+/;
if ($version < 5.18) {
die "\nERROR: Version of your interpreter ($current_interpreter) is v$version. Please upgrate to 5.18 or higher verion\n\n";
} else {
say STDOUT "Version of your perl interpreter ($current_interpreter) is v$version\n";
}
	
# check whether perl modules are properly installed in interpreter
my ($Bioseq, $ParallelForkManager, $SysInfo);
eval {
$Bioseq = system("$current_interpreter -e 'use Bio::Seq'");
$ParallelForkManager = system("$current_interpreter -e 'use Parallel::ForkManager'");
$SysInfo = system("$current_interpreter -e 'use Sys::Info'");
};

if ($Bioseq) {
die "\nERROR: Bio::Seq (in Bioperl) is not properly installed\n\n";
} elsif ($ParallelForkManager) {
die "\nERROR: Parallel::ForkManager is not properly installed\n\n";
} elsif ($SysInfo) {
die "\nERROR: Sys::Info is not properly installed\n\n";
} else {
say STDOUT "All modules are properly installed\n";
}
	
# check whether software are properly installed
`usearch` or die "\nERROR: Cannot find usearch in \$PATH ($!)\n\n";
`sga` or die "\nERROR: Cannot find sga in \$PATH ($!)\n\n";
`exonerate` or die "\nERROR: Cannot find exonerate in \$PATH ($!)\n\n";
say STDOUT "All softwares are properly installed\n";
	
# check whether wrapped scripts are placed in correct place
my $script_path = $opthash{"script_path"};

# if --script_path is not specified, we assume they are under $PATH
if (! $script_path) {
	my $rmdup = `which rmdup.pl`;
	my $ubxandp = `which ubxandp.pl`;
	my $sga_assemble = `which sga_assemble.pl`;
	my $exonerate_best = `which exonerate_best.pl`;
	my $merge = `which merge.pl`;
	my $reblast = `which reblast.pl`;
	
	if (! $rmdup) {
	die "\nERROR: Cannot find rmdup.pl in \$PATH. You can specify a interpreter and path to script to run the wrapper\n\n";
	}
	if (! $ubxandp) {
	die "\nERROR: Cannot find ubxandp.pl in \$PATH. You can specify a interpreter and path to script to run the wrapper\n\n";
	}
	if (! $sga_assemble) {
	die "\nERROR: Cannot find sga_assemble.pl in \$PATH. You can specify a interpreter and path to script to run the wrapper\n\n";
	}
	if (! $exonerate_best) {
	die "\nERROR: Cannot find exonerate_best.pl in \$PATH. You can specify a interpreter and path to script to run the wrapper\n\n";
	}
	if (! $merge) {
	die "\nERROR: Cannot find merge.pl in \$PATH. You can specify a interpreter and path to script to run the wrapper\n\n";
	}
	if (! $reblast) {
	die "\nERROR: Cannot find reblast.pl in \$PATH. You can specify a interpreter and path to script to run the wrapper\n\n";
	}
	say STDOUT "All scripts are found in \$PATH.\n";
} else { # if --script_path is specified, find whether all of them are under the specified path
	if (!(-e $script_path)) {
	die "\nERROR: Cannot find path to called scripts \"$script_path\", please check --script_path ($!)\n\n";
	}
	
	my $provided_script_path = $script_path;
	($script_path) = `cd $script_path; pwd` =~ /(\S+)/;
	$script_path .= "/";
	if (!(-e "${script_path}rmdup.pl")) {
	die "\nERROR: Cannot find rmdup.pl under \"$provided_script_path\", please specify a correct path ($!)\n\n";
	}
	if (!(-e "${script_path}ubxandp.pl")) {
	die "\nERROR: Cannot find ubxandp.pl under \"$provided_script_path\", please specify a correct path ($!)\n\n";
	}
	if (!(-e "${script_path}sga_assemble.pl")) {
	die "\nERROR: Cannot find sga_assemble.pl under \"$provided_script_path\", please specify a correct path ($!)\n\n";
	}
	if (!(-e "${script_path}exonerate_best.pl")) {
	die "\nERROR: Cannot find exonerate_best.pl under \"$provided_script_path\", please specify a correct path ($!)\n\n";
	}
	if (!(-e "${script_path}merge.pl")) {
	die "\nERROR: Cannot find merge.pl under \"$provided_script_path\", please specify a correct path ($!)\n\n";
	}
	if (!(-e "${script_path}reblast.pl")) {
	die "\nERROR: Cannot find reblast.pl under \"$provided_script_path\", please specify a correct path ($!)\n\n";
	}
	say STDOUT "All scripts are found under \"$provided_script_path\"\n";
}

exit;
}

# subroutine to get number of cpu
sub get_cpu_num {
use Sys::Info;
my $info = Sys::Info->new;
my $cpu = $info->device('CPU');
my @info = $cpu->identify =~ /(\S+)/g;
my $cpunum = $info[0];

say STDOUT "Number of avaliable CPU is $cpunum\n";

return $cpunum;
}

# subroutine to calculate run time
sub runtime {
my $duration = shift;

	my ($hour, $min, $sec) = (0, 0, 0);
	if ($duration < 60) {
	$sec = $duration;
	} elsif (($duration >= 60) && ($duration < 3600)) {
	$min = int($duration / 60);
	$sec = $duration % 60;
	} elsif ($duration >= 3600) {
	$hour = int($duration / 3600);
	my $timeleft = $duration % 3600;
		if ($timeleft < 60) {
		$sec = $timeleft;
		} else {
		$min = int($timeleft / 60);
		$sec = $timeleft % 60;
		}
	}

return "${hour}h ${min}m ${sec}s";
}

# subroutine to remove intermediate output
sub clean {
if (-e $run_dir) {
say STDOUT "$run_dir and all intermediate output under it will be removed";
`find $run_dir -delete`;
} else {
die "\nERROR: Cannot find $run_dir containing intermediate output ($!)\n\n";
}

exit;
}