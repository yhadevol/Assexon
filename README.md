# Assexon
Assexon: Assembling Exon Using Gene Capture Data

### INTRODUCTION

This pipeline is to recover targeted exons and their flanking sequences from exon capture data

Pipeline can be divided into 3 steps:

(1) Data preparation

(2) Assembling

(3) Further processing

This documentation includes a step-by-step tutorial using a small test dataset. It will
help you to familiarize with this pipeline.

### SYSTEM REQUIREMENTS

(1) Data preparation and assembling:

Softwares: (Please put them under `$PATH`)

* Perl v5.18 or higher
* trim_galore v0.4.1 or higher
* cutadapt v1.2.1 or higher
* USEARCH 10.0.240 or higher
* SGA v0.10.15 or higher
* Exonerate v2.2.0 or higher

Perl module:

* Bio::Seq (Included in Bioperl)
* Parallel::Forkmanager
* Sys::Info

(2) Further processing:

This step is optional, so system requirements for this step are not listed here. Please check
requirements for these scripts by `-h` or `--help` options

**NOTE:** All softwares need to be installed and can be found in `$PATH`.

# INSTALLATION
Clone the Assexon repostory from GitHub:

	$ git clone https://github.com/yhadevol/Assexon.git

Scripts are placed under "pipeline_scripts", and it has three subfolders:

	pipeline_scripts
		└── data_preparation
		└── assemble
		└── postprocess

Open `~/.bashrc` then add paths to those folders to `$PATH`, so that all scripts can be easily called:

	export PATH=$PATH:/path/to/pipeline_scripts/data_preparation:/path/to/pipeline_scripts/assemble:/path/to/pipeline_scripts/postprocess

We have successfully installed Assexon and all dependencies on Linux (CentOS 6).

# DOCUMENTATION
A manual page containing a description of all options can be accessed by option `-h` or
`--help`. For example, full options of `assemble.pl` can be accessed by:

	$ assemble.pl -h

or

	$ assemble.pl --help

### SCRIPTS FOR EACH STEP

(1) Data preparation (pipeline_scripts/data_preparation):

*	gunzip_Files.pl (Expand gunzipped raw data)
*	trim_adaptor.pl (Trim low quality bases and adaptors)
*	predict_frames.pl (predict frame for reference, then generate coding and AA sequences of reference)

(2) Assembling (pipeline_scripts/assemble):

Main script:

* assemble.pl (Wrapper around several scripts to recover orthologous exons and their flanking sequences from target enrichment data)

Called scripts:

*	rmdup.pl (Remove PCR duplication)
*	ubxandp.pl (Parse reads to homologous loci)
*	sga_assemble.pl (De novo assemble parsed reads)
*	exonerate_best.pl (Filter unqualified assemblies and find assemblies can be further assembled)
*	merge.pl (Assemble contigs further and retrieve best contigs for each locus)
*	reblast.pl (Remove potential paralogs)

(3) Further processing (pipeline_scripts/postprocess):

Manipulate dataset:
*	pick_taxa.pl (Pick out needed taxa or discard unneeded taxa)
*	merge_loci.pl (Merge sequences under several directories from the same loci)
*	get_orthologues.pl (Find sequences orthology to reference from existing genomes)

Align:
*	mafft_aln.pl (Align nucleotide sequences in codon or normally)

Filter:
*	filter.pl (Remove badly aligned sequences)
*	flank_filter.pl (Discard too variable flanking regions)
*	clocklikeness_test.pl (Pick out loci which follows molecular clock hypothesis)
*	monophyly_test.pl (Pick out loci which topology is not congurence with provided monophyly group)

Statistics:
*	statistics.pl (Summary statistics for each locus and sample)
*	map_statistics.pl (Summary statistics for each sample from duplication-marked bam file)
*	count_reads_bases.pl (Count number of base pairs and reads for fastq file)

SNP-based analysis
*	consensus.pl (Make majority consensus sequences)
*	gatk.sh (Wrapper to call SNPs by GATK)
*	vcftosnps.pl (Convert output from GATK into other format)

Phylogenetic analysis
*	concat_loci.pl (Concatenate all loci into a master gene)
*	construct_tree.pl (Construct constrained or not constrained ML trees in batch)

Others:
*	unixlb_unwarp.pl (Substitute line breaks and unwrap sequences of fasta file)

**NOTE:** Only part of scripts will be demonstrated in following tutorial. Please use "-h" or "--help" to
see the detailed usage and options for each script

# TUTORIAL

The purpose of this tutorial is to help familiarize you with the format of the input you need and output you should expect from running this pipeline. The tutorial uses a test dataset in `test_data` that is a subset of real Illumina data from an enriched library.

This tutorial assumes that you have some experience executing programs from the command line on a UNIX-like system. If you would like to learn more about the command line, or need a refresher, you can find a command line tutorial in `introToCmdLine.pdf`.


### TEST DATA

Change directory to test data:

	$ cd test_data

Under `test_data` there are:
* raw_reads: A folder containing 2 gzipped raw reads. Structure of folder looks like:

		raw_reads
		    └── test1
		    |     ├── test1_R1.fq.gz
		    |     └── test1_R2.fq.gz
		    └── test2
		          ├── test2_R1.fq.gz
		          └── test2_R2.fq.gz		

* Oreochromis_niloticus.fas: DNA sequences of reference in fasta format.

* Oreochromis_niloticus.pep.fas: Reference protein sequences mined from Ensembl

* Oreochromis_niloticus.onehitCDSmarkers.column1.txt: First coloumn of OnehitCDSmarker
generated from Evolmarker, which should be generated during baits designing

* Oreochromis_niloticus.genome.fas: Soft-masked genomic DNA sequences in fasta format.
(It's not a real genome, but we treat it as a genome in this tutorial for convenience.
In real case, well soft-masked genomic sequences can be downloaded from Ensembl).

* species1.genome.fas and species2.genome.fas: genome sequences of other species

##### Let's start to run Assexon now
### DATA PREPARATION

#### Step 1: Gunzip data
First, let's expand gzipped reads under `raw_reads`:

	$ gunzip_Files.pl \
	--gzip raw_reads \
	--gunzipped gunzipped_raw_reads

Involved options:

* --gzip: Directory containing gzipped raw data
* --gunzipped: Directory containing expanded raw data

Output:
* gunzipped_raw_reads: Directory containing expanded raw data

#### Step 2: Trim adaptor and low quality bases
Then, trim illumina adaptor and low quality bases.

	$ trim_adaptor.pl \
	--raw_reads gunzipped_raw_reads \
	--trimmed trimmed

Output:
* trimmed: Directory containing reads without adaptor and low quality bases
* trimmed_reads_bases_count.txt: file summarized number of reads and bases in raw and trimmed reads
* trimming_report: Directory containg trimming report for each sample

Involved options:

* --raw_reads: Directory containing raw data
* --trimmed: Output directory containing reads without adaptor and low quality bases

#### Step 3: Correct frame and prepare coding and AA sequences of reference

First, we need to correct frame and prepare coding and AA sequences of reference in fasta format:

  	$ predict_frames.pl \
  	--baits Oreochromis_niloticus.fas \
  	--cds Oreochromis_niloticus.onehitCDSmarkers.column1.txt \
  	--proteome Oreochromis_niloticus.pep.fas

Involved options:

*	--baits: Frame-uncorrected sequences of reference
*	--cds: OnehitCDSmarker generated from Evolmarker. Only info in first columns will be used,so just input file with only first column is fine
*	--proteome: Reference protein sequences mined from ENSMBL

Output:
* Oreochromis_niloticus.dna.fas: Coding DNA sequences of reference
* Oreochromis_niloticus.aa.fas: Amino acid sequences of reference
* frame_result.txt: Frame prediction result for each locus

### ASSEMBLE

All inputs for assembly has been prepared. Let's start assembling now. The main script
is assemble.pl. This script calls another 6 scripts to recover assemblies.

6 scripts represent 6 steps of assembly. They are called by main script in following procedure:
* rmdup.pl: Remove PCR duplicates
* ubxandp.pl: Parse reads to targeted loci
* sga_assemble.pl: Assemble reads for each locus
* exonerate_best.pl: Filter unqualified contigs and find contigs which might be furtherly assembled  
* merge.pl: Assemble contigs further and retrieve best contigs for each locus
* reblast.pl: Remove potential paralogs

#### Run whole pipeline

Normally, we run the whole assembling pipeline (input cleaned reads, output orthologous assemblies), which includes 3 steps:

##### Step 1: Check requirements of assembling
Before running the script, we need to check requirements which can be checked by `--check_depends`.

	$ assemble.pl --check_depends

Involved options:

*	--check_depends: Check all dependencies for assemble.pl

If all dependencies are properly installed, you will see the following text in STDOUT:

	Currently used interpreter is "/XXX/perl"

	Version of your perl interpreter (/XXX/perl) is v5.xx

	All modules are properly installed

	All softwares are properly installed

	All scripts are found under $PATH

##### Step 2: Check the existence of reference sequences in given genome:
Determination of orthology between reference and enriched sequence is based on whether they can be aligned to same position on the genome of reference species. Thus, existence of reference in genome must be verified first to avoid false negative detection resulting from missing targeted loci:

	$ assemble.pl \
	--check_query \
	--queryn Oreochromis_niloticus.dna.fas \
	--db Oreochromis_niloticus.genome.fas \
	--dbtype nucleo

Involved options:

*	--check_query: Check whether reference sequences existing in given database, and return list of missing loci, then exit
*	--queryn: DNA sequences of reference in fasta format
*	--db: Path to DNA or AA database, either in Fasta or UDB format
*	--dbtype: Database type either 'nucleo' for DNA or 'prot' for AA database

If input genome is in Fasta format, a corresponding UDB database will be generated. If input database is in udb format, then no file will be generated. You will see the following text in STDOUT in both case:

	Start constructing ublast database in parallel
	Something generated by usearch...
	Ublast database has been constructed
	Something generated by usearch...
	#### All genes are found in provided database ####

If some genes do not exist in given genome, STDOUT will be:

	#### 2 genes below are not found in provided database  ####
	Danio_rerio.1.46410167.46410317
	Danio_rerio.14.21871634.21871111

##### Step 3: Assemble:
Requirements and existence of target loci in given genome have been checked. Let's start assemble:

	$ assemble.pl \
	--trimmed trimmed \
	--queryp Oreochromis_niloticus.aa.fas \
	--queryn Oreochromis_niloticus.dna.fas \
	--db Oreochromis_niloticus.genome.fas \
	--dbtype nucleo \
	--ref_name Oreochromis_niloticus \
	--outdir assemble_result

Involved options:

*	--trimmed: Directory containing reads without adaptor and low quality bases
*	--queryp: Amino acid sequences of target loci in fasta format
*	--queryn: Nucleotide sequences of target loci in fasta format
*	--db: Path to DNA or amino acid database, either in fasta or udb format
*	--dbtype: Database type either 'nucleo' for DNA or 'prot' for amino acid database
*	--ref_name: Substitute name of target loci as --ref_name in the output of last step (reblast.pl), disabled in default
*	--outdir: Directory to pipeline output

Several folders and files will be generated during the execution:
* run_dir: All intermediate outputs will be generated under this folder.
* samplelist.txt: A list includes name of all samples
* rmdup_reads_bases_count.txt: A table records number of reads and bases before and after removing PCR duplicates
* enriched_loci.txt: A table records number of total loci, number of enriched loci and percentage of enriched loci for each sample
* Oreochromis_niloticus.genome.fas.udb: UDB of `Oreochromis_niloticus.genome.fas`. This can be used as input database.

Output will be placed under `assemble_result` including 3 folders:
* nf: folder containing coding nucleotide sequences
* f: folder containing coding sequences with flanking regions
* p: folder containing AA sequences

#### Clean intermediate output

Intermediate output under "run_dir" would occupy lot of memory. Remove "run_dir" and all files under it by:

	$ assemble.pl --clean


#### Run partial pipeline

If something goes wrong at the intermediate step, don't worry, assemble.pl is able to restart from intermediate step. It can also stop at the step you want.

To restart from intermediate step, 4 things are essentially needed:
* Intermediate output from previous one step
* Essential options for the following step
* Options for restart or stop
* sample list

1) Intermediate output from each step:
*	Step 1: rmdup.pl: `./run_dir/rmdup`
*	Step 2: ubxandp.pl: `./run_dir/parsed`
*	Step 3: sga_assemble.pl: `./run_dir/assembled`
*	Step 4: exonerate_best.pl: `./run_dir/filtered`
*	Step 5: merge.pl: `./run_dir/merged`
*	Step 6: reblast.pl: `./run_dir/reblastout`

2) Essential options for each step:
*	Step 1: rmdup.pl: `--trimmed`
*	Step 2: ubxandp.pl: `--queryp`
*	Step 3: sga_assemble.pl: nothing
*	Step 4: exonerate_best.pl: `--queryp`
*	Step 5: merge.pl: `--queryp` and `--queryn`
*	Step 6: reblast.pl: `--db`, `--dbtype`, `--ref_name` and `--outdir`

3) Option for restart or stop at a step
* To restart from a step: `--restart_from_xxx`
* To stop at a step: `--stop_after_xxx`

For example, I want restart from step 4 (`exonerate_best.pl`). The option is:

`--restart_from_exonerate_best`

I want stop at step 5 (`merge.pl`). The option is:

`	--stop_after_merge`

I want restart from step 4 and stop at step 5, then specify 2 options:

`--restart_from_exonerate_best` and `--stop_after_merge`

4) Sample list
Sample list is named as `samplelist.txt` in default. It contains the list of sample names, one sample name per line. It is automatically generated from first step. It looks like:

	test1
	test2

Here are 3 examples of running partial pipeline:

**EXAMPLE 1**: From an intermediate step to the end (exonerate_best.pl -> end):
* Output from previous one step `sga_assemble.pl` (`./run_dir/assembled`)
* Essential inputs of `exonerate_best.pl` (`--queryp`), `merge.pl` (`--queryn`, `--queryp`) and
`reblast.pl` (`--db`, `--dbtype`, `--ref_name`, `--outdir`)
* samplelist.txt
* option `--restart_from_exonerate_best`

So the command is:

	$ assemble.pl \
	--queryp Oreochromis_niloticus.aa.fas \
	--queryn Oreochromis_niloticus.dna.fas \
	--db Oreochromis_niloticus.genome.fas \
	--dbtype nucleo \
	--outdir assemble_result \
	--ref_name Oreochromis_niloticus \
	--samplelist samplelist.txt \
	--restart_from_exonerate_best

**EXAMPLE 2**: Restart from an intermediate step to another intermediate step (sga_assemble.pl -> merge.pl):

* output from previous one step `ubxandp.pl` (`./run_dir/parsed`)
* Essential inputs of `sga_assemble.pl` (nothing), exonerate_best.pl (`--queryp`) and `merge.pl`
(`--queryn`, `--queryp`)
* samplelist.txt
* 2 options `--restart_from_sga_assemble` as well as `--stop_after_merge`

Command:

	$ assemble.pl \
	--queryp Oreochromis_niloticus.aa.fas \
	--queryn Oreochromis_niloticus.dna.fas \
	--samplelist samplelist.txt \
	--restart_from_sga_assemble \
	--stop_after_merge

**EXAMPLE 3**: Stop at an intermediate step (start -> sga_assemble.pl)

We start from the first step, so there's no input from previous one step. We just need:

* Essential inputs of `rmdup.pl` (`--trimmed`), `ubxandp.pl` (`--queryp`), `sga_assemble.pl` (nothing)
* samplelist.txt
* option `--stop_after_sga_assemble`

Command:

	$ assemble.pl \
	--trimmed trimmed \
	--queryp Oreochromis_niloticus.aa.fas \
	--samplelist samplelist.txt \
	--stop_after_sga_assemble

#### Only assemble part of samples

Samples exist in `samplelist.txt` will be assembled. You can only write name of samples which need to be assembled. For example, I want to assemble test1 only, then the list is:

	test1

Then, specify the option `--samplelist` to input your list:

	$ assemble.pl \
	--trimmed trimmed \
	--queryp Oreochromis_niloticus.aa.fas \
	--queryn Oreochromis_niloticus.dna.fas \
	--db Oreochromis_niloticus.genome.fas \
	--dbtype nucleo \
	--outdir assemble_result \
	--ref_name Oreochromis_niloticus \
	--samplelist samplelist.txt

### FURTHER PROCESSING

Before downstream analysis, datasets probably need to be modified. Sequences are required to be aligned, and poorly aligned sequences should be discarded. We also need to access statistics of filtered alignments. So further processing mainly includes:
* Manipulate dataset
* Aligning
* Filtering
* Summary statistics

#### Manipulate dataset (optional):

Before aligning and filtering, Some users may want to add orthologue sequences from existing
genomes or delete poorly enriched sequences. Before introducing how to do it, we emphasize
once again that:

`SEQUENCES MUST BE ADDED OR DELETED BEFORE ALIGNING!!!`

#### Add orthologous sequences

Dependencies:
* USEARCH v10.0.240 or higher
* BioPerl v1.007001 or higher

##### Step 1: Extract orthologous sequences from existing genomes
First, we extract sequences orthology to loci in `Oreochromis_niloticus.dna.fas` from `species1.genome.fas` and `species2.genome.fas`

	$ get_orthologues.pl \
	--query Oreochromis_niloticus.dna.fas \
	--querydb Oreochromis_niloticus.genome.fas \
	--subdb "species1.genome.fas|species2.genome.fas" \
	--subname "species1 species2" \
	--outdir orthologs \
	--cpu 12

Output:
* orthologs: Directory includes sequences orthology to reference
* Oreochromis_niloticus.genome.fas.udb: udb of `Oreochromis_niloticus.genome.fas`.
* species1.genome.fas.udb: udb of `species1.genome.fas`
* species2.genome.fas.udb: udb of `species2.genome.fas`

Involved options:

*	--query: Coding DNA sequences of reference
*	--querydb: Space delimited list of one or more DNA databases of reference species in either in FASTA or UDB format
*	--subdb: List of DNA databases of subjects in FASTA or UDB format ONLY, but database format of the same subject need to consistent. Input database list of different subjects are delimited by '|', and database belonging to the same subject are delimited by space. e.g. \"sp1.genome.1.fas sp1.genome.2.fas|sp2.genome.1.udb sp2.genome.2.udb\"
*	--subname: Space delimited list of subject name in output, which is one-to-one match to the list of subject databases.
*	--outdir: Name of output directory, which has 2 subfolders including 'nf' for coding sequences and 'p' for AA sequences.
*	--cpu: Limit the number of CPUs, 1 in default

##### Step 2: Add them into datasets
Then, we add coding sequences from `species1.genome.fas` and `species2.genome.fas` to enriched datasets. Each locus should have at least 3 sequences:

	$ merge_loci.pl \
	--indir "assemble_result/nf orthologs/nf" \
	--outdir merged_nf \
	--min_seq 3

Output:
* merged_nf: Dir containing merged loci files

Involved options:

* --indir: List of dir containing sequences
* --outdir: Dir containing merged loci files
*	--min_seq: Minimum sequences required in merged file, 2 in default

#### Delete unneeded sequences

Dependencies:
* Nothing

Discard taxa `species2` by "--deselected_taxa":

	$ pick_taxa.pl \
	--indir merged_nf \
	--outdir merged_nf_deselected \
	--deselected_taxa "species2"

Output:
* merged_nf_deselected: Dir containing sequences of without discarded taxon

Involved options:

* --indir: Dir containing unaligned sequences
* --outdir: Dir containing sequences of selected taxon
*	--deselected_taxa: List of taxa want to be discarded, each taxon is delimited by space

#### Aligning

Then, let's align each loci.

Dependencies:
* BioPerl v1.007001 or higher
* Mafft v7.294b or higher (rename it as `mafft` and put it under `$PATH`)
* Parallel::ForkManager;

1) If input sequences are full-coding sequence:

	$ mafft_aln.pl \
	--dna_unaligned merged_nf \
	--dna_aligned merged_nf_aligned \
	--cpu 12

Output:
* merged_nf_aligned: Dir containing nucleotide sequences aligned in codon

2) If input sequences are not coding sequence or coding sequences with flanks:

	$ mafft_aln.pl \
	--dna_unaligned assemble_result/f \
	--dna_aligned f_aligned \
	--non_codon_aln \
	--cpu 12

Output:
* f_aligned: Dir containing aligned nucleotide sequences

Involved options:

*	--dna_unaligned: Dir containing unaligned DNA sequences
*	--dna_aligned: Dir containing aligned DNA sequences, named as "xx_aligned" if this option is not specified
*	--non_codon_aln: Do not align DNA sequences in codon. This option is turned off by default
*	--cpu: Limit the number of CPUs, 1 in default.

#### Filtering

##### Filter poorly aligned coding sequences:

Dependencies:
* BioPerl v1.007001 or higher
* Mafft v7.294b or higher (rename it as `mafft` and put it in `$PATH`)

Coding sequences are input. Only remove poorly aligned sequences, and `Oreochromis_niloticus` is reference:

	$ filter.pl \
	--indir merged_nf_aligned \
	--filtered merged_nf_filtered \
	--ref_taxa "Oreochromis_niloticus" \
	--cpu 12

Output:
* merged_nf_filtered: Dir containing filtered alignments which are aligned in codon

Involved options:

* --indir: Dir containing unfiltered alignments
*	--filtered: Dir containing filtered alignments
*	--ref_taxa: A space delimit list of reference taxa
* --cpu: Limit the number of CPUs, 1 in default

**NOTE:** Remember to check filtered alignments. Please tuning the parameters, if too much alignments are filtered. Use `-h` or `--help` to access detailed parameter

##### Filter poorly aligned flanking sequences:

Dependencies:
* BioPerl v1.007001 or higher
* Mafft v7.294b or higher (rename it as `mafft` and put it under `$PATH`)

Filter sequences with flanking regions in `f_aligned`, and only files and sequences co-exist in `f_aligned` and `nf_filtered` will be written to `f_filtered`. Sequence of `Oreochromis_niloticus` is reference. Run script with 4 process:

	$ flank_filter.pl \
	--flank f_aligned \
	--nonflank_filtered merged_nf_filtered \
	--flank_filtered f_filtered \
	--ref_taxa "Oreochromis_niloticus" \
	--cpu 4

Output:
* f_filtered: Dir containing well-aligned coding sequences with flanks

Involved options:

*	--flank: Dir containing unfiltered alignments with flanks
*	--nonflank_filtered: Dir containing well-aligned coding sequences
*	--flank_filtered: Dir containing well-aligned coding sequences with flanks
*	--ref_taxa: A space delimit list of reference taxa
*  --cpu: Limit the number of CPUs, 1 in default

#### Filter loci for different purposes:

1. Users can pick out loci which topologies are not congurence with provided monophyletic group using `monophyly_test.pl`. This analysis needs unconstrained ML tree and ML tree constrained by provided monophyletic group for each locus, which can be generated using `construct_tree.pl`.

**NOTE:** Format of input file after `--constrain` (option in `construct_tree.pl`) can be found in `pipeline_scripts/postprocess/monogroup.txt`

2. Some analysis need genes follow the molecular clock hypotheses (like construction of time-recalibrated tree). Users can filter loci which disobey the molecular clock hypotheses using `clocklikeness_test.pl`

#### Summary statistics

Finally, let's summary statistics of filter alignment.

Dependencies:
* Bio::AlignIO (included in Bioperl)
* Bio::Align::DNAStatistics (included in Bioperl)
* Parallel::ForkManager

Summarized statistics for coding region of each locus including:
* Number of enriched samples
* Alignment length
* GC content
* Percentage of Missing data
* Average pairwise distance among sequences

Summarized statistics for flanking region of each locus including:
* Alignment length of left flanking region
* Alignment length of right flanking region
* Average pairwise distance in left flanking region (only calculated for flanking region >= 20 bp)
* Average pairwise distance in right flanking region (only calculated for flanking region >= 20 bp)

Summarized statistcis for each sample including:
* Number of enriched loci
* GC content
* Average length of flanking region (left + right)

Command:

	$ statistics.pl \
	--coding_aligned merged_nf_filtered \
	--flanking_aligned f_aligned

Output:
* coding_summary.txt: Tab delimited table of summarized statistics for coding region of each locus
* flanking_summary.txt: Tab delimited table of summarized statistics for flanking region of each locus
* sample_summary.txt: Tab delimited table of summarized statistics for each sample

Involved options:
* --coding_aligned: Folder comprising aligned coding sequences
* --flanking_aligned: Folder comprising aligned coding sequences with flanking region

#### Phylogenetic analysis

Several scripts were provided to reformat filtered alignments as input of phylogenetic analysis. Alignment files can be rearranged and input for analysis in following procedure:

Concatenated tree:

	                 Filter alignments
	                         ↓
	  Codeoncatenate loci into master gene (concat_loci.pl)
	                         ↓
	            Build concatenated tree (RAxML)

Species tree:

	                 Filter alignments
	                         ↓
	    Construct trees for each loci (construct_tree.pl)
	                         ↓
	 Merge resulting gene trees into one file ("cat" command)
	                         ↓
	             Build species tree (ASTRAL)

SNP-based analysis:

	      Trimmed reads            Filter alignments
	            |                           ↓
	            |              Generate majority consensus
	            |               reference (consensus.pl)
	            |                           |
	            ————————————————————————————
	                         ↓
	           Call SNPs from reads (gatk.sh)
	                         ↓
	          Reformat vcf file (vcftosnps.pl)
	                         ↓
	            BEAST, STRUCTURE, dudi.pca
