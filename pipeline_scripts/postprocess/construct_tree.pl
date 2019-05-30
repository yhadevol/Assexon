#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long; 
use Parallel::ForkManager;

# input dir of aligned loci, input dir of constrained ML tree of each loci, input dir of ML tree of each loci, temp dir
my ($dir, $constrain, $mldir, $constraindir, $copy_dir, $help);
my $cpu = 1;

my $mlrundir = "mlrun"; # run dir for unconstrained ml tree
my $constrainrundir = "constrain_run"; # run dir for constrained ml tree
my $constraintreedir = "constrain_tree"; # dir for constrained tree for each gene
my $constrain_tree = "constrain_tree.tre"; # constrained tree containing all taxa provided in --constrain

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'indir:s', \$dir,
					  'constrain:s', \$constrain,
					  'ml_dir:s', \$mldir,
					  'constrained_ml_dir', \$constraindir,
					  'cpu:i', \$cpu,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";
                      
# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/indir/;
check_option(\@essential, \@ARGVoptions);

# absolute path of dir of aligned loci
$dir = absolute_path($dir);

# set name of ml tree dir and transform into absolute path
if (! $mldir) {
$mldir = "${dir}_ml"; # output dir for unconstrained ml tree
}
$mldir = absolute_path($mldir);

# set name of constrained tree dir and transform into absolute path
if (! $constraindir) {
$constraindir = "${dir}_constrained"; # output dir for constrained ml tree
}
$constraindir = absolute_path($constraindir);

$copy_dir = "${dir}_copy";

# check existence of raxml
`raxmlHPC-SSE3 -h` or die "\nERROR: Cannot find raxmlHPC-SSE3 in \$PATH\n\n";

# check existence and correctness of pre-defined monophyletic group
if ($constrain) {
my $count = 0;
open MONO, $constrain or die "\nERROR: Cannot find file of constrained group \"$constrain\", please check --constrain ($!)\n\n";
while (my $line = <MONO>) {
my @line = $line =~ /(\S+)/g;
	if (@line == 1) {
	my $line_num = $count+1;
	die "\nERROR: One taxon cannot be defined as a group at line $line_num\n\n";
	}
$count++;
}
close MONO;
}

# grep all fas file under input dir
opendir DIR, $dir or die "\nERROR: Cannot find input directory \"$dir\", please check --indir ($!)\n\n";
my @fas = grep {$_ !~ /^\./} readdir DIR;
@fas = grep {$_ =~ /fasta$|fas$|fa$/i} @fas;
closedir DIR; 

print "######### Start Construct ML tree ########\n";

# create output and run dir for unconstrained ml tree
mkdir $mldir;
mkdir $mlrundir;

# create copy dir
mkdir $copy_dir;

chdir($mlrundir); # change to run dir

# split dataset into $split parts
my $split = $cpu;
my $splited_array = split_array(\@fas, $split);

# initiate multi-process
my $pm = Parallel::ForkManager->new(int($cpu/2));

DATA_LOOP: while (my $array = shift @$splited_array) {
$pm->start and next DATA_LOOP;
	
	LOOP: while (my $fas = shift @$array) {
	
		my $path = "$dir/$fas";
		my ($gene) = $fas =~ /(\S+)\.fa\S*$/;
		
		# check input file is wrongly formatted or have any strange characters 
		my ($seqhash, $input_order, $aln_len) = readfas($path);
		
		# check whether input file is nucl or prot
		my $test_seq = $seqhash->{$input_order->[0]};
		my @nucleo = $test_seq =~ /[A|T|C|G]/g;

		if (@nucleo/length($test_seq) < 0.5) {
		print STDERR "WARNING: Input sequences must be nucleotide sequences, while amino acid sequence found in $fas\n";
		next LOOP;
		} 
		
		# number of taxa
		my $taxa_num = @$input_order;
		
		# skip if too few taxa in alignment
		if ($taxa_num < 4) {
		print STDERR "WARNING: At least 4 taxa are required. There are only $taxa_num in $fas\n";
		next LOOP;
		}
		
		# skip if sequences are not aligned
		if ($aln_len eq "na") {
		print STDERR "\nERROR: Number of characters among sequences are not the same in $path. It may be not aligned. Skip this file\n";
		next LOOP;
		}
		
		# copy alignment file to $copy_dir
		my $copy_path = "$copy_dir/$fas";
		open COPY_FAS, ">$copy_path";
		foreach my $taxon (@$input_order) {
		print COPY_FAS ">$taxon\n$seqhash->{$taxon}\n";
		}
		close COPY_FAS;
		
		# construct tree
		`raxmlHPC-SSE3 -p 12345 -m GTRGAMMA -s $copy_path -n $gene.tre`; # run raxml
		`mv RAxML_bestTree.$gene.tre $mldir/$gene.tre`; # mv tree to outdir
		unlink "RAxML_info.$gene.tre";
		unlink "RAxML_log.$gene.tre";
		unlink "RAxML_parsimonyTree.$gene.tre";
		unlink "RAxML_result.$gene.tre";
	}

$pm->finish;
}
$pm->wait_all_children;

chdir "../";

`rm -rf $mlrundir`; # remove run dir for unconstrained ml tree

print "######### ML tree construction Finished ########\n";

`find $copy_dir -name "*reduced" -delete`; # remove reduced dataset(no whole gap col)

if ($constrain) { # build constrained tree if constrained file is provided
print "######### Start Constrained ML tree ########\n";
	
	# hash of monogroup{taxa} = line number where it locate
	my (%monogroup, @constrain_tree);
	my $count = 0;
	open MONO, $constrain;
	while (my $line = <MONO>) {
	my @line = $line =~ /(\S+)/g;
		foreach (@line) {
		$monogroup{$_} = $count;
		}
	$count++;	
	
	# join each taxon together
	my $joint = join "," , @line;
	$joint = "($joint)";
	push @constrain_tree, $joint;
	}
	close MONO;
	
	# join all constrained group
	my $joint_constrain_tree = join ",", @constrain_tree;
	$joint_constrain_tree = do {
		if (@constrain_tree > 1) {
		"($joint_constrain_tree);"; 
		} else {
		"$joint_constrain_tree;"; 
		}
	};

	# write out constrained tree containing all taxa provided in --constrain
	open CONSTRAINTREE, ">$constrain_tree";
	print CONSTRAINTREE "$joint_constrain_tree\n";
	close CONSTRAINTREE;
	
	# grep all fas file under input dir
	opendir DIR, $copy_dir or die "\nERROR: Cannot find copied directory of --indir \"$copy_dir\" ($!)\n\n";
	my @copy_fas = grep {$_ !~ /^\./} readdir DIR;
	@copy_fas = grep {$_ =~ /fasta$|fas$|fa$/i} @copy_fas;
	closedir DIR; 
	
	$splited_array = split_array(\@copy_fas, $split);
	
	# make output, tree and run dir
	mkdir $constraindir;
	mkdir $constraintreedir;
	mkdir $constrainrundir;
	
	# change to run dir
	chdir($constrainrundir);
	
	DATA_LOOP: while (my $array = shift @$splited_array) {
	$pm->start and next DATA_LOOP;
	
		LOOP: while (my $fas = shift @$array) {
		my $path = "$copy_dir/$fas";
		my ($gene) = $fas =~ /(\S+)\.fa\S*$/;
		
			# record each taxon
			my @taxa;
			open FAS, $path;
			<FAS>;
			while (my $line = <FAS>) {
				if ($line =~ />(\S+)/) {
				my $taxa = $1;
				push @taxa, $taxa;
				}
			}
			close FAS;
		
			# extract sub-group from whole group generated from --constrained, and treat other taxon as single group
			my %genemono;
			my $key_num = scalar keys %monogroup;
			foreach my $taxa (@taxa) {
				if (exists $monogroup{$taxa}) { # if taxa exist in provided group
				push @{$genemono{$monogroup{$taxa}}}, $taxa; # push them together
				} else { # if taxa is not exist in provided group
				push @{$genemono{$key_num}}, $taxa; # each of them push to a single element
				$key_num++;
				}
			}
		
			# get number of all taxa in current gene
			my $alltaxa = 0;
			foreach my $cat_num (sort keys %genemono) {
			my $taxanum = scalar @{$genemono{$cat_num}};
			$alltaxa += $taxanum if ($taxanum == 1);
			}
		
			# get number of group
			$key_num = scalar keys %genemono;
			if ($alltaxa == $key_num) { # if number of all taxa == number of group means there's no constrain group
			`cp $mldir/$gene.tre $constraindir/$gene.tre`;
			next LOOP;
			}
		
			# generate constrained tree
			# get each group
			my @gene_constrain;
			foreach my $monogroup (sort keys %genemono) {
				if (@{$genemono{$monogroup}} > 1) {
				my $joint_mono = join ",", @{$genemono{$monogroup}};
				$joint_mono = "($joint_mono)";
				push @gene_constrain, $joint_mono;
				} else {
				push @gene_constrain, $genemono{$monogroup}->[0];
				}
			}
			# join all group together
			my $gene_constrain = join "," , @gene_constrain;
			$gene_constrain = do {
			if (@gene_constrain > 1) {
				"($gene_constrain);";
				} else {
				"$gene_constrain;";
				}
			};
		
			# write out sub constrained tree which only containing taxa in current gene
			my $constrain_tree_path = "../$constraintreedir/$gene.constrain.tre";
			open CONSTRAIN_GENETREE, ">$constrain_tree_path";
			print CONSTRAIN_GENETREE "$gene_constrain\n";
			close CONSTRAIN_GENETREE;
	
			# run raxml and mv result to output dir
			`raxmlHPC-SSE3 -p 12345 -m GTRGAMMA -s $path -g $constrain_tree_path -n $gene.tre`;
			`mv RAxML_bestTree.$gene.tre $constraindir/$gene.tre`;
			unlink "RAxML_info.$gene.tre";
			unlink "RAxML_log.$gene.tre";
			unlink "RAxML_parsimonyTree.$gene.tre";
			unlink "RAxML_result.$gene.tre";
			unlink $constrain_tree_path;
		}
		
	$pm->finish;
	}
	$pm->wait_all_children;
	
	chdir("../");

	# remove run dir, tree dir and all constrained tree
	`rm -rf $constrainrundir`;
	`rm -rf $constraintreedir`;
	unlink $constrain_tree; 

print "######### Constrained ML tree Construction Finished ########\n";
}

`rm -rf $copy_dir`; # remove reduced dataset(no whole gap col)

######################################################
# Subroutines
######################################################

# subroutine to translate relative path into absolute path
sub absolute_path {
my $path = shift; # relative path

my @str = $path =~ /([^\/]+)/g;
my $suffix = $str[-1];
$path =~ s/$suffix$//;
my $prefix = $path;

if (! $prefix) {
$prefix = "./";
}
my ($absolute_path) = `cd $prefix;pwd` =~ /(\S+)/;
$absolute_path .= "/$suffix";

# return absolute path
return $absolute_path;
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

# subroutine to reformat input fasta file
sub readfas {
my $file = shift;

open INFILE, $file or die "\nERROR: Cannot find input file \"$file\" ($!)\n\n";

my (%seq, %seqlength, @input_order, $seq, $taxon, $lasttaxon, %taxa_cnt);

while (my $line = <INFILE>) {

	# remove enter
	$line =~ s/\r//g;
	chomp $line;
	
	# if find >
	if ($line =~ /^>(\S+)/) {
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
$taxon_num++;
print STDERR "WARNING: Found $taxon_num \"$taxon\" in $file. Only first sequence named in \"$taxon\" will be kept\n";
}

my @seqlength = sort keys %seqlength;

my $aln_len;
if (@seqlength == 1) {
$aln_len = $seqlength[0];
} else {
$aln_len = "na";
}

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

# subroutine to print usage
sub usage {
print STDERR "
Script name: construct_tree.pl

This is a script to construct ML trees for each locus in fasta format. Topologies of trees can also be constrained based on user-definition.

Dependencies: 
(1) raxmlHPC-SSE3 (download from 'https://github.com/stamatak/standard-RAxML')
(2) Perl module:
	1. Parallel::ForkManager

Example usage:
(1) Construct ML trees for files under 'nf_aligned'. Write ML trees to 'nf_aligned_ml'. Run scripts with 4 process:

	perl construct_tree.pl --indir nf_aligned --cpu 4

(2) Construct ML trees and constrained trees for files under 'nf_aligned'. Write ML trees to 'nf_aligned_ml'. Write constrained ML trees to 'nf_aligned_constrained'. Run scripts with 4 process:

	perl construct_tree.pl --indir nf_aligned --constrain monogroup.txt --cpu 4

Input files:
(1) nf_aligned
(2) monogroup.txt (if --constrain is specified)

Output files:
(1) nf_aligned_ml
(2) nf_aligned_constrained (if --constrain is specified)

Options:
--indir
  Directory containing aligned nucleotide sequences
--constrain
  Text file with given group, one group per line, example please refer to monogroup.txt
--ml_dir
  Directory containing unconstrained ML trees, named as xxx_ml in default, xxx is the name of input directory
--constrained_ml_dir
  Directory containing constrained ML trees, named as xxx_constrained in default, xxx is the name of input directory
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