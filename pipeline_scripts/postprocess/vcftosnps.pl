#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;   # include the module for input

my ($infilename, $outfile, $selected_taxa, $help);
my $SNPs_pct_trsd = 80;
my $qual_trsd = 0;

my @ARGVoptions = @ARGV;

my $opt = GetOptions( 'vcf:s', \$infilename,
					  'outfile:s', \$outfile,
					  'selected_taxa:s', \$selected_taxa,
					  'SNPs_pct:f', \$SNPs_pct_trsd,
					  'qual:f', \$qual_trsd,
                      'help|h!', \$help) or die "\nERROR: Found unknown options. Use -h or --help for more information on usage\n\n";

# print help if no options or -h is specified
if ((@ARGVoptions == 0)||($help)) {
usage();
}

# check missing options
my @essential = qw/vcf/;
check_option(\@essential, \@ARGVoptions);

my ($nexusout, $structout, $beastout, $rout);
if (! $outfile) {
my ($file) = $infilename =~ /(\S+)\.vcf$/;
$nexusout = $file . ".nex";
$structout = $file . "_struct.txt";
$beastout = $file . "_beast.nex";
$rout = $file ."_rout.txt";
} else {
$nexusout = $outfile . ".nex";
$structout = $outfile . "_struct.txt";
$beastout = $outfile . "_beast.nex";
$rout = $outfile ."_rout.txt";
} 

my ($IN_FILE, $NEX_OUT, $STRUCT_OUT, $BEAST_OUT, $R_OUT);
open ($IN_FILE, $infilename) or die "\nERROR: Cannot find input vcf file \"$infilename\", please check --vcf ($!)\n\n";

#process the headline
my $line_cnt = 0;
my @head;
while (my $line = readline($IN_FILE)) {
	$line_cnt++;
    if ($line =~ /^#CHROM/){ #if find the head line
        chomp $line;
        @head = split/\s/, $line;
        last;
    }
}


my (@selected_taxa, %selected_taxa);
if ($selected_taxa) {
	# num of minimum selected taxa
	@selected_taxa = $selected_taxa =~ /\S+/g;
	die "\nERROR: Taxa list must have at least 2 taxa, please check --selected_taxa \n\n" if (@selected_taxa < 2);

	# check whether input taxa is right
	my @taxalist = @head;
	foreach (1..9) {
	shift @taxalist;
	}
	my %taxalist = map {$_ , ""} @taxalist;

	foreach my $selected_taxa (@selected_taxa) {
	die "\nERROR: Input taxa \"$selected_taxa\" do not exist, please check --selected_taxa \n\n" if (! exists $taxalist{$selected_taxa});
	}

	# get name of selected taxa
	%selected_taxa = map {$_ , ""} @selected_taxa;
}

#get the selected taxa name and order
my (@selected_taxalist, @selected_taxa_order);
my $arraylen = scalar(@head);
for (my $i = 9; $i < $arraylen; $i ++) {
    if ((exists $selected_taxa{$head[$i]})||(! $selected_taxa)) { # add taxon into list if it exist in --selected_taxa. add all taxa if --selected_taxa is not specified
    push @selected_taxa_order, $i;
    push @selected_taxalist, $head[$i];
    }
}

my $total_taxa_num = scalar @selected_taxa_order;

#process the snps
my (%snps, %gene_snps, $idlag);
impossible:while (my $line = readline($IN_FILE)) {
	$line_cnt++;
	
    chomp $line;
    my @tmp = split/\s/, $line;
    
    my @ref = $tmp[3] =~ /[^\,]+/g; # nucleotides for reference
    my @alt = $tmp[4] =~ /[^\,]+/g; # nucleotides for alternative SNPs
    
    # skip indel
    foreach (@ref) {
		if (length($_) > 1) {
			print STDERR "WARNING: Find indel in line $line_cnt\n";
			next;
		}
    }
    
    foreach (@alt) {
    	if (length($_) > 1) {
			print STDERR "WARNING: Find indel in line $line_cnt\n";
			next;
		}
    }
    
    foreach my $i (@selected_taxa_order) {
    	# get snps of both loci (./0/1)
        my ($gt1, $gt2) = $tmp[$i] =~ /^(\S)\S(\S)/; 
        # skip if found unrecognized char
        next impossible if ((($gt1 ne ".") && ($gt1 ne "0") && ($gt1 ne "1")) || (($gt2 ne ".") && ($gt2 ne "0") && ($gt2 ne "1")));  
    }
    
    if ($tmp[6] ne "PASS"){ # skip unqualified SNPs
        next;
    }
    
    # Skip SNP if quality is too low;
    my $qual = $tmp[5]; # SNPs quality
    next impossible if ($qual < $qual_trsd);
    
    my $id = $tmp[0]; # Gene where SNP lies
    
    # if both allele is ".", means missing data
    my $missingcnt = 0;
	foreach my $i (@selected_taxa_order) {
        my ($gt1, $gt2) = $tmp[$i] =~ /^(\S)\S(\S)/;
        $missingcnt ++ if (($gt1 eq ".") && ($gt2 eq "."));
    }
    
    my $missingcnt_pct = sprintf("%.2f", $missingcnt/$total_taxa_num*100);
    my $SNPs_pct = 100-$missingcnt_pct;
    next impossible if ($SNPs_pct < $SNPs_pct_trsd);
	
    # if all SNPs in previous gene has been read, randomly select a qualified gene
    if (($idlag)&&($id ne $idlag)) { 
    	my @cand_snps;
    	foreach my $entry (sort keys %gene_snps) { # find entry with snp
			my @tmp2 = split/\s/, $entry;
		
			my %gt_type;
			foreach my $i (@selected_taxa_order) {
				my ($gt1, $gt2) = $tmp2[$i] =~ /^(\S)\S(\S)/;
				
				if (($gt1 ne ".") && ($gt2 ne ".")) {
				$gt_type{$gt1} = "";
				$gt_type{$gt2} = "";
				}
			}
			
			my $type_num = scalar keys %gt_type;
			if ($type_num > 1) {
			push @cand_snps, $entry;
			}
    	}
    	
    	if (@cand_snps > 0) { # randomly select one
			my $gene_snps_num = scalar @cand_snps; 
			my $rand_select = int(rand($gene_snps_num));
			my $selected_entry = $cand_snps[$rand_select];
		
			my @tmp2 = split/\s/, $selected_entry;
			foreach my $i (@selected_taxa_order) {
				my ($gt1, $gt2) = $tmp2[$i] =~ /^(\S)\S(\S)/;
			
				$snps{$tmp2[0]}->{$head[$i]}->{1} = $gt1;
				$snps{$tmp2[0]}->{$head[$i]}->{2} = $gt2;
				
				if (($gt1 eq ".") || ($gt2 eq ".")) {
				$snps{$tmp2[0]}->{$head[$i]}->{3} = "-";
				} else {
				$snps{$tmp2[0]}->{$head[$i]}->{3} = $gt1+$gt2;
				}
            }
    	}
        
        %gene_snps = ();  
    }
	
	# save the name of current $gene
	$gene_snps{$line} = "";
    	
    $idlag = $id;
}

#need print the last entry
my @cand_snps;
foreach my $entry (sort keys %gene_snps) {
	my @tmp2 = split/\s/, $entry;

	my %gt_type;
	foreach my $i (@selected_taxa_order) {
		my ($gt1, $gt2) = $tmp2[$i] =~ /^(\S)\S(\S)/;
	
		if (($gt1 ne ".") && ($gt2 ne ".")) {
		$gt_type{$gt1} = "";
		$gt_type{$gt2} = "";
		}
	}

	my $type_num = scalar keys %gt_type;
	if ($type_num > 1) {
	push @cand_snps, $entry;
	}
}
	
if (@cand_snps > 0) {
	my $gene_snps_num = scalar @cand_snps; 
	my $rand_select = int(rand($gene_snps_num));
	my $selected_entry = $cand_snps[$rand_select];

	my @tmp2 = split/\s/, $selected_entry;
	foreach my $i (@selected_taxa_order) {
		my ($gt1, $gt2) = $tmp2[$i] =~ /^(\S)\S(\S)/;
	
		$snps{$tmp2[0]}->{$head[$i]}->{1} = $gt1;
		$snps{$tmp2[0]}->{$head[$i]}->{2} = $gt2;
		
		if (($gt1 eq ".") || ($gt2 eq ".")) {
		$snps{$tmp2[0]}->{$head[$i]}->{3} = "-";
		} else {
		$snps{$tmp2[0]}->{$head[$i]}->{3} = $gt1+$gt2;
		}
	}
}

close ($IN_FILE);

#generate output files
open ($NEX_OUT, ">$nexusout") or die "\nERROR: Cannot write output file \"$nexusout\", please check --outfile ($!)\n\n";
open ($STRUCT_OUT, ">$structout") or die "\nERROR: Cannot write output file \"$structout\", please check --outfile ($!)\n\n";
open ($BEAST_OUT, ">$beastout") or die "\nERROR: Cannot write output file \"$beastout\", please check --outfile ($!)\n\n";
open ($R_OUT, ">$rout") or die "\nERROR: Cannot write output file \"$rout\", please check --outfile ($!)\n\n";

# header of nexus output
my @loci = keys (%snps);
my @sortedloci = sort @loci;
print $NEX_OUT "#NEXUS\n";
print $NEX_OUT "Begin characters;\n";
print $BEAST_OUT "#NEXUS\n";
print $BEAST_OUT "Begin data;\n";

# allele num, char num, taxa num
my $ntax = 2 * scalar(@selected_taxa_order);
my $nchar = scalar(@loci);
my $stax = scalar(@selected_taxa_order);

# part of nexus output
print $NEX_OUT "        Dimensions ntax=$ntax nchar=$nchar;\n";
print $NEX_OUT "        Format symbols=\"012\" missing=?;\n";
print $NEX_OUT "        Matrix\n";
print $BEAST_OUT "        Dimensions ntax=$stax nchar=$nchar;\n";
print $BEAST_OUT "        Format datatype=integerdata symbols=\"012\" gap=-;\n";
print $BEAST_OUT "        Matrix\n";

foreach my $taxon (@selected_taxalist) {
	# count loci with missing data 
    my $count1 = 0;
    my $count2 = 0;
    
    # print taxa name of single line file
    print $BEAST_OUT "$taxon\t";
    print $R_OUT "$taxon";
	
	# print taxa name of allele 1
    print $NEX_OUT "$taxon" . "_1\t";
    print $STRUCT_OUT "$taxon";
    
    # print loci
    foreach my $loci (@sortedloci){
    	# print merged data (0, 1 or 2, missing data(-))
        print $BEAST_OUT "$snps{$loci}->{$taxon}->{3}";
        print $R_OUT "\t$snps{$loci}->{$taxon}->{3}";
        
        # print data of allele 1 (0, 1, missing data(?/-9))
        if ($snps{$loci}->{$taxon}->{1} eq "."){ # missing data
            print $NEX_OUT "?";
            print $STRUCT_OUT "\t-9";
            $count1 ++;
        }
        else{ # not missing data
            print $NEX_OUT "$snps{$loci}->{$taxon}->{1}";
            print $STRUCT_OUT "\t$snps{$loci}->{$taxon}->{1}";
        }
    }
    
    # print taxa name of allele 2
    print $NEX_OUT "\n$taxon" . "_2\t";
    print $STRUCT_OUT "\n$taxon";
    
    # print loci
    foreach my $loci (@sortedloci){
        if ($snps{$loci}->{$taxon}->{2} eq "."){ # missing data
            print $NEX_OUT "?";
            print $STRUCT_OUT "\t-9";
            $count2 ++;
        }
        else{ # not missing data
            print $NEX_OUT "$snps{$loci}->{$taxon}->{2}";
            print $STRUCT_OUT "\t$snps{$loci}->{$taxon}->{2}";     
        }
    }
    
    # last enter
    print $NEX_OUT "\n";
    print $STRUCT_OUT "\n";
    print $BEAST_OUT "\n";
    print $R_OUT "\n";
    
    # report number of loci with missing data
    print "Number of loci with missing data in $taxon is $count1 \n";

}

# print end of file
print $NEX_OUT ";\nEnd;\n";
print $BEAST_OUT "\t;\nEnd;\n";

# report the total number of loci
print "The total number of loci is $nchar \n";

close ($NEX_OUT);
close ($STRUCT_OUT);
close ($BEAST_OUT);
close ($R_OUT);

#####################################################
# Subroutines
#####################################################

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
Script name: vcftosnps.pl

Overview:
This script is to read vcf file (generated by GATK), and randomly select SNPs from each locus. Output 4 files for different kinds of downstream analysis:
(1) outfile.nex: nex file in which each sample has 2 lines, 1 allele per line (input for SNAPP analysis)
(2) outfile_beast.nex: nex file in which each sample has 1 line, 2 alleles combined into 1 line
(3) outfile_struct.txt: txt file in which each sample has 2 lines, 1 allele per line (input for STRUCTURE analysis)
(4) outfile_rout.txt: txt file in which each sample has 1 line, 2 alleles combined into 1 line (input for PCA analysis)

Example usage:
(1) Generate 4 kind of outputs from 'species.vcf':

	perl vcftosnps.pl --vcf species.vcf

Input files:
(1) species.vcf

Output files:
(1) species.nex
(2) species_beast.nex
(3) species_struct.txt
(4) species_rout.txt

Options:
--vcf
  Vcf file generated from GATK, which ONLY have SNPs and DO NOT have INDEL
--outfile
  Prefix of outfile. if --outfile is not specified, the prefix of outfile is the prefix of vcf file 
--selected_taxa
  Space delimited list of taxa need to be output. All taxa are output if this option is not specified
--SNPs_pct
  The minimum percentage of nucleotides required for a SNP, $SNPs_pct_trsd in default. Smaller value means more missing data in the SNPs are allowed
--qual
  The minimum quality required for for a SNP, $qual_trsd in default.
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