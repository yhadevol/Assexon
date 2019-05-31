#!/bin/bash

###############################################
#
# gatk.sh
#
# This shell is for automatic gatk pipeline
#
# Usage: ./gatk.sh consensus_ref num_of_thread sub1 sub2 sub3...
#
# Written by Hao Yuan
#
# 2017.1.1 at SHOU
#
###############################################

#1 index reference by bwa
bwa index -a is ${1}.fa
echo "step 1 processed"

#2 make index file using samtools
samtools faidx ${1}.fa
echo "step 2 processed"

#3 make reference dictionary using picard
java -jar CreateSequenceDictionary.jar OUTPUT=${1}.dict R=${1}.fa
echo "step 3 processed"

#4 mapping
for i in *_R1.fq;
do bwa mem -t $2 -R "@RG\tID:${i%_R1.fq}\tSM:${i%_R1.fq}\tPL:ILLUMINA" ${1}.fa ${i%_R1.fq}_R1.fq ${i%_R1.fq}_R2.fq > ${i%_R1.fq}.sam;
done
echo "step 4 processed"

#5 SAM to BAM conversion using picard
for i in *.sam;
do java -Xmx4g -Djava.io.tmpdir=/tmp \
-jar SortSam.jar \
SO=coordinate \
INPUT=${i%.sam}.sam \
OUTPUT=${i%.sam}.bam \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true;
done
echo "step 5 processed"

#6 Marking PCR duplicates using picard
for i in *.bam;
do java -Xmx4g -Djava.io.tmpdir=/tmp \
-jar MarkDuplicates.jar \
INPUT=${i%.bam}.bam \
OUTPUT=${i%.bam}.marked.bam \
METRICS_FILE=metrics \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT;
done
echo "step 6 processed"

#7 first make realign list
for i in *.marked.bam;
do java -Xmx4g -jar GenomeAnalysisTK.jar \
-nt $2 \
-T RealignerTargetCreator \
-R ${1}.fa \
-I ${i%.marked.bam}.marked.bam \
-o ${i%.marked.bam}.bam.list;
done
echo "step 7 processed"

#8 realign
for i in *.marked.bam;
do java -Xmx4g -Djava.io.tmpdir=/tmp \
-jar GenomeAnalysisTK.jar \
-T IndelRealigner \
-R ${1}.fa \
-I ${i%.marked.bam}.marked.bam \
-targetIntervals ${i%.marked.bam}.bam.list \
-o ${i%.marked.bam}.marked.realigned.bam;
done
echo "step 8 processed"

#9 fix pair end using picard
for i in *marked.realigned.bam;
do java -Djava.io.tmpdir=/tmp/flx-auswerter \
-jar FixMateInformation.jar \
INPUT=${i%.marked.realigned.bam}.marked.realigned.bam \
OUTPUT=${i%.marked.realigned.bam}.marked.realigned.fixed.bam \
SO=coordinate \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true;
done
echo "step 9 processed"

#10 Variant-only calling on DNAseq using HaplotypeCaller
for i in *marked.realigned.fixed.bam;
do java  -Xmx4g -Djava.io.tmpdir=/tmp/flx-auswerter \
-jar GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R ${1}.fa \
-I ${i%.marked.realigned.fixed.bam}.marked.realigned.fixed.bam \
-stand_call_conf 30 \
-stand_emit_conf 10 \
--minPruning 5 \
-o ${i%.marked.realigned.fixed.bam}.raw.snps.indels.vcf;
done
echo "step 10 processed"

#11 BaseRecalibrator
for i in *marked.realigned.fixed.bam;
do java -Xmx4g  -Djava.io.tmpdir=/tmp/flx-auswerter \
-jar GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-I ${i%.marked.realigned.fixed.bam}.marked.realigned.fixed.bam \
-R ${1}.fa \
-knownSites ${i%.marked.realigned.fixed.bam}.raw.snps.indels.vcf \
-o ${i%.marked.realigned.fixed.bam}.recal_data.table;
done
echo "step 11 processed"

#12 Creating a recalibrated BAM（用printreads对recal_data.table进行校正生成bam file）#
for i in *marked.realigned.fixed.bam;
do java -Xmx4g  -Djava.io.tmpdir=/tmp/flx-auswerter \
-jar GenomeAnalysisTK.jar \
-T PrintReads \
-R ${1}.fa \
-I ${i%.marked.realigned.fixed.bam}.marked.realigned.fixed.bam \
-BQSR ${i%.marked.realigned.fixed.bam}.recal_data.table \
-o ${i%.marked.realigned.fixed.bam}.marked.realigned.fixed.recal.bam
done
echo "step 12 processed"

#13 Single-sample all-sites calling on DNAseq using HaplotypeCaller (for GVCF-based cohort analysis workflow)＃
for i in *marked.realigned.fixed.bam;
do java -Djava.io.tmpdir=/tmp/flx-auswerter \
 -jar GenomeAnalysisTK.jar \
 -T HaplotypeCaller \
 -R ${1}.fa \
 -I ${i%.marked.realigned.fixed.bam}.marked.realigned.fixed.recal.bam \
 --emitRefConfidence GVCF \
 --variant_index_type LINEAR \
 --variant_index_parameter 128000 \
 -o ${i%.marked.realigned.fixed.bam}.raw.snps.indels.g.vcf;
done
echo "step 13 processed"

#14 GenotypeGVCFs merges gVCF records
var=""
for i in "${@:3}"
do
j="--variant $i.raw.snps.indels.g.vcf "
var="${var}${j} "
done

java -Xmx4g -Djava.io.tmpdir=/tmp/flx-auswerter \
-jar GenomeAnalysisTK.jar \
-nt $2 \
-R ${1}.fa \
-T GenotypeGVCFs \
$var \
-o ${1}_snp.all.vcf
echo "step 14 processed"

#15 Extract the SNPs from the call set
java -jar GenomeAnalysisTK.jar \
-nt $2 \
-T SelectVariants \
-R ${1}.fa \
-V ${1}_snp.all.vcf \
-selectType SNP \
-o ${1}_snp.all.raw_snps.vcf
echo "step 15 processed"

# 16 Apply the filter to the SNP call set
java -jar GenomeAnalysisTK.jar \
-T VariantFiltration \
-R ${1}.fa \
-V ${1}_snp.all.raw_snps.vcf \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filterName "my_snp_filter" \
-o ${1}_snp.all.filtered_snps.vcf
echo "step 16 processed"

echo "all processed"

exit 0