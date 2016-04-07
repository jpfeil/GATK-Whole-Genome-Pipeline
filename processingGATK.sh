#!/bin/bash

# This script was created to run the GATK germline variant calling best pratices pipeline.
# The pipeline takes unprocessed bam file(s) as input and outputs the biological genomic variance
# of the input to a given reference

# Run the GATKsetup.sh script first. See the associated README

# stderr redirected to descriptively named files.report

set -e
set -x
set -o pipefail

# Set $RAM to use for all java processes
RAM=-Xmx200g
# Set number of threads to the number of cores/machine for speed optimization
THREADS=32
# Set the dir for the reference files
dir="reference"
# Set the reference fasta
#ref=human_g1k_v37.fasta
#phase2 reference
ref="hg19.fa"
mills="Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
phase1_indels="1000G_phase1.indels.hg19.sites.vcf"
dbsnp="dbsnp_138.hg19.vcf"
#choose how to log time
Time=/usr/bin/time

# Create Variable for input file
INPUT=$1
# Root input name
RINPUT=${INPUT%.*}

# sort sample
samtools sort -@ $THREADS $INPUT $RINPUT.sorted

# index bam file
samtools index $RINPUT.sorted

#START PIPELINE
# calculate flag stats using samtools
$Time samtools flagstat \
    $RINPUT.sorted.bam \
    > flagstat.report 2>&1

# sort reads in picard
$Time java $RAM \
    -jar ~/picard/picard.jar \
    SortSam \
    INPUT=$RINPUT.sorted.bam \
    OUTPUT=$RINPUT.picardSorted.bam \
    SORT_ORDER=coordinate \
    > sortReads.report 2>&1

# mark duplicates reads in picard
$Time java $RAM \
    -jar ~/picard/picard.jar \
    MarkDuplicates \
    INPUT=$RINPUT.picardSorted.bam \
    OUTPUT=$RINPUT.mkdups.bam \
    METRICS_FILE=metrics.txt \
    ASSUME_SORTED=true \
    > markDups.report 2>&1

rm -r $RINPUT.picardSorted.bam

# GATK Indel Realignment
# There are 2 steps to the realignment process
# 1. Determining (small) suspicious intervals which are likely in need of realignment (RealignerTargetCreator)
# 2. Running the realigner over those intervals (IndelRealigner)

# index reference genome, There's no need for this, we have the .fai file
# samtools faidx ${dir}/${ref}

# create sequence dictionary for reference genome
# We could skip this step. I've download the b37.dict

#$Time java $RAM \
#    -jar ~/picard/picard.jar \
#    CreateSequenceDictionary \
#    REFERENCE=${dir}/${ref} \
#    OUTPUT=${dir}/$i{ref%.*}.dict


# Index the markdups.bam for realignment
$Time samtools index $RINPUT.mkdups.bam

#Here Next
# create target indels (find areas likely in need of realignment)
$Time java $RAM \
    -jar ~/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -known ${dir}/$mills \
    -known ${dir}/$phase1_indels \
    -R ${dir}/${ref} \
    -I $RINPUT.mkdups.bam \
    -o ${dir}/target_intervals.list \
    -nt $THREADS \
    > targetIndels.report 2>&1

# realign reads
$Time java $RAM \
    -jar ~/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -known ${dir}/$mills \
    -known ${dir}/$phase1_indels \
    -R ${dir}/${ref} \
    -targetIntervals ${dir}/target_intervals.list \
    -I $RINPUT.mkdups.bam \
    -o $RINPUT.realigned.bam \
    > realignIndels.report 2>&1

rm -r $RINPUT.mkdups.bam

# GATK BQSR

# create gatk recalibration table
$Time java $RAM \
    -jar ~/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R ${dir}/${ref} \
    -I $RINPUT.realigned.bam \
    -knownSites ${dir}/$dbsnp \
    -knownSites ${dir}/$mills \
    -knownSites ${dir}/$phase1_indels \
    -o ${dir}/recal_data.table \
    -nct $THREADS \
    > recalibrationTable.report 2>&1

# recalibrate reads
$Time java $RAM \
    -jar ~/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R ${dir}/${ref} \
    -I $RINPUT.realigned.bam \
    -BQSR ${dir}/recal_data.table \
    -o $RINPUT.bqsr.bam \
    -nct $THREADS \
    > bqsr.report 2>&1

echo "This is the body" | mail -s "Preprocessing is done for $INPUT" $2
