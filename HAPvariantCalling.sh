#! /bin/bash
# GATK Haplotype VARIANT CALLING

set -e
set -x
set -o pipefail

cd ~

# Set $RAM to use for all java processes
RAM=-Xmx200g
# Set number of threads to the number of cores/machine for speed optimization
THREADS=32
# Set the dir for the reference files and input bam
dir="reference"
# Set the reference fasta
#phase2 reference
ref="hg19.fa"
mills="Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
phase1_indels="1000G_phase1.indels.hg19.sites.vcf"
dbsnp="dbsnp_138.hg19.vcf"
hapmap="hapmap_3.3.hg19.sites.vcf"
omni="1000G_omni2.5.hg19.sites.vcf"

# Create Variable for input file
INPUT=$1

# Root input name
RINPUT=${INPUT%.*}
#choose how to log time
Time=/usr/bin/time

# index bam file
samtools index $INPUT

# For running on multipe samples (a cohort) see:
#http://gatkforums.broadinstitute.org/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode

#snps & indels together
$Time java $RAM \
-jar ~/GenomeAnalysisTK.jar \
-nct $THREADS \
-R ${dir}/${ref} \
-T HaplotypeCaller \
--genotyping_mode DISCOVERY \
--output_mode EMIT_VARIANTS_ONLY \
-I $INPUT \
-o $RINPUT.unified.raw.BOTH.gatk.vcf \
-stand_emit_conf 10.0 \
-stand_call_conf 30.0 \
> HapVars.report 2>&1
	#-I $INPUT2.recal.bam \
	#-I $INPUT3.recal.bam \
	#-I $INPUT4.recal.bam \
	#-metrics snps.metrics \
	#-L /b37/target_intervals.bed \ intervals to operate over. See docstring

## GATK VARIANT QUALITY SCORE RECALIBRATION
#SNP
$Time java $RAM \
  -Djava.io.tmpdir=/tmp \
  -jar ~/GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R ${dir}/${ref} \
  -input $RINPUT.unified.raw.BOTH.gatk.vcf \
  -nt $THREADS \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${dir}/$hapmap \
  -resource:omni,known=false,training=true,truth=false,prior=12.0 ${dir}/$omni \
  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dir}/$dbsnp \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${dir}/$phase1_indels \
  -an QD \
  -an DP \
  -an FS \
  -an ReadPosRankSum \
  -mode SNP \
  -minNumBad 1000 \
  -recalFile $RINPUT.HAPSNP.recal \
  -tranchesFile $RINPUT.HAPSNP.tranches \
  -rscriptFile $RINPUT.HAPSNP.plots.R \
  > HAPVariantRecalibrator_HAPSNP.report 2>&1
  #-percentBad 0.01 \ depreciated

#Apply Snp Recalibration
$Time java $RAM \
  -jar ~/GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -input $RINPUT.unified.raw.BOTH.gatk.vcf \
  -o $RINPUT.HAPSNP.vqsr.SNP.vcf \
  -R ${dir}/${ref} \
  -nt $THREADS \
  -ts_filter_level 99.0 \
  -tranchesFile $RINPUT.HAPSNP.tranches \
  -recalFile $RINPUT.HAPSNP.recal \
  -mode SNP \
  > HAPApplyRecalibration_HAPSNP.report 2>&1

#Indel Recalibration
$Time java $RAM -Djava.io.tmpdir=/tmp \
  -jar ~/GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R ${dir}/${ref} \
  -input $RINPUT.unified.raw.BOTH.gatk.vcf \
  -nt $THREADS \
  -resource:mills,known=true,training=true,truth=true,prior=12.0 ${dir}/$mills \
  -an DP \
  -an FS \
  -an ReadPosRankSum \
  -mode INDEL \
  -minNumBad 1000 \
  -recalFile $RINPUT.HAPINDEL.recal \
  -tranchesFile $RINPUT.HAPINDEL.tranches \
  -rscriptFile $RINPUT.HAPINDEL.plots.R \
  > HAPVariantRecalibrator_HAPINDEL.report 2>&1
  #-percentBad 0.01 \depreciated

#Apply Indel Recalibration
$Time java $RAM \
  -jar ~/GenomeAnalysisTK.jar \
  -T ApplyRecalibration \
  -input $RINPUT.unified.raw.BOTH.gatk.vcf \
  -o $RINPUT.HAPvqsr_INDEL.vcf \
  -R ${dir}/${ref} \
  -nt $THREADS \
  -ts_filter_level 99.0 \
  -tranchesFile $RINPUT.HAPINDEL.tranches \
  -recalFile $RINPUT.HAPINDEL.recal \
  -mode INDEL \
  > HAPApplyRecalibration_HAPINDEL.report 2>&1
