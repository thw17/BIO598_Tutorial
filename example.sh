#!/usr/bin/env bash

# Prepare reference
samtools faidx reference/human_g1k_v37_MT.fasta
picard CreateSequenceDictionary R=reference/human_g1k_v37_MT.fasta o=reference/human_g1k_v37_MT.dict
bwa index reference/human_g1k_v37_MT.fasta

# Process sample 1 (ind1)
bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samblaster -M | samtools fixmate - - | samtools sort -O bam -o bam/ind1.rmdup.sorted.bam -
samtools index bam/ind1.rmdup.sorted.bam
samtools stats bam/ind1.rmdup.sorted.bam | grep ^SN | cut -f 2- > stats/ind1.rmdup.sorted.bam.stats

# Process sample 2 (ind2)
bwa mem -M -R '@RG\tID:ind2\tSM:ind2\tLB:ind2\tPU:ind2\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind2_1.fastq.gz fastq/ind2_2.fastq.gz | samblaster -M | samtools fixmate - - | samtools sort -O bam -o bam/ind2.rmdup.sorted.bam -
samtools index bam/ind2.rmdup.sorted.bam
samtools stats bam/ind2.rmdup.sorted.bam | grep ^SN | cut -f 2- > stats/ind2.rmdup.sorted.bam.stats

# Jointly call variants for both samples
freebayes -f reference/human_g1k_v37_MT.fasta bam/ind1.rmdup.sorted.bam bam/ind2.rmdup.sorted.bam > vcf/joint.raw.vcf
