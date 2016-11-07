### Creating fastqs

These fastqs were created by subsetting data from 1000 genomes high coverage bams (HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam and HG00268.wgs.ILLUMINA.bwa.FIN.high_cov_pcr_free.20140203.bam) available from here: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/

The following commands use tools that I installed with bioconda:

```
conda install samtools bbmap seqtk
```


To extract reads mapped to MT, I used the commands:

```
samtools view -b HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam MT | samtools bam2fq -1 ind1_reads_1.fq -2 ind1_reads_2.fq -

samtools view -b HG00268.wgs.ILLUMINA.bwa.FIN.high_cov_pcr_free.20140203.bam MT | samtools bam2fq -1 ind2_reads_1.fq -2 ind2_reads_2.fq -

```

These reads need to be sorted and read pairing needs to be restored, so I used repair.sh available from bbmap:

```
repair.sh in1=ind1_reads_1.fq in2=ind1_reads_2.fq out1=ind1_reads_1.fastq out2=ind1_reads_2.fastq

repair.sh in1=ind2_reads_1.fq in2=ind2_reads_2.fq out1=ind2_reads_1.fastq out2=ind2_reads_2.fastq

gzip *.fastq
```

Finally, I subset a total of 10,000 paired reads from each sample with the commands:

```
seqtk sample -s100 ind1_reads_1.fastq.gz 5000 | gzip > ind1_1.fastq.gz

seqtk sample -s100 ind1_reads_2.fastq.gz 5000 | gzip > ind1_2.fastq.gz

seqtk sample -s100 ind2_reads_1.fastq.gz 5000 | gzip > ind2_1.fastq.gz

seqtk sample -s100 ind2_reads_2.fastq.gz 5000 | gzip > ind2_2.fastq.gz
```
