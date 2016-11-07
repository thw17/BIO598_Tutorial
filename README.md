# BIO598_Tutorial
A hands-on tutorial introducing users to reproducible reference-based genome assembly and variant calling.

This tutorial has been tested on Mac and Linux operating systems and will assume you're working with one of these operating systems.

## Setting up your environment
For today's tutorial, you'll need this repository and Anaconda.

#### Getting the repo
To get this repository (repo), move to the directory on your computer or cluster where you'd like to work and type the command:

```
git clone https://github.com/thw17/BIO598_Tutorial
```
This should create a directory called BIO598_Tutorial containing all of the contents of this repo.

#### Setting up Anaconda
Anaconda is an environment and package manager for Python and it makes installation, environment management, etc. simple without requiring root or administrator access.  Fortunately, its framework has been leveraged for a project called bioconda that extends these capabailities to external programs.  All of the packages and programs we're using today can easily be installed with Anaconda/bioconda with the following steps:

* First, install Miniconda version 3.5 [available here](http://conda.pydata.org/miniconda.html).  During installation, be sure to allow miniconda to append to your .bashrc or .bash_profile (this will add it and all programs it installs to your PATH).  If installation goes well, the command ``` which python ``` should result in something like ```/Users/<yourusername>/anaconda/bin/python ``` or ```/home/<yourusername>/anaconda/bin/python ```

* Add bioconda channels to conda with the following commands:
  ```
  conda config --add channels r
  conda config --add channels conda-forge
  conda config --add channels bioconda
  ```
* Create the environment we'll be working in and install required packages with the command:

  ```
  conda create --name BIO598 python=3.5 snakemake fastqc bwa samtools bamtools picard samblaster freebayes bcftools snpsift bioawk
  ```
This will create a working environment called BIO598 containing python 3.5 (python 3 is required for snakemake) and all of the tools listed in the command.  You can see the full list of programs available through bioconda [listed here](https://bioconda.github.io/) and the full list of python packages available through Anaconda [listed here](https://docs.continuum.io/anaconda/pkg-docs).

You can load this environment with the command:
```
source activate BIO598
```
and leave the environment with the command:
```
source deactivate
```

If you're in your environment, you can easily add additional programs and packages with the command:
```
conda install <program/package name>
```

For example, if we also want to take a look at Bowtie2, another read mapper (we'll use bwa today), we can easily add it by entering our environment ``` source activate BIO598 ``` and typing ```conda install bowtie2 ```

## Reference-based genome assembly
There are, in general, two main flavors of genome assembly.  _De novo_ assembly involves taking raw sequencing reads and piecing them together into a genome assembly using only the information contained in the reads and their metadata (e.g., the sequences themselves, insert sizes, etc.).  While a number of _de novo_ assemblers exist and there's a great deal of work being done to improve alogrithms, lengthen sequencing reads, develop methods to increase insert sizes, etc.), _de novo_ assembly remains challenging, expensive (usually 100x or greater sequencing depth), and computationally demanding.  Fortunately, if we have a reference genome available to us, we can make do with much less sequencing (often 30x coverage or less; low coverage - 1-5x - are not uncommon for other purposes), and use tools that require far less memory and storage.

In this tutorial we'll walk through the basics of reference-based genome assembly.  While the dataset we're working with is tiny (we're using the human mitochondrial genome and a tiny subset of reads from the 1000 genomes project), you should be able to use this as a starting point for working with larger datasets down the road.

### What you you'll need
In the simplest cases, you need:
* **a reference genome assembly (in fasta format)**
* **sequencing reads (in fastq format - many processes are quicker if files are gzipped as well)**
* **a computing environment with the storage, memory, and software required**

#### Reference genome
This is the pre-prepared reference genome assembly that you're going to use to map the sequencing reads from your project/experiment.  In a perfect world, this assembly is of a high-quality, has a good set of annotations available (e.g., genes, functional elements, repeats, etc.), and is relatively closely related to the species that you're studying.  This probably won't be a problem if you're working with a model organism, but in other situations you'll have to consider whether a good assembly is available for a taxon evolutionarily close enough for your purposes (if not, you might need to think about assembling a reference for your project _de novo_).  For today, we're working with example sequencing reads from human samples and we have the human reference genome available, so we'll be fine.

Fasta format is a file format designed for DNA or protein sequences.  It looks something like (the first 10 lines of the ``` human_g1k_v37_MT.fasta``` file in the ```references``` directory:
  ```
  >MT
  GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTT
  CGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTC
  GCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATT
  ACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATA
  ACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCA
  AACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAA
  ACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCAC
  TTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAAT
  CTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACACCGCTGCTAACCCCATA
  ```
Here, the name of the sequence, ```MT``` is given after ```>```.  The lines that follow contain the sequence itself.  Because most fasta files wrap lines every 50-80 characters (this isn't uniform across files, unfortunately), there will often be many lines composing each sequence.  This file should only contain a single sequence, the 1000 genomes reference MT sequence, that's 16-17kb in length.  We can quickly check to make sure using (the very, very powerful) [bioawk](https://github.com/lh3/bioawk), which we installed earlier:
  ```
  bioawk -c fastx '{print ($name); print length($seq)}' human_g1k_v37_MT.fasta
  MT
  16569
  ```
We see that we do indeed have a single sequence called "MT" that's 16,569 bases in length.  

#### Sequencing reads
There are a few types of sequencing technologies out there, but Illumina is still dominant, so that's what we'll focus on today.  We won't get into the nuts and bolts of sequencing itself (there are plenty of resources available online that describe it well like [this video](https://www.youtube.com/watch?v=fCd6B5HRaZ8) or [this review](http://www.nature.com/nrg/journal/v17/n6/abs/nrg.2016.49.html)).

If you've sent your samples to a core for sequencing, they'll likely return to you a series of FASTQ files.  If you used single-end sequencing, your files can be concatenated into a single FASTQ file per sample per lane. On the other hand, if you used paired-end sequencing, you'll end up with two FASTQ files per sample per lane - one for the forward read and one for the reverse read (see the video I linked to in the previous paragraph for more information about the latter).  It's generally very important that paired reads are in the same order in both files.  If you're getting reads directly from a sequencing center, they're probably already organized this way.  However, you might have to sort and re-pair reads if you have, for example, stripped reads from a bam file containing an alignment (the README in the ```fastq``` directory explains how to do this, if needed).

For today's tutorial, we have paired-end reads from two individuals: ind1 and ind2.  Forward reads are in files ending in "1.fastq.gz" and reverse reads are in the two files ending in "2.fastq.gz".  We can take a look at the first two reads in the file ```ind1_1.fastq.gz``` using ```zcat``` on Linux or ```gzcat``` on Mac/Unix:
  ```
zcat ind1_1.fastq.gz | head -n 8

@H7AGFADXX131213:1:2110:18963:43686
ATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACGTTGTAGCCCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCCATTTCACTATCATATTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCCTATCCGGAATGCCCCGAC
+
>>=>>>???>??=?<?=?>>><@>=@<?>>7>><?><?<7?<@<7?>=?><8>@>@>?>=>>=A@>?>=B?????=B???@A@??>B@AB@@@?@A@?>@?>A@?;;@?@A@@@@>@>?=@@@?@@?<>?A@@A>@?A??@AA=@=>@A??@>>@A@>@BA>4@@AAAA??>@????@=B??@A?>??@?>>5>?:??AA??A?@>AB@?@A@@?A?@99C?@@@=@@9<B?@>?@?34@A?>@@@=199
@H7AGFADXX131213:2:1212:12308:62676
CATATGAAGTCACCCTAGCCATCATTCTACTATCAACATTACTAATAAGTGGCTCCTTTAACCTCTCCACCCTTATCACAACACAAGAACACCTCTGATTACTCCTGCCATCATGACCCTTGGCCATAATATGATTTATCTCCACACTAGCAGAGACCAACCGAACCCCCTTCGACCTTGCCGAAGGGGAGTCCGAACTAGTCTCAGGCTTCAACATCGAATACGCCGCAGGCCCCTTCGCCCTATTCTT
+
>?<>>>@?>?>?<>=?>>>>>>>?>>?A>;??????<?>?><??>?>?>>>=>??>A@>?@<?A???>@=>>??>>?@=@@=@=AA>@@=A>AA?A?@?A@>A@>A?@A@>?@??@=??@A>=?>@?@A?@??@>@A???@??@=@<@???@>@>@=A@A>>3?A?A???A??:A>@BA;>@9?;><4;@?@=<5:A=A?=@?@?A?@@AA@@A=@?@:@A?@>8A?8@@?=?AAAA?@7>A?B@?=770
```
In these files, each sequencing read is listed in a series of four lines:
* the ID line, beginning with @
* the sequence line (consisting of A, T, C, G, or N)
* a comment line (here, the comment lines only contain +)
* a quality score line (ASCII characters)



For today's tutorial, the reference genome is in the ``` reference ``` directory of this repository, the sequencing reads are in two files for each sample in the ``` fastq ``` directory of this repository, and the "Setting Up Anaconda" section above should take care of the software requirements.  Because we're working with a small dataset today, we won't need too much in the way of memory/storage, but bigger projects will often require at minimum a high-memory computer, but more likely high-performance computing clusters, dedicated servers, or online services such as Amazon Web Services.





