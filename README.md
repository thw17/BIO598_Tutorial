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
* a reference genome assembly (in fasta format)
* sequencing reads (in fastq format - many processes are quicker if files are gzipped as well)
* a computing environment with the storage, memory, and software required.

For today's tutorial, the reference genome is in the ``` reference ``` directory of this repository, the sequencing reads are in two files for each sample in the ``` fastq ``` directory of this repository, and the "Setting Up Anaconda" section above should take care of the software requirements.  Because we're working with a small dataset today, we won't need too much in the way of memory/storage, but bigger projects will often require at minimum a high-memory computer, but more likely high-performance computing clusters, dedicated servers, or online services such as Amazon Web Services.

##### Reference genome assembly
The reference genome assembly



