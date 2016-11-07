# BIO598_Tutorial
A hands-on tutorial introducing users to reproducible genome reassembly and variant calling

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

