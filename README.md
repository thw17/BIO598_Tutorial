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

For today's tutorial, we have paired-end reads from two individuals: ind1 and ind2.  Forward reads are in files ending in "1.fastq.gz" and reverse reads are in the two files ending in "2.fastq.gz".  We can take a look at the first two reads in the file ```ind1_1.fastq.gz``` using ```zcat``` on Linux or ```gzcat``` on Mac/Unix combined with ```head```:
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
* the sequence identifier line, beginning with @
* the sequence line (consisting of A, T, C, G, or N)
* a comment line (here, the comment lines only contain +)
* a quality score line (ASCII characters)

The sequence identifier can contain a lot of information (see [Illumina's description for more information](http://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm)), the combination of which will identify individual reads uniquely.  For whatever reason, we've don't have records for the instrument or run number, and instead have IDs organized as: flowcellID:lane:tile:x_position:y_position.

The sequence line here is 250 bases long and contains the nucleotide sequence corresponding to each read.

The comment line usually either is a lone + or a + and then the sequence identifier repeated.  You'll generally ignore this line, but it's good to be aware of what's on it in case you're counting the number of reads belonging to a certain tile, for example.  In this case if the sequence identifier is repeated, your count will be double the actual number of reads coming from that tile.

The quality score line contains an ASCII character for every nucleotide in the sequence (i.e., it'll be 250 characters long for a 250 base read, 75 characters long for a 75 base read, etc.).  Each ASCII character can be converted into a integer PHRED quality score, ranging from 0 to 93, that indicates the probability that a particular base call is incorrect.  [See more details here](http://www.drive5.com/usearch/manual/quality_score.html).

If we want to quickly take a look at a summary of the reads in our fastq, we can use ```fastqc```. To analyze all of the files in the ```fastq``` directory, change to this directory and enter the command:
  ```
  fastqc *.fastq.gz
  ```
and it will quickly generate a series of reports for each fastq file that are summarized in the .html files (one per fastq).

We can also use ```bioawk``` to parse and analyze fastq files as well.  For example we can count all of the reads in each file:
  ```
  for i in $(ls *.fastq.gz); do echo $i; bioawk -c fastx 'END{print NR}' $i; done
  ind1_1.fastq.gz
  5000
  ind1_2.fastq.gz
  5000
  ind2_1.fastq.gz
  5000
  ind2_2.fastq.gz
  5000
  ```
Or count the number of reads from tile 2111 in ```ind1_1.fastq.gz```:
  ```
  bioawk -c fastx '{print $name}' ind1_1.fastq.gz | grep ':2111:' | wc -l
  ```
Or count the number of reads from each tile in ```ind1_1.fastq.gz```:
  ```
  bioawk -c fastx '{print $name}' ind1_1.fastq.gz | cut -d':' -f 3 | sort | uniq -c
  ```
#### Programs
Reference-based assembly requires a number of tools for steps including (but not limited to) read trimming, read mapping, bam processing (sorting, removing duplicates, adding read group information, etc.), variant calling, and variant filtering.  We installed everything we need for today's tutorial with the conda command above.  I'll discuss these steps in more detail below.

### Assembly: Step-by-step
For today's tutorial, the reference genome is in the ``` reference ``` directory of this repository, the sequencing reads are in two files for each sample in the ``` fastq ``` directory of this repository, and the "Setting Up Anaconda" section above should take care of the software requirements.  Because we're working with a small dataset today, we won't need too much in the way of memory/storage, but bigger projects will often require at minimum a high-memory computer, but more likely high-performance computing clusters, dedicated servers, or online services such as Amazon Web Services.

So now we'll walk through the major steps of reference-based genome assembly and build a pipeline along the way.

#### Preparing your reference
Most reference genomes are quite large, so it's very inefficient to linearly search through, say, 3 billion characters spread across 25 sequences.  So, many programs use various hashing/indexing strategies to more efficiently handle the data. We'll create two of the most commonly required index files: .dict and .fai.  We'll also create the required index for ```bwa```, our read mapper.  This is all we'll need for our purposes today, but check any additional tools you incorporate into your work down the line to see if they require additional indexing or processing of reference assemblies.

The three commands we'll use to prepare our reference, in the ```reference``` directory, are:
  ```
  samtools faidx human_g1k_v37_MT.fasta
  
  picard CreateSequenceDictionary R=human_g1k_v37_MT.fasta o=human_g1k_v37_MT.dict
  
  bwa index human_g1k_v37_MT.fasta
  ```
This will be quick on our small reference, but the bwa indexing in particular can take much longer on a full, human-sized reference.

And that's it.  All of our reference files and indices are now contained in our reference directory.

#### Fastq Quality Control
The first thing you should do when you get fastq files is get a sense of their quality.  The quickest way to do this is to run ```fastqc``` and take a look at the reports.  We can run fastqc on every fastq file in references with the command:
  ```
  fastqc *.fastq.gz
  ```
This should take less than a minute to complete and output a .zip and .html file for each fastq.  You can open the .html locally on your computer in your browser.

If your sequences are of an unexpectedly low quality it might be worth contacting your sequencing center.  Otherwise, lower quality towards the ends of the reads, some PCR duplication, and sequence content that's slightly off are all pretty typical in sequencing experiments.  If you're interested in trimming low-quality ends of reads or removing sequencing adapters, take a look at programs like [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) or [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).  Many read mappers like BWA handle adapters and low quality ends very well and downstream variant callers will build base quality into their genotype likelihoods, so it's often not necessary to do any preprocessing if your goal is variant calling.  For today's purposes, we'll take this approach and move onto our next step.

#### Read mapping
Our next step involves mapping our reads to our reference.  We'll use the [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml) algorithm to do this, as it's among the most popular mappers and works very well mapping reads to a closely related reference.  Other popular mappers include [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [Novoalign](http://www.novocraft.com/products/novoalign/), and [Stampy](http://www.well.ox.ac.uk/stampy) (Stampy, in particular, for mapping to a very diverged reference).

The command line for bwa mem is quite straightforward.  From the main directory (/path/to/BIO598_Tutorial), we can execute the following command:
  ```
  bwa mem -M reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz > bam/ind1.sam
  ```
This will map both sets of paired-end reads from ind1 to our reference genome and output the alignments in SAM format in our ```bam``` directory.  The ```-M``` flag ensures compatibility of the alignments with downstream tools.

If we take a look at the top of the file using ```head -n 8 ind1.sam```, we see the following:
  ```
  @SQ     SN:MT   LN:16569
  @PG     ID:bwa  PN:bwa  VN:0.7.15-r1140 CL:bwa mem -M human_g1k_v37_MT.fasta ind1_1.fastq.gz ind1_2.fastq.gz
  H7AGFADXX131213:1:2110:18963:43686      99      MT      6969    60      250M    =       7187    468     ATTGTATTAGCAAACTCATCACTAGACATCGTACTACACGACACGTACTACGTTGTAGCCCACTTCCACTATGTCCTATCAATAGGAGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCTACGCCAAAATCCATTTCACTATCATATTCATCGGCGTAAATCTAACTTTCTTCCCACAACACTTTCTCGGCCTATCCGGAATGCCCCGAC      >>=>>>???>??=?<?=?>>><@>=@<?>>7>><?><?<7?<@<7?>=?><8>@>@>?>=>>=A@>?>=B?????=B???@A@??>B@AB@@@?@A@?>@?>A@?;;@?@A@@@@>@>?=@@@?@@?<>?A@@A>@?A??@AA=@=>@A??@>>@A@>@BA>4@@AAAA??>@????@=B??@A?>??@?>>5>?:??AA??A?@>AB@?@A@@?A?@99C?@@@=@@9<B?@>?@?34@A?>@@@=199      NM:i:0  MD:Z:250        AS:i:250        XS:i:0
  H7AGFADXX131213:1:2110:18963:43686      147     MT      7187    60      250M    =       6969    -468    ACACTTCCTCGGCCTAACCGGAATGCCCCGACGTTACTCGGACTACCCCGATGCATACACCACATGAAACATCCTATCATCTGTAGGCTCATTCATTTCTCTAACAGCAGTAATATTAATAATTTTCATGATTTGAGAAGCCTTCGCTTCGAAGCGAAAAGTCCTAATAGTAGAAGAACCCCCCATAAACCTGGAGTGACTATATGGATGCCCCCCACCCTACCACACATTCGAAGACCCCGTATACATA      #################?6=,?<:4<821+08=@???B9<>==77?=<6+<8'?=;9:9<==?@=<5<>==<>?>,*<=?=>>>?=<>@>?A?>?@@@>?=>><>@>>?<A@@?<<;7<=<,??=@<=@<?@@?<>>?>>7<?A8><>;87>B@6=>@@@=@===<+83(=?9598&0;91%<;=2<=79=<809;>78;3;:>><98>80;9:/45==<69<16<8;=>64;);<=*>?'<=0<:;<;9      NM:i:4  MD:Z:6T9T164T55A12      AS:i:230        XS:i:0
  H7AGFADXX131213:2:1212:12308:62676      99      MT      3728    60      250M    =       3978    500     CATATGAAGTCACCCTAGCCATCATTCTACTATCAACATTACTAATAAGTGGCTCCTTTAACCTCTCCACCCTTATCACAACACAAGAACACCTCTGATTACTCCTGCCATCATGACCCTTGGCCATAATATGATTTATCTCCACACTAGCAGAGACCAACCGAACCCCCTTCGACCTTGCCGAAGGGGAGTCCGAACTAGTCTCAGGCTTCAACATCGAATACGCCGCAGGCCCCTTCGCCCTATTCTT      >?<>>>@?>?>?<>=?>>>>>>>?>>?A>;??????<?>?><??>?>?>>>=>??>A@>?@<?A???>@=>>??>>?@=@@=@=AA>@@=A>AA?A?@?A@>A@>A?@A@>?@??@=??@A>=?>@?@A?@??@>@A???@??@=@<@???@>@>@=A@A>>3?A?A???A??:A>@BA;>@9?;><4;@?@=<5:A=A?=@?@?A?@@AA@@A=@?@:@A?@>8A?8@@?=?AAAA?@7>A?B@?=770      NM:i:0  MD:Z:250        AS:i:250        XS:i:0
  H7AGFADXX131213:2:1212:12308:62676      147     MT      3978    60      250M    =       3728    -500    CATAGCCGAATACACAAACATTATTATAATAAACCCCCTCACCACTACAATCTTCCTAGGAACAACATATGACGCCCTCTCCCCTGAACTCTACACAACATATTTTGTCACCAAGACCCTACTTCTAACCTCCCTGTTCTTATGAATTCGAACAGCATACCCCCGATTCCGCTACGACCAACTCATACACCTCCTATGAAAAAACTTCCTACCACTCACCCTAGCATTACTTATATGATATGTCTCCATA      @>=B?>95@>?@>@?@?@>??>;??>+@@@7@?8&?@@@?>?@@><=@A@B?@@><?A?>@@<:;7=826?=6//(;1;>>=>;>=?==>@>?@@?A@@>@@@???<:?=>?@@=0==>>>?@?>>@>@>A@@?@=@@?@@?=@A@AA6?A@@A?@?>;??@?5??A??9?><:9?@=>>>?@>@?>>??>>=?>>>??????>>??=>>>><7>>>8>>>>@=>?>=>>?>>>>>>>>>?<@>@><;<<      NM:i:2  MD:Z:34A40A174  AS:i:240        XS:i:0
  H7AGFADXX131213:2:1114:9857:73128       73      MT      15987   60      250M    =       15987   0       CACCCAAAGCTAAGATTCTAATTTAAACTATTCTCTGTCCTTTCATGGGGAAGCAGATTTGGGTACCACCCAAGTATTGACTCACCCATCAACAACCGCTATTCATCCAGTACATTACTGCCAGCCACCATGAATATTTTACGGTACCATAAATACTTGACCACCTGTAGTACATCAAAACCCACTCCACATCAAAACCCCCTCCCCATGCTTACTACCTAGTACATCAATCACTCCACAACTATCACAC      ;?;>?0??>?@<@>+>7<?>1/?@=??<A7>?>?>?:8/-B<2?1=-*((:+9??=A>@8=;&5;<A0:>?@9-3;7<=A>A@@*?>@?8@@>@A###########################################################################################################################################################      NM:i:16 MD:Z:38T63G0T2T0T0C29G36A8A30A1G1A6G6A0C2T12    AS:i:170        XS:i:0
  H7AGFADXX131213:2:1114:9857:73128       133     MT      15987   0       *       =       15987   0       ATGGCCCAAAAGAGAGGAGTGCAGGAAGATTTTGCGGGATAATGAATACCCGAAGAAGGGTGGACAAGGGATCCCTATCTCCGGTGGAACATACATAGGGCCGAGAAAGGACTTAACTGTAATGAGCTATGCATAATAGATAACTGTACAGTTCATCAAATTGTGAGGATGAGTATGATTATTTGTGTCCTCGAGGGAGAGGCGAGGACATGGATAAGCCGCTGAGTTGGGGACGTTGAAGGTTAATTGC      ##########################################################################################################################################################################################################################################################      AS:i:0  XS:i:0
  ```
As you can see, our first line (@SQ) contains information about the reference genome (there would be more lines if there were more sequences in our reference), our second line (@PG) contains information about our bwa commands, and the remaining lines contain information about each mapped read.  You can find more information about what information is contained in each record in the [SAM/BAM specificaions](https://samtools.github.io/hts-specs/SAMv1.pdf).

While it's exciting that we've successfully mapped our first sample, there are two potential problems with our command above.  First, while SAM format is convenient in that it's human-readable, alignment files are ENORMOUS, so we'll want to compress to minimize our data footprint.  Further, we want to ensure that all read-pairing, etc. didn't get lost in the mapping process.  We can easily add these steps to our pipeline by piping our output to ```samtools``` which handles this easily.  So our complete command now looks like this:
  ```
  bwa mem -M reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samtools fixmate -O bam - bam/ind1.bam
  ```
#### Processing our bam file: adding read groups, removing duplicates, and sorting.
The bam file we just created contains all of our alignments, but it's currently unordered, unlabeled, and unfiltered.

##### Adding read groups
Read groups are very useful when you're working with multiple samples, sequencing lanes, flowcells, etc.  Importantly for us, it'll help our downstream variant caller label samples correctly and handle potential sequencing batch effects.  [Picard](https://broadinstitute.github.io/picard/command-line-overview.html) is very commonly used to add read groups to bam files, but ```bwa``` also has the ability to add read groups on the fly while mapping.  This latter option will save us time and space, so we'll add read groups to individual 1 with ```bwa``` by adding to our previous command:
  ```
  bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samtools fixmate -O bam - bam/ind1.bam
  ```
In a perfect world, we'd know more about our sample and could use these tags more appropriately.  ID is the name of a read group (containing a unique combination of the following tags).  SM is the sample name.  LB is the sequencing library. PU is the flowcell barcode.  PL is the sequencing technology.  A number of other options exist (see the [@RG section of the SAM/BAM specifications](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information).

*Note that the fastq files actually contain reads from multiple lanes, etc., as I grabbed them from high-coverage 1000 genomes bam files.  But for simplicity's sake in this tutorial, we'll ignore that.

##### Removing Duplicates
During library preparation for sequencing, amplification steps can lead to PCR duplicates of reads.  The inclusion of duplicate reads can negatively affect our downstream analyses, so we need to remove them (this is the case with DNA sequencing - there's debate whether or not to do this in RNA-seq).  Most available tools, such as [Picard](https://broadinstitute.github.io/picard/command-line-overview.html), take a sorted bam file as input and output a sorted bam with duplicates either flagged or removed.  [Samblaster](https://github.com/GregoryFaust/samblaster), on the other hand, can take streaming output from bwa, which speeds up duplicate removal and allows us to produce one less bam file (saving us space).  We'll opt for Samblaster here, and incorporate it into our bwa command like so:
  ```
  bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samblaster -M | samtools fixmate -O bam - bam/ind1.rmdup.bam
  ```
This will flag, but not remove, duplicates from our bam.

*Note that the fastq files actually contain reads from multiple lanes, etc., as I grabbed them from high-coverage 1000 genomes bam files.  But for simplicity's sake in this tutorial, we'll ignore that.

#### 





