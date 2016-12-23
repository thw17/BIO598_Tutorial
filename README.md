# BIO598_Tutorial
A hands-on tutorial introducing users to reproducible reference-based genome assembly and variant calling.

This tutorial has been tested on a Mac and Linux operating systems and will assume you're working with one of them.  

## Setting up your environment
For today's tutorial, you'll need this repository and Anaconda.

#### Getting the repo
To get this repository (repo), move to the directory on your computer or cluster where you'd like to work and type the command:

```
git clone https://github.com/thw17/BIO598_Tutorial
```
This should create a directory called BIO598_Tutorial containing all of the contents of this repo.

#### Setting up Anaconda
Anaconda is an environment and package manager for Python and it makes installation, environment management, etc. simple without requiring root or administrator privileges.  Fortunately, its framework has been leveraged for a project called bioconda that extends these capabailities to external programs as well.  All of the packages and programs we're using today can easily be installed with Anaconda/bioconda with the following steps:

* First, install Miniconda version 3.5 [available here](http://conda.pydata.org/miniconda.html).  During installation, be sure to allow miniconda to append to your .bashrc or .bash_profile (this will add it and all programs it installs to your PATH).  If installation goes well, the command ``` which python ``` should result in something like ```/Users/<yourusername>/anaconda/bin/python ``` or ```/home/<yourusername>/anaconda/bin/python ```

* Add bioconda channels to conda with the following commands:
  ```
  conda config --add channels r
  conda config --add channels conda-forge
  conda config --add channels bioconda
  ```
* Create the environment we'll be working in and install required packages with the command:

  ```
  conda create --name BIO598 python=3.5 snakemake fastqc bwa samtools picard freebayes bcftools snpsift bioawk
  ```
* Load the new environment and add the samblaster package

  ```
  source activate BIO598
  conda install -c biobuilds samblaster
  ```

This will create a working environment called BIO598 containing python 3.5 (python 3 is required for snakemake) and all of the tools listed in the command.  You can see the full list of programs available through bioconda [listed here](https://bioconda.github.io/) and the full list of python packages available through Anaconda [listed here](https://docs.continuum.io/anaconda/pkg-docs).

The reason for installing samblaster separately in the final command is that we have to get it from the biobuilds channel.  Bioconda does support samblaster, but only on Linux.  So to ensure that your environment works across both Mac and Linux operating systems and to avoid confusion over the source of other programs that might be available from biobuilds and bioconda, we're installing the from biobuilds only in this single instance.

If you want to load our new environment, you can do so with the command:
```
source activate BIO598
```
and leave the environment with the command:
```
source deactivate
```

If you're in your environment, you can easily add additional programs (like we did for samblaster) and packages with the command:
```
conda install <program/package name>
```

For example, if we also want to take a look at Bowtie2, another read mapper (we'll use bwa today), we can easily add it by entering our environment ``` source activate BIO598 ``` and typing ```conda install bowtie2 ```

## Reference-based genome assembly
There are, in general, two main flavors of genome assembly.  _De novo_ assembly involves taking raw sequencing reads and piecing them together into a genome assembly using only the information contained in the reads and their metadata (e.g., the sequences themselves, insert sizes, etc.).  While a number of _de novo_ assemblers exist and there's a great deal of work being done to improve alogrithms, lengthen sequencing reads, develop methods to increase insert sizes, etc.), _de novo_ assembly remains challenging, expensive (usually 100x or greater sequencing depth), and computationally demanding.  Fortunately, if we have a reference genome available to us, we can make do with much less sequencing (often 30x coverage or less; low coverage - 1-5x - are not uncommon for some purposes), and use tools that require far less memory and storage.

In this tutorial we'll walk through the basics of reference-based genome assembly.  While the dataset we're working with is tiny (we're using the human mitochondrial genome and a tiny subset of reads from the 1000 genomes project), you should be able to use this as a starting point for working with larger datasets down the road.

### What you you'll need
In the simplest cases, you need:
* **a reference genome assembly (in fasta format)**
* **sequencing reads (in fastq format - many processes are quicker if files are gzipped as well)**
* **a computing environment with adequate storage and memory, and required software installed**

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
Here, the name of the sequence, ```MT``` is given after ```>```.  The lines that follow contain the sequence itself.  Because most fasta files wrap lines every 50-80 characters (this isn't uniform across files, unfortunately), there will often be many lines composing each sequence.  Today's (```human_g1k_v37_MT.fasta```) file should only contain a single sequence, the 1000 genomes reference MT sequence, that's 16-17kb in length.  We can quickly check to make sure using (the very, very powerful) [bioawk](https://github.com/lh3/bioawk), which we installed earlier:
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

The comment line usually either is a lone + or a + and then the sequence identifier repeated.  You'll generally ignore this line, but it's good to be aware of what's on it in case you're counting the number of reads belonging to a certain tile, for example.  In this case if the sequence identifier is repeated and you're counting the number of times a tile number is present in a file, your count will be double the actual number of reads coming from that tile.

The quality score line contains an ASCII character for every nucleotide in the sequence (i.e., it'll be 250 characters long for a 250 base read, 75 characters long for a 75 base read, etc.).  Each ASCII character can be converted into a integer PHRED quality score, ranging from 0 to 93, that indicates the probability that a particular base call is incorrect.  [See more details here](http://www.drive5.com/usearch/manual/quality_score.html).

If we want to quickly take a look at a summary of the reads in our fastq, we can use ```fastqc```. To analyze all of the files in the ```fastq``` directory and output the results to the ```stats``` directory, enter the following command from our main directory:
  ```
  fastqc -o stats fastq/*.fastq.gz
  ```
and it will quickly generate a series of reports for each fastq file that are summarized in the .html files (one per fastq) that you can view using your browser (e.g., Firefox, Chrome, etc.).

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
The first thing you should do when you get fastq (sequencing read) files is get a sense of their quality.  The quickest way to do this is to run ```fastqc``` and take a look at the reports.  We can run fastqc on every fastq file in ```fastq``` and output reports into ```stats``` with the command (from our main directory):
  ```
  fastqc -o stats fastq/*.fastq.gz
  ```
This should take less than a minute to complete and output a .zip and .html file for each example fastq (it'll take longer for full-sized files).  You can open the .html locally on your computer in your browser.

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

While it's exciting that we've successfully mapped our first sample, there are two potential problems with our command above.  First, while SAM format is convenient in that it's human-readable, alignment files are ENORMOUS, so we'll want to compress (BAM format) to minimize our data footprint.  Further, we want to ensure that all read-pairing, etc. didn't get lost in the mapping process.  We can easily add these steps to our pipeline by piping our output to ```samtools``` which handles this easily.  So our complete command now looks like this:
  ```
  bwa mem -M reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samtools fixmate -O bam - bam/ind1.bam
  ```
#### Processing our bam file: adding read groups, removing duplicates, sorting, and indexing.
The bam file we just created contains all of our alignments, but it's currently unordered, unlabeled, unfiltered, and unindexed.

##### Adding read groups
Read groups are very useful when you're working with multiple samples, sequencing lanes, flowcells, etc.  Importantly for us, it'll help our downstream variant caller label samples correctly and handle potential sequencing batch effects.  [Picard](https://broadinstitute.github.io/picard/command-line-overview.html) is very commonly used to add read groups to bam files, but ```bwa``` also has the ability to add read groups on the fly while mapping.  This latter option will save us time and space, so we'll add read groups to individual 1 with ```bwa``` by adding to our previous command:
  ```
  bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samtools fixmate -O bam - bam/ind1.bam
  ```
In a perfect world, we'd know more about our sample and could use these tags more appropriately.  ID is the name of a read group (containing a unique combination of the following tags).  SM is the sample name.  LB is the sequencing library. PU is the flowcell barcode.  PL is the sequencing technology.  A number of other options exist (see the [@RG section of the SAM/BAM specifications](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information).

*Note that the fastq files we're using today actually contain reads from multiple lanes, etc., as I randomly grabbed them from high-coverage 1000 genomes bam files.  But for simplicity's sake in this tutorial, we'll ignore that.

##### Removing Duplicates
During library preparation for sequencing, amplification steps can lead to PCR duplicates of reads.  The inclusion of duplicate reads can negatively affect our downstream analyses, so we need to remove them (this is the case with DNA sequencing - there's debate whether or not to do this in RNA-seq).  Most available tools, such as [Picard](https://broadinstitute.github.io/picard/command-line-overview.html), take a sorted bam file as input and output a sorted bam with duplicates either flagged or removed.  [Samblaster](https://github.com/GregoryFaust/samblaster), on the other hand, can take streaming output from bwa, which speeds up duplicate removal and allows us to produce one less bam file (saving us space).  We'll opt for Samblaster here, and incorporate it into our bwa command like so:
  ```
  bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samblaster -M | samtools fixmate -O bam - bam/ind1.rmdup.bam
  ```
This will flag, but not remove, duplicates in our bam.

*Note, again, that the fastq files we're using today actually contain reads from multiple lanes, etc., as I randomly grabbed them from high-coverage 1000 genomes bam files.  But for simplicity's sake in this tutorial, we'll ignore that.

##### Sorting bam files
If we weren't piping data directly to Samblaster, and instead using a different tool for duplicate removal, we'd have to sort our bam (in genome coordinate order) first.  But because we were able to build duplicate removal and adding our read groups into our pipeline, sorting will come last.

Sorting doesn't require too much explanation.  Most genomic datasets are huge, so it's inefficient to move across unsorted bam files (we need access to all reads covering a given base for variant calling, for example).  ```samtools sort ``` is widely used, and that's what we'll employ here.  Like our previous tools, it handles streaming input, so we can simply add to our previous command to save space:
  ```
  bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samblaster -M | samtools fixmate - - | samtools sort -O bam -o bam/ind1.rmdup.sorted.bam -
  ```


So, now, after running our command, we can view the first 10 lines of ```ind1.rmdup.sorted.bam``` with the command:
  ```
  samtools view -h ind1.rmdup.sorted.bam | head
  ```
```samtools``` is required to view the bam file because of how it is compressed.  The ```-h``` flag ensures that the bam header is printed along with the records.

Anyway, that command will give the following output:
  ```
  @HD     VN:1.3  SO:coordinate
  @SQ     SN:MT   LN:16569
  @RG     ID:ind1 SM:ind1 LB:ind1 PU:ind1 PL:Illumina
  @PG     ID:bwa  PN:bwa  VN:0.7.15-r1140 CL:bwa mem -M -R @RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina human_g1k_v37_MT.fasta ind1_1.fastq.gz ind1_2.fastq.gz
  @PG     ID:SAMBLASTER   VN:0.1.23       CL:samblaster -i stdin -o stdout -M
  H7AGFADXX131213:1:1216:6248:58707       353     MT      1       60      216H34M =       1       206     GATCACAGGTCTATCACCCTATTAACCACTCACG      ;@8<B=A>>17@???A=@@B@@A@A>?A?A<95)      NM:i:0  MD:Z:34 AS:i:34 XS:i:0  RG:Z:ind1       SA:Z:MT,16354,+,216M34S,60,2;
  H7AGFADXX131213:1:1115:9742:18153       99      MT      1       60      122S128M        =       1       193     TCGCTCCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATC      >>5?@>?7=>?>>?>>><?<?>>>>==>>>@?>??>>>?><A>=>>>=8@<?>>@=>?<=>@?=AA?@?>?@@@@AAA??@BAAA@???@>A>A>9??@@@@A??@A@AA??=@?@A>9B?@?A?@@>@?@??@@>@A=??BA?@>?>?@=@@@>:?@B@AA?A@@A@??A@@@>=???B@B@9??@@>>@@<&==>@A>9?:B?@?@=>A>@9B???7=@<<B?@A9>B?AA=A?A@??>?7B@?9968      NM:i:0  MD:Z:128        AS:i:128        XS:i:0  RG:Z:ind1       SA:Z:MT,16448,+,122M128S,60,1;  MQ:i:60
  H7AGFADXX131213:1:1203:19236:9314       99      MT      1       60      93S157M =       1       240     GCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGAATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCCATTAT      >>=?@@=@>??<@6:>>><6@<??=@:6=>==A><@>=><3>>?@?=>@>>@A??@>>?@=?>>@=9?A>==@@?@@@=>=<@=B??@<'@5=;@>?@=@==@=@=?=?;A?B??@?@==?>B?A588<A>A<>@<?@?==><?@=);<>@@B;&?;A>4;3/1&3<;5>;'27B<53:<=?/:4?&?=2=@=>@??A9@??=A???C@?=1=,:&9@,=?@>@>AAA@A@?==@>?AA@>?=>?1:77;      NM:i:2  MD:Z:71T79T5    AS:i:147        XS:i:0  RG:Z:ind1       SA:Z:MT,16477,+,93M157S,60,1;   MQ:i:60
  H7AGFADXX131213:1:2110:2821:77416       163     MT      1       60      89S161M =       1       225     AAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCCATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATGGATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATCCCATTATTTAT      <=;=><8<A>@>>>=71<?>>>>>>?>>?><??>>>=>>>>>>?>>>>?>?>>>>@?>?<@=7??????A@?@?@?@?@=@??@=9B===A@@A>@>@@@A@@@A>?AB??@A@=>@=@?@=5@<B@AB?@??A?@@@@@@=?@@@B@A@:?@B@??8;;8<?@@A>:@:B?@@@A@@>@:B@A<7AA=@@=B?:@<?AA>A?B@@@@@7>>;???@A>@@AAA?A@@@?19=><@A@@?@A@??@?>==      NM:i:1  MD:Z:151T9      AS:i:156        XS:i:0  RG:Z:ind1       SA:Z:MT,16481,+,89M161S,60,1;   MQ:i:60
  H7AGFADXX131213:1:1211:16620:80940      353     MT      1       60      213H37M =       1       221     GATCACAGGTCTATCACCCTATTAACCACTCACGGGA   @A??@=@??@@A@>@A=??A@?A>@=?A>A@@=251/   NM:i:0  MD:Z:37 AS:i:37 XS:i:0  RG:Z:ind1       SA:Z:MT,16357,+,213M37S,60,1;
  ```
You can see that there are a few changes.  Most notably, there is a @HD line in the header indicating that the file is sorted in coordinate order, we have an additional @PG flag with information about the SAMBLASTER command, and we now have a @RG line in the header with our read group information and a RG tag on each record with its corresponding read group ID (here they're all ind1).

##### Indexing
Again, bam files can get pretty big and they're compressed, so we need to index them for other tools to use them.  Here's a simple command for indexing our bam (from our main directory):
  ```
  samtools index bam/ind1.rmdup.sorted.bam
  ```
This creates a file with a .bai extension that needs to remain in the same directory as its corresponding bam.

##### Summarizing the contents of your bam file
Now that we have our final bam file, we might be interested in knowing what it contains.  More specifically, we should get a sense of how many duplicates there were, how successful mapping was, and other measures along those lines.  Fortunately, [samtools offers tools to calculate summary statistics](http://www.htslib.org/doc/samtools.html).  We'll calculate stats using ```samtools stats``` because it provides a bit more detail, but you should have a look at ```samtools flagstat``` as well.  From our main directory, enter the command:
  ```
  samtools stats bam/ind1.rmdup.sorted.bam | grep ^SN | cut -f 2- > stats/ind1.rmdup.sorted.bam.stats
  ```
Because samtools stats offers a huge range of statistics including a number of very big tables in our output, we'll just grab the summary statistics.  This is what ```grep ^SN | cut -f 2-``` does.

We can print the contents of each file to screen using the ```cat``` command.  Here's the command and its result for the samtools stats output:
  ```
  cat stats/ind1.rmdup.sorted.bam.stats
  
  raw total sequences:	10000
  filtered sequences:	0
  sequences:	10000
  is sorted:	1
  1st fragments:	5000
  last fragments:	5000
  reads mapped:	9882
  reads mapped and paired:	9764	# paired-end technology bit set + both mates mapped
  reads unmapped:	118
  reads properly paired:	9616	# proper-pair bit set
  reads paired:	10000	# paired-end technology bit set
  reads duplicated:	12	# PCR or optical duplicate bit set
  reads MQ0:	0	# mapped and MQ=0
  reads QC failed:	0
  non-primary alignments:	117
  total length:	2490364	# ignores clipping
  bases mapped:	2460864	# ignores clipping
  bases mapped (cigar):	2411834	# more accurate
  bases trimmed:	0
  bases duplicated:	3000
  mismatches:	21935	# from NM fields
  error rate:	9.094738e-03	# mismatches / bases mapped (cigar)
  average length:	249
  maximum length:	250
  average quality:	25.9
  insert size average:	566.3
  insert size standard deviation:	902.9
  inward oriented pairs:	4669
  outward oriented pairs:	80
  pairs with other orientation:	2
  pairs on different chromosomes:	0
  ```
As you can see, it gives us a lot of information about number of reads, mapping, pairing, etc.  A quick glance shows us that our mapping was quite successful (9882 out of 10000 reads mapped).  We also had very few PCR duplicates (12 reads - note that this can be calculated because we flagged, but didn't remove duplicates using samblaster earlier), but this is probably because I randomly sampled read pairs to create our fastq files.

#### Variant calling
Now that our bam is processed and indexed, it's time to call variants!!  There are a few popular variant callers in use, but today we'll be using [Freebayes](https://github.com/ekg/freebayes) because it works well on both Mac and Linux and can be easily installed via conda.  Freebayes is extremely flexible (it allows a great deal of performance tuning, handles atypical ploidy, can take tumor/normal pairings, etc.), but it's simplest case is perfect for us today.  From our main directory, enter the command:
  ```
  freebayes -f reference/human_g1k_v37_MT.fasta bam/ind1.rmdup.sorted.bam > vcf/ind1.raw.vcf
  ```
The resulting list of variants is quite small.  Because the mitochondrial genome is haploid, we really only expect to see variant records for sites where individual 1 differs from the reference. However, because there are so many copies of mitochondria sequenced during normal genome sequencing, we might expect a bit of heteroplasmy (variation among different mitochondrial genomes within an individual) and the possibility of heterozygous calls.

We can print all heterozygous calls with a site quality greater than 30 (less than a 1 in 1000 chance of being incorrectly identified as a polymorphic site) to the screen using SnpSift with the following command (from the main directory):
  ```
  cat vcf/ind1.raw.vcf | SnpSift filter "'(isHet(GEN[0])) & (QUAL >= 30)'" 
  ```
* Note that we have to use double quotes around the single quotes because of how conda wraps SnpSift.  SnpSift is actually a .jar file that has to be called using java (e.g., ```java -Xmx2g -jar SnpSift.jar ...```).  Conda wraps the command line in a bash script, and to pass your filter expression ```'(isHet(GEN[0])) & (QUAL >= 30)'``` to a bash script from a command line, you have to enclose it in double quotes.  If you were to use a version of SnpSift that you downloaded directly (as a .jar), you'd have to add the java arguments before calling SnpSift and remove the double quotes (leaving only the single quotes).

You should see four heterozygous calls passing filters, along with the full VCF header.  While each caller, unfortunately, produces its own version of a VCF, the header usually contains a very detailed description of how to interpret the file.

That's all we're going to do in terms of filtering vcf files today. In case you're interested in doing more filtering, SnpSift is a very flexible and powerful vcf parser.  You can find out more about what it can do [here](http://snpeff.sourceforge.net/SnpSift.html).  We also downloaded bcftools, which is another useful tool for working with vcf files ([check out the manual here](https://samtools.github.io/bcftools/bcftools.html))

### Putting it all together
To run our pipeline on our two samples, we can simply run the following commands:

```
bwa mem -M -R '@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind1_1.fastq.gz fastq/ind1_2.fastq.gz | samblaster -M | samtools fixmate - - | samtools sort -O bam -o bam/ind1.rmdup.sorted.bam -
```
```
samtools index bam/ind1.rmdup.sorted.bam
```
```
samtools stats bam/ind1.rmdup.sorted.bam | grep ^SN | cut -f 2- > stats/ind1.rmdup.sorted.bam.stats
```
```
bwa mem -M -R '@RG\tID:ind2\tSM:ind2\tLB:ind2\tPU:ind2\tPL:Illumina' reference/human_g1k_v37_MT.fasta fastq/ind2_1.fastq.gz fastq/ind2_2.fastq.gz | samblaster -M | samtools fixmate - - | samtools sort -O bam -o bam/ind2.rmdup.sorted.bam -
```
```
samtools index bam/ind2.rmdup.sorted.bam
```
```
samtools stats bam/ind2.rmdup.sorted.bam | grep ^SN | cut -f 2- > stats/ind2.rmdup.sorted.bam.stats
```
And while we ran Freebayes on a single bam file before, it will just as easily take two files for joint calling:

```
freebayes -f reference/human_g1k_v37_MT.fasta bam/ind1.rmdup.sorted.bam bam/ind2.rmdup.sorted.bam > vcf/joint.raw.vcf
```

## Making your pipeline reproducible
So, now that we have our pipeline, it's time to make it as reproducible as possible.  There are a few reasons for this, including (but not limited to), allowing us to keep track of and control which versions of programs we're using, making it easy to share with collaborators, ensuring you know exactly what you did down the line, and making your research open so that reviewers and other researchers can reproduce your analyses and adapt them to their own research.

### Sharing your environment
Luckily, by working with Anaconda, we can easily share our environment (including exact versions of tools).  To do this, we can enter the command:
  ```
  conda env export > environment.yml
  ```
and we'll get a file that looks something like:
  ```
  name: BIO598
  channels: !!python/tuple
  - !!python/unicode
    'bioconda'
  - !!python/unicode
    'defaults'
  dependencies:
  - biobuilds::samblaster=0.1.23=0
  - bioconda::bcftools=1.3.1=1
  - bioconda::bioawk=1.0=0
  - bioconda::bwa=0.7.15=0
  - bioconda::curl=7.45.0=2
  - bioconda::dropbox=5.2.1=py35_0
  - bioconda::fastqc=0.11.5=1
  - bioconda::filechunkio=1.6=py35_0
  - bioconda::freebayes=1.0.2.29=py35_2
  - bioconda::ftputil=3.2=py35_0
  - bioconda::java-jdk=8.0.92=1
  - bioconda::picard=2.5.0=1
  - bioconda::pysftp=0.2.8=py35_0
  - bioconda::samtools=1.3.1=4
  - bioconda::snakemake=3.8.2=py35_0
  - bioconda::snpsift=4.3=1
  - bioconda::urllib3=1.12=py35_0
  - cffi=1.8.3=py35_0
  - cryptography=1.5.3=py35_0
  - docutils=0.12=py35_2
  - idna=2.1=py35_0
  - libgcc=4.8.5=1
  - ncurses=5.9=8
  - openssl=1.0.2j=0
  - paramiko=2.0.2=py35_0
  - pip=9.0.1=py35_0
  - pyasn1=0.1.9=py35_0
  - pycparser=2.16=py35_0
  - python=3.5.2=0
  - pyyaml=3.12=py35_0
  - readline=6.2=2
  - requests=2.11.1=py35_0
  - setuptools=27.2.0=py35_0
  - six=1.10.0=py35_0
  - sqlite=3.13.0=0
  - tk=8.5.18=0
  - wheel=0.29.0=py35_0
  - wrapt=1.10.8=py35_0
  - xz=5.2.2=0
  - yaml=0.1.6=0
  - zlib=1.2.8=3
  prefix: /Users/thw/anaconda/envs/BIO598
  
  ```
This is perfect, except for the prefix line at the end (which will get in the way for future users that don't have the same prefix).  So, we'll need to delete that line.  You can do this in a text editor like ```nano``` or ```vi``` (both come standard in most Linux/UNIX environments).  If you haven't used ```vi``` or something similar before, you're probably going to have huge issues writing, saving, and exiting, so I would probably avoid it for now.  ```nano``` can be used by simply typing ```nano <filename>```, using your keyboard arrows to move down to the bottom of the file, manually deleting the last two lines, and following the instructions on the bottom of the screen to exit.

We can also delete the last line of the file easily with ```sed```.  The command for doing so is ```sed '$d' <filename>```. This will print the file, without its last line, to the screen.  If your file is like mine, with one blank line at the end, and the prefix statement on the second to last line, we can delete both and print to a new file with:
  ```
  sed '$d' environment.yml | sed '$d' > BIO598.yml
  ```
Now we can easily share our environment file, ```BIO598.yml```, with anyone.  They can then create an identical environment with the command:
  ```
  conda env create -f BIO598.yml
  ```
This will create an environment called ```BIO598``` on their computer (because of the "name" row in the .yml file).

### Sharing your pipeline

#### Bash script
Now that you're able to share your environment, what about your pipeline?  One option is that you can simply share your list of commands, and ask users to run them on their own.  You could make this easy for them by creating a short shell script.  I've included one in the main directory, ```example.sh```.  If we look inside (try ```cat example.sh```), we can see that it's very similar to our list of commands:
```
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
```
The one difference is the first line, ```#!/usr/bin/env bash```, which is called the shebang (the other lines beginning with ```#``` are just comments).  It allows us to run the script like a program from the command line.  For example, to run ```example.sh```, we could type (from the main directory):
```
./example.sh
```
and you should see that our pipeline is runs exactly as it did earlier, only now using a single command (you might have gotten an "Exception in thread "main" picard.PicardException..." error along the way - this will happen if you did not delete the .dict file before starting, as it won't overwrite).  Note that the ```./``` is required for your shell to recognize your script as a program.

One quick note: if we want others to be able to use our script (either to run it, or adapt it for their own purposes), we need to ensure that the script gives everyone permission to read, write, or execute it.  You can do this with the ```chmod``` command:
```
chmod 777 example.sh
```
This will allow _anyone_ to read, write, or run your script.  If you need to limit some aspects of reading/writing/executing, see [this page](http://ss64.com/bash/chmod.html) for more information on the different codes you can use.

#### Snakemake
So, for the general purpose of sharing our analyses, a bash script works well.  But what if we want to allow users to flexibly add more samples or change the reference genome?  What about if we need to run processes in parallel on a high performance cluster (this isn't a worry for our tiny example, but is often critical when working with real genomic datasets)?

One solution is [Snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home), a program that allows users to construct workflows in Python.  I will only provide a very, very brief introduction in this tutorial.  For more information, I highly recommend checking out the [official tutorial](http://snakemake.bitbucket.org/snakemake-tutorial.html), [reading the documentation](https://bitbucket.org/snakemake/snakemake/wiki/Documentation), and browsing the [Google group](https://groups.google.com/forum/#!forum/snakemake).

Snakemake is a Python 3 package (for Python users: note that it is NOT compatible with Python 2).  If you completed Anaconda/Miniconda installation as directed under the section "Setting Up Anaconda" above, you should have installed Snakemake into your ```BIO598``` environment.  

Snakemake style and syntax is based on Make, so each snake workflow is composed of a series of rules, each with a minimum setup similar to:
```
rule <rule_name>:
  input:
    <input_file_name>
  output:
    <output_file_name>
  shell:
    "command to be run"
```
A more contrete example, such as our very first ```bwa``` command from above, would look like:
```
rule bwa_mem_mapping:
  input:
    ref="reference/human_g1k_v37_MT.fasta",
    fastq1="fastq/ind1_1.fastq.gz",
    fastq2="fastq/ind1_2.fastq.gz"
  output:
    "bam/ind1.sam"
  shell:
    "bwa mem -M {input.ref} {input.fastq1} {input.fastq2} > {output}"
```
Breaking down each part, we first have to declare our "rule" with ```rule bwa_mem_mapping:```.  Note the lack of indentation (tabs, spaces, and whitespace are very important in python) and the colon.  We then indent and declare our input with the keyword ```input:```.  Again, note the indentation (a single tab) and the colon.  Next, we indent and declare variables that are required for this rule as input.  Here, we declare ```ref``` (a reference genome), ```fastq1``` (the first of the two paired-end fastq files), and ```fastq2``` (the second of the two paired-end fastq files). Note the indentation (two tabs), the quotes around the file paths, and the commas.




