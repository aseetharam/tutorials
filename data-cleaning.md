# Data cleaning

Based on fastqc report, we will have to remove the based with low quality scores, adapters from the reads before proceeding further.

First open an interactive session on HPC-class

```
sallco -N 1 -n 16 -t 4:00:00
```

Load the required modules:

```
module load trimmomatic
module load parallel
```

Since the reads we have are paired-end, we should process them together. If not, it will result in removal of pair (if they are composed of low quality bases) and cause problem while mapping (the reads cannot be mixed type: paired and single end).

Let's see what syntax `trimmomatic` command needs in order to run the trimming on our files. You can either do it online or by trying some commands to see the help page  (-h, or --help or just running the command itself). Online page can be found [here](http://www.usadellab.org/cms/?page=trimmomatic)

```
trimmomatic
```

gave us a brief help on how to use this command. Glacning at the command it looks like we will need the options `PE` (for paired end), `-threads` (to run faster), `-phred33` (quality scoring scale) ` <inputFile1>` and `<inputFile2>` (input files), `-baseout` (prefix of output files). However, the trimming parameters seems to be missing. Opening up online page, we see that by default it does:
```
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
This will perform the following:

1. Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
2. Remove leading low quality or N bases (below quality 3) (LEADING:3)
3. Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
4. Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
6. Drop reads below the 36 bases long (MINLEN:36)

Let's go with the default for now and see what it does for our files.

based on these information, our command should look something like this:

```
trimmomatic \
    PE \
    -phred33 \
    -trimlog SRR6826996_trim.log
    -threads 16 \
    SRR6826996_1.fastq.gz SRR6826996_2.fastq.gz \
    -baseout SRR6826996 \
    ILLUMINACLIP:${TRIMMOMATIC_HOME}/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
```

Run this for all files (hint you can make a slurm submission script and submit them)

or run them all using parallel command. Note that you only have 16 procs in a compute node. Running them all at a time, each using 16 processors will fail, so you have to run one job at a time. With the option `-j 1` you can make parallel to run serially

```
cat srr.ids | parallel -j 1 "trimmomatic PE -phred33 -trimlog {}_trim.log -threads 16 {}_1.fastq.gz {}_2.fastq.gz -baseout {} ILLUMINACLIP:${TRIMMOMATIC_HOME}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
```

Now, let's check the quality score again. before that, let's rename the files.

```
rename _2P _trimmed_2.fq *_2P
rename _1P _trimmed_1.fq *_1P
```


Let's run fastqc again.

```
module load fastqc
cat srr.ids | parallel "fastqc {}_trimmed_1.fq {}_trimmed_2.fq"
```

Download them to your local machine and inspect them to see if they are good.
