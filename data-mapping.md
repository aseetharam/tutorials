# Mapping reads to the reference genome


Let's find the latest reference genome along with its feature file (GFF3 or GTF files describing genes and other features on the reference genome). Popular choices to download the genomes include [UCSC](http://hgdownload.cse.ucsc.edu/downloads.html), [NCBI](https://www.ncbi.nlm.nih.gov/genome/), [ensembl](http://ensemblgenomes.org/) etc.

Let's use the `[ensemblgenomes](http://useast.ensembl.org/Mus_musculus/Info/Index)` for our project. You're welcome to try other genomes as well, and there is no rational behind selecting one source over the other.

```
wget ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/Mus_musculus.GRCm38.92.gtf.gz
```

Many different aligners are available to map the RNAseq reads to the reference genome. We will use STAR program as it is very accurate and fast than other programs such as HiSat2, GSNAP [see ref](https://www.nature.com/articles/nmeth.2722).

Before mapping, you need to create index files for your genomes. Open up the manual page and see how it is done. You can download a pdf copy of manual from here: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf


For indexing files:

```
module load star
mkdir mm10_star
STAR \
  --runMode genomeGenerate \
  --runThreadN 16 \
  --genomeDir mm10_star \
  --genomeFastaFiles Mus_musculus.GRCm38.dna.toplevel.fa \
  --sjdbGTFfile Mus_musculus.GRCm38.92.gtf \
  --sjdbOverhang 74
  --limitGenomeGenerateRAM=33524399488
```

Once this is done, we can map each pair of fastq files (trimmed) to the refernece genome (indexed)

Again, look up the manual and see the command that you need to run for mapping fastq files in pairs.


```
STAR \
 --runMode alignReads \
 --runThreadN 16 \
 --genomeDir mm10_star \
 --readFilesCommand zcat \
 --outFileNamePrefix OUTPURNAME_star \
 --readFilesIn INPUT_1.fq INPUT_2.fq
```

Since we need to run this for all 12 files, we will create a submission script and run them all separately.

First, let's create a run script:

```
nano runSTAR.sh
```

paste the following contents in the file you're editing:

```
#!/bin/bash
R1="$1"
R2="$2"
OUT=$(basename ${R1} |cut -f 1 -d "_");
GFF="Mus_musculus.GRCm38.92.gtf"
DB="mm10_star"
STAR \
 --runMode alignReads \
 --runThreadN 16 \
 --genomeDir ${DB} \
 --readFilesCommand zcat \
 --outFileNamePrefix ${OUT}_star \
--readFilesIn ${R1} ${R2}
```

Now, if you run this as:
```
./runSTAR.sh SRR6827002_trimmed_1.fq SRR6827002_trimmed_2.fq
```
it will execute  the STAR aligner with all the options that you put in that run script. It is easier to run like this on many files (12 in this case)
