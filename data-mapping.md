# Mapping reads to the reference genome


Let's find the latest reference genome along with its feature file (GFF3 or GTF files describing genes and other features on the reference genome). Popular choices to download the genomes include [UCSC](http://hgdownload.cse.ucsc.edu/downloads.html), [NCBI](https://www.ncbi.nlm.nih.gov/genome/), [ensembl](http://ensemblgenomes.org/) etc.

Let's use the [ensemblgenomes](http://useast.ensembl.org/Mus_musculus/Info/Index) for our project. You're welcome to try other genomes as well, and there is no rational behind selecting one source over the other.

```
wget ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/Mus_musculus.GRCm38.92.gtf.gz
```

Many different aligners are available to map the RNAseq reads to the reference genome. We will use STAR program as it is very accurate and fast than other programs such as HiSat2, GSNAP [see ref](https://www.nature.com/articles/nmeth.2722).

Before mapping, you need to create index files for your genomes. Open up the manual page and see how it is done. You can download a pdf copy of manual from here: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf


For indexing files:

```
salloc -N 16 -n 1 -t 8:00:00
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

Once this is done, we can map each pair of `fastq` files (trimmed) to the reference genome (indexed)

Again, look up the [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) and see the command that you need to run for mapping `fastq` files in pairs.

An example command you need to run, looks like this
```
STAR \
 --runMode alignReads \
 --runThreadN 16 \
 --genomeDir mm10_star \
 --readFilesCommand zcat \
 --outFileNamePrefix OUTPURNAME_star \
 --readFilesIn INPUT_1.fq INPUT_2.fq
```
`runMode` tells that we are running STAR for aligning the reads, `runThreadN` is to make the run faster (use more CPUs), `genomeDir` is where the index files are located (from the previous command), `readFilesCommand` tells that the files we are providing as input are compressed (gz) and it needs to be uncompressed before aligning,  `outFileNamePrefix` prefix for the output files, and finally, the input files (2, forward and reverse).

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
chmod +x runSTAR.sh
./runSTAR.sh SRR6827002_trimmed_1.fq SRR6827002_trimmed_2.fq
```
it will execute  the STAR aligner with all the options that you put in that run script. It is easier to run like this on many files (12 in this case)


Copy the slurm script from your home directory:

```
cp ~/template.sub ./star.sub
```
edit this file to add the runSTAR commands.

```
nano star.sub
```

```
#!/bin/bash
#SBATCH -J star
#SBATCH -N 1
#SBATCH  -n 16
#SBATCH -t 24:00:00
#SBATCH -o star.stdout
#SBATCH -e star.stderr
#SBATCH --mail-user=username@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
module load star
./runSTAR.sh SRR6827002_trimmed_1.fq SRR6827002_trimmed_2.fq
..
..
..
(12 such commands)
```

run the sbatch script:
```
sbatch star.sub
```
