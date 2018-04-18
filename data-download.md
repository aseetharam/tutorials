## Downloading Data from NCBI:

Since we decided to use the  bioproject id [PRJNA437979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111730), we have to download 12 files. We will use the SRA ids to download. You should already have the SRA ids by now. If not, go the the SRA database,search for the bioproject id, go to run selector and download the table.

Next, copy the SRR files to the HPC class cluster (not in the home dir, but it should go to project temp directory).
 
1. Request a allocation (compute node)

```
salloc -N 1 -n 16 -t 8:00:00
```

Once alloted, load the module `sratoolkit` for prefetching and converting the SRA files.

2. Downlaod data

Load the module
```
module load sratoolkit
```

For downloading the file from NCBI, use the sra id to do so. The command:

```
prefetch <SRA_ID>

```
replace the `<SRA_ID>` with each of the 12 SRA ids in your list. This will download the SRA file to the `ncbi` directory you created. This is the directory that we created in the class called `ncbi`. It should be located `/ptmp/net_id/rnaseq/ncbi` (double check if this path is correct). 


Next step is to convert each of the 12 downloaded sra files to fastq files (splitting them to forward and reverse). Again, we will use the `sratoolkit` to do so (the command `fastq-dump`)

you can do this using:

```
fastq-dump --origfmt --gzip --split-files ~/ncbi/public/sra/<SRA_ID>.sra
```

replace the `<SRA_ID>` with SRA ids.


After finishing this, you will have 12 `fastq` files that are gzipped (files ending with `fastq.gz`). You should be good to run `fastqc` on them.


3. Running fastqc:

make sure you are on compute node. If you are not, run the salloc command to get one. Load the module needed for running quality check. 

```
module load fastqc
```

run it as follows:

```
fastqc <SRA_ID>_1.fastq.gz <SRA_ID>_2.fastq.gz
```

(2 files at a time, replacing <SRA_ID> with the acutal sra number, 12 times)

Try to undersntad the results and come up with conclusions.



