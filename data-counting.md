# Estimating gene counts using `featureCounts`

`featureCounts` is a program in the `subread` package, developed for counting reads to genomic features such as genes, exons, promoters and genomic bins. We now have the reads mapped to the genome, but we don't know what are the regions the reads are mapped and how many.

The `featureCounts` as a nice tutorial page [here](http://bioinf.wehi.edu.au/featureCounts/). You will use the GTF/GFF file that we used in the previous section to do the counting.

First, open the help page to see/read the options that are available.

```
featureCounts -h
```

_What are the parameters that you think will affect the total number of reads counted per feature?_

1. Do you include or exclude multi-mapped reads (reads mapped more than once)?
2. Do you count chimeric reads (reads split aligned to multiple genes)?
3. Do you count strand specific reads (reads can map to both strands of the DNA, if the sequencing was stranded and the gene is in one strand, reads on the other strand might not be transcripts for that gene)?
4. What do you do with the reads with only one end of the read mapped?
5. What feature do you want to count the reads for? Gene? exon? Transcript?
6. How much overlap should be there to be included them in the count?
7. Can you make it run faster by running it on multiple CPUs?


Next, open a interactive session:

```
salloc -N 1 -n 16 -t 8:00:00
module load subread
featureCounts -t exon -g gene_id -p -B -T 16 -a Mus_musculus.GRCm38.92.gtf  -o count-exon.txt SRR6826996_starAligned.out.sam
featureCounts -t gene -g gene_id -p -B -T 16 -a Mus_musculus.GRCm38.92.gtf  -o count-gene.txt SRR6826996_starAligned.out.sam
```

What is the difference between the two count files?


Run the feature counts on all the sam files
```
featureCounts -t gene -g gene_id -p -B -T 16 -a Mus_musculus.GRCm38.92.gtf  -o count_full.txt *.sam
```

# Quality plots

Once we have the counts data, we can proceed to run DGE analysis. But, before that, we need to make sure all our steps had reasonable results. We can quickly do so by checkign the results we got after each step. Since it is hard to do manually we will use a program for that called `multiqc`

`multiqc` can be found [here](http://multiqc.info/) it compiles results from multiple sources to a easy readable format that lets you check the quality of your results after each step.

```
module load python
multiqc .
```
Now transfer the file to local computer and view the results.
