---
knit: bookdown::preview_chapter
---

# Construction of expression matrix

## Reads QC

The output from a scRNA-seq experiment is a large collection of short
cDNA reads. The first step is to ensure that the reads are of high
quality. In application to scRNA-seq, this can be performed by using standard bulk RNA-seq QC tools, such as [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) or [Kraken](http://www.ebi.ac.uk/research/enright/software/kraken). Currently, there are no specialized tools for scRNA-seq
available, so we use
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

Assuming that our reads are in experiment.bam, we run FastQC as
```
$<path_to_fastQC>/fastQC experiment.bam
```

Below is an example of the output from FastQC for a dataset of 125 bp
reads. Here, the quality of the reads is overall high, so we can
proceed the analysis with confidence.

![](figures/per_base_quality.png)

Additionally, the data can be visualized using the [Integrative Genomics Browser (IGV)](https://www.broadinstitute.org/igv/).

## Reads alignment

Once the low-quality reads have been removed, the remaining reads can
be mapped to a reference genome. Again, there are no special purpose
methods for this, so we can use the
[STAR](https://github.com/alexdobin/STAR) or the [TopHat](https://ccb.jhu.edu/software/tophat/index.shtml) aligner.

An example of how to map reads.bam to using STAR hg38 is

```
$<path_to_STAR>/STAR --runThreadN 1 --runMode alignReads
--readFilesIn reads1.fq.gz reads2.fq.gz --readFilesCommand zcat --genomeDir <path>
--parametersFiles FileOfMoreParameters.txt --outFileNamePrefix <outpath>/output
```

__Note__, if the _spike-ins_ are used, the reference sequence should be augmented with the DNA sequence of the _spike-in_ molecules before mapping.

__Note__, when UMIs are used, their barcodes should be removed from every read.

Once we have mapped the reads for each cell to the reference genome,
we need to make sure that a sufficient number of reads from each cell
could be mapped to the reference genome. In our experience, the
fraction of mappable reads is 60-70%. However, this result may vary
depending on protocol, read length and settings for the read
alignment. As a general rule, we expect all cells to have a similar
fraction of mapped reads, so any outliers should be inspected and
possibly removed.

## Alignment example

The histogram below shows the total number of reads mapped to each
cell for an scRNA-seq experiment. Each bar represents one cell, and
they have been sorted in ascending order by the total number of reads
per cell. The three red arrows indicate cells that are outliers in
terms of their coverage and they should be removed from further
analysis. The two yellow arrows point to cells with a surprisingly
large number of unmapped reads. However, we deem the discrepancy as
small and we retain the cells for now.

![](figures/Bergiers_exp1_mapping_by_cell.png)

## Mapping QC

After mapping the raw sequencing to the genome we should evaluate the quality of the mapping. There are many ways to measure this including: amount of reads mapping to rRNA/tRNAs, proportion of uniquely mapping reads, reads mapping across splice junctions, read depth along the transcripts. Methods developed for bulk RNA-seq, such as [RSeQC](http://rseqc.sourceforge.net/), are applicable to single-cell data:

```
python <RSeQCpath>/geneBody_coverage.py -i input.bam -r genome.bed -o output.txt
python <RSeQCpath>/bam_stat.py -i input.bam -r genome.bed -o output.txt
python <RSeQCpath>/split_bam.py -i input.bam -r rRNAmask.bed -o output.txt
```

However the expected results will depend on the experimental protocol, e.g. many scRNA-seq methods use poly-A selection to avoid sequencing rRNAs which results in a 3' bias in the read coverage across the genes (aka gene body coverage). The figure below shows this 3' bias as well as three cells which were outliers and removed from the dataset:

![](figures/Exp1_RSEQC_geneBodyCoverage_plot_Combined.png)

## Reads quantification

The next step is to quantify the expression level of each gene for
each cell. For mRNA data, we can use one of the tools which has been
developed for bulk RNA-seq data, e.g. [HT-seq](http://www-huber.embl.de/users/anders/HTSeq/) or [FeatureCounts](http://subread.sourceforge.net/)

```
# include multimapping
<featureCounts_path>/featureCounts -O -M -Q 30 -p -a genome.gtf -o outputfile input.bam
# exclude multimapping
<featureCounts_path>/featureCounts -Q 30 -p -a genome.gtf -o outputfile input.bam
```

__Note__, when UMIs are used, the expression counts can be collapsed by summing the number of unique barcodes associated with all reads mapped to a given gene.
