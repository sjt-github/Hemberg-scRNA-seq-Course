---
knit: bookdown::preview_chapter
---

# Introduction to single-cell RNA-seq

## Bulk RNA-seq

* A major breakthrough (replaced microarrays) in the late 00's and has been widely used since
* Provides an __average expression level__ for each gene from a large population of input cells
* Useful for comparative transcriptomics, e.g. samples of the same tissue from different species
* Useful for quantifying expression signatures from ensembles, e.g. in disease studies
* __Insufficient__ for studying heterogeneous systems, e.g. early development studies, complex tissues (brain)
* Does __not__ provide insights into the stochastic nature of gene expression
    
## scRNA-seq

* A __new__ technology, first publication by [Tang et al](http://www.nature.com/nmeth/journal/v6/n5/abs/nmeth.1315.html) in 2009
* Instead of providing an average of expression of a population of cells, scRNA-seq provides a __distribution of expression levels__ from a population of cells
* Allows to study new biological questions in which __cell-specific changes in transcriptome are important__, e.g. cell type identification, inference of gene regulatory networks across the cells, stochastic component of transcription
* Datasets range __from $10^2$ to $10^4$ cells__ and increase in size every year
* Currently there are several different protocols in use, e.g. [SMART-seq2](http://www.nature.com/nmeth/journal/v10/n11/full/nmeth.2639.html), [CELL-seq](http://www.cell.com/cell-reports/abstract/S2211-1247%2812%2900228-8) and [Drop-seq](http://mccarrolllab.com/dropseq/)
* Several computational analysis methods from bulk RNA-seq __can__ be used
* __In most cases__ requires adaptation of the existing methods or development of new ones

## Protocol

![](figures/RNA-Seq_workflow-5.pdf.jpg)

image from [Wikipedia - Single cell sequencing](https://en.wikipedia.org/wiki/Single_cell_sequencing)

Overall, experimental scRNA-seq protocols are similar to the methods used for bulk RNA-seq. For a discussion on experimental methods, please see reviews by [Saliba et al](http://nar.oxfordjournals.org/content/42/14/8845), [Handley et al](http://www.sciencedirect.com/science/article/pii/S1097276515003068)  or [Kolodziejczyk et al](http://www.sciencedirect.com/science/article/pii/S1097276515002610).

## Analysis

This course is concerned with the computational analysis of the data
obtained from scRNA-seq experiments. Even though the format of the
data is the same as for bulk RNA-seq, the fact that individual cells
are sampled means that the data __should__ be analyzed differently. We
will provide an overview of how to process a scRNA-seq sample and how
to analyze the data to provide biological insights.

![](figures/flowchart.png)

Flowchart for analyzing scRNA-seq data.

## Challenges

The main difference between bulk and single cell RNA-seq is that each sequencing library represents a single cell, instead of a population of cells. Therefore, significant attention has to be paid to comparison of the results from different cells (sequencing libraries). The main sources of discrepancy between the libraries are:

* __Amplification__ (up to 1 million fold)
* __Gene 'dropouts'__ in which a gene is observed at a moderate expression level in one cell but is not detected in another cell ([Kharchenko et al](http://www.nature.com/nmeth/journal/v11/n7/full/nmeth.2967.html)).

In both cases the discrepancies are introduced due to low starting transcript amount (from one cell only) and have yet to be addressed by improving the experimental protocol.

## Controls

To address challenges in technical variation between scRNA sequencing libraries two quantitative standards were introduced. They allow to normalize gene expression levels across different cells.

### ERCC _spike-ins_

One standard, strongly recommended for all scRNA-seq experiments is to use extrinsic _spike-in_ molecules. These molecules are added to the lysate before the reverse transcription. The most popular and widely used artificial spike-ins are from [External RNA Control Consortium (ERCC)](https://www.thermofisher.com/order/catalog/product/4456740). It contains 92 synthetic spikes based on bacterial sequences. Normalization using _spike-ins_ is based on the fact that the number of molecules of each _spike-in_ RNA species should be the same across all single-cell libraries.

### UMIs

Another method of standartisation is to use [Unique Molecular Identifiers (UMIs)](http://www.nature.com/nmeth/journal/v9/n1/full/nmeth.1778.html). Instead of sequencing the amplified reads from a cell library, it allows for sequencing reads derived solely from 3' or 5' end of the amplified transcript. UMIs are added as barcodes to the individual RNA molecules. This approach provides an estimate of the number of transcripts that is independent of amplification biases.
