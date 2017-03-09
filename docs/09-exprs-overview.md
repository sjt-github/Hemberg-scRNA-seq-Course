---
output: html_document
---

# Data visualization

## Introduction

In this chapter we will continue to work with the filtered __blischak__ dataset produced in the previous chapter. We will explore different ways of visualizing the data to allow you to asses what happened to the expression matrix after the quality control step. [scater](https://github.com/davismcc/scater) package provides several very useful functions to simplify visualisation. 

One important aspect of single-cell RNA-seq is to control for batch effects. Batch effects are technical artefacts that are added to the samples during handling. For example, if two sets of samples were prepared in different labs or even on different days in the same lab, then we may observe greater similarities between the samples that were handled together. In the worst case scenario, batch effects may be [mistaken](http://f1000research.com/articles/4-121/v1) for true biological variation. The Blischak data allows us to explore these issues in a controlled manner since some of the salient aspects of how the samples were handled have been recorded. Ideally, we expect to see batches from the same individual grouping together and distinct groups corresponding to each individual. 




```r
library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
umi <- readRDS("blischak/umi.rds")
umi.qc <- umi[fData(umi)$use, pData(umi)$use]
endog_genes <- !fData(umi.qc)$is_feature_control
```

## PCA plot {#visual-pca}

The easiest thing to overview the data is to transform it using the principal component analysis and then visualize the first two principal components.

[Principal component analysis (PCA)](https://en.wikipedia.org/wiki/Principal_component_analysis) is a statistical procedure that uses a transformation to convert a set of observations into a set of values of linearly uncorrelated variables called principal components (PCs). The number of principal components is less than or equal to the number of original variables.

Mathematically, the PCs correspond to the eigenvectors of the covariance matrix. Typically, the eigenvectors are sorted by eigenvalue so that the first principal component accounts for as much of the variability in the data as possible, and each succeeding component in turn has the highest variance possible under the constraint that it is orthogonal to the preceding components (the figure below is taken from [here](http://www.nlpca.org/pca_principal_component_analysis.html)).

<div class="figure" style="text-align: center">
<img src="figures/pca.png" alt="Schematic representation of PCA dimensionality reduction" width="100%" />
<p class="caption">(\#fig:clust-pca)Schematic representation of PCA dimensionality reduction</p>
</div>

### Before QC


```r
scater::plotPCA(umi[endog_genes, ],
                ntop = 500,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "counts")
```

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-pca-before-qc-1.png" alt="PCA plot of the blischak data" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-before-qc)PCA plot of the blischak data</p>
</div>

### After QC


```r
scater::plotPCA(umi.qc[endog_genes, ],
                ntop = 500,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "counts")
```

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-pca-after-qc-1.png" alt="PCA plot of the blischak data" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-after-qc)PCA plot of the blischak data</p>
</div>

Comparing Figure \@ref(fig:expr-overview-pca-before-qc) and Figure \@ref(fig:expr-overview-pca-after-qc), it is clear that after quality control the NA19098.r2 cells no longer form a group of outliers.

By default only the top 500 most variable genes are used by scater to calculate the PCA. This can be adjusted by changing the `ntop` argument. 

__Exercise 1__
How do the PCA plots change if when all 14,214 genes are used? Or when only top 50 genes are used?

__Our answer__

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-pca-after-qc-exercise1-1-1.png" alt="PCA plot of the blischak data (14214 genes)" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-after-qc-exercise1-1)PCA plot of the blischak data (14214 genes)</p>
</div>

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-pca-after-qc-exercise1-2-1.png" alt="PCA plot of the blischak data (50 genes)" width="90%" />
<p class="caption">(\#fig:expr-overview-pca-after-qc-exercise1-2)PCA plot of the blischak data (50 genes)</p>
</div>

If your answers are different please compare your code with [ours](https://github.com/hemberg-lab/scRNA.seq.course/blob/master/07-exprs-overview.Rmd) (you need to search for this exercise in the opened file).

## tSNE map {#visual-tsne}

An alternative to PCA for visualizing scRNASeq data is a tSNE plot. [tSNE](https://lvdmaaten.github.io/tsne/) (t-Distributed Stochastic Neighbor Embedding) combines dimensionality reduction (e.g. PCA) with random walks on the nearest-neighbour network to map high dimensional data (i.e. our 14,214 dimensional expression matrix) to a 2-dimensional space while preserving local distances between cells. In contrast with PCA, tSNE is a stochastic algorithm which means running the method multiple times on the same dataset will result in different plots. Due to the non-linear and stochastic nature of the algorithm, tSNE is more difficult to intuitively interpret tSNE. To ensure reproducibility, we fix the "seed" of the random-number generator in the code below so that we always get the same plot. 


### Before QC


```r
scater::plotTSNE(umi[endog_genes, ],
                 ntop = 500,
                 perplexity = 130,
                 colour_by = "batch",
                 size_by = "total_features",
                 shape_by = "individual",
                 exprs_values = "counts",
                 rand_seed = 123456)
```

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-tsne-before-qc-1.png" alt="tSNE map of the blischak data" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-before-qc)tSNE map of the blischak data</p>
</div>

### After QC


```r
scater::plotTSNE(umi.qc[endog_genes, ],
                 ntop = 500,
                 perplexity = 130,
                 colour_by = "batch",
                 size_by = "total_features",
                 shape_by = "individual",
                 exprs_values = "counts",
                 rand_seed = 123456)
```

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-tsne-after-qc-1.png" alt="tSNE map of the blischak data" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-after-qc)tSNE map of the blischak data</p>
</div>

Interpreting PCA and tSNE plots is often challenging and due to their stochastic and non-linear nature, they are less intuitive. However, in this case it is clear that they provide a similar picture of the data. Comparing Figure \@ref(fig:expr-overview-tsne-before-qc) and \@ref(fig:expr-overview-tsne-after-qc), it is again clear that the samples from NA19098.r2 are no longer outliers after the QC filtering.

Furthermore tSNE requires you to provide a value of "perplexity" which reflects the number of neighbours used to build the nearest-neighbour network; a high value creates a dense network which clumps cells together while a low value makes the network more sparse allowing groups of cells to separate from each other. __scater__ uses a default perplexity of the total number of cells divided by five (rounded down).

You can read more about the pitfalls of using tSNE [here](http://distill.pub/2016/misread-tsne/).

__Exercise 2__
How do the tSNE plots change when a perplexity of 10 or 200 is used?

__Our answer__

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-tsne-after-qc-exercise2-1-1.png" alt="tSNE map of the blischak data (perplexity = 10)" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-after-qc-exercise2-1)tSNE map of the blischak data (perplexity = 10)</p>
</div>

<div class="figure" style="text-align: center">
<img src="09-exprs-overview_files/figure-html/expr-overview-tsne-after-qc-exercise2-2-1.png" alt="tSNE map of the blischak data (perplexity = 200)" width="90%" />
<p class="caption">(\#fig:expr-overview-tsne-after-qc-exercise2-2)tSNE map of the blischak data (perplexity = 200)</p>
</div>

If your answers are different please compare your code with [ours](https://github.com/hemberg-lab/scRNA.seq.course/blob/master/07-exprs-overview.Rmd) (you need to search for this exercise in the opened file).

## Big Exercise

Perform the same analysis with read counts of the Blischak data. Use `blischak/reads.rds` file to load the reads SCESet object. Once you have finished please compare your results to ours (next chapter).
