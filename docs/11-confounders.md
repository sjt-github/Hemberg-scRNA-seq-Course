---
knit: bookdown::preview_chapter
---

# Confounding factors

## Introduction

There is a large number of potential confounders, artifacts and biases in sc-RNA-seq data. One of the main challenges in analyzing scRNA-seq data stems from the fact that it is difficult to carry out a true technical replicate (why?) to distinguish biological and technical variability. In the previous chapters we considered batch effects and in this chapter we will continue to explore how experimental artifacts can be identified and removed. We will continue using the [scater](https://github.com/davismcc/scater) package since it provides a set of methods specifically for quality control of experimental and explanatory variables. Moreover, we will continue to work with the Blischak data that was used in the previous chapter.




```r
library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
umi <- readRDS("blischak/umi.rds")
umi.qc <- umi[fData(umi)$use, pData(umi)$use]
endog_genes <- !fData(umi.qc)$is_feature_control
```

The `umi.qc` dataset contains filtered cells and genes. Our next step is to explore technical drivers of variability in the data to inform data normalisation before downstream analysis.

## Correlations with PCs

Let's first look again at the PCA plot of the QCed dataset:

```r
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                exprs_values = "counts")
```

<div class="figure" style="text-align: center">
<img src="11-confounders_files/figure-html/confound-pca-1.png" alt="PCA plot of the blischak data" width="90%" />
<p class="caption">(\#fig:confound-pca)PCA plot of the blischak data</p>
</div>

scater allows one to identify principle components that correlate with experimental and QC variables of interest (it ranks principle components by $R^2$ from a linear model regressing PC value against the variable of interest).

Let's test whether some of the variables correlate with any of the PCs.

### Detected genes


```r
scater::plotQC(umi.qc[endog_genes, ],
               type = "find-pcs",
               variable = "total_features",
               exprs_values = "counts")
```

<div class="figure" style="text-align: center">
<img src="11-confounders_files/figure-html/confound-find-pcs-total-features-1.png" alt="PC correlation with the number of detected genes" width="90%" />
<p class="caption">(\#fig:confound-find-pcs-total-features)PC correlation with the number of detected genes</p>
</div>

Indeed, we can see that PC1 can be completely explained by the number of the detected genes. In fact, it was also visible on the PCA plot above. This is a well-known issue in scRNA-seq and was described [here](http://biorxiv.org/content/early/2015/12/27/025528).

## Explanatory variables

scater can also compute the marginal $R^2$ for each variable when fitting a linear model regressing expression values for each gene against just that variable, and display a density plot of the gene-wise marginal $R^2$ values for the variables.


```r
scater::plotQC(umi.qc[endog_genes, ],
               type = "expl",
               exprs_values = "counts",
               variables = c("total_features",
                             "total_counts",
                             "batch",
                             "individual",
                             "pct_counts_feature_controls_ERCC",
                             "pct_counts_feature_controls_MT"))
```

<div class="figure" style="text-align: center">
<img src="11-confounders_files/figure-html/confound-find-expl-vars-1.png" alt="Explanatory variables" width="90%" />
<p class="caption">(\#fig:confound-find-expl-vars)Explanatory variables</p>
</div>

This analysis indicates that the number of detected genes (again) and also the sequencing depth (number of counts) have substantial explanatory power for many genes, so these variables are good candidates for conditioning out in a normalisation step, or including in downstream statistical models. Expression of ERCCs also appears to be an important explanatory variable.

## Other confounders

In addition to correcting for batch, there are other factors that one
may want to compensate for. As with batch correction, these
adjustments require extrinsic information. One popular method is
[scLVM](https://github.com/PMBio/scLVM) which allows you to identify
and subtract the effect from processes such as cell-cycle or
apoptosis.

In addition, protocols may differ in terms of their coverage of each transcript, 
their bias based on the average content of __A/T__ nucleotides, or their ability to capture short transcripts.
Ideally, we would like to compensate for all of these differences and biases.

## Exercise

Perform the same analysis with read counts of the Blischak data. Use `blischak/reads.rds` file to load the reads SCESet object. Once you have finished please compare your results to ours (next chapter).
