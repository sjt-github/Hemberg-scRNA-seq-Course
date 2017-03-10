---
output: html_document
---

# Dealing with confounders

## Introduction

In the previous chapter we normalized for library size, effectively removing it as a confounder. Now we will consider removing other less well defined confounders from our data. Technical confounders (aka batch effects) can arise from difference in reagents, isolation methods, the lab/experimenter who performed the experiment, even which day/time the experiment was performed. Accounting for technical confounders, and batch effects particularly, is a large topic that also involves principles of experimental design. Here we address approaches that can be taken to account for confounders when the experimental design is appropriate.

Fundamentally, accounting for technical confounders involves identifying and, ideally, removing sources of variation in the expression data that are not related to (i.e. are confounding) the biological signal of interest. Various approaches exist, some of which use spike-in or housekeeping genes, and some of which use endogenous genes.

The use of spike-ins as control genes is appealing, since the same amount of ERCC (or other) spike-in was added to each cell in our experiment. In principel, all the variablity we observe for these genes is due to technical noise; whereas endogenous genes are affected by both technical noise and biological variability. Technical noise can be removed by fitting a model to the spike-ins and "substracting" this from the endogenous genes. There are several methods available based on this premise (eg. [BASiCS](https://github.com/catavallejos/BASiCS), [scLVM](https://github.com/PMBio/scLVM), [RUVg](http://bioconductor.org/packages/release/bioc/html/RUVSeq.html)); each using different noise models and different fitting procedures. Alternatively, one can identify genes which exhibit significant variation beyond technical noise (eg. Distance to median, [Highly variable genes](http://www.nature.com/nmeth/journal/v10/n11/full/nmeth.2645.html)). However, there are issues with the use of spike-ins for normalisation (particularly ERCCs, derived from bacterial sequences), including that their variability can, for various reasons, actually be *higher* than that of endogenous genes.

Given the issues with using spike-ins, better results can often be obtained by using endogenous genes instead. Where we have a large number of endogenous genes that, on average, do not vary systematically between cells and where we expect technical effects to affect a large number of genes (a very common and reasonable assumption), then such methods (for example, the RUVs method) can perform well. 

We explore both general approaches below.





```r
library(scRNA.seq.funcs)
library(RUVSeq)
library(scater, quietly = TRUE)
library(scran)
library(edgeR)
options(stringsAsFactors = FALSE)
umi <- readRDS("blischak/umi.rds")
umi.qc <- umi[fData(umi)$use, pData(umi)$use]
endog_genes <- !fData(umi.qc)$is_feature_control
erccs <- fData(umi.qc)$is_feature_control
```

## Remove Unwanted Variation

Factors contributing to technical noise frequently appear as "batch
effects" where cells processed on different days or by different
technicians systematically vary from one another. Removing technical
noise and correcting for batch effects can frequently be performed
using the same tool or slight variants on it. We will be considering
the [Remove Unwanted Variation (RUVSeq)](http://bioconductor.org/packages/RUVSeq). Briefly, RUVSeq works as follows. For $n$ samples and $J$ genes, consider the following generalized linear model (GLM), where the RNA-Seq read counts are regressed on both the known covariates of interest and unknown factors of unwanted variation:
\[\log E[Y|W,X,O] = W\alpha + X\beta + O\]
Here, $Y$ is the $n \times J$ matrix of observed gene-level read counts, $W$ is an $n \times k$ matrix corresponding to the factors of “unwanted variation” and $O$ is an $n \times J$ matrix of offsets that can either be set to zero or estimated with some other normalization procedure (such as upper-quartile normalization). The simultaneous estimation of $W$, $\alpha$, $\beta$, and $k$ is infeasible. For a given $k$, instead the following three
approaches to estimate the factors of unwanted variation $W$ are used:

* _RUVg_ uses negative control genes (e.g. ERCCs), assumed to have constant expression across samples;
* _RUVs_ uses centered (technical) replicate/negative control samples for which the covariates of interest are
constant;
* _RUVr_ uses residuals, e.g., from a first-pass GLM regression of the counts on the covariates of interest.

We will concentrate on the first two approaches.

### RUVg


```r
ruvg <- RUVg(counts(umi.qc), erccs, k = 1)
set_exprs(umi.qc, "ruvg1") <- ruvg$normalizedCounts
ruvg <- RUVg(counts(umi.qc), erccs, k = 2)
set_exprs(umi.qc, "ruvg2") <- ruvg$normalizedCounts
set_exprs(umi.qc, "ruvg2_logcpm") <- log2(t(t(ruvg$normalizedCounts) / 
                                           colSums(ruvg$normalizedCounts) * 1e6) + 1)
```

### RUVs


```r
scIdx <- matrix(-1, ncol = max(table(umi.qc$individual)), nrow = 3)
tmp <- which(umi.qc$individual == "NA19098")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(umi.qc$individual == "NA19101")
scIdx[2, 1:length(tmp)] <- tmp
tmp <- which(umi.qc$individual == "NA19239")
scIdx[3, 1:length(tmp)] <- tmp
cIdx <- rownames(umi.qc)
ruvs <- RUVs(counts(umi.qc), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
set_exprs(umi.qc, "ruvs1") <- ruvs$normalizedCounts
ruvs <- RUVs(counts(umi.qc), cIdx, k = 2, scIdx = scIdx, isLog = FALSE)
set_exprs(umi.qc, "ruvs2") <- ruvs$normalizedCounts
set_exprs(umi.qc, "ruvs2_logcpm") <- log2(t(t(ruvs$normalizedCounts) / 
                                           colSums(ruvs$normalizedCounts) * 1e6) + 1)
```

## Effectiveness 1

We evaluate the effectiveness of the normalization by inspecting the
PCA plot where colour corresponds the technical replicates and shape
corresponds to different biological samples (individuals). Separation of biological samples and
interspersed batches indicates that technical variation has been
removed. 


```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvg1") +
    ggtitle("PCA - RUVg normalisation: k = 1")
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvg2") +
    ggtitle("PCA - RUVg normalisation: k = 2")
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-5-2.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs1") +
    ggtitle("PCA - RUVs normalisation: k = 1")
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-5-3.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2") +
    ggtitle("PCA - RUVs normalisation: k = 2")
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-5-4.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2_logcpm") +
    ggtitle("PCA - RUVs normalisation log2-cpm: k = 2")
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-5-5.png" width="672" style="display: block; margin: auto;" />

Plotting log2-normalized CPM from RUVs with k = 2 looks to give the best separation of cells by individual.

## Effectiveness 2

We can also examine the effectiveness of correction using the relative log expression (RLE) across cells to confirm technical noise has been removed from the dataset.


```r
boxplot(
    list(
        "Raw counts" = calc_cell_RLE(counts(umi.qc), erccs),
        "RUVg (k = 1)" = calc_cell_RLE(assayData(umi.qc)$ruvg1, erccs),
        "RUVg (k = 2)" = calc_cell_RLE(assayData(umi.qc)$ruvg2, erccs),
        "RUVs (k = 1)" = calc_cell_RLE(assayData(umi.qc)$ruvs1, erccs),
        "RUVs (k = 2)" = calc_cell_RLE(assayData(umi.qc)$ruvs2, erccs)
    )
)
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-6-1.png" width="672" style="display: block; margin: auto;" />

## Effectiveness 3

Another way of evaluating the effectiveness of correction is to look at the differentially expressed (DE) genes among the batches of the same individual Theoretically, these batches should not differ from each other. Let's take the most promising individual (__NA19101__, whose batches are the closest to each other) and check whether it is true.

For demonstration purposes we will only use a subset of cells. You should not do that with your real dataset, though.

```r
keep <- c(
    sample(which(umi.qc$batch == "NA19101.r1"), 20), 
    sample(which(umi.qc$batch == "NA19101.r2"), 20),
    sample(which(umi.qc$batch == "NA19101.r3"), 20)
)
design <- model.matrix(~umi.qc[, keep]$batch)
```

We will use the [edgeR](http://bioconductor.org/packages/edgeR) package to calculate DE genes between plates for this particular individual. Recall that the input data for edgeR (and similar methods like DESeq2) must always be raw counts.

The particular coefficient that we test for DE in each case below tests to for genes that show a difference in expression between replicate plate 3 and replicate plate 1.

### DE (raw counts)

```r
dge1 <- DGEList(
    counts = counts(umi.qc[, keep]), 
    norm.factors = rep(1, length(keep)),
    group = umi.qc[, keep]$batch
)
dge1 <- estimateDisp(dge1, design = design, trend.method = "none")
plotBCV(dge1)
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />

```r
fit1 <- glmFit(dge1, design)
res1 <- glmLRT(fit1)
topTags(res1)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                     logFC   logCPM       LR       PValue          FDR
## ENSG00000163106 -1.536007 7.143216 74.72736 5.404195e-18 7.599919e-14
## ENSG00000187193 -1.896129 7.700405 73.12248 1.218490e-17 8.567813e-14
## ENSG00000185885 -1.167476 8.911075 60.69379 6.668004e-15 2.399627e-11
## ENSG00000125144 -3.158568 5.729870 60.64788 6.825364e-15 2.399627e-11
## ENSG00000008311 -1.264010 7.572071 48.89941 2.694318e-12 6.954757e-09
## ENSG00000131969 -1.217598 7.332086 48.71016 2.967257e-12 6.954757e-09
## ENSG00000205358 -3.262638 5.370962 43.45007 4.349150e-11 8.737442e-08
## ENSG00000125148 -2.869212 5.509874 42.93426 5.661065e-11 9.951445e-08
## ENSG00000198417 -2.447159 5.631385 41.94370 9.393969e-11 1.467860e-07
## ENSG00000170542 -2.362922 5.438883 30.15752 3.983395e-08 5.601849e-05
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   116
## 0  13875
## 1     72
```

```r
plotSmear(
    res1, lowess = TRUE,
    de.tags = rownames(topTags(res1, n = sum(abs(decideTestsDGE(res1))))$table)
)
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-8-2.png" width="672" style="display: block; margin: auto;" />

### DE (RUVg, k = 2)

```r
design_ruvg <- model.matrix(~ruvg$W[keep,] + umi.qc[, keep]$batch)
head(design_ruvg)
```

```
##   (Intercept) ruvg$W[keep, ]W_1 ruvg$W[keep, ]W_2
## 1           1        0.04086903       0.006991291
## 2           1        0.02055774       0.011285430
## 3           1        0.08859852      -0.003545212
## 4           1        0.10783009       0.018594563
## 5           1        0.03359971       0.028625429
## 6           1        0.02600381      -0.007498933
##   umi.qc[, keep]$batchNA19101.r2 umi.qc[, keep]$batchNA19101.r3
## 1                              0                              0
## 2                              0                              0
## 3                              0                              0
## 4                              0                              0
## 5                              0                              0
## 6                              0                              0
```

```r
dge_ruvg <- estimateDisp(dge1, design = design_ruvg, trend.method = "none")
plotBCV(dge_ruvg)
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

```r
fit2 <- glmFit(dge_ruvg, design_ruvg)
res2 <- glmLRT(fit2)
topTags(res2)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                     logFC   logCPM       LR       PValue          FDR
## ENSG00000187193 -2.051919 7.700185 57.88242 2.782633e-14 3.913217e-10
## ENSG00000163106 -1.498600 7.142889 49.64311 1.844170e-12 1.296728e-08
## ENSG00000185885 -1.097866 8.911132 37.76592 7.976295e-10 3.739021e-06
## ENSG00000125144 -3.020386 5.729542 34.15402 5.091785e-09 1.790144e-05
## ENSG00000131969 -1.045967 7.331250 29.53486 5.491997e-08 1.452963e-04
## ENSG00000008311 -1.153158 7.571210 29.30015 6.199087e-08 1.452963e-04
## ENSG00000205358 -3.240439 5.370816 28.31482 1.031033e-07 2.071346e-04
## ENSG00000198417 -2.595526 5.631366 27.37434 1.676411e-07 2.917155e-04
## ENSG00000177105 -3.091291 6.794804 27.16623 1.866913e-07 2.917155e-04
## ENSG00000129749  3.935091 5.143148 25.98877 3.434092e-07 4.829364e-04
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1    49
## 0  13978
## 1     36
```

```r
plotSmear(
    res2, lowess = TRUE,
    de.tags = rownames(topTags(res2, n = sum(abs(decideTestsDGE(res2))))$table)
)
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-9-2.png" width="672" style="display: block; margin: auto;" />

### DE (RUVs, k = 2)

```r
design_ruvs <- model.matrix(~ruvs$W[keep,] + umi.qc[, keep]$batch)
head(design_ruvs)
```

```
##   (Intercept) ruvs$W[keep, ]W_1 ruvs$W[keep, ]W_2
## 1           1         0.2354309       -0.08767774
## 2           1         0.2672065       -0.10296152
## 3           1         0.2177944       -0.08524802
## 4           1         0.1944238       -0.07758570
## 5           1         0.2744731       -0.07996901
## 6           1         0.2625222       -0.07040784
##   umi.qc[, keep]$batchNA19101.r2 umi.qc[, keep]$batchNA19101.r3
## 1                              0                              0
## 2                              0                              0
## 3                              0                              0
## 4                              0                              0
## 5                              0                              0
## 6                              0                              0
```

```r
dge_ruvs <- estimateDisp(dge1, design = design_ruvs, trend.method = "none")
plotBCV(dge_ruvs)
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" />

```r
fit3 <- glmFit(dge_ruvs, design_ruvs)
res3 <- glmLRT(fit3)
topTags(res3)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC    logCPM       LR       PValue          FDR
## ENSG00000187193 -2.1689109  7.699086 52.61349 4.060909e-13 5.710857e-09
## ENSG00000125144 -4.3278004  5.728014 47.34518 5.952423e-12 4.185446e-08
## ENSG00000137818  0.3210379 11.621421 44.43501 2.629429e-11 1.232589e-07
## ENSG00000198417 -3.7262197  5.631261 40.21592 2.273871e-10 7.325280e-07
## ENSG00000144713  0.5190602 10.591252 39.95077 2.604451e-10 7.325280e-07
## ENSG00000198918  0.8203125  9.516352 39.32345 3.591023e-10 7.675058e-07
## ENSG00000137154  0.3775859 11.178331 39.20259 3.820337e-10 7.675058e-07
## ENSG00000105372  0.3753038 11.587063 34.23058 4.895370e-09 8.605449e-06
## ENSG00000181163  0.2913025 11.519308 30.22039 3.856331e-08 6.025731e-05
## ENSG00000229117  0.5268648 10.138517 29.23712 6.404042e-08 7.826899e-05
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1    62
## 0  13948
## 1     53
```

```r
plotSmear(
    res3, lowess = TRUE,
    de.tags = rownames(topTags(res3, n = sum(abs(decideTestsDGE(res3))))$table)
)
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-10-2.png" width="672" style="display: block; margin: auto;" />

In the above analyses, we have ignored size factors between cells. A typical edgeR analysis would always include these.


```r
umi.qc <- scran::computeSumFactors(umi.qc, sizes = 15)
dge_ruvs$samples$norm.factors <- sizeFactors(umi.qc)[keep]
dge_ruvs_sf <- estimateDisp(dge_ruvs, design = design_ruvs, trend.method = "none")
plotBCV(dge_ruvs_sf)
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

```r
fit4 <- glmFit(dge_ruvs_sf, design_ruvs)
res4 <- glmLRT(fit4)
topTags(res4)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC    logCPM       LR       PValue          FDR
## ENSG00000187193 -2.0355550  7.781858 50.10075 1.460513e-12 1.994689e-08
## ENSG00000198918  0.9588144  9.692246 48.79835 2.836791e-12 1.994689e-08
## ENSG00000125144 -4.1034864  5.710389 46.51964 9.070472e-12 4.251935e-08
## ENSG00000144713  0.6538171 10.788988 36.73516 1.353164e-09 4.757388e-06
## ENSG00000198417 -3.5189479  5.599718 36.20493 1.776208e-09 4.995764e-06
## ENSG00000177954  0.6213097 11.281223 33.66262 6.554836e-09 1.491619e-05
## ENSG00000229117  0.6555321 10.332680 33.42025 7.424681e-09 1.491619e-05
## ENSG00000240972  0.6325741 10.550878 30.73824 2.952887e-08 4.831587e-05
## ENSG00000162244  0.6708129 10.030581 30.64885 3.092106e-08 4.831587e-05
## ENSG00000145425  0.6450358 10.534072 30.02201 4.271698e-08 6.007289e-05
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1    26
## 0  13960
## 1     77
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```

<img src="15-remove-conf_files/figure-html/unnamed-chunk-11-2.png" width="672" style="display: block; margin: auto;" />


## Exercise

Perform the same analysis with read counts of the Blischak data. Use `blischak/reads.rds` file to load the reads SCESet object. Once you have finished please compare your results to ours (next chapter). Additionally, experiment with other combinations of normalizations and compare the results.
