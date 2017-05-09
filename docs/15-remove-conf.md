---
output: html_document
---

# Dealing with confounders

## Introduction

In the previous chapter we normalized for library size, effectively removing it as a confounder. Now we will consider removing other less well defined confounders from our data. Technical confounders (aka batch effects) can arise from difference in reagents, isolation methods, the lab/experimenter who performed the experiment, even which day/time the experiment was performed. Accounting for technical confounders, and batch effects particularly, is a large topic that also involves principles of experimental design. Here we address approaches that can be taken to account for confounders when the experimental design is appropriate.

Fundamentally, accounting for technical confounders involves identifying and, ideally, removing sources of variation in the expression data that are not related to (i.e. are confounding) the biological signal of interest. Various approaches exist, some of which use spike-in or housekeeping genes, and some of which use endogenous genes.

### Advantages and disadvantages of using spike-ins to remove confounders

The use of spike-ins as control genes is appealing, since the same amount of ERCC (or other) spike-in was added to each cell in our experiment. In principle, all the variablity we observe for these genes is due to technical noise; whereas endogenous genes are affected by both technical noise and biological variability. Technical noise can be removed by fitting a model to the spike-ins and "substracting" this from the endogenous genes. There are several methods available based on this premise (eg. [BASiCS](https://github.com/catavallejos/BASiCS), [scLVM](https://github.com/PMBio/scLVM), [RUVg](http://bioconductor.org/packages/release/bioc/html/RUVSeq.html)); each using different noise models and different fitting procedures. Alternatively, one can identify genes which exhibit significant variation beyond technical noise (eg. Distance to median, [Highly variable genes](http://www.nature.com/nmeth/journal/v10/n11/full/nmeth.2645.html)). However, there are issues with the use of spike-ins for normalisation (particularly ERCCs, derived from bacterial sequences), including that their variability can, for various reasons, actually be *higher* than that of endogenous genes.

Given the issues with using spike-ins, better results can often be obtained by using endogenous genes instead. Where we have a large number of endogenous genes that, on average, do not vary systematically between cells and where we expect technical effects to affect a large number of genes (a very common and reasonable assumption), then such methods (for example, the RUVs method) can perform well. 

We explore both general approaches below.

### How to evaluate and compare confounder removal strategies

A key question when considering the different methods for removing confounders is how to quantitatively determine which one is the most effective. The main reason why comparisons are challenging is because it is often difficult to know what corresponds to technical counfounders and what is interesting biological variability. Here, we consider three different metrics which are all reasonable based on our knowledge of the experimental design. Depending on the biological question that you wish to address, it is important to choose a metric that allows you to remove the confounders that are likely to be the biggest concern for the given situation.





```r
library(scRNA.seq.funcs)
library(RUVSeq)
library(scater, quietly = TRUE)
library(scran)
library(edgeR)
set.seed(1234567)
options(stringsAsFactors = FALSE)
umi <- readRDS("tung/umi.rds")
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



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-5-1} \end{center}

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvg2") +
    ggtitle("PCA - RUVg normalisation: k = 2")
```



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-5-2} \end{center}

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs1") +
    ggtitle("PCA - RUVs normalisation: k = 1")
```



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-5-3} \end{center}

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2") +
    ggtitle("PCA - RUVs normalisation: k = 2")
```



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-5-4} \end{center}

```r
plotPCA(
    umi.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2_logcpm") +
    ggtitle("PCA - RUVs normalisation log2-cpm: k = 2")
```



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-5-5} \end{center}

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



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-6-1} \end{center}

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



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-8-1} \end{center}

```r
fit1 <- glmFit(dge1, design)
res1 <- glmLRT(fit1)
topTags(res1)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC   logCPM       LR       PValue          FDR
## ENSG00000187193 -1.8179825 7.436706 75.64199 3.400520e-18 4.782151e-14
## ENSG00000185885 -1.1909683 8.812688 62.22578 3.062525e-15 2.153414e-11
## ENSG00000125144 -3.3823146 5.590605 58.16741 2.407340e-14 1.128481e-10
## ENSG00000150459  0.7135016 8.414351 40.20293 2.289043e-10 6.467710e-07
## ENSG00000182463 -3.6331080 5.160311 40.19399 2.299548e-10 6.467710e-07
## ENSG00000008311 -1.1891979 7.398917 37.66873 8.383786e-10 1.965020e-06
## ENSG00000164265 -2.1530458 5.547251 37.26738 1.029932e-09 2.069134e-06
## ENSG00000186439 -3.1023383 5.209191 36.68009 1.391937e-09 2.371266e-06
## ENSG00000134369  1.2492359 7.112593 36.51166 1.517556e-09 2.371266e-06
## ENSG00000198417 -2.6086556 5.456231 34.46438 4.341236e-09 6.105081e-06
```

```r
summary(decideTestsDGE(res1))
```

```
##    umi.qc[, keep]$batchNA19101.r3
## -1                            125
## 0                           13875
## 1                              63
```

```r
plotSmear(
    res1, lowess = TRUE,
    de.tags = rownames(topTags(res1, n = sum(abs(decideTestsDGE(res1))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-8-2} \end{center}

### DE (RUVg, k = 2)

```r
design_ruvg <- model.matrix(~ruvg$W[keep,] + umi.qc[, keep]$batch)
head(design_ruvg)
```

```
##   (Intercept) ruvg$W[keep, ]W_1 ruvg$W[keep, ]W_2
## 1           1       0.031566414       0.028446254
## 2           1       0.008323015       0.031336158
## 3           1       0.010699708      -0.011631828
## 4           1       0.019683744       0.006157921
## 5           1       0.033731033       0.020439301
## 6           1       0.031992504       0.057216310
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



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-9-1} \end{center}

```r
fit2 <- glmFit(dge_ruvg, design_ruvg)
res2 <- glmLRT(fit2)
topTags(res2)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC   logCPM       LR       PValue          FDR
## ENSG00000187193 -1.6237017 7.434075 66.71822 3.132259e-16 4.404895e-12
## ENSG00000185885 -1.1420992 8.811712 50.87408 9.848549e-13 6.925008e-09
## ENSG00000125144 -3.3692260 5.589868 49.87533 1.638318e-12 7.679889e-09
## ENSG00000186439 -3.4569861 5.208899 35.54571 2.491356e-09 8.758984e-06
## ENSG00000134369  1.2360744 7.112423 32.43121 1.234871e-08 3.473198e-05
## ENSG00000196683 -0.4292702 9.485778 30.67591 3.049281e-08 7.147005e-05
## ENSG00000198417 -2.6320892 5.455931 29.86782 4.625256e-08 9.292139e-05
## ENSG00000196591 -0.5111855 8.779698 29.39948 5.889344e-08 1.035273e-04
## ENSG00000143570 -0.7983990 7.670748 28.97874 7.317717e-08 1.143434e-04
## ENSG00000150459  0.5770140 8.417628 28.77421 8.132704e-08 1.143702e-04
```

```r
summary(decideTestsDGE(res2))
```

```
##    umi.qc[, keep]$batchNA19101.r3
## -1                             94
## 0                           13920
## 1                              49
```

```r
plotSmear(
    res2, lowess = TRUE,
    de.tags = rownames(topTags(res2, n = sum(abs(decideTestsDGE(res2))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-9-2} \end{center}

### DE (RUVs, k = 2)

```r
design_ruvs <- model.matrix(~ruvs$W[keep,] + umi.qc[, keep]$batch)
head(design_ruvs)
```

```
##   (Intercept) ruvs$W[keep, ]W_1 ruvs$W[keep, ]W_2
## 1           1         0.2786262       -0.07496082
## 2           1         0.2825106       -0.09287973
## 3           1         0.2688981       -0.08130043
## 4           1         0.2152265       -0.09724690
## 5           1         0.2727856       -0.08117065
## 6           1         0.2351730       -0.07046587
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



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-10-1} \end{center}

```r
fit3 <- glmFit(dge_ruvs, design_ruvs)
res3 <- glmLRT(fit3)
topTags(res3)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC    logCPM       LR       PValue          FDR
## ENSG00000137818  0.3806669 11.625728 45.53881 1.496440e-11 1.554820e-07
## ENSG00000144713  0.5006206 10.588483 44.77416 2.211221e-11 1.554820e-07
## ENSG00000105372  0.3973437 11.574721 37.81839 7.764625e-10 3.639797e-06
## ENSG00000197958  0.3942581 10.207469 33.58036 6.837956e-09 2.404054e-05
## ENSG00000137154  0.3291951 11.209969 31.28073 2.232847e-08 5.913156e-05
## ENSG00000187193 -1.3710013  7.429212 31.04369 2.522857e-08 5.913156e-05
## ENSG00000181163  0.3355133 11.486981 28.79788 8.033888e-08 1.614008e-04
## ENSG00000089157  0.3079715 11.155625 27.89563 1.280385e-07 2.076654e-04
## ENSG00000150459  0.6754420  8.423888 27.82351 1.329011e-07 2.076654e-04
## ENSG00000114391  0.5009652 10.200433 27.52536 1.550481e-07 2.180441e-04
```

```r
summary(decideTestsDGE(res3))
```

```
##    umi.qc[, keep]$batchNA19101.r3
## -1                             83
## 0                           13911
## 1                              69
```

```r
plotSmear(
    res3, lowess = TRUE,
    de.tags = rownames(topTags(res3, n = sum(abs(decideTestsDGE(res3))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-10-2} \end{center}

In the above analyses, we have ignored size factors between cells. A typical edgeR analysis would always include these.


```r
umi.qc <- scran::computeSumFactors(umi.qc, sizes = 15)
dge_ruvs$samples$norm.factors <- sizeFactors(umi.qc)[keep]
dge_ruvs_sf <- estimateDisp(dge_ruvs, design = design_ruvs, trend.method = "none")
plotBCV(dge_ruvs_sf)
```



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-11-1} \end{center}

```r
fit4 <- glmFit(dge_ruvs_sf, design_ruvs)
res4 <- glmLRT(fit4)
topTags(res4)
```

```
## Coefficient:  umi.qc[, keep]$batchNA19101.r3 
##                      logFC    logCPM       LR       PValue          FDR
## ENSG00000187193 -1.3230540  7.418234 31.81968 1.691697e-08 0.0002379033
## ENSG00000144713  0.5465584 10.690938 27.50251 1.568906e-07 0.0011031766
## ENSG00000153246  4.1765109  4.644252 23.73794 1.103848e-06 0.0039359534
## ENSG00000150459  0.6934232  8.365003 23.71081 1.119520e-06 0.0039359534
## ENSG00000125144 -2.6309945  5.400342 22.94655 1.665689e-06 0.0046849159
## ENSG00000162244  0.5584180  9.918262 22.01551 2.704563e-06 0.0063390446
## ENSG00000198918  0.5993234  9.603605 21.15631 4.233047e-06 0.0084570542
## ENSG00000174748  0.4802203 10.706546 20.52616 5.882171e-06 0.0084570542
## ENSG00000181163  0.3668285 11.635978 20.39748 6.291268e-06 0.0084570542
## ENSG00000114391  0.5512702 10.275265 20.07975 7.427918e-06 0.0084570542
```

```r
summary(decideTestsDGE(res4))
```

```
##    umi.qc[, keep]$batchNA19101.r3
## -1                             43
## 0                           13973
## 1                              47
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```



\begin{center}\includegraphics{15-remove-conf_files/figure-latex/unnamed-chunk-11-2} \end{center}


## Exercise

Perform the same analysis with read counts of the `tung` data. Use `tung/reads.rds` file to load the reads SCESet object. Once you have finished please compare your results to ours (next chapter). Additionally, experiment with other combinations of normalizations and compare the results.
