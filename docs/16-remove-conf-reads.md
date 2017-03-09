---
output: html_document
---

# Dealing with confounders (Reads)




```r
library(scRNA.seq.funcs)
library(RUVSeq)
library(scater, quietly = TRUE)
library(scran)
library(edgeR)
options(stringsAsFactors = FALSE)
reads <- readRDS("blischak/reads.rds")
reads.qc <- reads[fData(reads)$use, pData(reads)$use]
endog_genes <- !fData(reads.qc)$is_feature_control
erccs <- fData(reads.qc)$is_feature_control
```

## Remove Unwanted Variation

### RUVg


```r
ruvg <- RUVg(counts(reads.qc), erccs, k = 1)
set_exprs(reads.qc, "ruvg1") <- ruvg$normalizedCounts
ruvg <- RUVg(counts(reads.qc), erccs, k = 2)
set_exprs(reads.qc, "ruvg2") <- ruvg$normalizedCounts
set_exprs(reads.qc, "ruvg2_logcpm") <- log2(t(t(ruvg$normalizedCounts) / 
                                           colSums(ruvg$normalizedCounts)) + 1)
```

### RUVs


```r
scIdx <- matrix(-1, ncol = max(table(reads.qc$individual)), nrow = 3)
tmp <- which(reads.qc$individual == "NA19098")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(reads.qc$individual == "NA19101")
scIdx[2, 1:length(tmp)] <- tmp
tmp <- which(reads.qc$individual == "NA19239")
scIdx[3, 1:length(tmp)] <- tmp
cIdx <- rownames(reads.qc)
ruvs <- RUVs(counts(reads.qc), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
set_exprs(reads.qc, "ruvs1") <- ruvs$normalizedCounts
ruvs <- RUVs(counts(reads.qc), cIdx, k = 2, scIdx = scIdx, isLog = FALSE)
set_exprs(reads.qc, "ruvs2") <- ruvs$normalizedCounts
set_exprs(reads.qc, "ruvs2_logcpm") <- log2(t(t(ruvs$normalizedCounts) / 
                                           colSums(ruvs$normalizedCounts)) + 1)
```

## Effectiveness 1


```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvg1") +
    ggtitle("PCA - RUVg normalisation: k = 1")
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvg2") +
    ggtitle("PCA - RUVg normalisation: k = 2")
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-5-2.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs1") +
    ggtitle("PCA - RUVs normalisation: k = 1")
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-5-3.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2") +
    ggtitle("PCA - RUVs normalisation: k = 2")
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-5-4.png" width="672" style="display: block; margin: auto;" />

```r
plotPCA(
    reads.qc[endog_genes, ],
    colour_by = "batch",
    size_by = "total_features",
    shape_by = "individual",
    exprs_values = "ruvs2_logcpm") +
    ggtitle("PCA - RUVs normalisation log2-cpm: k = 2")
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-5-5.png" width="672" style="display: block; margin: auto;" />

## Effectiveness 2


```r
boxplot(
    list(
        "Raw counts" = calc_cell_RLE(counts(reads.qc), erccs),
        "RUVg (k = 1)" = calc_cell_RLE(assayData(reads.qc)$ruvg1, erccs),
        "RUVg (k = 2)" = calc_cell_RLE(assayData(reads.qc)$ruvg2, erccs),
        "RUVs (k = 1)" = calc_cell_RLE(assayData(reads.qc)$ruvs1, erccs),
        "RUVs (k = 2)" = calc_cell_RLE(assayData(reads.qc)$ruvs2, erccs)
    )
)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-6-1.png" width="672" style="display: block; margin: auto;" />

## Effectiveness 3


```r
keep <- c(
    sample(which(reads.qc$batch == "NA19101.r1"), 20), 
    sample(which(reads.qc$batch == "NA19101.r2"), 20),
    sample(which(reads.qc$batch == "NA19101.r3"), 20)
)
design <- model.matrix(~reads.qc[, keep]$batch)
```

### DE (raw counts)

```r
dge1 <- DGEList(
    counts = counts(reads.qc[, keep]), 
    norm.factors = rep(1, length(keep)),
    group = reads.qc[, keep]$batch
)
dge1 <- estimateDisp(dge1, design = design, trend.method = "none")
plotBCV(dge1)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />

```r
fit1 <- glmFit(dge1, design)
res1 <- glmLRT(fit1)
topTags(res1)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000182463 -9.573354 3.6640827 32.23511 1.366001e-08 0.0002193934
## ENSG00000049246  6.401912 3.2783678 27.10222 1.929762e-07 0.0011677345
## ENSG00000174652  7.081533 3.0280497 26.63709 2.454856e-07 0.0011677345
## ENSG00000109743  5.937359 0.7640719 25.92383 3.551585e-07 0.0011677345
## ENSG00000117480  7.555217 2.7930678 25.87885 3.635311e-07 0.0011677345
## ENSG00000131378 -6.525665 1.8117480 24.56998 7.165935e-07 0.0016621738
## ENSG00000187398  6.698311 1.2553328 24.50748 7.402190e-07 0.0016621738
## ENSG00000135905  6.504905 1.1138480 24.29175 8.279304e-07 0.0016621738
## ENSG00000182240  7.489052 2.5682191 23.80076 1.068403e-06 0.0017394207
## ENSG00000168016  6.619726 2.0662337 23.77462 1.083009e-06 0.0017394207
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   436
## 0  14988
## 1    637
```

```r
plotSmear(
    res1, lowess = TRUE,
    de.tags = rownames(topTags(res1, n = sum(abs(decideTestsDGE(res1))))$table)
)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-8-2.png" width="672" style="display: block; margin: auto;" />

### DE (RUVg, k = 2)

```r
design_ruvg <- model.matrix(~ruvg$W[keep,] + reads.qc[, keep]$batch)
head(design_ruvg)
```

```
##   (Intercept) ruvg$W[keep, ]W_1 ruvg$W[keep, ]W_2
## 1           1       -0.05126269       0.022687305
## 2           1        0.03417440      -0.016840208
## 3           1        0.04848977      -0.007775879
## 4           1        0.02696928       0.014887328
## 5           1        0.03000176       0.049402874
## 6           1        0.01814940      -0.063270368
##   reads.qc[, keep]$batchNA19101.r2 reads.qc[, keep]$batchNA19101.r3
## 1                                0                                0
## 2                                0                                0
## 3                                0                                0
## 4                                0                                0
## 5                                0                                0
## 6                                0                                0
```

```r
dge_ruvg <- estimateDisp(dge1, design = design_ruvg, trend.method = "none")
plotBCV(dge_ruvg)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

```r
fit2 <- glmFit(dge_ruvg, design_ruvg)
res2 <- glmLRT(fit2)
topTags(res2)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000182463 -7.219321 3.6640107 30.57915 3.205196e-08 0.0005147865
## ENSG00000049246  6.279866 3.2782361 25.74985 3.886569e-07 0.0017113716
## ENSG00000174652  6.963708 3.0279039 25.06331 5.547838e-07 0.0017113716
## ENSG00000178860 -6.461578 2.0367363 24.97424 5.810149e-07 0.0017113716
## ENSG00000117480  7.460653 2.7927804 24.85675 6.175256e-07 0.0017113716
## ENSG00000131378 -6.684604 1.8118141 24.78771 6.400423e-07 0.0017113716
## ENSG00000101951 -6.200653 1.0357779 24.47627 7.523062e-07 0.0017113716
## ENSG00000165140  7.399741 1.5105880 24.16201 8.856202e-07 0.0017113716
## ENSG00000262601 -5.726100 0.8111904 23.83838 1.047723e-06 0.0017113716
## ENSG00000137478  8.080566 2.4664688 23.55534 1.213722e-06 0.0017113716
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   374
## 0  15465
## 1    222
```

```r
plotSmear(
    res2, lowess = TRUE,
    de.tags = rownames(topTags(res2, n = sum(abs(decideTestsDGE(res2))))$table)
)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-9-2.png" width="672" style="display: block; margin: auto;" />

### DE (RUVs, k = 2)

```r
design_ruvs <- model.matrix(~ruvs$W[keep,] + reads.qc[, keep]$batch)
head(design_ruvs)
```

```
##   (Intercept) ruvs$W[keep, ]W_1 ruvs$W[keep, ]W_2
## 1           1         0.3411501         0.1898753
## 2           1         0.2864326         0.1671973
## 3           1         0.3147846         0.1785769
## 4           1         0.3482541         0.1820908
## 5           1         0.3425186         0.2093401
## 6           1         0.2702873         0.2169784
##   reads.qc[, keep]$batchNA19101.r2 reads.qc[, keep]$batchNA19101.r3
## 1                                0                                0
## 2                                0                                0
## 3                                0                                0
## 4                                0                                0
## 5                                0                                0
## 6                                0                                0
```

```r
dge_ruvs <- estimateDisp(dge1, design = design_ruvs, trend.method = "none")
plotBCV(dge_ruvs)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" />

```r
fit3 <- glmFit(dge_ruvs, design_ruvs)
res3 <- glmLRT(fit3)
topTags(res3)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000182240  8.699558 2.5680302 29.89979 4.549630e-08 0.0007307161
## ENSG00000167964  5.148967 0.7827290 27.26449 1.774400e-07 0.0007643465
## ENSG00000133110 -6.214904 0.9607716 26.16169 3.139897e-07 0.0007643465
## ENSG00000183690  9.112379 2.3911197 26.13672 3.180771e-07 0.0007643465
## ENSG00000107742  4.375768 0.1625397 25.86890 3.654095e-07 0.0007643465
## ENSG00000213965 -7.175267 2.0430530 25.46287 4.509801e-07 0.0007643465
## ENSG00000108231  6.519664 1.4512848 25.20939 5.143123e-07 0.0007643465
## ENSG00000131378 -7.005807 1.8118380 25.05468 5.572745e-07 0.0007643465
## ENSG00000168016  6.965939 2.0663181 25.05088 5.583734e-07 0.0007643465
## ENSG00000115468 -7.326474 1.7147165 24.88247 6.093389e-07 0.0007643465
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   322
## 0  15345
## 1    394
```

```r
plotSmear(
    res3, lowess = TRUE,
    de.tags = rownames(topTags(res3, n = sum(abs(decideTestsDGE(res3))))$table)
)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-10-2.png" width="672" style="display: block; margin: auto;" />


```r
reads.qc <- scran::computeSumFactors(reads.qc, sizes = 15)
dge_ruvs$samples$norm.factors <- sizeFactors(reads.qc)[keep]
dge_ruvs_sf <- estimateDisp(dge_ruvs, design = design_ruvs, trend.method = "none")
plotBCV(dge_ruvs_sf)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

```r
fit4 <- glmFit(dge_ruvs_sf, design_ruvs)
res4 <- glmLRT(fit4)
topTags(res4)
```

```
## Coefficient:  reads.qc[, keep]$batchNA19101.r3 
##                     logFC       logCPM       LR       PValue          FDR
## ENSG00000182240  9.124483 2.7335107311 29.40509 5.872314e-08 0.0009431523
## ENSG00000167964  5.203390 0.7595403429 26.46423 2.684640e-07 0.0010687770
## ENSG00000107742  4.561050 0.0002506229 26.30103 2.921331e-07 0.0010687770
## ENSG00000049246  6.263140 3.5143251970 25.91424 3.569270e-07 0.0010687770
## ENSG00000183690  9.134500 2.4689087459 25.16456 5.264079e-07 0.0010687770
## ENSG00000172828  5.635551 0.6286632236 25.12629 5.369589e-07 0.0010687770
## ENSG00000133110 -6.177699 0.7491222321 25.01096 5.700532e-07 0.0010687770
## ENSG00000168016  7.150925 2.0408255132 24.83088 6.258664e-07 0.0010687770
## ENSG00000108231  6.675345 1.4507250717 24.82778 6.268752e-07 0.0010687770
## ENSG00000156509  6.710909 1.2009230966 24.69732 6.707722e-07 0.0010687770
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   265
## 0  15407
## 1    389
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-11-2.png" width="672" style="display: block; margin: auto;" />
