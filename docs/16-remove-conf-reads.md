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
## ENSG00000166848  6.810057 3.6187405 35.08835 3.150796e-09 5.060493e-05
## ENSG00000183036 -7.525802 1.8885497 27.85308 1.308853e-07 7.620103e-04
## ENSG00000196597  5.640177 1.3306809 27.38280 1.669092e-07 7.620103e-04
## ENSG00000196981  6.760194 1.3219345 26.91059 2.130878e-07 7.620103e-04
## ENSG00000007516 -7.164159 1.6148472 26.70324 2.372238e-07 7.620103e-04
## ENSG00000185019  6.591430 1.2005208 26.15535 3.150226e-07 8.212195e-04
## ENSG00000166016  6.782962 2.7234501 25.90888 3.579189e-07 8.212195e-04
## ENSG00000206530 -6.318844 1.0335869 25.35382 4.772063e-07 9.580512e-04
## ENSG00000183826  6.591499 1.2340493 25.04430 5.602814e-07 9.998533e-04
## ENSG00000175646  5.952694 0.8025416 24.76954 6.461053e-07 1.022525e-03
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   365
## 0  14853
## 1    843
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
## 1           1      -0.020147219       0.118496638
## 2           1       0.024597369      -0.001815257
## 3           1      -0.007980392       0.018899126
## 4           1       0.028875292      -0.021003497
## 5           1      -0.011047886      -0.053796508
## 6           1      -0.020946537      -0.020047056
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
## ENSG00000166848  6.883979 3.6187198 34.03516 5.412508e-09 0.0000869303
## ENSG00000160097 -6.570645 1.0250226 31.59050 1.903582e-08 0.0001038310
## ENSG00000117586 -6.354790 1.1574387 31.55427 1.939438e-08 0.0001038310
## ENSG00000178385 -6.688025 1.4716597 30.23779 3.821897e-08 0.0001534587
## ENSG00000251201 -6.065053 0.9278176 29.50366 5.581112e-08 0.0001792765
## ENSG00000185019  6.302586 1.2004061 28.87139 7.734683e-08 0.0002070446
## ENSG00000183036 -7.642065 1.8884464 27.56470 1.519263e-07 0.0003485841
## ENSG00000171695 -5.305964 0.4506830 27.19101 1.843140e-07 0.0003504195
## ENSG00000188483 -5.939006 1.2956822 27.06859 1.963624e-07 0.0003504195
## ENSG00000114638  6.112550 1.0086310 26.85588 2.192067e-07 0.0003520678
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   382
## 0  15332
## 1    347
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
## 1           1         0.2961611         0.1507961
## 2           1         0.3257709         0.1987312
## 3           1         0.3521279         0.2034467
## 4           1         0.3689274         0.2143278
## 5           1         0.2576765         0.1992972
## 6           1         0.3439868         0.1636682
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
## ENSG00000088386 -6.721947 0.8992321 45.79918 1.310182e-11 1.367919e-07
## ENSG00000103888 -6.635174 0.8484026 45.28510 1.703404e-11 1.367919e-07
## ENSG00000007516 -7.677010 1.6150320 42.14232 8.486670e-11 4.145172e-07
## ENSG00000127533 -6.061633 0.5457230 41.75919 1.032357e-10 4.145172e-07
## ENSG00000225556 -5.878064 0.4619444 40.48819 1.978077e-10 6.353978e-07
## ENSG00000203780 -5.806261 0.4307903 39.98647 2.557287e-10 6.845430e-07
## ENSG00000131126 -6.767955 0.9393334 39.07228 4.084024e-10 9.370501e-07
## ENSG00000104432 -6.106352 0.5137518 38.20535 6.367766e-10 1.278409e-06
## ENSG00000183423 -6.656951 0.9306423 37.30989 1.007724e-09 1.788404e-06
## ENSG00000137252 -6.081055 0.5755341 37.11523 1.113507e-09 1.788404e-06
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   354
## 0  15332
## 1    375
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
##                     logFC    logCPM       LR       PValue          FDR
## ENSG00000088386 -6.668977 0.7464544 43.70711 3.813843e-11 4.027362e-07
## ENSG00000103888 -6.580943 0.6961582 43.17130 5.015083e-11 4.027362e-07
## ENSG00000127533 -6.003116 0.3966470 39.65297 3.033427e-10 1.337849e-06
## ENSG00000007516 -7.575707 1.5275829 39.46969 3.331918e-10 1.337849e-06
## ENSG00000225556 -5.818520 0.3137224 38.38878 5.796459e-10 1.861938e-06
## ENSG00000203780 -5.746332 0.2828807 37.88908 7.488314e-10 2.004497e-06
## ENSG00000131126 -6.716682 0.7863491 37.06479 1.142683e-09 2.621804e-06
## ENSG00000104432 -6.037617 0.3688667 35.98346 1.989991e-09 3.956199e-06
## ENSG00000166848  7.120329 3.6527279 35.77307 2.216910e-09 3.956199e-06
## ENSG00000137252 -6.018899 0.4293032 35.02451 3.255807e-09 5.229151e-06
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   257
## 0  15413
## 1    391
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-11-2.png" width="672" style="display: block; margin: auto;" />
