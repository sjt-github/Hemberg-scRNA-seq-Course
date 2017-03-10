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
##                     logFC   logCPM       LR       PValue          FDR
## ENSG00000118855  8.621408 3.329487 36.92639 1.226744e-09 1.970274e-05
## ENSG00000050767 -6.830721 1.272193 28.89991 7.621657e-08 6.120572e-04
## ENSG00000203814 -7.790881 3.366541 27.81471 1.335072e-07 7.147528e-04
## ENSG00000171174  7.227995 2.883476 26.48678 2.653480e-07 1.065439e-03
## ENSG00000164430 -7.879346 2.147507 25.68488 4.019650e-07 1.137170e-03
## ENSG00000146856  7.561196 1.822867 25.57817 4.248193e-07 1.137170e-03
## ENSG00000145284  7.573230 1.825023 25.06133 5.553560e-07 1.274225e-03
## ENSG00000130005  5.692851 3.028824 24.60747 7.027866e-07 1.376189e-03
## ENSG00000100505 -6.917710 2.287720 24.12320 9.036509e-07 1.376189e-03
## ENSG00000146350  5.883349 2.755230 24.07409 9.269924e-07 1.376189e-03
```

```r
summary(decideTestsDGE(res1))
```

```
##    [,1] 
## -1   516
## 0  14872
## 1    673
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
## 1           1       0.030955368       0.054997999
## 2           1       0.039005432       0.008739589
## 3           1      -0.020946537      -0.020047056
## 4           1      -0.009896783       0.019871106
## 5           1       0.026969278       0.014887328
## 6           1      -0.002485334       0.039227992
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
##                     logFC      logCPM       LR       PValue          FDR
## ENSG00000118855  8.463866  3.32947862 36.34741 1.650982e-09 2.651642e-05
## ENSG00000105383 -4.442921  0.23356542 32.76664 1.039130e-08 4.332399e-05
## ENSG00000106384 -4.442921  0.23356542 32.76664 1.039130e-08 4.332399e-05
## ENSG00000261603 -5.100550  0.58341361 32.69348 1.078986e-08 4.332399e-05
## ENSG00000125337 -5.045237  0.53583923 31.37855 2.123132e-08 6.819924e-05
## ENSG00000183208 -4.610764  0.27118319 30.91130 2.700950e-08 7.229993e-05
## ENSG00000149201 -2.657566  0.39110994 30.22098 3.855169e-08 8.845410e-05
## ENSG00000197140 -6.628971  1.12370449 29.44053 5.765918e-08 1.050819e-04
## ENSG00000050767 -6.818698  1.27245990 29.39979 5.888410e-08 1.050819e-04
## ENSG00000023445 -3.746487 -0.04975573 28.30273 1.037498e-07 1.659159e-04
```

```r
summary(decideTestsDGE(res2))
```

```
##    [,1] 
## -1   550
## 0  15285
## 1    226
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
## 1           1         0.2701635         0.1806164
## 2           1         0.3437096         0.2013995
## 3           1         0.3439868         0.1636682
## 4           1         0.3176779         0.2226902
## 5           1         0.3482541         0.1820908
## 6           1         0.3760740         0.1978700
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
## ENSG00000103888 -6.285342 0.8025081 39.90423 2.667244e-10 4.283860e-06
## ENSG00000137252 -5.788134 0.5136255 37.23517 1.047084e-09 6.365277e-06
## ENSG00000127533 -5.743047 0.4897261 36.98739 1.188957e-09 6.365277e-06
## ENSG00000203780 -5.505964 0.3704316 35.66849 2.339174e-09 7.751590e-06
## ENSG00000267941 -5.411256 0.2983855 35.60782 2.413172e-09 7.751590e-06
## ENSG00000104432 -5.679853 0.4404412 34.67016 3.905735e-09 1.045500e-05
## ENSG00000118855  8.247953 3.3294767 33.57198 6.867504e-09 1.506759e-05
## ENSG00000147174 -5.210265 0.2636056 33.39928 7.505181e-09 1.506759e-05
## ENSG00000000005  8.077882 1.6975357 32.72750 1.060262e-08 1.892097e-05
## ENSG00000144649 -8.599711 2.0479901 31.47450 2.020762e-08 3.245546e-05
```

```r
summary(decideTestsDGE(res3))
```

```
##    [,1] 
## -1   317
## 0  15479
## 1    265
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
##                     logFC     logCPM       LR       PValue          FDR
## ENSG00000103888 -6.386901 0.58792252 38.70312 4.934191e-10 7.924803e-06
## ENSG00000137252 -5.889957 0.28880830 36.05362 1.919616e-09 1.165796e-05
## ENSG00000127533 -5.844999 0.26394192 35.80796 2.177565e-09 1.165796e-05
## ENSG00000203780 -5.609078 0.13950779 34.50022 4.262030e-09 1.224810e-05
## ENSG00000267941 -5.518057 0.06448738 34.42550 4.428829e-09 1.224810e-05
## ENSG00000118855  8.392195 3.31140228 34.36205 4.575592e-09 1.224810e-05
## ENSG00000104432 -5.785359 0.21290934 33.31585 7.834163e-09 1.578081e-05
## ENSG00000000005  8.487394 1.66527082 33.30934 7.860437e-09 1.578081e-05
## ENSG00000147174 -5.320627 0.02348538 31.98118 1.556733e-08 2.778076e-05
## ENSG00000144649 -8.676705 1.80776695 29.76526 4.876515e-08 7.832171e-05
```

```r
summary(decideTestsDGE(res4))
```

```
##    [,1] 
## -1   249
## 0  15539
## 1    273
```

```r
plotSmear(
    res4, lowess = TRUE,
    de.tags = rownames(topTags(res4, n = sum(abs(decideTestsDGE(res4))))$table)
)
```

<img src="16-remove-conf-reads_files/figure-html/unnamed-chunk-11-2.png" width="672" style="display: block; margin: auto;" />
