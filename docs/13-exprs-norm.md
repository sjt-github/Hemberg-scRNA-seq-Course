---
knit: bookdown::preview_chapter
---

# Normalization for library size

## Introduction



In the previous chapter we identified important confounding factors and explanatory variables. scater allows one to account for these variables in subsequent statistical models or to condition them out using `normaliseExprs()`, if so desired. This can be done by providing a design matrix to `normaliseExprs()`. We are not covering this topic here, but you can try to do it yourself as an exercise.

Instead we will explore how simple size-factor normalisations correcting for library size can remove the effects of some of the confounders and explanatory variables.

## Library size

Library sizes vary because scRNA-seq data is often sequenced on highly multiplexed platforms the total reads which are derived from each cell may differ substantially. Some quantification methods
(eg. [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/), [RSEM](http://deweylab.github.io/RSEM/)) incorporated library size when determining gene expression estimates thus do not require this normalization.

However, if another quantification method was used then library size must be corrected for by multiplying or dividing each column of the expression matrix by a "normalization factor" which is an estimate of the library size relative to the other cells. Many methods to correct for library size have been developped for bulk RNA-seq and can be equally applied to scRNA-seq (eg. UQ, SF, CPM, RPKM, FPKM, TPM). 


## Normalisations

The simplest way to normalize this data is to convert it to counts per
million (__CPM__) by dividing each column by its total then multiplying by
1,000,000. Note that spike-ins should be excluded from the
calculation of total expression in order to correct for total cell RNA
content, therefore we will only use endogenous genes. 


```r
calc_cpm <-
function (expr_mat, spikes = NULL) 
{
    norm_factor <- colSums(expr_mat[-spikes, ])
    return(t(t(expr_mat)/norm_factor)) * 10^6
}
```

One potential drawback of __CPM__ is if your sample contains genes that are both very highly expressed and differentially expressed across the cells. In this case, the total molecules in the cell may depend of whether such genes are on/off in the cell and normalizing by total molecules may hide the differential expression of those genes and/or falsely create differential expression for the remaining genes. 

__Note__: RPKM, FPKM and TPM are variants on CPM which further adjust counts by the length of the respective gene/transcript.

To deal with this potentiality several other measures were devised:

The __size factor (SF)__ was proposed and popularized by DESeq ([Anders and Huber (2010)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106)). First the geometric mean of each gene across all cells is calculated. The size factor for each cell is the median across genes of the ratio of the expression to the gene's geometric mean. A drawback to this method is that since it uses the geometric mean only genes with non-zero expression values across all cells can be used in its calculation, making it unadvisable for large low-depth scRNASeq experiments. edgeR & scater call this method __RLE__ for "relative log expression".


```r
calc_sf <-
function (expr_mat, spikes = NULL) 
{
    geomeans <- exp(rowMeans(log(expr_mat[-spikes, ])))
    SF <- function(cnts) {
        median((cnts/geomeans)[(is.finite(geomeans) & geomeans > 
            0)])
    }
    norm_factor <- apply(expr_mat[-spikes, ], 2, SF)
    return(t(t(expr_mat)/norm_factor))
}
```

The __upperquartile (UQ)__ was proposed by [Bullard et al (2010)](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-94). Here each column is divided by the 75% quantile of the counts for each library. Often the calculated quantile is scaled by the median across cells to keep the absolute level of expression relatively consistent. A drawback to this method is that for low-depth scRNASeq experiments the large number of undetected genes may result in the 75% quantile being zero (or close to it). This limitation can be overcome by generalizing the idea and using a higher quantile (eg. the 99% quantile is the default in scater) or by excluding zeros prior to calculating the 75% quantile.


```r
calc_uq <-
function (expr_mat, spikes = NULL) 
{
    UQ <- function(x) {
        quantile(x[x > 0], 0.75)
    }
    uq <- unlist(apply(expr_mat[-spikes, ], 2, UQ))
    norm_factor <- uq/median(uq)
    return(t(t(expr_mat)/norm_factor))
}
```

Another method is called __TMM__ is the weighted trimmed mean of M-values (to the reference) proposed by [Robinson and Oshlack (2010)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25). The M-values in question are the gene-wise log2 fold changes between individual cells. One cell is used as the reference then the M-values for each other cell is calculated compared  to this reference. These values are then trimmed by removing the top and bottom ~30%, and the average of the remaining values is calculated by weighting them to account for the effect of the log scale on variance. Each non-reference cell is multiplied by the calculated factor. Two potential issues with this method are insufficient non-zero genes left after trimming, and the assumption that most genes are not differentially expressed.

Finally the `scran` package implements a variant on CPM specialized for single-cell data ([Lun et al 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7)). Briefly this method deals with the problem of vary large numbers of zero values per cell by pooling cells together calculating a normalization factor (similar to CPM) for the sum of each pool. Since each cell is found in many different pools, cell-specific factors can be deconvoluted from the collection of pool-specific factors using linear algebra. 


We will use visual inspection of PCA plots and calculation of cell-wise relative log expression (`calc_cell_RLE()`) to compare the efficiency of different normalization methods. Cells with many[few] reads have higher[lower] than median expression for most genes resulting in a positive[negative] RLE across the cell, whereas normalized cells have an RLE close to zero.


```r
calc_cell_RLE <-
function (expr_mat, spikes = NULL) 
{
    RLE_gene <- function(x) {
        if (median(unlist(x)) > 0) {
            log((x + 1)/(median(unlist(x)) + 1))/log(2)
        }
        else {
            rep(NA, times = length(x))
        }
    }
    if (!is.null(spikes)) {
        RLE_matrix <- t(apply(expr_mat[-spikes, ], 1, RLE_gene))
    }
    else {
        RLE_matrix <- t(apply(expr_mat, 1, RLE_gene))
    }
    cell_RLE <- apply(RLE_matrix, 2, median, na.rm = T)
    return(cell_RLE)
}
```

The __RLE__, __TMM__, and __UQ__ size-factor methods were developed for bulk RNA-seq data and, depending on the experimental context, may not be appropriate for single-cell RNA-seq data, as their underlying assumptions may be problematically violated. [Lun et al (2016)](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7) recently published a size-factor normalisation method specifically designed for scRNA-seq data and accounting for single-cell biases, which we will call __LSF__ (Lun Sum Factors). Briefly, expression values are summed across pools of cells and the summed values are used to compute normalization size-factors per pool. The pool-based size factors can then be deconvolved into cell-specific size factors, which can be used for normalization in the same way as other size factors. 

 __scater__ acts as a wrapper for the `calcNormFactors` function from [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) which implements several library size normalization methods making it easy to apply any of these methods to our data. The __LSF__ method is implementated in the Bioconductor package [scran](https://bioconductor.org/packages/release/bioc/html/scran.html), which allows seamless integrationinto the `scater` workflow. (The `scran` package itself depends on `scater`).

__Note:__ edgeR makes extra adjustments to some of the normalization methods which may result in somewhat different results than if the original methods are followed exactly, e.g. edgeR's and scater's "RLE" method which is based on the "size factor" used by [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html) may give different results to the `estimateSizeFactorsForMatrix` method in the DESeq/DESeq2 packages. In addition, some versions of edgeR will not calculate the normalization factors correctly unless `lib.size` is set at 1 for all cells.

We will continue to work with the `tung` data that was used in the previous chapter.


```r
library(scRNA.seq.funcs)
library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
set.seed(1234567)
umi <- readRDS("tung/umi.rds")
umi.qc <- umi[fData(umi)$use, pData(umi)$use]
endog_genes <- !fData(umi.qc)$is_feature_control
```

### Raw

```r
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-pca-raw-1} 

}

\caption{PCA plot of the tung data}(\#fig:norm-pca-raw)
\end{figure}


```r
boxplot(calc_cell_RLE(counts(umi.qc[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim=c(-1,1))
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-ours-rle-raw-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-raw)
\end{figure}

### CPM
scater performs this normalisation by default, you can control it by changing `exprs_values` parameter to `"exprs"`.

```r
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "exprs")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-pca-cpm-1} 

}

\caption{PCA plot of the tung data after CPM normalisation}(\#fig:norm-pca-cpm)
\end{figure}

```r
boxplot(calc_cell_RLE(exprs(umi.qc[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim = c(-1,1))
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-ours-rle-cpm-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-cpm)
\end{figure}


### TMM

```r
umi.qc <- 
    scater::normaliseExprs(umi.qc,
                           method = "TMM",
                           feature_set = endog_genes)
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "norm_counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-pca-tmm-1} 

}

\caption{PCA plot of the tung data after TMM normalisation}(\#fig:norm-pca-tmm)
\end{figure}

```r
boxplot(calc_cell_RLE(norm_counts(umi.qc[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim=c(-1,1))
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-ours-rle-tmm-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-tmm)
\end{figure}

### scran

```r
qclust <- scran::quickCluster(umi.qc, min.size = 30)
umi.qc <- scran::computeSumFactors(umi.qc, sizes = 15, clusters = qclust)
umi.qc <- scater::normalize(umi.qc)
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "exprs")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-pca-lsf-1} 

}

\caption{PCA plot of the tung data after LSF normalisation}(\#fig:norm-pca-lsf)
\end{figure}

```r
boxplot(calc_cell_RLE(exprs(umi.qc[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim=c(-1,1))
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-ours-rle-scran-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-scran)
\end{figure}

### Size-factor (RLE)

```r
umi.qc <- 
    scater::normaliseExprs(umi.qc,
                           method = "RLE", 
                           feature_set = endog_genes)
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "norm_counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-pca-rle-1} 

}

\caption{PCA plot of the tung data after RLE normalisation}(\#fig:norm-pca-rle)
\end{figure}

```r
boxplot(calc_cell_RLE(norm_counts(umi.qc[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim=c(-1,1))
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-ours-rle-rle-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-rle)
\end{figure}


### Upperquantile

```r
umi.qc <- 
    scater::normaliseExprs(umi.qc,
                           method = "upperquartile", 
                           feature_set = endog_genes,
                           p = 0.99)
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "norm_counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-pca-uq-1} 

}

\caption{PCA plot of the tung data after UQ normalisation}(\#fig:norm-pca-uq)
\end{figure}

```r
boxplot(calc_cell_RLE(norm_counts(umi.qc[endog_genes, ])),
        col = "grey50",
        ylab = "RLE",
        main = "", ylim = c(-1, 1))
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-ours-rle-uq-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-uq)
\end{figure}
## Downsampling 

A final way to correct for library size is to downsample the expression matrix so that each cell has approximately the same total number of molecules. The benefit of this method is that zero values will be introduced by the down sampling thus eliminating any biases due to differing numbers of detected genes. However, the major drawback is that the process is not deterministic so each time the downsampling is run the resulting expression matrix is slightly different. Thus, often analyses must be run on multiple downsamplings to ensure results are robust.


```r
Down_Sample_Matrix <-
function (expr_mat) 
{
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
        prob <- min_lib_size/sum(x)
        return(unlist(lapply(x, function(y) {
            rbinom(1, y, prob)
        })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
}
```


```r
norm_counts(umi.qc) <- 
    scRNA.seq.funcs::Down_Sample_Matrix(counts(umi.qc))
scater::plotPCA(umi.qc[endog_genes, ],
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "norm_counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-pca-downsample-1} 

}

\caption{PCA plot of the tung data after downsampling}(\#fig:norm-pca-downsample)
\end{figure}

```r
tmp <- norm_counts(umi.qc[endog_genes, ])
# ignore genes which are not detected in any cells following downsampling
boxplot(calc_cell_RLE(tmp[rowMeans(tmp) > 0, ]), 
        col = "grey50",
        ylab = "RLE",
        main = "", ylim = c(-1, 1))
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-ours-rle-downsample-1} 

}

\caption{Cell-wise RLE of the tung data}(\#fig:norm-ours-rle-downsample)
\end{figure}

## Normalizing for gene/transcript length

Some methods combine library size and fragment/gene length normalization such as:

* __RPKM__ - Reads Per Kilobase Million (for single-end sequencing)
* __FPKM__ - Fragments Per Kilobase Million (same as __RPKM__ but for paired-end sequencing, makes sure that paired ends mapped to the same fragment are not counted twice)
* __TPM__ - Transcripts Per Kilobase Million (same as __RPKM__, but the order of normalizations is reversed - length first and sequencing depth second)

These methods are not applicable to our dataset since the end
of the transcript which contains the UMI was preferentially
sequenced. Furthermore in general these should only be calculated
using appropriate quantification software from aligned BAM files not
from read counts since often only a portion of the entire
gene/transcript is sequenced, not the entire length. If in doubt check 
for a relationship between gene/transcript length and expression level.

However, here we show how these normalisations can be calculated using scater. First, we need to find the effective transcript length in Kilobases. However, our dataset containes only gene IDs, therefore we will be using the gene lengths instead of transcripts. scater uses the [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) package, which allows one to annotate genes by other attributes:

```r
umi.qc <-
    scater::getBMFeatureAnnos(umi.qc,
                              filters = "ensembl_gene_id", 
                              attributes = c("ensembl_gene_id",
                                             "hgnc_symbol",
                                             "chromosome_name",
                                             "start_position",
                                             "end_position"), 
                              feature_symbol = "hgnc_symbol",
                              feature_id = "ensembl_gene_id",
                              biomart = "ENSEMBL_MART_ENSEMBL", 
                              dataset = "hsapiens_gene_ensembl",
                              host = "www.ensembl.org")

# If you have mouse data, change the arguments based on this example:
# scater::getBMFeatureAnnos(object,
#                           filters = "ensembl_transcript_id", 
#                           attributes = c("ensembl_transcript_id", 
#                                        "ensembl_gene_id", "mgi_symbol", 
#                                        "chromosome_name",
#                                        "transcript_biotype",
#                                        "transcript_start",
#                                        "transcript_end", 
#                                        "transcript_count"), 
#                           feature_symbol = "mgi_symbol",
#                           feature_id = "ensembl_gene_id",
#                           biomart = "ENSEMBL_MART_ENSEMBL", 
#                           dataset = "mmusculus_gene_ensembl",
#                           host = "www.ensembl.org") 
```

Some of the genes were not annotated, therefore we filter them out:

```r
umi.qc.ann <-
    umi.qc[!is.na(fData(umi.qc)$ensembl_gene_id), ]
```

Now we compute the total gene length in Kilobases by using the `end_position` and `start_position` fields:

```r
eff_length <- abs(fData(umi.qc.ann)$end_position -
                      fData(umi.qc.ann)$start_position)/1000
```


```r
plot(eff_length, rowMeans(counts(umi.qc.ann)))
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/length-vs-mean-1} 

}

\caption{Gene length vs Mean Expression for the raw data}(\#fig:length-vs-mean)
\end{figure}
There is no relationship between gene length and mean expression so FPKMs & TPMs are inappropriate for this dataset. 
But we will demonstrate them anyway.

Now we are ready to perform the normalisations:

```r
tpm(umi.qc.ann) <-
    calculateTPM(
        umi.qc.ann,
        eff_length
    )
fpkm(umi.qc.ann) <-
    calculateFPKM(
        umi.qc.ann,
        eff_length
    )
```

Plot the results as a PCA plot:

```r
scater::plotPCA(umi.qc.ann,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "fpkm")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-pca-fpkm-1} 

}

\caption{PCA plot of the tung data after FPKM normalisation}(\#fig:norm-pca-fpkm)
\end{figure}


```r
scater::plotPCA(umi.qc.ann,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "tpm")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{13-exprs-norm_files/figure-latex/norm-pca-tpm-1} 

}

\caption{PCA plot of the tung data after TPM normalisation}(\#fig:norm-pca-tpm)
\end{figure}

__Note:__ The PCA looks for differences between cells. Gene length is the same across cells for each gene thus FPKM is almost identical to the CPM plot (it is just rotated) since it performs CPM first then normalizes gene length. Whereas, TPM is different because it weights genes by their length before performing CPM. 



## Exercise

Perform the same analysis with read counts of the `tung` data. Use `tung/reads.rds` file to load the reads SCESet object. Once you have finished please compare your results to ours (next chapter).
