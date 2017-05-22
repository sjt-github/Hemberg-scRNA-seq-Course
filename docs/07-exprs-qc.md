---
output: html_document
---

# Expression QC (UMI) {#exprs-qc}

## Introduction

Once gene expression has been quantified it is summarized as an __expression matrix__ where each row corresponds to a gene (or transcript) and each column corresponds to a single cell. This matrix should be examined to remove poor quality cells which were not detected in either read QC or mapping QC steps. Failure to remove low quality cells at this
stage may add technical noise which has the potential to obscure
the biological signals of interest in the downstream analysis. 

Since there is currently no standard method for performing scRNASeq the expected values for the various QC measures that will be presented here can vary substantially from experiment to experiment. Thus, to perform QC we will be looking for cells which are outliers with respect to the rest of the dataset rather than comparing to independent quality standards. Consequently, care should be taken when comparing quality metrics across datasets collected using different protocols.


## Tung dataset

To illustrate cell QC, we consider a
[dataset](http://jdblischak.github.io/singleCellSeq/analysis/) of
 induced pluripotent stem cells generated from three different individuals [@Tung2017-ba] in [Yoav Gilad](http://giladlab.uchicago.edu/)'s lab at the
University of Chicago. The experiments were carried out on the
Fluidigm C1 platform and to facilitate the quantification both unique
molecular identifiers (UMIs) and ERCC _spike-ins_ were used. The data files are located in the `tung` folder in your working directory. These files are the copies of the original files made on the 15/03/16. We will use these copies for reproducibility purposes.




```r
library(scater, quietly = TRUE)
library(knitr)
options(stringsAsFactors = FALSE)
```

Load the data and annotations:

```r
molecules <- read.table("tung/molecules.txt", sep = "\t")
anno <- read.table("tung/annotation.txt", sep = "\t", header = TRUE)
```

Inspect a small portion of the expression matrix

```r
knitr::kable(
    head(molecules[ , 1:3]), booktabs = TRUE,
    caption = 'A table of the first 6 rows and 3 columns of the molecules table.'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-4)A table of the first 6 rows and 3 columns of the molecules table.}
\centering
\begin{tabular}[t]{lrrr}
\toprule
  & NA19098.r1.A01 & NA19098.r1.A02 & NA19098.r1.A03\\
\midrule
ENSG00000237683 & 0 & 0 & 0\\
ENSG00000187634 & 0 & 0 & 0\\
ENSG00000188976 & 3 & 6 & 1\\
ENSG00000187961 & 0 & 0 & 0\\
ENSG00000187583 & 0 & 0 & 0\\
ENSG00000187642 & 0 & 0 & 0\\
\bottomrule
\end{tabular}
\end{table}

```r
knitr::kable(
    head(anno), booktabs = TRUE,
    caption = 'A table of the first 6 rows of the anno table.'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-4)A table of the first 6 rows of the anno table.}
\centering
\begin{tabular}[t]{lllll}
\toprule
individual & replicate & well & batch & sample\_id\\
\midrule
NA19098 & r1 & A01 & NA19098.r1 & NA19098.r1.A01\\
NA19098 & r1 & A02 & NA19098.r1 & NA19098.r1.A02\\
NA19098 & r1 & A03 & NA19098.r1 & NA19098.r1.A03\\
NA19098 & r1 & A04 & NA19098.r1 & NA19098.r1.A04\\
NA19098 & r1 & A05 & NA19098.r1 & NA19098.r1.A05\\
NA19098 & r1 & A06 & NA19098.r1 & NA19098.r1.A06\\
\bottomrule
\end{tabular}
\end{table}

The data consists of 3 individuals and 3 replicates and therefore has 9 batches in total.

We standardize the analysis by using the scater package. First, create the scater SCESet classes:

```r
pheno_data <- new("AnnotatedDataFrame", anno)
rownames(pheno_data) <- pheno_data$sample_id
umi <- scater::newSCESet(
    countData = molecules,
    phenoData = pheno_data
)
```

Remove genes that are not expressed in any cell:

```r
keep_feature <- rowSums(counts(umi) > 0) > 0
umi <- umi[keep_feature, ]
```

Define control features (genes) - ERCC spike-ins and mitochondrial genes ([provided](http://jdblischak.github.io/singleCellSeq/analysis/qc-filter-ipsc.html) by the authors):

```r
ercc <- featureNames(umi)[grepl("ERCC-", featureNames(umi))]
mt <- c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
        "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
        "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
        "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
        "ENSG00000198840")
```

Calculate the quality metrics:

```r
umi <- scater::calculateQCMetrics(
    umi,
    feature_controls = list(ERCC = ercc, MT = mt)
)
```


## Cell QC

### Library size

Next we consider the total number of RNA molecules detected per
sample (if we were using read counts rather than UMI counts this would
be the total number of reads). Wells with few reads/molecules are likely to have
been broken or failed to capture a cell, and should thus be removed.


```r
hist(
    umi$total_counts,
    breaks = 100
)
abline(v = 25000, col = "red")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-qc_files/figure-latex/total-counts-hist-1} 

}

\caption{Histogram of library sizes for all cells}(\#fig:total-counts-hist)
\end{figure}

__Exercise 1__

1. How many cells does our filter remove?

2. What distribution do you expect that the
total number of molecules for each cell should follow?

__Our answer__

\begin{table}

\caption{(\#tab:unnamed-chunk-9)The number of cells removed by total counts filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_by\_total\_counts & Freq\\
\midrule
FALSE & 46\\
TRUE & 818\\
\bottomrule
\end{tabular}
\end{table}

### Detected genes (1)

In addition to ensuring sufficient sequencing depth for each sample, we also want to make sure that the reads are distributed across the transcriptome. Thus, we count the total number of unique genes detected in each sample.


```r
hist(
    umi$total_features,
    breaks = 100
)
abline(v = 7000, col = "red")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-qc_files/figure-latex/total-features-hist-1} 

}

\caption{Histogram of the number of detected genes in all cells}(\#fig:total-features-hist)
\end{figure}

From the plot we conclude that most cells have between 7,000-10,000 detected genes,
which is normal for high-depth scRNA-seq. However, this varies by
experimental protocol and sequencing depth. For example, droplet-based methods
or samples with lower sequencing-depth typically detect fewer genes per cell. The most notable feature in the above plot is the __"heavy tail"__ on the left hand side of the
distribution. If detection rates were equal across the cells then the
distribution should be approximately normal. Thus we remove those
cells in the tail of the distribution (fewer than 7,000 detected genes).

__Exercise 2__

How many cells does our filter remove?

__Our answer__

\begin{table}

\caption{(\#tab:unnamed-chunk-10)The number of cells removed by total features filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_by\_expr\_features & Freq\\
\midrule
FALSE & 120\\
TRUE & 744\\
\bottomrule
\end{tabular}
\end{table}

### ERCCs and MTs

Another measure of cell quality is the ratio between ERCC _spike-in_
RNAs and endogenous RNAs. This ratio can be used to estimate the total amount
of RNA in the captured cells. Cells with a high level of _spike-in_ RNAs
had low starting amounts of RNA, likely due to the cell being
dead or stressed which may result in the RNA being degraded.


```r
scater::plotPhenoData(
    umi,
    aes_string(x = "total_features",
               y = "pct_counts_feature_controls_MT",
               colour = "batch")
)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-qc_files/figure-latex/mt-vs-counts-1} 

}

\caption{Percentage of counts in MT genes}(\#fig:mt-vs-counts)
\end{figure}


```r
scater::plotPhenoData(
    umi,
    aes_string(x = "total_features",
               y = "pct_counts_feature_controls_ERCC",
               colour = "batch")
)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-qc_files/figure-latex/ercc-vs-counts-1} 

}

\caption{Percentage of counts in ERCCs}(\#fig:ercc-vs-counts)
\end{figure}

The above analysis shows that majority of the cells from NA19098.r2 batch have a very high ERCC/Endo ratio. Indeed, it has been shown by the authors that this batch contains cells of smaller size. 

__Exercise 3__

Create filters for removing batch NA19098.r2 and cells with high expression of mitochondrial genes (>10% of total counts in a cell).

__Our answer__

\begin{table}

\caption{(\#tab:unnamed-chunk-11)The number of cells removed by ERCC filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_by\_ERCC & Freq\\
\midrule
FALSE & 96\\
TRUE & 768\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}

\caption{(\#tab:unnamed-chunk-11)The number of cells removed by MT filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_by\_MT & Freq\\
\midrule
FALSE & 31\\
TRUE & 833\\
\bottomrule
\end{tabular}
\end{table}

__Exercise 4__

What would you expect to see in the ERCC vs counts plot if you were examining a dataset containing cells of different sizes (eg. normal & senescent cells)?

__Answer__

You would expect to see a group corresponding to the smaller cells (normal) with a higher fraction of ERCC reads than a separate group corresponding to the larger cells (senescent).

## Cell filtering

### Manual

Now we can define a cell filter based on our previous analysis:


```r
umi$use <- (
    # sufficient features (genes)
    filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # sufficient endogenous RNA
    filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT
)
```


```r
knitr::kable(
  as.data.frame(table(umi$use)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by manual filter (FALSE)'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-13)The number of cells removed by manual filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
Var1 & Freq\\
\midrule
FALSE & 210\\
TRUE & 654\\
\bottomrule
\end{tabular}
\end{table}

### Default Thresholds

Results from the biological analysis of single-cell RNA-seq data are often strongly influenced by outliers. Thus, it is important to include multiple filters. A robust way of detecting outliers is through the median absolute difference which is defined as $d_i = |r_i - m|$, where $r_i$ is the number of reads in cell $i$ and $m$ is the median number of reads across the all cells in the sample. By default, scater removes all cells where $d_i>5*$median($d_i$). A similar filter is used for the number of detected genes. Furthermore, scater removes all cells where >80% of counts were assigned to control genes and any cells that have been marked as controls.


```r
umi$use_default <- (
    # remove cells with unusual numbers of genes
    !umi$filter_on_total_features &
    # remove cells with unusual numbers molecules counted
    !umi$filter_on_total_counts &
    # < 80% ERCC spike-in
    !umi$filter_on_pct_counts_feature_controls_ERCC &
    # < 80% mitochondrial
    !umi$filter_on_pct_counts_feature_controls_MT &
    # controls shouldn't be used in downstream analysis
    !umi$is_cell_control
)
```


```r
knitr::kable(
  as.data.frame(table(umi$use_default)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by default filter (FALSE)'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-15)The number of cells removed by default filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
Var1 & Freq\\
\midrule
FALSE & 6\\
TRUE & 858\\
\bottomrule
\end{tabular}
\end{table}

### Automatic

Another option available in __scater__ is to conduct PCA on a set of QC metrics and then use automatic outlier detection to identify potentially problematic cells. 

By default, the following metrics are used for PCA-based outlier detection:

* __pct_counts_top_100_features__
* __total_features__
* __pct_counts_feature_controls__
* __n_detected_feature_controls__
* __log10_counts_endogenous_features__
* __log10_counts_feature_controls__

scater first creates a matrix where the rows represent cells and the columns represent the different QC metrics. Here, the PCA plot provides a 2D representation of cells ordered by their quality metrics. The outliers are then detected using methods from the mvoutlier package. 


```r
umi <-
scater::plotPCA(umi,
                size_by = "total_features", 
                shape_by = "use",
                pca_data_input = "pdata",
                detect_outliers = TRUE,
                return_SCESet = TRUE)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-qc_files/figure-latex/auto-cell-filt-1} 

}

\caption{PCA plot used for automatic detection of cell outliers}(\#fig:auto-cell-filt)
\end{figure}


```r
knitr::kable(
  as.data.frame(table(umi$outlier)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by automatic filter (FALSE)'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-16)The number of cells removed by automatic filter (FALSE)}
\centering
\begin{tabular}[t]{r}
\toprule
Freq\\


\bottomrule
\end{tabular}
\end{table}

## Compare filterings

__Exercise 5__

Compare the default, automatic and manual cell filters. Plot a Venn diagram of the outlier cells from these filterings.

__Hint__: Use `limma::vennCounts` and `limma::vennDiagram` functions from the [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) package to make a Venn diagram.

__Answer__

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-qc_files/figure-latex/cell-filt-comp-1} 

}

\caption{Comparison of the default, automatic and manual cell filters}(\#fig:cell-filt-comp)
\end{figure}

## Gene analysis

### Gene expression

In addition to removing cells with poor quality, it is usually a good idea to exclude genes where we suspect that technical artefacts may have skewed the results. Moreover, inspection of the gene expression profiles may provide insights about how the experimental procedures could be improved.

It is often instructive to consider the number of reads consumed by the top 50 expressed genes.


```r
scater::plotQC(umi, type = "highest-expression")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{07-exprs-qc_files/figure-latex/top50-gene-expr-1} 

}

\caption{Number of total counts consumed by the top 50 expressed genes}(\#fig:top50-gene-expr)
\end{figure}

The distributions are relatively flat indicating (but not guaranteeing!) good coverage of the full transcriptome of these cells. However, there are several spike-ins in the top 15 genes which suggests a greater dilution of the spike-ins may be preferrable if the experiment is to be repeated.


### Gene filtering

It is typically a good idea to remove genes whose expression level is considered __"undetectable"__. We define a gene as  detectable if at least two cells contain more than 1 transcript from the gene. If we were considering read counts rather than UMI counts a reasonable threshold is to require at least five reads in at least two cells. However, in both cases the threshold strongly depends on the sequencing depth. It is important to keep in mind that genes must be filtered after cell filtering since some genes may only be detected in poor quality cells (__note__ `pData(umi)$use` filter applied to the `umi` dataset).


```r
filter_genes <- apply(counts(umi[ , pData(umi)$use]), 1, 
                      function(x) length(x[x > 1]) >= 2)
fData(umi)$use <- filter_genes
```


```r
knitr::kable(
    as.data.frame(table(filter_genes)),
    booktabs = TRUE,
    row.names = FALSE,
    caption = 'The number of genes removed by gene filter (FALSE)'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-18)The number of genes removed by gene filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_genes & Freq\\
\midrule
FALSE & 4663\\
TRUE & 14063\\
\bottomrule
\end{tabular}
\end{table}

Depending on the cell-type, protocol and sequencing depth, other cut-offs may be appropriate.


## Save the data

Dimensions of the QCed dataset (do not forget about the gene filter we defined above):

```r
dim(umi[fData(umi)$use, pData(umi)$use])
```

```
## Features  Samples 
##    14063      654
```

Let's create an additional slot with log-transformed counts (we will need it in the next chapters):

```r
set_exprs(umi, "log2_counts") <- log2(counts(umi) + 1)
```

Save the data:

```r
saveRDS(umi, file = "tung/umi.rds")
```

## Big Exercise

Perform exactly the same QC analysis with read counts of the same Blischak data. Use `tung/reads.txt` file to load the reads. Once you have finished please compare your results to ours (next chapter).
