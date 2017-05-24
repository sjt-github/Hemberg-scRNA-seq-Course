---
output: html_document
---

# Expression QC (Reads)

This chapter contains the summary plots and tables for the QC exercise based on the reads for the Bischak data discussed in the previous chapter.


```r
library(scater, quietly = TRUE)
library(knitr)
options(stringsAsFactors = FALSE)
```




```r
reads <- read.table("tung/reads.txt", sep = "\t")
anno <- read.table("tung/annotation.txt", sep = "\t", header = TRUE)
```


```r
knitr::kable(
    head(reads[ , 1:3]), booktabs = TRUE,
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
ENSG00000188976 & 57 & 140 & 1\\
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


```r
pheno_data <- new("AnnotatedDataFrame", anno)
rownames(pheno_data) <- pheno_data$sample_id
reads <- scater::newSCESet(
    countData = reads,
    phenoData = pheno_data
)
```


```r
keep_feature <- rowSums(counts(reads) > 0) > 0
reads <- reads[keep_feature, ]
```


```r
ercc <- featureNames(reads)[grepl("ERCC-", featureNames(reads))]
mt <- c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
        "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
        "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
        "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
        "ENSG00000198840")
```


```r
reads <- scater::calculateQCMetrics(
    reads,
    feature_controls = list(ERCC = ercc, MT = mt)
)
```


```r
hist(
    reads$total_counts,
    breaks = 100
)
abline(v = 1.3e6, col = "red")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{08-exprs-qc-reads_files/figure-latex/total-counts-hist-reads-1} 

}

\caption{Histogram of library sizes for all cells}(\#fig:total-counts-hist-reads)
\end{figure}


```r
filter_by_total_counts <- (reads$total_counts > 1.3e6)
```


```r
knitr::kable(
    as.data.frame(table(filter_by_total_counts)),
    booktabs = TRUE,
    row.names = FALSE,
    caption = 'The number of cells removed by total counts filter (FALSE)'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-10)The number of cells removed by total counts filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_by\_total\_counts & Freq\\
\midrule
FALSE & 180\\
TRUE & 684\\
\bottomrule
\end{tabular}
\end{table}


```r
hist(
    reads$total_features,
    breaks = 100
)
abline(v = 7000, col = "red")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{08-exprs-qc-reads_files/figure-latex/total-features-hist-reads-1} 

}

\caption{Histogram of the number of detected genes in all cells}(\#fig:total-features-hist-reads)
\end{figure}


```r
filter_by_expr_features <- (reads$total_features > 7000)
```


```r
knitr::kable(
    as.data.frame(table(filter_by_expr_features)),
    booktabs = TRUE,
    row.names = FALSE,
    caption = 'The number of cells removed by total features filter (FALSE)'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-12)The number of cells removed by total features filter (FALSE)}
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


```r
scater::plotPhenoData(
    reads,
    aes_string(x = "total_features",
               y = "pct_counts_feature_controls_MT",
               colour = "batch")
)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{08-exprs-qc-reads_files/figure-latex/mt-vs-counts-reads-1} 

}

\caption{Percentage of counts in MT genes}(\#fig:mt-vs-counts-reads)
\end{figure}


```r
scater::plotPhenoData(
    reads,
    aes_string(x = "total_features",
               y = "pct_counts_feature_controls_ERCC",
               colour = "batch")
)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{08-exprs-qc-reads_files/figure-latex/ercc-vs-counts-reads-1} 

}

\caption{Percentage of counts in ERCCs}(\#fig:ercc-vs-counts-reads)
\end{figure}


```r
filter_by_ERCC <- reads$batch != "NA19098.r2" &
    reads$pct_counts_feature_controls_ERCC < 25
```


```r
knitr::kable(
  as.data.frame(table(filter_by_ERCC)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by ERCC filter (FALSE)'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-14)The number of cells removed by ERCC filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_by\_ERCC & Freq\\
\midrule
FALSE & 103\\
TRUE & 761\\
\bottomrule
\end{tabular}
\end{table}


```r
filter_by_MT <- reads$pct_counts_feature_controls_MT < 30
```


```r
knitr::kable(
  as.data.frame(table(filter_by_MT)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by MT filter (FALSE)'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-16)The number of cells removed by MT filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_by\_MT & Freq\\
\midrule
FALSE & 18\\
TRUE & 846\\
\bottomrule
\end{tabular}
\end{table}


```r
reads$use <- (
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
  as.data.frame(table(reads$use)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by manual filter (FALSE)'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-18)The number of cells removed by manual filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
Var1 & Freq\\
\midrule
FALSE & 259\\
TRUE & 605\\
\bottomrule
\end{tabular}
\end{table}


```r
reads$use_default <- (
    # remove cells with unusual numbers of genes
    !reads$filter_on_total_features &
    # sufficient molecules counted
    !reads$filter_on_total_counts &
    # sufficient endogenous RNA
    !reads$filter_on_pct_counts_feature_controls_ERCC &
    # remove cells with unusual number of reads in MT genes
    !reads$filter_on_pct_counts_feature_controls_MT &
    # controls shouldn't be used in downstream analysis
    !reads$is_cell_control
)
```


```r
knitr::kable(
  as.data.frame(table(reads$use_default)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by default filter (FALSE)'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-20)The number of cells removed by default filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
Var1 & Freq\\
\midrule
FALSE & 37\\
TRUE & 827\\
\bottomrule
\end{tabular}
\end{table}


```r
reads <-
scater::plotPCA(reads,
                size_by = "total_features", 
                shape_by = "use",
                pca_data_input = "pdata",
                detect_outliers = TRUE,
                return_SCESet = TRUE)
```

```
## The following cells/samples are detected as outliers:
## NA19098.r1.B10
## NA19098.r1.D07
## NA19098.r1.E04
## NA19098.r1.F06
## NA19098.r1.H08
## NA19098.r1.H09
## NA19098.r2.A01
## NA19098.r2.A06
## NA19098.r2.A09
## NA19098.r2.A12
## NA19098.r2.B01
## NA19098.r2.B11
## NA19098.r2.B12
## NA19098.r2.C04
## NA19098.r2.C09
## NA19098.r2.D02
## NA19098.r2.D03
## NA19098.r2.D09
## NA19098.r2.E04
## NA19098.r2.E07
## NA19098.r2.F01
## NA19098.r2.F11
## NA19098.r2.G01
## NA19098.r2.G05
## NA19098.r2.G10
## NA19098.r2.H01
## NA19098.r2.H07
## NA19098.r2.H08
## NA19098.r2.H12
## NA19098.r3.A05
## NA19098.r3.A07
## NA19098.r3.B02
## NA19098.r3.C07
## NA19098.r3.E05
## NA19098.r3.E08
## NA19098.r3.E09
## NA19098.r3.F11
## NA19098.r3.F12
## NA19098.r3.G02
## NA19098.r3.G03
## NA19098.r3.G04
## NA19098.r3.G11
## NA19098.r3.G12
## NA19098.r3.H08
## NA19101.r1.A01
## NA19101.r1.A12
## NA19101.r1.B01
## NA19101.r1.B06
## NA19101.r1.E09
## NA19101.r1.E11
## NA19101.r1.F05
## NA19101.r1.F10
## NA19101.r1.G01
## NA19101.r1.G06
## NA19101.r1.H04
## NA19101.r1.H09
## NA19101.r2.A03
## NA19101.r2.C10
## NA19101.r2.E05
## NA19101.r2.F02
## NA19101.r2.H04
## NA19101.r2.H10
## NA19101.r3.A02
## NA19101.r3.A03
## NA19101.r3.A05
## NA19101.r3.A09
## NA19101.r3.B05
## NA19101.r3.C01
## NA19101.r3.C09
## NA19101.r3.C12
## NA19101.r3.D01
## NA19101.r3.D04
## NA19101.r3.D07
## NA19101.r3.D09
## NA19101.r3.E08
## NA19101.r3.F09
## NA19101.r3.G09
## NA19101.r3.H01
## NA19101.r3.H03
## NA19101.r3.H07
## NA19101.r3.H09
## NA19239.r1.F05
## NA19239.r1.G05
## NA19239.r2.B01
## NA19239.r2.B03
## NA19239.r2.B10
## NA19239.r2.B11
## NA19239.r2.C03
## NA19239.r2.C06
## NA19239.r2.C08
## NA19239.r2.D07
## NA19239.r2.D09
## NA19239.r2.E09
## NA19239.r2.F04
## NA19239.r2.F06
## NA19239.r2.F07
## NA19239.r2.F12
## NA19239.r2.G03
## NA19239.r2.G08
## NA19239.r2.H02
## NA19239.r2.H03
## NA19239.r2.H07
## NA19239.r3.A01
## NA19239.r3.B09
## NA19239.r3.C04
## NA19239.r3.C07
## NA19239.r3.E01
## NA19239.r3.E03
## NA19239.r3.E12
## NA19239.r3.H02
## NA19239.r3.H10
## Variables with highest loadings for PC1 and PC2:
## \begin{tabular}{l|r|r}
## \hline
##   & PC1 & PC2\\
## \hline
## pct\_counts\_feature\_controls & 0.5057646 & -0.2473134\\
## \hline
## pct\_counts\_top\_100\_features & 0.4888852 & -0.2277068\\
## \hline
## n\_detected\_feature\_controls & 0.0231277 & -0.6235516\\
## \hline
## log10\_counts\_feature\_controls & -0.1226860 & -0.6576822\\
## \hline
## total\_features & -0.4655518 & -0.2219694\\
## \hline
## log10\_counts\_endogenous\_features & -0.5223679 & -0.1278782\\
## \hline
## \end{tabular}
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{08-exprs-qc-reads_files/figure-latex/auto-cell-filt-reads-1} 

}

\caption{PCA plot used for automatic detection of cell outliers}(\#fig:auto-cell-filt-reads)
\end{figure}


```r
knitr::kable(
  as.data.frame(table(reads$outlier)),
  booktabs = TRUE,
  row.names = FALSE,
  caption = 'The number of cells removed by automatic filter (FALSE)'
)
```

\begin{table}

\caption{(\#tab:unnamed-chunk-21)The number of cells removed by automatic filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
Var1 & Freq\\
\midrule
FALSE & 753\\
TRUE & 111\\
\bottomrule
\end{tabular}
\end{table}


```r
def <- colnames(reads)[!reads$use_default]
auto <- colnames(reads)[reads$outlier]
man <- colnames(reads)[!reads$use]
venn.diag <- limma::vennCounts(cbind(colnames(reads) %in% def,
                                     colnames(reads) %in% auto,
                                     colnames(reads) %in% man))
limma::vennDiagram(venn.diag,
                   names = c("Default", "Automatic", "Manual"),
                   circle.col = c("magenta", "blue", "green"))
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{08-exprs-qc-reads_files/figure-latex/cell-filt-comp-reads-1} 

}

\caption{Comparison of the default, automatic and manual cell filters}(\#fig:cell-filt-comp-reads)
\end{figure}


```r
scater::plotQC(reads, type = "highest-expression")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{08-exprs-qc-reads_files/figure-latex/top50-gene-expr-reads-1} 

}

\caption{Number of total counts consumed by the top 50 expressed genes}(\#fig:top50-gene-expr-reads)
\end{figure}


```r
filter_genes <- apply(counts(reads[, pData(reads)$use]), 1, 
                      function(x) length(x[x > 1]) >= 2)
fData(reads)$use <- filter_genes
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

\caption{(\#tab:unnamed-chunk-23)The number of genes removed by gene filter (FALSE)}
\centering
\begin{tabular}[t]{lr}
\toprule
filter\_genes & Freq\\
\midrule
FALSE & 2665\\
TRUE & 16061\\
\bottomrule
\end{tabular}
\end{table}


```r
dim(reads[fData(reads)$use, pData(reads)$use])
```

```
## Features  Samples 
##    16061      605
```


```r
set_exprs(reads, "log2_counts") <- log2(counts(reads) + 1)
```


```r
saveRDS(reads, file = "tung/reads.rds")
```

By comparing Figure \@ref(fig:cell-filt-comp) and Figure \@ref(fig:cell-filt-comp-reads), it is clear that the reads based filtering removed 49 more cells than the UMI based analysis. If you go back and compare the results you should be able to conclude that the ERCC and MT filters are more strict for the reads-based analysis.
