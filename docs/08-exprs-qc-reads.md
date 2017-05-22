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
\begin{tabular}[t]{r}
\toprule
Freq\\


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
