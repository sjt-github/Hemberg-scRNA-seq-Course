---
output: html_document
---

# Data visualization (Reads)


```r
library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
reads <- readRDS("tung/reads.rds")
reads.qc <- reads[fData(reads)$use, pData(reads)$use]
endog_genes <- !fData(reads.qc)$is_feature_control
```




```r
scater::plotPCA(reads[endog_genes, ],
                ntop = 500,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{10-exprs-overview-reads_files/figure-latex/expr-overview-pca-before-qc-reads1-1} 

}

\caption{PCA plot of the tung data}(\#fig:expr-overview-pca-before-qc-reads1)
\end{figure}


```r
scater::plotPCA(reads[endog_genes, ],
                ntop = 500,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "log2_counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{10-exprs-overview-reads_files/figure-latex/expr-overview-pca-before-qc-reads2-1} 

}

\caption{PCA plot of the tung data}(\#fig:expr-overview-pca-before-qc-reads2)
\end{figure}


```r
scater::plotPCA(reads.qc[endog_genes, ],
                ntop = 500,
                colour_by = "batch",
                size_by = "total_features",
                shape_by = "individual",
                exprs_values = "log2_counts")
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{10-exprs-overview-reads_files/figure-latex/expr-overview-pca-after-qc-reads-1} 

}

\caption{PCA plot of the tung data}(\#fig:expr-overview-pca-after-qc-reads)
\end{figure}


```r
scater::plotTSNE(reads[endog_genes, ],
                 ntop = 500,
                 perplexity = 130,
                 colour_by = "batch",
                 size_by = "total_features",
                 shape_by = "individual",
                 exprs_values = "log2_counts",
                 rand_seed = 123456)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{10-exprs-overview-reads_files/figure-latex/expr-overview-tsne-before-qc-reads-1} 

}

\caption{tSNE map of the tung data}(\#fig:expr-overview-tsne-before-qc-reads)
\end{figure}


```r
scater::plotTSNE(reads.qc[endog_genes, ],
                 ntop = 500,
                 perplexity = 130,
                 colour_by = "batch",
                 size_by = "total_features",
                 shape_by = "individual",
                 exprs_values = "log2_counts",
                 rand_seed = 123456)
```

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{10-exprs-overview-reads_files/figure-latex/expr-overview-tsne-after-qc-reads-1} 

}

\caption{tSNE map of the tung data}(\#fig:expr-overview-tsne-after-qc-reads)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{10-exprs-overview-reads_files/figure-latex/expr-overview-tsne-after-qc-exercise2-1-1} 

}

\caption{tSNE map of the tung data (perplexity = 10)}(\#fig:expr-overview-tsne-after-qc-exercise2-1)
\end{figure}

\begin{figure}

{\centering \includegraphics[width=0.9\linewidth]{10-exprs-overview-reads_files/figure-latex/expr-overview-tsne-after-qc-exercise2-2-1} 

}

\caption{tSNE map of the tung data (perplexity = 200)}(\#fig:expr-overview-tsne-after-qc-exercise2-2)
\end{figure}
