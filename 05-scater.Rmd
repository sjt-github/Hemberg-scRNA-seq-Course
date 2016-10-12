---
knit: bookdown::preview_chapter
---

# scater package

## Introduction

[scater](https://github.com/davismcc/scater) is a R package single-cell RNA-seq analysis. The package contains several useful methods for quality control, visualisation and pre-processing of data prior to further downstream analysis.

[scater](https://github.com/davismcc/scater) features the following functionality:

* Automated computation of QC metrics
* Transcript quantification from read data with pseudo-alignment
* Data format standardisation
* Rich visualizations for exploratory analysis
* Seamless integration into the Bioconductor universe
* Simple normalisation methods

We highly recommend to use [scater](https://github.com/davismcc/scater) for all single-cell RNA-seq analyses and [scater](https://github.com/davismcc/scater) is the basis of the first part of the course.

## scater workflow

![](figures/scater_qc_workflow.png)

## Terminology

(this chapter is taken from the [scater vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/vignette.html))

* The capabilities of scater are built on top of Bioconductor’s Biobase.
* In Bioconductor terminology we assay numerous __“features”__ for a number of __“samples”__.
* Features, in the context of [scater](https://github.com/davismcc/scater), correspond most commonly to genes or transcripts, but could be any general genomic or transcriptomic regions (e.g. exon) of interest for which we take measurements.
* In the following chapters it may be more intuitive to mentally replace __“feature”__ with __“gene”__ or __“transcript”__ (depending on the context of the study) wherever __“feature”__ appears.
* In the scater context, __“samples”__ refer to individual cells that we have assayed.

## SCESet class

(this chapter is taken from the [scater vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/vignette.html))

In scater we organise single-cell expression data in objects of the __SCESet__ class. The class inherits the Bioconductor ExpressionSet class, which provides a common interface familiar to those who have analyzed microarray experiments with Bioconductor. The class requires three input files:

* __exprs__, a numeric matrix of expression values, where rows are features, and columns are cells
* __phenoData__, an _AnnotatedDataFrame_ object, where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
* __featureData__, an _AnnotatedDataFrame_ object, where rows are features (e.g. genes), and columns are feature attributes, such as biotype, gc content, etc.

For more details about other features inherited from Bioconductor’s __ExpressionSet__ class, type `?ExpressionSet` at the R prompt.

When the data are encapsulated in the __SCESet__ class, [scater](https://github.com/davismcc/scater) will automatically calculate several different properties. This will be demonstrated in the subsequent chapters.
