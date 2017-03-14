---
output: html_document
---

# Clustering example {#clust-methods}




```r
library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(pheatmap)
set.seed(1234567)
```

To illustrate clustering of scRNA-seq data, we consider the `Pollen` dataset of cells from 
different human tissues [@Pollen2014-cu]. We have preprocessed the dataset and created a 
scater object in advance. We have also annotated the cells with the cell type information 
(it is the `cell_type1` column in the `phenoData` slot).

## Pollen dataset

Let's load the data and look at it:

```r
pollen <- readRDS("clustering/pollen.rds")
pollen
```

```
## SCESet (storageMode: lockedEnvironment)
## assayData: 23730 features, 301 samples 
##   element names: exprs, is_exprs, tpm 
## protocolData: none
## phenoData
##   rowNames: Hi_2338_1 Hi_2338_2 ... Hi_GW16_26 (301 total)
##   varLabels: cell_type1 cell_type2 ... is_cell_control (33 total)
##   varMetadata: labelDescription
## featureData
##   featureNames: A1BG A1BG-AS1 ... ZZZ3 (23730 total)
##   fvarLabels: mean_exprs exprs_rank ... feature_symbol (11 total)
##   fvarMetadata: labelDescription
## experimentData: use 'experimentData(object)'
## Annotation:
```

Let's look at the cell type annotation:

```r
table(pData(pollen)$cell_type1)
```

```
## 
##   2338   2339     BJ   GW16   GW21 GW21+3  hiPSC   HL60   K562   Kera 
##     22     17     37     26      7     17     24     54     42     40 
##    NPC 
##     15
```

A simple PCA analysis already separates some strong cell types and provides some insights in the data structure:

```r
plotPCA(pollen, colour_by = "cell_type1")
```

<img src="18-clustering_files/figure-html/unnamed-chunk-5-1.png" width="672" style="display: block; margin: auto;" />

## SC3

Let's run `SC3` clustering on the Pollen data. The advantage of the `SC3` is that it can directly take a [scater](http://bioconductor.org/packages/scater/) object (see previous chapters) as an input.

Now let's image we do not know the number of clusters _k_ (cell types). `SC3` can estimate a number of clusters for you:

```r
pollen <- sc3_prepare(pollen, ks = 2:5)
```

```
## Setting SC3 parameters...
```

```
## Setting a range of k...
```

```r
pollen <- sc3_estimate_k(pollen)
```

```
## Estimating k...
```

```r
pollen@sc3$k_estimation
```

```
## [1] 11
```

Interestingly, the number of cell types predicted by `SC3` is the same as the number of cell types in the Pollen data annotation.

Now we are ready to run `SC3` (we also ask it to calculate biological properties of the clusters): 

```r
pollen <- sc3(pollen, ks = 11, biology = TRUE)
```

```
## Setting SC3 parameters...
```

```
## Setting a range of k...
```

```
## Calculating distances between the cells...
```

```
## Performing transformations and calculating eigenvectors...
```

```
## Performing k-means clustering...
```

```
## Calculating consensus matrix...
```

```
## Calculating biology...
```

`SC3` result consists of several different outputs (please look in [@Kiselev2016-bq] and [SC3 vignette](http://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/my-vignette.html) for more details). Here we show some of them:

Consensus matrix:

```r
sc3_plot_consensus(pollen, k = 11, show_pdata = "cell_type1")
```

<img src="18-clustering_files/figure-html/unnamed-chunk-8-1.png" width="672" style="display: block; margin: auto;" />

Silhouette plot:

```r
sc3_plot_silhouette(pollen, k = 11)
```

<img src="18-clustering_files/figure-html/unnamed-chunk-9-1.png" width="672" style="display: block; margin: auto;" />

Heatmap of the expression matrix:

```r
sc3_plot_expression(pollen, k = 11, show_pdata = "cell_type1")
```

<img src="18-clustering_files/figure-html/unnamed-chunk-10-1.png" width="672" style="display: block; margin: auto;" />

Identified marker genes:

```r
sc3_plot_markers(pollen, k = 11, show_pdata = "cell_type1")
```

<img src="18-clustering_files/figure-html/unnamed-chunk-11-1.png" width="672" style="display: block; margin: auto;" />

PCA plot with highlighted `SC3` clusters:

```r
plotPCA(pollen, colour_by = "sc3_11_clusters")
```

<img src="18-clustering_files/figure-html/unnamed-chunk-12-1.png" width="672" style="display: block; margin: auto;" />

Note, that one can also run `SC3` in an interactive `Shiny` session:

```r
sc3_interactive(pollen)
```

This command will open `SC3` in a web browser.

* __Exercise 1__: Run `SC3` for $k$ from 9 to 13 and explore different clustering solutions in your web browser.

* __Exercise 2__: Which clusters are the most stable when $k$ is changed from 9 to 13? (Look at the "Stability" tab)

* __Exercise 3__: Check out differentially expressed genes and marker genes for the obtained clusterings. Please use $k=11$.

* __Exercise 4__: Change the marker genes threshold (the default is 0.85). Does __SC3__ find more marker genes?

## pcaReduce

`pcaReduce` operates directly on the expression matrix. It is recommended to use a gene filter and log transformation before running `pcaReduce`. We will use the default `SC3` gene filter (note that the `exprs` slot of a `scater` object is log-transformed by default).


```r
# use the same gene filter as in SC3
input <- exprs(pollen[fData(pollen)$sc3_gene_filter, ])
```

There are several parameters used by `pcaReduce`:
* `nbt` defines a number of `pcaReduce` runs (it is stochastic and may have different solutions after different runs)
* `q` defines number of dimensions to start clustering with. The output will contain partitions for all $k$ from 2 to q+1.
* `method` defines a method used for clustering. `S` - to perform sampling based merging, `M` - to perform merging based on largest probability.

We will run `pcaReduce` 1 time:

```r
# run pcaReduce 1 time creating hierarchies from 1 to 30 clusters
pca.red <- PCAreduce(t(input), nbt = 1, q = 30, method = 'S')[[1]]
```


```r
pData(pollen)$pcaReduce <- as.character(pca.red[,32 - 11])
plotPCA(pollen, colour_by = "pcaReduce")
```

<img src="18-clustering_files/figure-html/unnamed-chunk-16-1.png" width="672" style="display: block; margin: auto;" />

__Exercise 5__: Run pcaReduce for $k=2$ and plot a similar PCA plot. Does it look good?

__Hint__: When running pcaReduce for different $k$s you do not need to rerun PCAreduce function, just use already calculated `pca.red` object.

__Our solution__:
<div class="figure" style="text-align: center">
<img src="18-clustering_files/figure-html/clust-pca-reduce2-1.png" alt="Clustering solutions of pcaReduce method for $k=2$." width="672" />
<p class="caption">(\#fig:clust-pca-reduce2)Clustering solutions of pcaReduce method for $k=2$.</p>
</div>

__Exercise 6__: Compare the results between `SC3` and `pcaReduce` for $k=11$. What is
the main difference between the solutions provided by the two
different methods?

__Our solution__:
<img src="18-clustering_files/figure-html/unnamed-chunk-17-1.png" width="672" style="display: block; margin: auto;" />


## tSNE + kmeans

[tSNE](https://lvdmaaten.github.io/tsne/) plots that we saw before (\@ref(visual-tsne)) when used the __scater__ package are made by using the [Rtsne](https://cran.r-project.org/web/packages/Rtsne/index.html) and [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) packages. Here we will do the same:

```r
pollen <- plotTSNE(pollen, rand_seed = 1, return_SCESet = TRUE)
```

<div class="figure" style="text-align: center">
<img src="18-clustering_files/figure-html/clust-tsne-1.png" alt="tSNE map of the patient data" width="672" />
<p class="caption">(\#fig:clust-tsne)tSNE map of the patient data</p>
</div>

Note that all points on the plot above are black. This is different from what we saw before, when the cells were coloured based on the annotation. Here we do not have any annotation and all cells come from the same batch, therefore all dots are black.

Now we are going to apply _k_-means clustering algorithm to the cloud of points on the tSNE map. How many groups do you see in the cloud?

We will start with $k=8$:

```r
pData(pollen)$tSNE_kmeans <- as.character(kmeans(pollen@reducedDimension, centers = 8)$clust)
plotTSNE(pollen, rand_seed = 1, colour_by = "tSNE_kmeans")
```

<div class="figure" style="text-align: center">
<img src="18-clustering_files/figure-html/clust-tsne-kmeans2-1.png" alt="tSNE map of the patient data with 8 colored clusters, identified by the k-means clustering algorithm" width="672" />
<p class="caption">(\#fig:clust-tsne-kmeans2)tSNE map of the patient data with 8 colored clusters, identified by the k-means clustering algorithm</p>
</div>

__Exercise 7__: Make the same plot for $k=11$.

__Exercise 8__: Compare the results between `SC3` and `tSNE+kmeans`. Can the
results be improved by changing the `perplexity` parameter?

__Our solution__:
<img src="18-clustering_files/figure-html/unnamed-chunk-18-1.png" width="672" style="display: block; margin: auto;" />

As you may have noticed, both `pcaReduce` and `tSNE+kmeans` are stochastic
and give different results every time they are run. To get a better
overview of the solutions, we need to run the methods multiple times. `SC3` is also stochastic, but thanks to the consensus step, it is more robust and less likely to produce different outcomes.

## SNN-Cliq

Here we run SNN-cliq with te default parameters provided in the author's example:


```r
distan <- "euclidean"
par.k <- 3
par.r <- 0.7
par.m <- 0.5
# construct a graph
scRNA.seq.funcs::SNN(
    data = t(input),
    outfile = "snn-cliq.txt",
    k = par.k,
    distance = distan
)
# find clusters in the graph
snn.res <- 
    system(
        paste0(
            "python snn-cliq/Cliq.py ", 
            "-i snn-cliq.txt ",
            "-o res-snn-cliq.txt ",
            "-r ", par.r,
            " -m ", par.m
        ),
        intern = TRUE
    )
cat(paste(snn.res, collapse = "\n"))
```

```
## input file snn-cliq.txt
## find 65 quasi-cliques
## merged into 15 clusters
## unique assign done
```

```r
snn.res <- read.table("res-snn-cliq.txt")
# remove files that were created during the analysis
system("rm snn-cliq.txt res-snn-cliq.txt")

pData(pollen)$SNNCliq <- as.character(snn.res[,1])
plotPCA(pollen, colour_by = "SNNCliq")
```

<img src="18-clustering_files/figure-html/unnamed-chunk-19-1.png" width="672" style="display: block; margin: auto;" />

__Exercise 9__: Compare the results between `SC3` and `SNN-Cliq`.

__Our solution__:
<img src="18-clustering_files/figure-html/unnamed-chunk-20-1.png" width="672" style="display: block; margin: auto;" />

## SINCERA

As mentioned in the previous chapter [SINCERA](https://research.cchmc.org/pbge/sincera.html) is based on hierarchical clustering. One important thing to keep in mind is that it performs a gene-level z-score transformation before doing clustering:


```r
# perform gene-by-gene per-sample z-score transformation
dat <- apply(input, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
# hierarchical clustering
dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
hc <- hclust(dd, method = "average")
```

If the number of cluster is not known [SINCERA](https://research.cchmc.org/pbge/sincera.html) can identify __k__ as the minimum height of the hierarchical tree that generates no more than a specified number of singleton clusters (clusters containing only 1 cell)

```r
num.singleton <- 0
kk <- 1
for (i in 2:dim(dat)[2]) {
    clusters <- cutree(hc, k = i)
    clustersizes <- as.data.frame(table(clusters))
    singleton.clusters <- which(clustersizes$Freq < 2)
    if (length(singleton.clusters) <= num.singleton) {
        kk <- i
    } else {
        break;
    }
}
cat(kk)
```

```
## 14
```

Let's now visualize the SINCERA results as a heatmap:

```r
pheatmap(
    t(dat),
    cluster_cols = hc,
    cutree_cols = 14,
    kmeans_k = 100,
    show_rownames = FALSE
)
```

<div class="figure" style="text-align: center">
<img src="18-clustering_files/figure-html/clust-sincera-1.png" alt="Clustering solutions of SINCERA method using $k=3$" width="672" />
<p class="caption">(\#fig:clust-sincera)Clustering solutions of SINCERA method using $k=3$</p>
</div>

__Exercise 10__: Compare the results between `SC3` and `SNN-Cliq`.

__Our solution__:
<img src="18-clustering_files/figure-html/unnamed-chunk-23-1.png" width="672" style="display: block; margin: auto;" />

__Exercise 11__: Is using the singleton cluster criteria for finding __k__ a good idea?

## SEURAT

Here we follow an [example](http://satijalab.org/seurat/get_started.html) created by the authors of `SEURAT` (8,500 Pancreas cells). We mostly use default values in various function calls, for more details please consult the documentation and the authors:


```r
library(Seurat)
library(Matrix)
pollen_seurat <- new("seurat", raw.data = get_exprs(pollen, exprs_values = "tpm"))
pollen_seurat <- Setup(pollen_seurat, project = "Pollen")
pollen_seurat <- MeanVarPlot(pollen_seurat)
```

<img src="18-clustering_files/figure-html/unnamed-chunk-24-1.png" width="672" style="display: block; margin: auto;" />

```r
pollen_seurat <- RegressOut(pollen_seurat, latent.vars = c("nUMI"), 
                            genes.regress = pollen_seurat@var.genes)
```

```
## [1] "Regressing out nUMI"
```

```r
pollen_seurat <- PCAFast(pollen_seurat)
```

```
## [1] "PC1"
##  [1] "MAP1B"     "TUBA1A"    "SOX11"     "DPYSL2"    "DCX"      
##  [6] "STMN2"     "RTN1"      "GPM6A"     "DPYSL3"    "MLLT11"   
## [11] "SOX4"      "MAP2"      "CRMP1"     "LOC150568" "NREP"     
## [16] "CALM1"     "FAM110B"   "NFIB"      "MAPT"      "TUBB2B"   
## [21] "RUFY3"     "NNAT"      "CXADR"     "FXYD6"     "MIR100HG" 
## [26] "FOXG1"     "POU3F2"    "KIF5C"     "MN1"       "NCAM1"    
## [1] ""
##  [1] "MYL12A"    "ARHGDIB"   "S100A11"   "IFITM1"    "CKS2"     
##  [6] "SRGN"      "IFITM3"    "ANXA1"     "IFI30"     "GTSF1"    
## [11] "MT2A"      "ISG15"     "KRT8"      "CD53"      "ANXA2"    
## [16] "HIST1H1C"  "HIST1H2BK" "LAPTM5"    "LGALS1"    "NOB1"     
## [21] "MPO"       "PRG2"      "KRT18"     "AIF1"      "NUCB2"    
## [26] "PRAME"     "PSMB9"     "SFXN4"     "TRAP1"     "RNF114"   
## [31] "CD52"     
## [1] ""
## [1] ""
## [1] "PC2"
##  [1] "S100A6"   "CD44"     "TPM1"     "ANXA2"    "IFITM3"   "ANXA1"   
##  [7] "S100A11"  "RND3"     "LGALS1"   "IGFBP3"   "DKK1"     "MT2A"    
## [13] "THBS1"    "TGFBI"    "PLS3"     "TMSB10"   "SERPINE1" "KRT7"    
## [19] "FN1"      "S100A10"  "TIMP1"    "CTGF"     "SAT1"     "PLAUR"   
## [25] "CAV1"     "DDAH1"    "TPM2"     "MYL6"     "EMP1"     "CCND1"   
## [1] ""
##  [1] "MPO"          "SRGN"         "PRG2"         "AIF1"        
##  [5] "MS4A3"        "GTSF1"        "CTSG"         "TESC"        
##  [9] "MZB1"         "LOC100190986" "CD53"         "LAPTM5"      
## [13] "ARHGDIB"      "TRAP1"        "SPNS3"        "RNASE2"      
## [17] "PRTN3"        "BTK"          "APOC2"        "CD33"        
## [21] "LOC100272216" "HSH2D"        "CLEC5A"       "SLC43A1"     
## [25] "PRAME"        "DOK3"         "SERPINB10"    "PLEK"        
## [29] "DOCK2"        "CD52"         "NDUFV1"      
## [1] ""
## [1] ""
## [1] "PC3"
##  [1] "KRT81"    "STEAP4"   "LCN2"     "S100A9"   "KRT15"    "ELF3"    
##  [7] "CEACAM6"  "CLDN4"    "RARRES1"  "SLPI"     "KLK5"     "GRB7"    
## [13] "DHRS3"    "CXorf61"  "RARRES3"  "KLK6"     "IFI27"    "AZGP1"   
## [19] "S100A8"   "MYEOV"    "CXCL17"   "KLK8"     "PDZK1IP1" "BMP3"    
## [25] "MUC1"     "FOLR1"    "TNFSF10"  "MIEN1"    "KRT23"    "VGLL1"   
## [1] ""
##  [1] "COL1A2"    "LUM"       "DCN"       "GREM1"     "PSG5"     
##  [6] "TNFRSF11B" "SERPINE1"  "TPM2"      "COL1A1"    "LOX"      
## [11] "TIMP3"     "VIM"       "TAGLN"     "CRYAB"     "S100A4"   
## [16] "SULF1"     "GLIPR1"    "DKK1"      "COX7A1"    "RGS4"     
## [21] "CCDC80"    "SERPINE2"  "FN1"       "THBS1"     "KRT34"    
## [26] "ALDH1A1"   "KIAA1199"  "FGF7"      "TIMP1"     "WBP5"     
## [31] "MT1E"     
## [1] ""
## [1] ""
## [1] "PC4"
##  [1] "PPME1"    "S100A11"  "KIAA1598" "NUCB2"    "SEMA3C"   "MRPS10"  
##  [7] "MAPT"     "STMN2"    "GTSF1"    "TRIM16"   "DLX6-AS1" "MYT1L"   
## [13] "KRT18"    "PQBP1"    "MT2A"     "ANXA1"    "CXADR"    "POLR2C"  
## [19] "CDK5R1"   "MPO"      "DLX5"     "INA"      "BCL11B"   "PRG2"    
## [25] "KRT8"     "ATAT1"    "PRAME"    "CBWD3"    "S100A10"  "TAF11"   
## [1] ""
##  [1] "CD74"     "HLA-DPA1" "HLA-DRA"  "MS4A1"    "HLA-DQA1" "HLA-DRB5"
##  [7] "HLA-DQB1" "HLA-DQA2" "BLNK"     "HLA-DMA"  "HLA-DRB1" "IRF4"    
## [13] "HLA-DPB1" "HLA-DMB"  "ELK2AP"   "MIR155HG" "CHI3L2"   "TCL1A"   
## [19] "CD48"     "LRMP"     "SLAMF1"   "BCL2A1"   "LY86"     "CLECL1"  
## [25] "HLA-A"    "CRIP1"    "CD27"     "CCL3"     "CYTIP"    "JUN"     
## [31] "PTN"     
## [1] ""
## [1] ""
## [1] "PC5"
##  [1] "GLI3"    "PTN"     "SLC1A3"  "PTPRZ1"  "MDK"     "GPX3"    "CDO1"   
##  [8] "CLU"     "ITGB8"   "JUN"     "AIF1L"   "FAM107A" "ATP1B2"  "FOS"    
## [15] "TSPAN6"  "HEPN1"   "PAX6"    "HOPX"    "FABP7"   "ID4"     "PON2"   
## [22] "SHISA2"  "GPM6B"   "PELI2"   "PMP2"    "ABAT"    "FBXO32"  "NFIA"   
## [29] "C1orf61" "CNN3"   
## [1] ""
##  [1] "HLA-DPA1" "CD74"     "HLA-DRA"  "HLA-DQA1" "MS4A1"    "HLA-DRB5"
##  [7] "HLA-DQB1" "HLA-DQA2" "BLNK"     "HLA-DMA"  "HLA-DRB1" "IRF4"    
## [13] "LPXN"     "MIR155HG" "ELK2AP"   "HLA-DPB1" "HLA-DMB"  "LRMP"    
## [19] "CHI3L2"   "CD48"     "BCL2A1"   "TCL1A"    "CLECL1"   "SLAMF1"  
## [25] "LY86"     "IFI30"    "CRIP1"    "CD27"     "CCL3"     "IGJ"     
## [31] "LAPTM5"  
## [1] ""
## [1] ""
```

```r
pollen_seurat <- RunTSNE(pollen_seurat)
pollen_seurat <- FindClusters(pollen_seurat)
TSNEPlot(pollen_seurat, do.label = T)
```

<img src="18-clustering_files/figure-html/unnamed-chunk-24-2.png" width="672" style="display: block; margin: auto;" />

__Exercise 12__: Compare the results between `SC3` and `SEURAT`.

__Our solution__:
<img src="18-clustering_files/figure-html/unnamed-chunk-25-1.png" width="672" style="display: block; margin: auto;" />


Seurat can also find marker genes, e.g. marker genes for cluster 2:

```r
markers <- FindMarkers(pollen_seurat, 2)
FeaturePlot(pollen_seurat, 
            head(rownames(markers)), 
            cols.use = c("lightgrey", "blue"), 
            nCol = 3)
```

<img src="18-clustering_files/figure-html/unnamed-chunk-26-1.png" width="672" style="display: block; margin: auto;" />

__Exercise 13__: Compare marker genes provided by `SEURAT` and `SC3`.
