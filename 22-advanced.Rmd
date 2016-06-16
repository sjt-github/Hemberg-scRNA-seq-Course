---
output: html_document
---

# Advanced exercises

For the final part of the course we provide you with a few open-ended
question. The goal of these questions is to stimulate new ideas about
what types of analyses that you can carry out with scRNA-seq data.

We have collected published scRNA-seq datasets and stored the
expression matrices as .rds files which you can download

* [Treutlein et al](https://www.dropbox.com/s/2964wgz27vn0jv1/treutlein.rds?dl=0)
* [Deng et al](https://www.dropbox.com/s/l6ycsod5qb4ryuz/deng.rds?dl=0)
* [Pollen et al](https://www.dropbox.com/s/fxpjno6sl3ui644/pollen1.rds?dl=0), [alt Pollen](https://www.dropbox.com/s/fxpjno6sl3ui644/pollen1.rds?dl=0)
* [Usoskin et al](https://www.dropbox.com/s/fxpjno6sl3ui644/pollen1.rds?dl=0), [alt Usoskin](https://www.dropbox.com/s/fxpjno6sl3ui644/pollen1.rds?dl=0), [another Usoskin](https://www.dropbox.com/s/fxpjno6sl3ui644/pollen1.rds?dl=0)
* [Klein et al](https://www.dropbox.com/s/pzj5mt8w8q2nl8p/klein.rds?dl=0)

Here are some suggestions for questions that you can explore:

* SC3 uses a different method compared to SCDE and DESeq2 to compare
  several groups at once to identify marker genes. Can you use SCDE or
  DESeq2 to compare the clusters from the above datasets to find
  marker genes? How does this set of genes differ from the one
  obtained by SC3?

* Both the [Pollen et
  al](http://www.nature.com/nbt/journal/v32/n10/abs/nbt.2967.html) and
  the [Usoskin et
  al](http://www.nature.com/neuro/journal/v18/n1/abs/nn.3881.html)
  come in more than one version. The only difference between the files
  are the cell labels. The different datasets have been clustered at
  different granularities. Use SC3 and other clustering tools to
  explore the different granularities. Do you think that there is
  stronger support for one choice of granularity over another? How can
  you best motivate one choice over another? Does the interpretation
  of the data change depending on the choice of $k$? Can you identify
  some other choice of $k$ that the authors did not suggest?

* In addition to identifying differences in expression between
  different cell types, we are also interested in finding differences
  in regulatory interactions. One way of characterizing regulatory
  interactions is through correlation coefficients. Can you identify
  differences in regulation between different clusters of cells in any
  of the above datasets?
