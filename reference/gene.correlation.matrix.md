# Compute Pairwise Gene Co-expression Correlation Matrix

Calculates pairwise correlations between a set of genes across all
samples, producing a symmetric correlation matrix that quantifies
co-expression relationships.

## Usage

``` r
gene.correlation.matrix(expression.matrix, probe.ids, method = "pearson")
```

## Arguments

- expression.matrix:

  Numeric matrix of gene expression values as returned by
  `extract.expression()$expression`. Rows are probes, columns are
  samples.

- probe.ids:

  Integer vector of probe IDs to include, as returned by
  [`find.probe.by.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/find.probe.by.gene.md).

- method:

  Character. Correlation method: `"pearson"` (default), `"spearman"`, or
  `"kendall"`. Spearman is recommended when distributions are skewed or
  outliers are a concern.

## Value

A symmetric numeric matrix of dimensions `length(probe.ids)` x
`length(probe.ids)`, where each cell contains the pairwise correlation
coefficient across all samples. Diagonal values are 1. Row and column
names correspond to probe IDs.

## Details

Co-expression correlations capture whether two genes tend to be
simultaneously up- or down-regulated across samples, which can suggest
shared regulatory control or pathway membership. The resulting matrix is
the direct input for `plot.correlation.heatmap()`.

## Examples

``` r
# \donttest{
geo <- extract.expression(load.geo.soft(accession = "GDS507", log.transform = TRUE))
#> GDS507 not found locally, downloading from NCBI GEO...
#> Using locally cached version of GDS507 found here:
#> /tmp/RtmpYU9fuS/GDS507.soft.gz 
#> Using locally cached version of GPL97 found here:
#> /tmp/RtmpYU9fuS/GPL97.annot.gz 
probe.ids <- find.probe.by.gene(geo$gene, c("BRCA1", "TP53", "MYC"))
cor.mat <- gene.correlation.matrix(geo$expression, probe.ids, method = "spearman")
plot.correlation.heatmap(cor.mat)
#> Error in plot.correlation.heatmap(cor.mat): could not find function "plot.correlation.heatmap"
# }
```
