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
set.seed(42)
expr.mat <- matrix(rnorm(400), nrow = 4, ncol = 100)
rownames(expr.mat) <- c(101, 102, 103, 104)
probe.ids <- c(101, 102, 103, 104)
cor.mat <- gene.correlation.matrix(expr.mat, probe.ids)
correlation.heatmap.plot(cor.mat, gene.names = c("BRCA1", "TP53", "MYC", "EGFR"))

# }
```
