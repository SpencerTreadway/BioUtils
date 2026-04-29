# Build Analysis Data Frame

Combines a multi-gene expression matrix with sample metadata into a
long-format, analysis-ready data frame. Each row represents one
gene-sample observation, with human-readable gene names resolved from
the annotation data frame.

## Usage

``` r
build.analysis.df(expr.matrix, phenotype, genes, group.col = "disease.state")
```

## Arguments

- expr.matrix:

  Numeric matrix of gene expression values as returned by
  [`get.gene.expression()`](https://github.com/SpencerTreadway/BioUtils/reference/get.gene.expression.md).
  Rows are probe IDs, columns are sample IDs.

- phenotype:

  Data frame of sample metadata as returned by
  `extract.expression()$phenotype`. Row names must correspond to the
  column names of `expr.matrix`.

- genes:

  Data frame of gene annotations as returned by
  `extract.expression()$gene`. Row indices must correspond to the probe
  IDs (row names) of `expr.matrix`. Used to resolve probe IDs to the
  human-readable gene names stored in the `"Gene symbol"` column.

- group.col:

  Character. Name of the column in `phenotype` to use as the grouping
  variable. Default is `"disease.state"`.

## Value

A long-format data frame with one row per gene-sample pair and three
columns:

- gene:

  Character. Human-readable gene name resolved from the gene annotation
  data frame.

- expression:

  Numeric. Expression value for that gene-sample pair.

- group:

  Character or factor. Group label for each sample, renamed from the
  `group.col` column in `phenotype` for consistency with downstream
  functions.

Rows with `NaN` or `NA` expression values are removed.

## Details

The function pivots `expr.matrix` from wide format (probes x samples) to
long format, merges sample metadata by sample ID, resolves probe IDs to
gene names using the annotation data frame, and selects the three
columns needed for analysis. The output is the standard input for
[`analyze.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/analyze.gene.md),
[`gene.analysis.plot()`](https://github.com/SpencerTreadway/BioUtils/reference/gene.analysis.plot.md),
and
[`fit.lasso()`](https://github.com/SpencerTreadway/BioUtils/reference/fit.lasso.md).

Requires `tidyr` for the pivot step. Ensure `tidyr` is listed under
`Imports` in the package DESCRIPTION.

## Examples

``` r
if (FALSE) { # \dontrun{
geo <- extract.expression(load.geo.soft("GDS3268.soft"))
probe <- find.probe.by.gene(geo$gene, c("MUC20", "ADH1A"))
expr <- get.gene.expression(geo$expression, probe)
df <- build.analysis.df(expr, geo$phenotype, geo$gene)
head(df)
} # }
```
