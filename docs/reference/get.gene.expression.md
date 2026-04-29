# Extract Gene Expression by Probe ID

Retrieves expression values for one or more genes from the expression
matrix using their integer probe IDs.

## Usage

``` r
get.gene.expression(expression, probe.id)
```

## Arguments

- expression:

  Numeric matrix of gene expression values as returned by
  `extract.expression()$expression`. Rows are probe IDs, columns are
  sample IDs.

- probe.id:

  Integer vector of one or more probe IDs as returned by
  [`find.probe.by.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/find.probe.by.gene.md).

## Value

A numeric matrix with one row per probe and one column per sample. Row
names are the probe IDs; column names are the sample IDs. When a single
probe ID is supplied the result is still a matrix (not a vector),
ensuring that
[`build.analysis.df()`](https://github.com/SpencerTreadway/BioUtils/reference/build.analysis.df.md)
receives a consistent input type.

## Details

The returned matrix is the direct input for
[`build.analysis.df()`](https://github.com/SpencerTreadway/BioUtils/reference/build.analysis.df.md).
Using [`as.matrix()`](https://rdrr.io/r/base/matrix.html) ensures the
single-probe case returns a one-row matrix rather than a named vector,
which keeps the `pivot_longer` step in
[`build.analysis.df()`](https://github.com/SpencerTreadway/BioUtils/reference/build.analysis.df.md)
consistent regardless of how many genes are queried.

## Examples

``` r
if (FALSE) { # \dontrun{
geo <- extract.expression(load.geo.soft("GDS3268.soft"))
probes <- find.probe.by.gene(geo$gene, c(
  "mucin 20, cell surface associated",
  "alcohol dehydrogenase 1A (class I), alpha polypeptide"
))
expr <- get.gene.expression(geo$expression, probes)
dim(expr)  # probes x samples
} # }
```
