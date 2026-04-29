# Extract Expression and Metadata from GEO Dataset

Decomposes an `ExpressionSet` object into its core components, returning
a named list that is used as the primary data structure throughout
BioUtils.

## Usage

``` r
extract.expression(eset)
```

## Arguments

- eset:

  An `ExpressionSet` object as returned by
  [`load.geo.soft()`](https://github.com/SpencerTreadway/BioUtils/reference/load.geo.soft.md).

## Value

A named list with three elements:

- expression:

  Numeric matrix of gene expression values. Rows are probe IDs, columns
  are sample IDs.

- phenotype:

  Data frame of sample metadata (e.g., disease state, age, tissue). Row
  names are sample IDs matching the column names of `expression`.

- gene:

  Data frame of probe-level gene annotations (e.g., gene symbol, gene
  title, chromosomal location). Row names are probe IDs matching the row
  names of `expression`.

## Details

The three returned components form the basis of all downstream BioUtils
workflows. The `gene` data frame is passed to
[`find.probe.by.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/find.probe.by.gene.md)
for probe lookup, `expression` is passed to
[`get.gene.expression()`](https://github.com/SpencerTreadway/BioUtils/reference/get.gene.expression.md),
and `phenotype` is passed to
[`build.analysis.df()`](https://github.com/SpencerTreadway/BioUtils/reference/build.analysis.df.md)
for group labeling.

## Examples

``` r
if (FALSE) { # \dontrun{
eset <- load.geo.soft("GDS3268.soft")
geo <- extract.expression(eset)
head(geo$phenotype)
dim(geo$expression)
} # }
```
