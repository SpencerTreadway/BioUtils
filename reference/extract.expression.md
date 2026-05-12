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
# \donttest{
eset <- load.geo.soft(accession = "GDS507", log.transform = TRUE)
#> GDS507 not found locally, downloading from NCBI GEO...
#> Using locally cached version of GDS507 found here:
#> /tmp/RtmpYU9fuS/GDS507.soft.gz 
#> Using locally cached version of GPL97 found here:
#> /tmp/RtmpYU9fuS/GPL97.annot.gz 
geo <- extract.expression(eset)
head(geo$phenotype)
#>            sample disease.state individual
#> GSM11815 GSM11815           RCC        035
#> GSM11832 GSM11832           RCC        023
#> GSM12069 GSM12069           RCC        001
#> GSM12083 GSM12083           RCC        005
#> GSM12101 GSM12101           RCC        011
#> GSM12106 GSM12106           RCC        032
#>                                                                                                                                   description
#> GSM11815 Value for GSM11815: C035 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
#> GSM11832 Value for GSM11832: C023 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
#> GSM12069 Value for GSM12069: C001 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
#> GSM12083 Value for GSM12083: C005 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
#> GSM12101 Value for GSM12101: C011 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
#> GSM12106 Value for GSM12106: C032 Renal Clear Cell Carcinoma U133B; src: Trizol isolation of total RNA from Renal Clear Cell Carcinoma tissue
dim(geo$expression)
#> [1] 22645    17
# }
```
