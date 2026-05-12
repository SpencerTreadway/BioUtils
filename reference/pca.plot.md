# Plot PCA of Sample Expression Profiles

Reduces the dimensionality of the full expression matrix using Principal
Component Analysis and plots samples in the first two principal
component axes, colored by a phenotype variable. Useful for quality
control, batch effect detection, and exploratory assessment of sample
clustering.

## Usage

``` r
pca.plot(
  expression.matrix,
  phenotype,
  color.by = "disease.state",
  scale = TRUE
)
```

## Arguments

- expression.matrix:

  Matrix of gene expression values. Rows are probes, columns are
  samples.

- phenotype:

  Data frame of sample metadata, as returned by
  `extract.expression()$phenotype`.

- color.by:

  Character. Column name in `phenotype` to use for point coloring.
  Default is `"disease.state"`.

- scale:

  Logical. Whether to scale probes to unit variance before PCA. Default
  is `TRUE`.

## Value

A ggplot object showing samples as points in PC1 vs PC2 space, colored
by the specified phenotype variable. The percentage of variance
explained is shown on each axis label.

## Details

PCA is typically applied at the whole-dataset level before any gene
filtering, making it a useful first diagnostic after
[`load.geo.soft()`](https://github.com/SpencerTreadway/BioUtils/reference/load.geo.soft.md)
and
[`extract.expression()`](https://github.com/SpencerTreadway/BioUtils/reference/extract.expression.md).
Tight clustering of samples by disease state in PC1/PC2 space suggests a
strong and detectable expression signal. Unexpected clustering by batch,
sex, or other covariates can indicate the need for covariate adjustment
in downstream
[`run.limma.de()`](https://github.com/SpencerTreadway/BioUtils/reference/run.limma.de.md)
models.

## Examples

``` r
# \donttest{
data <- extract.expression(load.geo.soft(accession = "GDS3268", log.transform = TRUE))
#> GDS3268 not found locally, downloading from NCBI GEO...
#> Using locally cached version of GDS3268 found here:
#> /tmp/RtmpYU9fuS/GDS3268.soft.gz 
#> Warning: NaNs produced
#> Using locally cached version of GPL1708 found here:
#> /tmp/RtmpYU9fuS/GPL1708.annot.gz 
pca.plot(data$expression, data$phenotype, color.by = "disease.state")
#> Error in prcomp.default(t(expression.matrix), scale. = scale): cannot rescale a constant/zero column to unit variance
# }
```
