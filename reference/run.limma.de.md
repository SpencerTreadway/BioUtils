# Run Differential Expression Analysis Using limma

Performs genome-wide differential expression analysis across all probes
using the limma linear modeling framework. This is the recommended
approach for microarray data from GEO, as it borrows variance
information across genes to produce more stable estimates than
gene-by-gene t-tests.

## Usage

``` r
run.limma.de(geo, condition.col = "disease.state", adjust.method = "BH")
```

## Arguments

- geo:

  Named list as returned by
  [`extract.expression()`](https://github.com/SpencerTreadway/BioUtils/reference/extract.expression.md),
  containing at minimum `$expression` (probe x sample matrix) and
  `$phenotype` (sample metadata data frame).

- condition.col:

  Character. Name of the column in `geo$phenotype` to use as the
  grouping variable. Default is `"disease.state"`.

- adjust.method:

  Character. P-value correction method passed to
  [`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html).
  Default is `"BH"` (Benjamini-Hochberg FDR).

## Value

A data frame (TopTable) with one row per probe, sorted by evidence of
differential expression, containing:

- logFC:

  Log2 fold-change between conditions.

- AveExpr:

  Average log2 expression across all samples.

- t:

  Moderated t-statistic.

- P.Value:

  Raw p-value.

- adj.P.Val:

  FDR-adjusted p-value (BH by default).

- B:

  Log-odds that the gene is differentially expressed.

## Details

Unlike running
[`analyze.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/analyze.gene.md)
on each probe individually, limma fits a linear model to all probes
simultaneously and applies empirical Bayes shrinkage to the variance
estimates. This stabilizes results for genes with few observations and
is the standard method for microarray DE analysis.

The returned TopTable is the primary input for downstream functions such
as
[`volcano.plot()`](https://github.com/SpencerTreadway/BioUtils/reference/volcano.plot.md)
and
[`run.gsea()`](https://github.com/SpencerTreadway/BioUtils/reference/run.gsea.md),
and can be filtered by `adj.P.Val` and `logFC` thresholds to define a
gene signature.

## Examples

``` r
# \donttest{
geo <- extract.expression(load.geo.soft(accession = "GDS3268", log.transform = TRUE))
#> GDS3268 not found locally, downloading from NCBI GEO...
#> Using locally cached version of GDS3268 found here:
#> /tmp/RtmpxRZSjV/GDS3268.soft.gz 
#> Warning: NaNs produced
#> Using locally cached version of GPL1708 found here:
#> /tmp/RtmpxRZSjV/GPL1708.annot.gz 
de.results <- run.limma.de(geo, condition.col = "disease.state")
#> Warning: Partial NA coefficients for 2141 probe(s)
#> Removing intercept from test coefficients
head(de.results)
#>            logFC   AveExpr         t      P.Value    adj.P.Val         B
#> 16982  0.5350351 -1.152244  6.134758 4.381867e-09 0.0001027521 10.435447
#> 13087 -1.5601595 -3.766461 -6.318251 5.868863e-09 0.0001027521 10.059365
#> 5334   0.8320762 -1.946530  6.013908 1.160081e-08 0.0001354047  9.512088
#> 13277  0.9150156 -1.580899  5.582388 7.648712e-08 0.0006695683  7.737907
#> 6416   0.9928177 -2.532741  5.493257 1.281664e-07 0.0007491880  7.258928
#> 14638  1.4593609 -3.346868  5.530253 1.283735e-07 0.0007491880  7.253982
# }
```
