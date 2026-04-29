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
if (FALSE) { # \dontrun{
geo <- extract.expression(load.geo.soft("GDS3268.soft"))
de.results <- run.limma.de(geo, condition.col = "disease.state")
head(de.results)
} # }
```
