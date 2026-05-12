# Analyze Gene Expression Between Groups

Performs a comprehensive statistical analysis comparing gene expression
between two groups. Integrates parametric and nonparametric testing,
effect size estimation, bootstrapped confidence intervals, and a
heuristic biological relevance assessment.

## Usage

``` r
analyze.gene(df, alpha = 0.05, n.boot = 1000)
```

## Arguments

- df:

  A data frame as produced by
  [`build.analysis.df()`](https://github.com/SpencerTreadway/BioUtils/reference/build.analysis.df.md),
  containing at least two columns:

  expression

  :   Numeric vector of gene expression values.

  group

  :   Character or factor vector with exactly two group labels.

  When the data frame contains multiple genes (i.e., a `gene` column
  with more than one unique value), expression values from all genes are
  pooled together. Subset the data frame to a single gene before calling
  this function for per-gene results.

- alpha:

  Numeric. Significance threshold. Default is `0.05`.

- n.boot:

  Integer. Number of bootstrap resamples for the confidence interval.
  Default is `1000`.

## Value

A named list with nine elements:

- p.value:

  Numeric. P-value from the adaptive t-test.

- effect.size:

  Numeric. Cohen's d.

- effect.size.class:

  Character. Qualitative magnitude label from
  [`classify.effect.size()`](https://github.com/SpencerTreadway/BioUtils/reference/classify.effect.size.md).

- confidence.interval:

  Named list with `lower` and `upper` bounds for the bootstrapped CI
  around Cohen's d.

- nonparametric.p:

  Numeric. P-value from the Wilcoxon rank-sum test.

- robustness:

  Character. Agreement assessment between parametric and nonparametric
  tests.

- biological.relevance:

  Character. Heuristic relevance label from
  [`flag.biological.relevance()`](https://github.com/SpencerTreadway/BioUtils/reference/flag.biological.relevance.md).

- interpretation:

  Character. Human-readable summary combining all of the above.

- raw:

  List. Raw output from the parametric and nonparametric tests for
  further inspection.

## Details

This function integrates multiple statistical perspectives to provide a
more complete picture of group differences than a p-value alone. The
result is designed to feed directly into `plot.gene.analysis()` for
visualization.

## Examples

``` r
# \donttest{
geo <- extract.expression(load.geo.soft(accession = "GDS3268", log.transform = TRUE))
#> GDS3268 not found locally, downloading from NCBI GEO...
#> Warning: NaNs produced
probe <- find.probe.by.gene(geo$gene, "mucin 20, cell surface associated")
expr <- get.gene.expression(geo$expression, probe)
df <- build.analysis.df(expr, geo$phenotype, geo$gene)
result <- analyze.gene(df)
cat(result$interpretation)
#> The difference is not statistically significant with a negligible effect size.
#> Result is no consistent evidence of difference.
#> No strong evidence of biological relevance.
# }
```
