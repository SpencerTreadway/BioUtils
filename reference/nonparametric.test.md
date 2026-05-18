# Perform Nonparametric Group Comparison

Conducts a Wilcoxon rank-sum test (Mann-Whitney U test) to compare
expression values between two groups without assuming normality.

## Usage

``` r
nonparametric.test(df)
```

## Arguments

- df:

  A data frame containing at least two columns:

  expression

  :   Numeric vector of gene expression values.

  group

  :   Character or factor vector with exactly two group labels.

## Value

A named list with two elements:

- p.value:

  Numeric. P-value from the Wilcoxon rank-sum test.

- test:

  Character. Name of the test used (`"Wilcoxon rank-sum"`).

## Details

This test is useful as a robustness check against parametric methods
such as the t-test, especially when data are skewed or contain outliers.
The result is used by
[`analyze.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/analyze.gene.md)
to populate the `robustness` field of its return value.

## Examples

``` r
# \donttest{
analysis.df <- data.frame(
  expression = c(1.2, 2.3, 1.8, 2.1, 3.4, 2.9, 3.1, 2.7,
                 1.5, 2.0, 1.9, 2.4),
  group      = rep(c("normal", "RCC"), each = 6)
)
res <- nonparametric.test(analysis.df)
res$p.value
#> [1] 1
# }
```
