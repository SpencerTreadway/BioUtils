# Compute Confidence Interval for Effect Size

Estimates a confidence interval for Cohen's d using bootstrap
resampling. This provides a measure of uncertainty around the effect
size estimate, which is especially useful when sample sizes are small or
distributions deviate from normality.

## Usage

``` r
compute.ci(df, alpha = 0.05, n.boot = 1000)
```

## Arguments

- df:

  A data frame containing at least two columns:

  expression

  :   Numeric vector of gene expression values.

  group

  :   Character or factor vector indicating group membership.

- alpha:

  Numeric. Significance level used to compute the confidence interval
  bounds. Default is `0.05`, yielding a 95 percent CI.

- n.boot:

  Integer. Number of bootstrap resamples to perform. Default is `1000`.

## Value

A named list with two elements:

- lower:

  Lower bound of the bootstrapped confidence interval.

- upper:

  Upper bound of the bootstrapped confidence interval.

## Details

The function repeatedly resamples the input data with replacement and
recomputes Cohen's d for each resample. The confidence interval is
derived from the empirical distribution of bootstrapped effect sizes
using the percentile method.

## Examples

``` r
# \donttest{
# Build a minimal example data frame
analysis.df <- data.frame(
  expression = c(1.2, 2.3, 1.8, 2.1, 3.4, 2.9, 3.1, 2.7),
  group = rep(c("normal", "RCC"), each = 4)
)
ci <- compute.ci(analysis.df, alpha = 0.05, n.boot = 100)
# }
```
