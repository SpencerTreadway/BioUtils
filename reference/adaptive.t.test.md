# Perform Adaptive T-Test

Automatically selects between Welch's and Student's t-test based on the
result of a variance equality test, then returns a unified result
structure regardless of which branch was taken.

## Usage

``` r
adaptive.t.test(df, alpha = 0.05)
```

## Arguments

- df:

  A data frame containing at least two columns:

  expression

  :   Numeric vector of gene expression values.

  group

  :   Character or factor vector with exactly two group labels.

- alpha:

  Numeric. Significance level for the variance equality test
  ([`var.test()`](https://rdrr.io/r/stats/var.test.html)). Default is
  `0.05`.

## Value

A named list with three elements:

- variance.test:

  The result object from
  [`var.test()`](https://rdrr.io/r/stats/var.test.html) or
  [`oneway.test()`](https://rdrr.io/r/stats/oneway.test.html), depending
  on the branch taken.

- t.test:

  The result object from
  [`t.test()`](https://rdrr.io/r/stats/t.test.html).

- p.value:

  Numeric. The p-value from the t-test.

## Details

If [`var.test()`](https://rdrr.io/r/stats/var.test.html) returns a
p-value below `alpha`, variances are considered unequal and Welch's
t-test (`var.equal = FALSE`) is used. Otherwise, Student's t-test
(`var.equal = TRUE`) is applied alongside a one-way ANOVA as a
confirmatory check. The returned list is consumed internally by
[`analyze.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/analyze.gene.md).

## Examples

``` r
# \donttest{
result <- adaptive.t.test(df, alpha = 0.05)
#> Error in model.frame.default(formula = expression ~ group, data = df): 'data' must be a data.frame, environment, or list
result$p.value
#> Error: object 'result' not found
# }
```
