# Compute Effect Size (Cohen's d)

Calculates the standardized mean difference between two groups.

## Usage

``` r
compute.effect.size(df)
```

## Arguments

- df:

  A data frame containing at least two columns:

  expression

  :   Numeric vector of gene expression values.

  group

  :   Character or factor vector with exactly two group labels, as
      produced by
      [`build.analysis.df()`](https://github.com/SpencerTreadway/BioUtils/reference/build.analysis.df.md).

## Value

Numeric. Cohen's d value. Positive values indicate the first group
(alphabetically) has higher mean expression; negative values indicate
the second group is higher.
