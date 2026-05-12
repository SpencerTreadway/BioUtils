# Assess Biological Relevance of Gene Expression Differences

Provides a heuristic interpretation of whether a statistically
significant difference in gene expression is likely to be biologically
meaningful, based on effect size magnitude and p-value.

## Usage

``` r
flag.biological.relevance(effect.size, p.value, alpha = 0.05)
```

## Arguments

- effect.size:

  Numeric. Cohen's d or other standardized effect size.

- p.value:

  Numeric. P-value from a statistical test.

- alpha:

  Numeric. Significance threshold used for the p-value comparison.
  Default is `0.05`.

## Value

Character string. One of three heuristic labels:

- "Potentially biologically meaningful":

  p \< alpha and \|d\| \> 0.5

- "Statistically significant but small effect":

  p \< alpha and \|d\| \<= 0.5

- "No strong evidence of biological relevance":

  p \>= alpha

## Details

These rules are heuristic and intended as guidance rather than
definitive biological conclusions. The thresholds for effect size (0.5)
and significance are conventional starting points; domain knowledge
should inform interpretation. This function is called internally by
[`analyze.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/analyze.gene.md)
using the same `alpha` passed to that function.

## Examples

``` r
# \donttest{
flag.biological.relevance(effect.size = 0.6, p.value = 0.01)
#> [1] "Potentially biologically meaningful"
# }
```
