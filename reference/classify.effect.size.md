# Classify Effect Size

Categorizes an absolute Cohen's d value into a qualitative magnitude
label using conventional thresholds.

## Usage

``` r
classify.effect.size(d)
```

## Arguments

- d:

  Numeric. Cohen's d effect size value (signed or unsigned).

## Value

Character string. One of `"negligible"`, `"small"`, `"moderate"`, or
`"large"` based on the absolute value of `d`, using the thresholds:
negligible \< 0.2 \<= small \< 0.5 \<= moderate \< 0.8 \<= large.

## Examples

``` r
# \donttest{
classify.effect.size(0.3)   # "small"
#> [1] "small"
classify.effect.size(-0.9)  # "large"
#> [1] "large"
# }
```
