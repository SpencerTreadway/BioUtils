# Resolve Probe IDs to Gene Names

Converts one or more integer probe IDs back to their human-readable gene
title strings or short gene symbols using the gene annotation data
frame. This is the inverse operation of
[`find.probe.by.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/find.probe.by.gene.md).

## Usage

``` r
get.gene.name(genes, probe.id, use.symbols = FALSE)
```

## Arguments

- genes:

  Data frame of gene annotations as returned by
  `extract.expression()$gene`. Must contain columns `"ID"`,
  `"Gene title"`, and `"Gene symbol"`.

- probe.id:

  Integer or character vector of one or more probe IDs to resolve, as
  returned by
  [`find.probe.by.gene()`](https://github.com/SpencerTreadway/BioUtils/reference/find.probe.by.gene.md).

- use.symbols:

  Logical. If `FALSE` (default), returns full descriptive gene titles
  from the `"Gene title"` column, suitable for interpretation and
  reporting. If `TRUE`, returns short gene symbols from the
  `"Gene symbol"` column, suitable for plot axis labels and any
  downstream tool that expects standard gene symbols such as
  [`run.gsea()`](https://github.com/SpencerTreadway/BioUtils/reference/run.gsea.md).

## Value

Character vector of gene names corresponding to the supplied probe IDs.
Returns empty strings `""` for probes with no annotation, which should
be filtered with `which(result != "")`.

## Examples

``` r
# \donttest{
geo <- extract.expression(load.geo.soft(accession = "GDS507",
                                        log.transform = TRUE))
#> GDS507 not found locally, downloading from NCBI GEO...
#> Using locally cached version of GDS507 found here:
#> /tmp/RtmpYU9fuS/GDS507.soft.gz 
#> Using locally cached version of GPL97 found here:
#> /tmp/RtmpYU9fuS/GPL97.annot.gz 
de.results <- run.limma.de(geo)
#> Removing intercept from test coefficients
top.probes <- rownames(head(de.results, 10))

# Full titles for reporting
top.titles  <- get.gene.name(geo$gene, top.probes)

# Short symbols for plot labels and GSEA
top.symbols <- get.gene.name(geo$gene, top.probes, use.symbols = TRUE)
# }
```
