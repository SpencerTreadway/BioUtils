# Find Probe IDs by Gene Name

Searches the gene annotation data frame to find the probe IDs
corresponding to one or more gene names or descriptions.

## Usage

``` r
find.probe.by.gene(genes, gene.names)
```

## Arguments

- genes:

  Data frame of gene annotations as returned by
  `extract.expression()$gene`. The second column (index 2) is expected
  to contain gene titles or descriptions.

- gene.names:

  Character vector of one or more gene names or descriptions to search
  for. Matching is exact and case-sensitive.

## Value

Integer vector of probe IDs corresponding to the matched genes. The
length of the vector equals the number of matches found. Returns an
empty integer vector if no matches are found.

## Details

Probe IDs are the integer row names of the `genes` annotation data
frame. These IDs are used directly to index rows of the expression
matrix in
[`get.gene.expression()`](https://github.com/SpencerTreadway/BioUtils/reference/get.gene.expression.md).
If multiple gene names are supplied, all matching probe IDs are returned
as a vector, which is then passed as-is to
[`get.gene.expression()`](https://github.com/SpencerTreadway/BioUtils/reference/get.gene.expression.md)
to retrieve a multi-row expression matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
geo <- extract.expression(load.geo.soft("GDS3268.soft"))

# Single gene
probe <- find.probe.by.gene(geo$gene, "mucin 20, cell surface associated")

# Multiple genes
probes <- find.probe.by.gene(geo$gene, c(
  "mucin 20, cell surface associated",
  "alcohol dehydrogenase 1A (class I), alpha polypeptide"
))
} # }
```
