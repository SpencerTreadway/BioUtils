# Load GEO Dataset from SOFT File

Imports a GEO dataset and converts it into a structured ExpressionSet
object. If a local SOFT file is found at `file.path` it is loaded
directly. Otherwise the dataset is downloaded from NCBI GEO using
`accession` and cached in
[`tempdir()`](https://rdrr.io/r/base/tempfile.html).

## Usage

``` r
load.geo.soft(file.path = NULL, accession = NULL, log.transform = FALSE)
```

## Arguments

- file.path:

  Character. Path to a local GEO SOFT file. If the file does not exist
  and `accession` is provided, the dataset is downloaded automatically.

- accession:

  Character. GEO accession ID (e.g., `"GDS507"`). Only used when
  `file.path` does not exist. Default is `NULL`.

- log.transform:

  Logical. Whether to apply a log2 transformation to the expression
  matrix during import. Default is `FALSE`.

## Value

An `ExpressionSet` object for use with
[`extract.expression()`](https://github.com/SpencerTreadway/BioUtils/reference/extract.expression.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Load from a local file
eset <- load.geo.soft("GDS507.soft")

# Download automatically if not found locally
eset <- load.geo.soft("GDS507.soft", accession = "GDS507", log.transform = TRUE)

# Download without a local file at all
eset <- load.geo.soft(NULL, accession = "GDS507", log.transform = TRUE)
} # }
```
