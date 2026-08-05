# Generate sample information table from raw data files

Read raw MS data files from a directory and generate a sample
information data frame.

## Usage

``` r
get_MS_sampleinfo(raw.data.dir, rawDataFormat = ".raw", verbose = T)
```

## Arguments

- raw.data.dir:

  Path to directory containing raw data files

- rawDataFormat:

  File extension of raw data files (default `".raw"`; also `".mzXML"`,
  `".mzML"`, `".wiff"`, `".lcd"`; matching is case-insensitive)

- verbose:

  Logical indicating whether to print messages

## Value

data.frame with sample information
