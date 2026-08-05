# featureCandidate

match feature with database by mz, return all spectra matched in a list
splited by feature

## Usage

``` r
featureCandidate(
  object,
  mz.ppm = 20,
  spectraDatabase =
    "C:\\Users\\91879\\OneDrive\\Documents\\Code\\R\\Projecct\\2022.1.17_Compounds.database\\Spectra.integrated.database.integration.2022_02_12.Rdata"
)
```

## Arguments

- object:

  a `MSdev` object

- mz.ppm:

  mz error

- spectraDatabase:

  databse path

## Value

a `MSdev` object
