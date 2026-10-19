# Create an MSdev object

Create a new MSdev object for mass spectrometry data analysis. This
function initializes the project structure, detects raw data format, and
creates sample information from the raw data directory.

## Usage

``` r
MSdev(rawDataDir, projectDir = dirname(rawDataDir), experimentInfo = MS_Exp())
```

## Arguments

- rawDataDir:

  path to directory containing raw mass spectrometry data files.
  Required. Supported extensions: `.wiff`, `.raw`, `.mzXML`, `.mzML`,
  `.lcd`; matching is case-insensitive.

- projectDir:

  project directory for storing processed data, defaults to parent of
  rawDataDir. This is the single relocatable project root. Other project
  paths are also stored relative to it in `projectInfo$paths`.

- experimentInfo:

  MS_Exp object containing experiment metadata

## Value

MSdev object
