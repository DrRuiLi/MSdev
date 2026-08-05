# Create an MSdev object

Create a new MSdev object for mass spectrometry data analysis. This
function initializes the project structure, detects raw data format, and
creates sample information from the raw data directory.

## Usage

``` r
MSdev(
  rawDataDir = "C:/Users/91879/OneDrive/Code/R/Projecct/2022.1.8_MS.demo/Demo/raw.data",
  projectDir = dirname(rawDataDir),
  experimentInfo = MS_Exp()
)
```

## Arguments

- rawDataDir:

  path to directory containing raw mass spectrometry data files
  (supported extensions: `.wiff`, `.raw`, `.mzXML`, `.mzML`, `.lcd`;
  matching is case-insensitive)

- projectDir:

  project directory for storing processed data, defaults to parent of
  rawDataDir

- experimentInfo:

  MS_Exp object containing experiment metadata

## Value

MSdev object
