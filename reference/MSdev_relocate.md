# Point an MSdev object at a new project directory

Rewrite absolute paths that sit under the saved `projectDir` so they sit
under `projectDir`. Sample files, `onDiskData` locations, Spectra
`dataOrigin` / `dataStorage`, and xcms experiment file paths are
updated. Files are not moved; the project tree is assumed to already be
at the new location.

`rawDataDir` overrides the raw-data directory when it is not inside the
project directory. `compoundDb` sets `projectInfo$CompoundDB_path`
without prefix replacement. After rewriting, `projectInfo$paths` holds
paths relative to `projectDir`.

## Usage

``` r
MSdev_relocate(object, projectDir, rawDataDir = NULL, compoundDb = NULL)
```

## Arguments

- object:

  MSdev object

- projectDir:

  New project directory.

- rawDataDir:

  Optional new raw-data directory. Used when raw data does not live
  under `projectDir`.

- compoundDb:

  Optional CompoundDb SQLite file for this machine.

## Value

The same MSdev object with updated paths.
