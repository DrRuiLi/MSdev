# Set RStudio Files Pane Directory

Navigates the RStudio files pane to the specified directory or the
directory of the currently active file.

## Usage

``` r
setRStudioDir(path = rstudioapi::getSourceEditorContext()$path)
```

## Arguments

- path:

  Path to a file or directory. Defaults to the path of the currently
  active file in the RStudio editor.

## Value

NULL (invisibly)
