# Changelog

## MSdev 3.1.0

- Project paths can move with the project directory. `projectDir` is the
  only machine-specific root. Other project paths are also stored
  relative to it in `projectInfo$paths`.
- [`MSdev()`](https://drruili.github.io/MSdev/reference/MSdev.md) now
  requires `rawDataDir`. The previous personal default path was removed.
- [`MSdev_relocate()`](https://drruili.github.io/MSdev/reference/MSdev_relocate.md)
  rewrites absolute paths inside an existing object (sample files,
  on-disk spectra and chromatograms, and xcms experiment file paths)
  after the project folder has been moved.
  [`MSdev_load()`](https://drruili.github.io/MSdev/reference/MSdev_IO.md)
  does this automatically when the saved `projectDir` no longer exists.
- CompoundDb is resolved from the function argument, a stored path that
  still exists, `options(MSdev.compoundDb)`, then the `MSDEV_COMPOUNDDB`
  environment variable. Annotation functions no longer default to a
  hardcoded database path.
- Objects saved before 3.1.0 still load on the original machine. Call
  [`MSdev_relocate()`](https://drruili.github.io/MSdev/reference/MSdev_relocate.md)
  after copying a project while the old directory is still present.

## MSdev 3.0.0

- Initial CRAN submission.
