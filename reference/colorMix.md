# Mix Colors with Alpha Blending

Mixes multiple colors together using RGB values weighted by their alpha
(transparency) channels. Useful for creating composite colors from
overlapping transparent elements.

## Usage

``` r
colorMix(...)
```

## Arguments

- ...:

  One or more color specifications (character strings or color names).
  NA values are treated as transparent.

## Value

A single color string in "#RRGGBBAA" format representing the mixed
color.
