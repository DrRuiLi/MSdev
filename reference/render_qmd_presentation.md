# Render a Quarto presentation and serve it on a local port

Starts
[`quarto::quarto_preview()`](https://quarto-dev.github.io/quarto-r/reference/quarto_preview.html)
so a revealjs (or html) slide deck is available locally. With
`frp = TRUE` the preview is also published through FRP using
`.cursor/frp.toml`. The R session stays usable. Call
`stop_qmd_presentation()` to stop the server.

Stops the FRP client (if running) and the Quarto preview, then kills any
leftover Quarto/Deno process still listening on `port`.

## Usage

``` r
render_qmd_presentation(
  qmd,
  port = NULL,
  host = "0.0.0.0",
  render = "auto",
  watch = TRUE,
  browse = TRUE,
  quiet = FALSE,
  frp = TRUE
)

stop_qmd_presentation(port = NULL, quiet = FALSE)
```

## Arguments

- qmd:

  Path to a `.qmd` file, or a directory that contains one. If a
  directory has several `.qmd` files, a name matching `presentation` is
  preferred.

- port:

  Local TCP port to free. If `NULL` (default), uses `remotePort` from
  `.cursor/frp.toml` when that file exists.

- host:

  Hostname to bind. Default `"0.0.0.0"` so other PCs on the LAN can
  connect. Use `"127.0.0.1"` for local-only access.

- render:

  Passed to
  [`quarto::quarto_preview()`](https://quarto-dev.github.io/quarto-r/reference/quarto_preview.html).
  For a document, `TRUE` / `"auto"` renders before serving; `FALSE`
  skips the initial render.

- watch:

  Re-render when the source changes (default `TRUE`).

- browse:

  Open a browser (default `TRUE`).

- quiet:

  Suppress messages.

- frp:

  If `TRUE` (default), start `frpc` using `.cursor/frp.toml`.

## Value

Invisibly, the preview URL (character).

## See also

[`quarto_preview`](https://quarto-dev.github.io/quarto-r/reference/quarto_preview.html)

## Examples

``` r
if (FALSE) { # \dontrun{
render_qmd_presentation("path/to/slides.qmd")
stop_qmd_presentation()
} # }
```
