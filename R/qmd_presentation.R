#' Render a Quarto presentation and serve it on a local port
#'
#' Starts \code{quarto::quarto_preview()} so a revealjs (or html) slide
#' deck is available locally. With \code{frp = TRUE} the preview is
#' also published through FRP using \code{.cursor/frp.toml}. The R
#' session stays usable. Call \code{\link{qmd_stop_presentation}()}
#' to stop the server.
#'
#' @param qmd Path to a \code{.qmd} file, or a directory that contains one.
#'   If a directory has several \code{.qmd} files, a name matching
#'   \code{presentation} is preferred.
#' @param port Local TCP port. If \code{NULL} (default), uses
#'   \code{remotePort} from \code{.cursor/frp.toml}. Use \code{"auto"}
#'   to let Quarto pick a free port.
#' @param host Hostname to bind. Default \code{"0.0.0.0"} so other PCs
#'   on the LAN can connect. Use \code{"127.0.0.1"} for local-only
#'   access.
#' @param render Passed to \code{quarto::quarto_preview()}. For a document,
#'   \code{TRUE} / \code{"auto"} renders before serving; \code{FALSE}
#'   skips the initial render.
#' @param watch Re-render when the source changes (default \code{TRUE}).
#' @param browse Open a browser (default \code{TRUE}).
#' @param quiet Suppress messages.
#' @param frp If \code{TRUE} (default), start \code{frpc} using
#'   \code{.cursor/frp.toml}.
#'
#' @return Invisibly, the preview URL (character).
#'
#' @examples
#' \dontrun{
#' qmd_render_presentation("path/to/slides.qmd")
#' qmd_stop_presentation()
#' }
#'
#' @seealso \code{\link[quarto]{quarto_preview}}
#' @export
qmd_render_presentation <- function(qmd,
                                    port = NULL,
                                    host = "0.0.0.0",
                                    render = "auto",
                                    watch = TRUE,
                                    browse = TRUE,
                                    quiet = FALSE,
                                    frp = TRUE) {
  if (!requireNamespace("quarto", quietly = TRUE)) {
    stop("Package 'quarto' is required for qmd_render_presentation()")
  }
  qmd <- .msdev_resolve_qmd(qmd)
  if (is.null(port) && isTRUE(frp)) {
    port <- .msdev_frp_settings()$remote_port
  }
  port <- .msdev_preview_port(port)
  host <- as.character(host)[[1]]
  if (!nzchar(host)) {
    stop("`host` must be a non-empty hostname")
  }
  if (isTRUE(render)) {
    render <- "auto"
  }

  qmd_stop_presentation(port, quiet = TRUE)

  lan_bind <- host %in% c("0.0.0.0", "::")
  url <- quarto::quarto_preview(
    file = qmd,
    render = render,
    port = port,
    host = host,
    browse = if (lan_bind) FALSE else browse,
    watch = isTRUE(watch),
    navigate = if (lan_bind) FALSE else TRUE,
    quiet = isTRUE(quiet)
  )
  local_url <- if (lan_bind) {
    sprintf("http://127.0.0.1:%s/", port)
  } else {
    as.character(url)
  }
  if (isTRUE(browse) && lan_bind) {
    utils::browseURL(local_url)
  }
  public_url <- NULL
  if (isTRUE(frp) && !identical(port, "auto")) {
    public_url <- .msdev_frp_try_start(port, host)
  }
  if (!quiet) {
    message("Local: ", local_url)
    lan <- .msdev_lan_ipv4()
    if (lan_bind && length(lan)) {
      message("LAN:   ", paste(sprintf("http://%s:%s/", lan, port), collapse = "  "))
    } else if (!identical(as.character(url), local_url)) {
      message("Preview: ", url)
    }
    if (!is.null(public_url)) {
      message("Public: ", public_url)
    }
    message("Stop with qmd_stop_presentation()")
  }
  invisible(local_url)
}

#' Stop the preview server
#'
#' Stops the FRP client (if running) and the Quarto preview, then kills
#' any leftover Quarto/Deno process still listening on \code{port}.
#'
#' @param port Local TCP port to free. If \code{NULL} (default), uses
#'   \code{remotePort} from \code{.cursor/frp.toml} when that file exists.
#' @param quiet Suppress messages.
#'
#' @rdname qmd_render_presentation
#' @export
qmd_stop_presentation <- function(port = NULL, quiet = FALSE) {
  .msdev_frp_stop(quiet = TRUE)
  if (is.null(port)) {
    port <- tryCatch(.msdev_frp_settings()$remote_port, error = function(e) NULL)
  }
  if (requireNamespace("quarto", quietly = TRUE)) {
    try(quarto::quarto_preview_stop(), silent = TRUE)
  }
  if (!is.null(port) && !identical(port, "auto")) {
    .msdev_stop_listen_port(port, quiet = quiet)
  }
  invisible(TRUE)
}

#' @keywords internal
.msdev_resolve_qmd <- function(qmd) {
  if (length(qmd) != 1L || is.na(qmd) || !nzchar(as.character(qmd)[[1]])) {
    stop("`qmd` must be a path to a .qmd file or a directory")
  }
  qmd <- path.expand(as.character(qmd)[[1]])
  if (dir.exists(qmd)) {
    files <- list.files(qmd, pattern = "\\.[Qq]md$", full.names = TRUE)
    if (!length(files)) {
      stop("No .qmd file found in: ", qmd)
    }
    pref <- files[grepl("presentation", basename(files), ignore.case = TRUE)]
    if (length(pref)) {
      files <- pref
    }
    if (length(files) > 1L) {
      message("Multiple .qmd files; using ", basename(files[[1]]))
    }
    qmd <- files[[1]]
  }
  if (!file.exists(qmd)) {
    stop("File not found: ", qmd)
  }
  normalizePath(qmd, winslash = "/", mustWork = TRUE)
}

#' @keywords internal
.msdev_lan_ipv4 <- function() {
  if (.Platform$OS.type == "windows") {
    out <- suppressWarnings(system(
      paste(
        "powershell -NoProfile -Command",
        "\"Get-NetIPAddress -AddressFamily IPv4 |",
        "Where-Object { $_.PrefixOrigin -eq 'Dhcp' -and $_.AddressState -eq 'Preferred' } |",
        "Select-Object -ExpandProperty IPAddress\""
      ),
      intern = TRUE
    ))
    ip <- trimws(as.character(out))
    return(ip[nzchar(ip) & grepl("^[0-9.]+$", ip)])
  }
  raw <- tryCatch(system("ipconfig", intern = TRUE), error = function(e) character(0))
  lines <- grep("IPv4 Address", raw, value = TRUE, ignore.case = TRUE)
  ip <- sub(".*: *", "", lines)
  ip <- trimws(ip)
  ip[grepl("^(10\\.|192\\.168\\.|172\\.(1[6-9]|2[0-9]|3[0-1])\\.)", ip)]
}

#' @keywords internal
.msdev_pids_listening <- function(port) {
  port <- as.integer(port)[[1]]
  lines <- suppressWarnings(system("netstat -ano -p TCP", intern = TRUE))
  if (!length(lines)) {
    return(integer(0))
  }
  listen <- grep("LISTENING", lines, value = TRUE, ignore.case = TRUE)
  re <- sprintf("[:.]%d\\s+\\S+\\s+LISTENING\\s+(\\d+)", port)
  pids <- integer(0)
  for (ln in listen) {
    m <- regexec(re, ln, perl = TRUE, ignore.case = TRUE)
    cap <- regmatches(ln, m)[[1]]
    if (length(cap) >= 2L) {
      pids <- c(pids, as.integer(cap[[2]]))
    }
  }
  unique(pids[!is.na(pids) & pids > 0L])
}

#' @keywords internal
.msdev_win_process <- function(pid) {
  pid <- as.integer(pid)[[1]]
  ps <- tempfile(fileext = ".ps1")
  on.exit(unlink(ps), add = TRUE)
  writeLines(
    sprintf(
      "$p = Get-CimInstance Win32_Process -Filter \"ProcessId=%d\" -ErrorAction SilentlyContinue; if (-not $p) { exit 0 }; $par = Get-CimInstance Win32_Process -Filter \"ProcessId=$($p.ParentProcessId)\" -ErrorAction SilentlyContinue; Write-Output ($p.Name + \"`t\" + $p.ProcessId + \"`t\" + $p.ParentProcessId + \"`t\" + $par.Name + \"`t\" + $p.CommandLine)",
      pid
    ),
    ps
  )
  out <- suppressWarnings(system2(
    "powershell",
    c("-NoProfile", "-File", ps),
    stdout = TRUE,
    stderr = TRUE
  ))
  out <- out[nzchar(trimws(out))]
  if (!length(out)) {
    return(NULL)
  }
  parts <- strsplit(out[[1]], "\t", fixed = TRUE)[[1]]
  list(
    name = parts[[1]],
    pid = as.integer(parts[[2]]),
    parent_id = suppressWarnings(as.integer(parts[[3]])),
    parent_name = if (length(parts) >= 4L) parts[[4]] else NA_character_,
    cmd = if (length(parts) >= 5L) parts[[5]] else ""
  )
}

#' @keywords internal
.msdev_stop_listen_port <- function(port, quiet = FALSE) {
  port <- as.integer(port)[[1]]
  pids <- .msdev_pids_listening(port)
  if (!length(pids)) {
    return(invisible(integer(0)))
  }
  killed <- integer(0)
  for (pid in pids) {
    info <- .msdev_win_process(pid)
    blob <- if (is.null(info)) {
      ""
    } else {
      paste(info$name, info$parent_name, info$cmd)
    }
    if (!grepl("quarto|deno", blob, ignore.case = TRUE)) {
      if (!quiet) {
        warning(
          "Port ", port, " is in use by PID ", pid,
          if (!is.null(info)) paste0(" (", info$name, ")"),
          "; not killing",
          call. = FALSE
        )
      }
      next
    }
    targets <- pid
    if (!is.null(info) && grepl("quarto", info$parent_name, ignore.case = TRUE)) {
      targets <- c(info$parent_id, pid)
    }
    for (t in unique(targets[!is.na(targets) & targets > 0L])) {
      suppressWarnings(system(
        sprintf("taskkill /PID %d /T /F", t),
        intern = TRUE,
        ignore.stderr = TRUE
      ))
      killed <- c(killed, t)
    }
  }
  Sys.sleep(0.5)
  if (!quiet && length(killed)) {
    message("Stopped preview on port ", port, " (PID ", paste(unique(killed), collapse = ", "), ")")
  }
  invisible(unique(killed))
}

#' @keywords internal
.msdev_preview_port <- function(port) {
  if (is.null(port) || (is.character(port) && identical(tolower(port[[1]]), "auto"))) {
    return("auto")
  }
  port <- suppressWarnings(as.integer(port)[[1]])
  if (!is.finite(port) || port < 1L || port > 65535L) {
    stop("`port` must be an integer in 1..65535, or 'auto'")
  }
  port
}
