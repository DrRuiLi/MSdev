.msdev_frp_env <- new.env(parent = emptyenv())

.msdev_frp_public_url <- function(settings = NULL) {
  if (is.null(settings)) {
    settings <- .msdev_frp_settings()
  }
  host <- settings$public_host
  if (is.null(host) || !nzchar(host)) {
    host <- settings$server_addr
  }
  port <- as.integer(settings$public_port)
  path <- settings$public_path
  if (is.null(path) || !nzchar(path)) {
    path <- "/"
  }
  if (!startsWith(path, "/")) {
    path <- paste0("/", path)
  }
  scheme <- if (identical(port, 443L)) "https" else "http"
  if (identical(port, 80L) || identical(port, 443L)) {
    paste0(scheme, "://", host, path)
  } else {
    paste0(scheme, "://", host, ":", port, path)
  }
}

.msdev_pkg_root <- function() {
  ns <- tryCatch(getNamespaceInfo("MSdev", "path"), error = function(e) NULL)
  if (!is.null(ns) && dir.exists(ns)) {
    return(ns)
  }
  getwd()
}

.msdev_frp_read_config_file <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    return(list())
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- sub("#.*$", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  out <- list()
  for (ln in lines) {
    if (!grepl("=", ln, fixed = TRUE)) {
      next
    }
    parts <- strsplit(ln, "=", fixed = TRUE)[[1L]]
    key <- trimws(parts[[1L]])
    val <- trimws(paste(parts[-1], collapse = "="))
    val <- sub("^['\"]", "", val)
    val <- sub("['\"]$", "", val)
    out[[key]] <- val
  }
  out
}

.msdev_frp_local_config_path <- function() {
  candidates <- unique(c(
    file.path(.msdev_pkg_root(), ".cursor", "frp.toml"),
    file.path(getwd(), ".cursor", "frp.toml")
  ))
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) {
    existing[[1]]
  } else {
    ""
  }
}

.msdev_frp_require <- function(cfg, key, path) {
  val <- cfg[[key]]
  if (is.null(val) || (is.character(val) && !nzchar(val))) {
    stop("Missing `", key, "` in ", path, call. = FALSE)
  }
  val
}

.msdev_frp_settings <- function() {
  path <- .msdev_frp_local_config_path()
  if (!nzchar(path)) {
    stop(
      "FRP config not found. Create `.cursor/frp.toml` under the MSdev package root.",
      call. = FALSE
    )
  }
  cfg <- .msdev_frp_read_config_file(path)
  list(
    server_addr = as.character(.msdev_frp_require(cfg, "serverAddr", path))[[1]],
    server_port = as.integer(.msdev_frp_require(cfg, "serverPort", path))[[1]],
    remote_port = as.integer(.msdev_frp_require(cfg, "remotePort", path))[[1]],
    token = as.character(.msdev_frp_require(cfg, "token", path))[[1]],
    public_host = as.character(.msdev_frp_require(cfg, "publicHost", path))[[1]],
    public_port = as.integer(.msdev_frp_require(cfg, "publicPort", path))[[1]],
    public_path = as.character(.msdev_frp_require(cfg, "publicPath", path))[[1]],
    version = as.character(.msdev_frp_require(cfg, "frpcVersion", path))[[1]]
  )
}

.msdev_frp_dir <- function() {
  if (.Platform$OS.type == "windows") {
    root <- Sys.getenv("LOCALAPPDATA", unset = path.expand("~"))
  } else {
    root <- file.path(path.expand("~"), ".cache")
  }
  dir <- file.path(root, "MSdev", "frp")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

.msdev_frpc_bin_name <- function() {
  if (.Platform$OS.type == "windows") "frpc.exe" else "frpc"
}

.msdev_ensure_frpc <- function() {
  bin <- .msdev_frpc_bin_name()
  candidates <- c(
    file.path(.msdev_frp_dir(), bin),
    file.path(Sys.getenv("LOCALAPPDATA"), "QuantLR", "frp", bin),
    unname(Sys.which("frpc"))
  )
  existing <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(existing)) {
    return(existing[[1]])
  }
  sys <- Sys.info()[["sysname"]]
  mach <- Sys.info()[["machine"]]
  is_arm <- mach %in% c("arm64", "aarch64")
  if (identical(sys, "Windows")) {
    os <- if (is_arm) "windows_arm64" else "windows_amd64"
    ext <- "zip"
  } else if (identical(sys, "Darwin")) {
    os <- if (is_arm) "darwin_arm64" else "darwin_amd64"
    ext <- "tar.gz"
  } else {
    os <- if (is_arm) "linux_arm64" else "linux_amd64"
    ext <- "tar.gz"
  }
  ver <- .msdev_frp_settings()$version
  url <- sprintf(
    "https://github.com/fatedier/frp/releases/download/v%s/frp_%s_%s.%s",
    ver, ver, os, ext
  )
  message("Downloading frpc ", ver, " ...")
  tmp <- tempfile("msdev-frp-")
  dir.create(tmp)
  archive <- file.path(tmp, paste0("frp.", ext))
  utils::download.file(url, archive, mode = "wb", quiet = TRUE)
  if (identical(ext, "zip")) {
    utils::unzip(archive, exdir = tmp)
  } else {
    utils::untar(archive, exdir = tmp)
  }
  found <- list.files(tmp, pattern = paste0("^", bin, "$"), recursive = TRUE, full.names = TRUE)
  if (!length(found)) {
    stop("Downloaded FRP archive did not contain ", bin, ".", call. = FALSE)
  }
  dest <- file.path(.msdev_frp_dir(), bin)
  ok <- file.copy(found[[1]], dest, overwrite = TRUE)
  if (!isTRUE(ok)) {
    stop("Failed to install frpc to ", dest, ".", call. = FALSE)
  }
  if (.Platform$OS.type != "windows") {
    Sys.chmod(dest, mode = "0755")
  }
  dest
}

.msdev_frp_config_path <- function() {
  file.path(.msdev_frp_dir(), "msdev-qmd-frpc.toml")
}

.msdev_frp_write_client_config <- function(local_host, local_port, settings) {
  if (local_host %in% c("0.0.0.0", "::", "[::]")) {
    local_host <- "127.0.0.1"
  }
  path <- .msdev_frp_config_path()
  lines <- c(
    paste0("serverAddr = \"", settings$server_addr, "\""),
    paste0("serverPort = ", settings$server_port),
    "auth.method = \"token\"",
    paste0("auth.token = \"", settings$token, "\""),
    "loginFailExit = false",
    "",
    "[[proxies]]",
    "name = \"msdev-qmd\"",
    "type = \"tcp\"",
    paste0("localIP = \"", local_host, "\""),
    paste0("localPort = ", as.integer(local_port)),
    paste0("remotePort = ", as.integer(settings$remote_port)),
    ""
  )
  writeLines(lines, path, useBytes = TRUE)
  path
}

.msdev_frp_is_running <- function() {
  proc <- .msdev_frp_env$proc
  isTRUE(tryCatch(proc$is_alive(), error = function(e) FALSE))
}

.msdev_frp_child_env <- function() {
  env <- Sys.getenv()
  drop <- c(
    "HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY",
    "http_proxy", "https_proxy", "all_proxy"
  )
  env[setdiff(names(env), drop)]
}

.msdev_frp_read_log <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    return("")
  }
  raw <- paste(readLines(path, warn = FALSE), collapse = "\n")
  gsub("\033\\[[0-9;]*m", "", raw)
}

.msdev_frp_kill_orphans <- function() {
  marker <- "msdev-qmd-frpc.toml"
  if (.Platform$OS.type == "windows") {
    ps <- tempfile(fileext = ".ps1")
    on.exit(unlink(ps), add = TRUE)
    writeLines(
      paste0(
        "$ErrorActionPreference='SilentlyContinue'; ",
        "Get-CimInstance Win32_Process | Where-Object { ",
        "$_.CommandLine -like '*", marker, "*' ",
        "} | ForEach-Object { Stop-Process -Id $_.ProcessId -Force }"
      ),
      ps
    )
    suppressWarnings(system2("powershell", c("-NoProfile", "-File", ps), stdout = TRUE, stderr = TRUE))
  } else {
    suppressWarnings(system(sprintf("pkill -f %s || true", marker), intern = TRUE))
  }
  Sys.sleep(0.4)
  invisible(NULL)
}

.msdev_frp_stop <- function(quiet = TRUE) {
  proc <- .msdev_frp_env$proc
  if (!is.null(proc)) {
    try(proc$kill(), silent = TRUE)
  }
  .msdev_frp_kill_orphans()
  .msdev_frp_env$proc <- NULL
  .msdev_frp_env$log <- NULL
  .msdev_frp_env$config <- NULL
  .msdev_frp_env$public_url <- NULL
  invisible(NULL)
}

.msdev_frp_start <- function(local_port, local_host = "127.0.0.1") {
  if (!requireNamespace("processx", quietly = TRUE)) {
    stop("Package 'processx' is required to project the presentation via FRP")
  }
  settings <- .msdev_frp_settings()
  if (!nzchar(settings$token)) {
    stop("FRP token missing in `.cursor/frp.toml`.", call. = FALSE)
  }
  if (local_host %in% c("0.0.0.0", "::", "[::]")) {
    local_host <- "127.0.0.1"
  }
  if (.msdev_frp_is_running()) {
    return(.msdev_frp_env$public_url)
  }
  .msdev_frp_kill_orphans()
  bin <- .msdev_ensure_frpc()
  cfg <- .msdev_frp_write_client_config(local_host, local_port, settings)
  logf <- tempfile("msdev-frpc-", fileext = ".log")
  proc <- processx::process$new(
    command = bin,
    args = c("-c", cfg),
    stdout = logf,
    stderr = "2>&1",
    env = .msdev_frp_child_env()
  )
  .msdev_frp_env$proc <- proc
  .msdev_frp_env$log <- logf
  .msdev_frp_env$config <- cfg
  .msdev_frp_env$public_url <- .msdev_frp_public_url(settings)

  deadline <- Sys.time() + 15
  ready <- FALSE
  failed <- FALSE
  last <- ""
  while (Sys.time() < deadline) {
    if (!isTRUE(proc$is_alive())) {
      last <- .msdev_frp_read_log(logf)
      break
    }
    last <- .msdev_frp_read_log(logf)
    if (grepl("start error|already exists|port already used|bind error", last, ignore.case = TRUE)) {
      failed <- TRUE
      break
    }
    if (grepl("start proxy success", last, ignore.case = TRUE)) {
      ready <- TRUE
      break
    }
    Sys.sleep(0.2)
  }
  if (!ready) {
    .msdev_frp_stop()
    hint <- ""
    if (isTRUE(failed) || grepl("already exists|port already used", last, ignore.case = TRUE)) {
      hint <- paste0(
        "\nRemote port ", settings$remote_port,
        " is already held by another FRP client."
      )
    }
    stop(
      "FRP client failed to publish at ",
      .msdev_frp_public_url(settings),
      hint,
      if (nzchar(last)) paste0(":\n", last) else ".",
      call. = FALSE
    )
  }
  .msdev_frp_env$public_url
}

.msdev_frp_try_start <- function(local_port, local_host = "127.0.0.1") {
  tryCatch(
    .msdev_frp_start(local_port, local_host),
    error = function(e) {
      warning("FRP not started: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
}
