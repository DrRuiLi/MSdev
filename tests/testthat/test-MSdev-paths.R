new_msdev <- function(project_info,
                     sample_info = data.frame(),
                     spectra = list(),
                     xcms = list()) {
  methods::new(
    "MSdev",
    projectInfo = project_info,
    processingInfo = list(),
    sampleInfo = sample_info,
    experimentInfo = MS_Exp(),
    xcmsData = xcms,
    spectra = spectra,
    annotation = list(),
    advancedAna = list()
  )
}

slash <- function(x) MSdev:::.msdev_slash(x)

test_that("prefix rewrite covers project paths and leaves outside paths", {
  maps <- list(list(old = "C:/old/proj", new = "C:/new/proj"))
  expect_equal(
    slash(MSdev:::.msdev_rewrite_paths("C:\\old\\proj\\raw.data\\a.raw", maps)),
    "C:/new/proj/raw.data/a.raw"
  )
  if (.Platform$OS.type == "windows") {
    expect_equal(
      slash(MSdev:::.msdev_rewrite_paths("c:/old/proj/raw.data/msData/a.mzML", maps)),
      "C:/new/proj/raw.data/msData/a.mzML"
    )
  }
  outside <- MSdev:::.msdev_rewrite_paths("E:/db/CompoundDB.sqlite", maps)
  expect_equal(slash(outside), slash(MSdev:::.msdev_path_prefix("E:/db/CompoundDB.sqlite")))
})

test_that("MSdev_path derives relatives when paths is missing", {
  info <- list(
    projectDir = "C:/old/proj",
    rawDataDir = "C:/old/proj/raw.data",
    msDataDir = "C:/old/proj/raw.data/msData",
    MSdevFile = "C:/old/proj/MSdev.Rdata"
  )
  obj <- new_msdev(info)
  expect_equal(slash(MSdev_path(obj, "rawData")), slash(MSdev:::.msdev_path_prefix(info$rawDataDir)))
  expect_equal(slash(MSdev_path(obj, "ms1")), slash(file.path(MSdev:::.msdev_path_prefix(info$projectDir), "MS1_Spectra.rds")))
  expect_error(MSdev_path(obj, "nope"), "Unknown MSdev path")
})

test_that("MSdev_relocate rewrites projectInfo, sampleInfo, and on-disk Spectra", {
  tmp <- tempfile("msdev-reloc-")
  old_root <- file.path(tmp, "old")
  new_root <- file.path(tmp, "new")
  dir.create(file.path(old_root, "raw.data", "msData"), recursive = TRUE)
  dir.create(new_root, recursive = TRUE)
  old_root <- MSdev:::.msdev_path_prefix(old_root)
  new_root <- MSdev:::.msdev_path_prefix(new_root)

  old_raw <- file.path(old_root, "raw.data", "a.raw")
  old_mzml <- file.path(old_root, "raw.data", "msData", "a.mzML")
  spd <- S4Vectors::DataFrame(msLevel = 1L, rtime = 1)
  sp <- Spectra::Spectra(spd, backend = Spectra::MsBackendMemory())
  sp$dataOrigin <- old_mzml
  sp$dataStorage <- old_mzml
  sp_rewritten <- MSdev:::.msdev_rewrite_spectra_paths(
    sp,
    list(list(old = old_root, new = new_root))
  )
  expect_equal(
    slash(Spectra::dataOrigin(sp_rewritten)),
    slash(file.path(new_root, "raw.data", "msData", "a.mzML"))
  )

  qs_old <- file.path(old_root, "MS1_Spectra.rds")
  od <- onDiskData(sp, path = qs_old)
  file.copy(od@path, file.path(new_root, "MS1_Spectra.rds"), overwrite = TRUE)

  obj <- new_msdev(
    project_info = list(
      projectDir = old_root,
      rawDataDir = file.path(old_root, "raw.data"),
      msDataDir = file.path(old_root, "raw.data", "msData"),
      MSdevFile = file.path(old_root, "MSdev.Rdata"),
      CompoundDB_path = "E:/db/CompoundDB.sqlite"
    ),
    sample_info = data.frame(
      raw.files = old_raw,
      msData.files = old_mzml,
      stringsAsFactors = FALSE
    ),
    spectra = list(MS1_Spectra = od)
  )

  obj2 <- MSdev_relocate(obj, projectDir = new_root)
  expect_equal(obj2@projectInfo$projectDir, new_root)
  expect_equal(slash(obj2@projectInfo$rawDataDir), slash(file.path(new_root, "raw.data")))
  expect_equal(slash(obj2@projectInfo$msDataDir), slash(file.path(new_root, "raw.data", "msData")))
  expect_equal(slash(obj2@projectInfo$MSdevFile), slash(file.path(new_root, "MSdev.Rdata")))
  expect_identical(obj2@projectInfo$CompoundDB_path, "E:/db/CompoundDB.sqlite")
  expect_equal(obj2@projectInfo$paths[["rawData"]], "raw.data")
  expect_equal(obj2@projectInfo$paths[["msData"]], "raw.data/msData")
  expect_equal(obj2@projectInfo$paths[["ms1"]], "MS1_Spectra.rds")
  expect_equal(slash(obj2@sampleInfo$raw.files), slash(file.path(new_root, "raw.data", "a.raw")))
  expect_equal(slash(obj2@sampleInfo$msData.files), slash(file.path(new_root, "raw.data", "msData", "a.mzML")))
  expect_equal(slash(obj2@spectra$MS1_Spectra@path), slash(file.path(new_root, "MS1_Spectra.rds")))

  sp2 <- onDiskData_retrieve(obj2@spectra$MS1_Spectra)
  expect_equal(
    slash(Spectra::dataOrigin(sp2)),
    slash(file.path(new_root, "raw.data", "msData", "a.mzML"))
  )
  expect_equal(
    slash(Spectra::dataStorage(sp2)),
    slash(file.path(new_root, "raw.data", "msData", "a.mzML"))
  )
  expect_equal(slash(MSdev_path(obj2, "ms1")), slash(obj2@spectra$MS1_Spectra@path))
})

test_that("MSdev_relocate can point raw data outside the project", {
  tmp <- tempfile("msdev-raw-")
  old_proj <- MSdev:::.msdev_path_prefix(file.path(tmp, "oldproj"))
  new_proj <- MSdev:::.msdev_path_prefix(file.path(tmp, "newproj"))
  old_raw <- MSdev:::.msdev_path_prefix(file.path(tmp, "raw1"))
  new_raw <- MSdev:::.msdev_path_prefix(file.path(tmp, "raw2"))
  dir.create(new_proj, recursive = TRUE, showWarnings = FALSE)
  dir.create(new_raw, recursive = TRUE, showWarnings = FALSE)

  obj <- new_msdev(
    project_info = list(
      projectDir = old_proj,
      rawDataDir = old_raw,
      msDataDir = file.path(old_raw, "msData"),
      MSdevFile = file.path(old_proj, "MSdev.Rdata")
    ),
    sample_info = data.frame(
      raw.files = file.path(old_raw, "a.raw"),
      msData.files = file.path(old_raw, "msData", "a.mzML"),
      stringsAsFactors = FALSE
    )
  )
  obj2 <- MSdev_relocate(obj, projectDir = new_proj, rawDataDir = new_raw)
  expect_equal(obj2@projectInfo$projectDir, new_proj)
  expect_equal(slash(obj2@projectInfo$rawDataDir), new_raw)
  expect_equal(slash(obj2@projectInfo$msDataDir), slash(file.path(new_raw, "msData")))
  expect_equal(slash(obj2@sampleInfo$msData.files), slash(file.path(new_raw, "msData", "a.mzML")))
  expect_true(MSdev:::.msdev_is_absolute(obj2@projectInfo$paths[["rawData"]]))
})

test_that("MSdev_relocate updates MsExperiment file paths", {
  tmp <- tempfile("msdev-mse-")
  old_root <- MSdev:::.msdev_path_prefix(file.path(tmp, "old"))
  new_root <- MSdev:::.msdev_path_prefix(file.path(tmp, "new"))
  dir.create(new_root, recursive = TRUE, showWarnings = FALSE)
  old_mzml <- file.path(old_root, "raw.data", "msData", "a.mzML")
  old_raw <- file.path(old_root, "raw.data", "a.raw")

  spd <- S4Vectors::DataFrame(msLevel = 1L, rtime = 1)
  sp <- Spectra::Spectra(spd, backend = Spectra::MsBackendMemory())
  sp$dataOrigin <- old_mzml
  sp$dataStorage <- old_mzml
  mse <- MsExperiment::MsExperiment()
  mse <- ProtGenerics::`spectra<-`(mse, value = sp)
  mse <- MsExperiment::`sampleData<-`(
    mse,
    value = S4Vectors::DataFrame(
      raw.files = old_raw,
      msData.files = old_mzml,
      dataOrigin = old_mzml
    )
  )

  obj <- new_msdev(
    project_info = list(
      projectDir = old_root,
      rawDataDir = file.path(old_root, "raw.data"),
      msDataDir = file.path(old_root, "raw.data", "msData"),
      MSdevFile = file.path(old_root, "MSdev.Rdata")
    ),
    xcms = list(PositiveMS1 = mse)
  )
  obj2 <- MSdev_relocate(obj, projectDir = new_root)
  mse2 <- obj2@xcmsData$PositiveMS1
  expect_equal(
    slash(Spectra::dataOrigin(ProtGenerics::spectra(mse2))),
    slash(file.path(new_root, "raw.data", "msData", "a.mzML"))
  )
  sd <- MsExperiment::sampleData(mse2)
  expect_equal(
    slash(as.character(sd$msData.files)),
    slash(file.path(new_root, "raw.data", "msData", "a.mzML"))
  )
  expect_equal(
    slash(as.character(sd$raw.files)),
    slash(file.path(new_root, "raw.data", "a.raw"))
  )
})

test_that("MSdev_load relocates only when the saved projectDir is missing", {
  tmp <- tempfile("msdev-load-")
  kept <- file.path(tmp, "kept")
  other <- file.path(tmp, "other")
  moved <- file.path(tmp, "moved")
  dir.create(kept, recursive = TRUE)
  dir.create(other, recursive = TRUE)
  dir.create(moved, recursive = TRUE)
  kept <- MSdev:::.msdev_path_prefix(kept)
  other <- MSdev:::.msdev_path_prefix(other)
  moved <- MSdev:::.msdev_path_prefix(moved)

  stay <- new_msdev(list(
    projectDir = kept,
    rawDataDir = file.path(kept, "raw.data"),
    msDataDir = file.path(kept, "raw.data", "msData"),
    MSdevFile = file.path(kept, "MSdev.Rdata")
  ), sample_info = data.frame(
    raw.files = file.path(kept, "raw.data", "a.raw"),
    msData.files = file.path(kept, "raw.data", "msData", "a.mzML"),
    stringsAsFactors = FALSE
  ))
  stay_file <- file.path(other, "MSdev.qs")
  qs2::qs_save(stay, file = stay_file)
  loaded_stay <- MSdev_load(stay_file)
  expect_equal(loaded_stay@projectInfo$projectDir, kept)

  missing_root <- "C:/msdev_missing_project/old"
  gone <- new_msdev(list(
    projectDir = missing_root,
    rawDataDir = file.path(missing_root, "raw.data"),
    msDataDir = file.path(missing_root, "raw.data", "msData"),
    MSdevFile = file.path(missing_root, "MSdev.Rdata")
  ), sample_info = data.frame(
    raw.files = file.path(missing_root, "raw.data", "a.raw"),
    msData.files = file.path(missing_root, "raw.data", "msData", "a.mzML"),
    stringsAsFactors = FALSE
  ))
  moved_file <- file.path(moved, "MSdev.qs")
  qs2::qs_save(gone, file = moved_file)
  loaded <- MSdev_load(moved_file)
  expect_equal(loaded@projectInfo$projectDir, moved)
  expect_equal(slash(loaded@sampleInfo$raw.files), slash(file.path(moved, "raw.data", "a.raw")))
  expect_equal(loaded@projectInfo$paths[["rawData"]], "raw.data")
})

test_that("CompoundDb resolution uses argument, stored path, option, then env", {
  obj <- new_msdev(list(
    projectDir = "C:/old/proj",
    CompoundDB_path = "C:/missing/CompoundDB.sqlite"
  ))
  old_opt <- getOption("MSdev.compoundDb")
  old_env <- Sys.getenv("MSDEV_COMPOUNDDB", unset = NA_character_)
  on.exit({
    options(MSdev.compoundDb = old_opt)
    if (is.na(old_env)) {
      Sys.unsetenv("MSDEV_COMPOUNDDB")
    } else {
      Sys.setenv(MSDEV_COMPOUNDDB = old_env)
    }
  }, add = TRUE)

  options(MSdev.compoundDb = NULL)
  Sys.unsetenv("MSDEV_COMPOUNDDB")
  expect_error(MSdev:::.msdev_require_cpdb_path(obj), "CompoundDb path not found")

  opt_file <- tempfile(fileext = ".sqlite")
  env_file <- tempfile(fileext = ".sqlite")
  arg_file <- tempfile(fileext = ".sqlite")
  file.create(opt_file, env_file, arg_file)
  opt_file <- MSdev:::.msdev_path_prefix(opt_file)
  env_file <- MSdev:::.msdev_path_prefix(env_file)
  arg_file <- MSdev:::.msdev_path_prefix(arg_file)

  options(MSdev.compoundDb = opt_file)
  Sys.setenv(MSDEV_COMPOUNDDB = env_file)
  expect_equal(MSdev:::.msdev_resolve_cpdb_path(obj, NULL), opt_file)

  options(MSdev.compoundDb = NULL)
  expect_equal(MSdev:::.msdev_resolve_cpdb_path(obj, NULL), env_file)

  expect_equal(MSdev:::.msdev_resolve_cpdb_path(obj, arg_file), arg_file)

  stored <- tempfile(fileext = ".sqlite")
  file.create(stored)
  stored <- MSdev:::.msdev_path_prefix(stored)
  obj@projectInfo$CompoundDB_path <- stored
  options(MSdev.compoundDb = opt_file)
  expect_equal(MSdev:::.msdev_resolve_cpdb_path(obj, NULL), stored)
})

test_that("MSdev requires rawDataDir", {
  expect_error(MSdev(), "rawDataDir")
})
