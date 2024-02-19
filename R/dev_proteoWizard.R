#' msconvert_raw2mzML
#'
#' @param raw.files file path
#' @param mzML.files file path
#' @param BPPARAM biocparall
#'
#' @return null
#' @importFrom BiocParallel bplapply SnowParam

msconvert_raw2mzML <- function(raw.files ,
                                mzML.files,
                                BPPARAM = SnowParam(workers = parallel::detectCores()-1)){


  raw.files <- gsub(pattern = "\\",x = raw.files,replacement = "/",fixed = T)%>%
    na.omit()
  mzML.files <- gsub(pattern = "\\",x = mzML.files,replacement = "/",fixed = T)%>%
    na.omit()

  ###check msconvert
  {
    msconvert_return <- try(system("msconvert",
                                   intern = T))
    if(!any(grepl(pattern = "Usage: msconvert", x = msconvert_return))){
      stop("Command msconvert error, please check environment variables")
    }

    }
  ###check file and directory
  {
    if(!any(file.exists(raw.files))){
      stop(paste0("File not found : ",sum(!file.exists(raw.files)),"/", length(raw.files)))
    }
    if(length(raw.files) != length(mzML.files)){
      stop("raw files and mzml files not match")
    }
    sapply(unique(dirname(mzML.files)),dir.create,recursive =T,showWarnings =F)

  }

  ###msconvert
  {

    shell.commomd <- paste0("msconvert  --filter \"peakPicking true 1-\" --mzML ",
                            raw.files,
                            " -o ",
                            dirname(mzML.files),
                            " --outfile ",
                            mzML.files)
    #system(shell.commomd,intern = T)


    bplapply(shell.commomd,
             FUN = function(x){ system(x,intern = T)},
             BPPARAM = BPPARAM)
    return(0)

  }


}



msconvert_wiff2mzML <- function(wiff.files ,
                                mzML.files,
                                BPPARAM = SnowParam(workers = parallel::detectCores()-1)){


  wiff.files <- gsub(pattern = "\\",x = wiff.files,replacement = "/",fixed = T)%>%
    na.omit()
  mzML.files <- gsub(pattern = "\\",x = mzML.files,replacement = "/",fixed = T)%>%
    na.omit()

  ###check msconvert
  {
    msconvert_return <- try(system("msconvert",
                                   intern = T))
    if(!any(grepl(pattern = "Usage: msconvert", x = msconvert_return))){
      stop("Command msconvert error, please check environment variables")
    }

    }
  ###check file and directory
  {
    if(!any(file.exists(wiff.files))){
      stop(paste0("File not found : ",sum(!file.exists(wiff.files)),"/", length(wiff.files)))
    }
    if(length(wiff.files) != length(mzML.files)){
      stop("wiff files and mzml files not match")
    }
    sapply(unique(dirname(mzML.files)),dir.create,recursive =T,showWarnings =F)

  }

  ###msconvert
  {

    shell.commomd <- paste0("msconvert  --filter \"peakPicking true 1-\" --mzML ",
                            wiff.files,
                            " -o ",
                            dirname(mzML.files),
                            " --outfile ",
                            mzML.files)
    #system(shell.commomd,intern = T)


    bplapply(shell.commomd,
             FUN = function(x){ system(x,intern = T)},
             BPPARAM = BPPARAM)
    return(0)

  }


}


msConvertDir<- function(raw.path,format.to = "mzML"){

  dir.create(paste0(raw.path,"/mzML"),recursive = T)
  raw.files <- data.frame(raw.file = dir(path = raw.path,full.names = T))%>%
    dplyr::mutate(format = case_when(grepl(pattern = ".raw$",x = raw.file)~".raw",
                                     grepl(pattern = ".wiff$",x = raw.file)~".wiff",
                                     T~"unknow"
                                     ))%>%
    dplyr::filter(format %in% c(".raw",".wiff"))%>%
    dplyr::group_by(raw.file)%>%
    dplyr::mutate(mzML = paste0(dirname(raw.file),
                                        "/mzML/",
                                        gsub(x = basename(raw.file) ,
                                                 replacement = ".mzML",
                                                 pattern = paste0(format,"$"))),
                  file.exist = file.exists(mzML))%>%
    dplyr::filter(!file.exist)
  msConvert2mzML(raw.files$raw.file,raw.files$mzML)
  return(raw.files$mzML)


}


msConvert2mzML <- function(raw.files ,
                           mzML.files,
                           BPPARAM = SnowParam(workers = parallel::detectCores()-1)){


 raw.files <- gsub(pattern = "\\",x = raw.files,replacement = "/",fixed = T)%>%
    na.omit()
  mzML.files <- gsub(pattern = "\\",x = mzML.files,replacement = "/",fixed = T)%>%
    na.omit()

     ###check msconvert
  {
    msconvert_return <- try(system("msconvert",
                                   intern = T))
    if(!any(grepl(pattern = "Usage: msconvert", x = msconvert_return))){
      stop("Command msconvert error, please check environment variables")
    }

    }
  ###check file and directory
  {
    if(!any(file.exists(raw.files))){
      stop(paste0("File not found : ",sum(!file.exists(raw.files)),"/", length(raw.files)))
    }
    if(length(raw.files) != length(mzML.files)){
      stop("raw files and mzml files not match")
    }
    sapply(unique(dirname(mzML.files)),dir.create,recursive =T,showWarnings =F)

  }

  ###msconvert
  {

    shell.commomd <- paste0("msconvert  --filter \"peakPicking true 1-\" --mzML ",
                            raw.files,
                            " -o ",
                            dirname(mzML.files),
                            " --outfile ",
                            mzML.files)
    #system(shell.commomd,intern = T)


    BiocParallel::bplapply(shell.commomd,
             FUN = function(x){ system(x,intern = T)},
             BPPARAM = BPPARAM)
    return(0)

  }


}
