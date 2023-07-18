msconvert_wiff2mzML <- function(wiff.files ,
                                mzML.files,
                                BPPARAM = SnowParam(workers = 10)){


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
