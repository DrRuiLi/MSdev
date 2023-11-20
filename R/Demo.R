load_demo <- function(demo = c("MSdev","XCMSnExp","xcms","SummarizedExperiment","data.se")){

  demo <-match.arg(demo)
  file_path <- switch(demo,
         "MSdev" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/MSdev_2023_11_17.Rdata",
         "XCMSnExp" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/XCMSnExp_2023_11_17.rda",
         "xcms" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/XCMSnExp_2023_11_17.rda",
          "SummarizedExperiment" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/SummarizedExperiment_2023_11_03.rda",
         "data.se" = "C:\\Users\\91879\\OneDrive\\Code\\R\\Projecct\\2022.1.8_MS.demo\\Demo/SummarizedExperiment_2023_11_03.rda"
         )
  demo <- load_as_var(file_path)

  return(demo)
}
