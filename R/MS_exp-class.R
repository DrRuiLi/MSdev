setClass(
  "MS_Exp",
  slots = list(
    "General" = "tbl",
    "Pre_process" = "list",
    "Moblie_phase_A" = "list",
    "Moblie_phase_B" = "list",
    "Chroma_column" = "list",
    "Chroma_gradient" = "list",
    "Mass_Spectrum" = "list"
  )

)

MS_Exp <- function() {
  new("MS_Exp")
}


setMethod("initialize" , "MS_Exp",
          function(.Object) {
            General <- tibble(
              "MSE_id" =  "MSE00001",
              "Name" =  "Metabolomics",
              "Creat_time" = paste0(Sys.time()),
              "Link" = "",
              "PubMed_ID" = "",
              "Note" = ""
            )
            Pre_process <- tibble("Step" = 1:5,
                                  "Process" = "")%>%list()

            Moblie_phase_A <- tibble(
              "Compound" = c("H2O", "Formic acid"),
              "Type" = c("solvent", "solute"),
              "Concentration" = c("100%", "0.001mM")
            )%>%list
            Moblie_phase_B <- tibble(
              "Compound" = c("ACN", "Formic acid"),
              "Type" = c("solvent", "solute"),
              "Concentration" = c("100%", "0.001mM")
            )%>%list
            Chroma_column <- tibble(
              "Column_name" = "Kinetex C18",
              "Manufacturer" = "Phenomenex",
              "Paricle_size" = "1.7 μm",
              "Length" = "50 mm",
              "Diameter" = "2.1 mm")%>%list
            Chroma_gradient <- tibble(
              time = c(0, 60, 100, 500),
              Contentration_B = c(0, 50, 100, 100))%>%list
            Mass_Spectrum <- tibble(
              "Instrument" = "SCIEX TripleTOF 6600",
              "MS_type" = "Q-TOF",
              "Data_aquisition" = "DDA TOP10",
            )%>%list

            .Object@General <- General
            .Object@Pre_process <-Pre_process
            .Object@Moblie_phase_A <- Moblie_phase_A
            .Object@Moblie_phase_B <- Moblie_phase_B
            .Object@Chroma_column <- Chroma_column
            .Object@Chroma_gradient <- Chroma_gradient
            .Object@Mass_Spectrum <- Mass_Spectrum
            .Object


          })

setMethod(
  "show" ,
  "MS_Exp",
  definition = function(object) {
    general_data <- object@General

    message(nrow(object@General)," Record in MS_Exp database")
    show(general_data)
  }
)


setMethod("[","MS_Exp",
          definition = function(x , i){
          new.object <- MS_Exp()
          new.object@General <- x@General[i,]
          new.object@Pre_process <- x@Pre_process[i]
          new.object@Moblie_phase_A <- x@Moblie_phase_A[i]
          new.object@Moblie_phase_B <- x@Moblie_phase_B[i]
          new.object@Chroma_column <- x@Chroma_column[i]
          new.object@Chroma_gradient <- x@Chroma_gradient[i]
          new.object@Mass_Spectrum <- x@Mass_Spectrum[i]
          new.object

          })




object <- MS_Exp()
object



