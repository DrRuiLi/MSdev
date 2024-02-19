#' @title MS_Exp-class
#'
#' @slot General tbl.
#' @slot Pre_process list.
#' @slot Moblie_phase_A list.
#' @slot Moblie_phase_B list.
#' @slot Chroma_column list.
#' @slot Chroma_gradient list.
#' @slot Mass_Spectrum list.
#'
#' @export
#'

setClass(
  "MS_Exp",
  slots = list(
    "General" = "tbl",
    "Pre_process" = "list",
    "Moblie_phase_A" = "list",
    "Moblie_phase_B" = "list",
    "Chroma_column" = "list",
    "Chroma_gradient" = "list",
    "Mass_Spectrum" = "list",
    "Internal_Standard" = "list"
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
              "Source" = "YHY Lab",
              "Link" = "",
              "PubMed_ID" = "",
              "Note" = ""
            )
            Pre_process <- tibble("Step" = 1:5,
                                  "Process" = "")%>%list()

            Moblie_phase_A <- tibble(
              "Compound" = c( "Formic acid","H2O"),
              "Type" = c( "solute","solvent"),
              "Concentration" = c("0.1 mM","100 %" )
            )%>%list
            Moblie_phase_B <- tibble(
              "Compound" = c( "Formic acid","ACN"),
              "Type" = c("solute","solvent"),
              "Concentration" = c( "0.1 mM","100 %")
            )%>%list
            Chroma_column <- tibble(
              "Column_name" = "Kinetex C18",
              "Manufacturer" = "Phenomenex",
              "Paricle_size" = "1.7 μm",
              "Length" = "50 mm",
              "Diameter" = "2.1 mm",
              "Item_No" = "",
            "Link" = "")%>%list
            Chroma_gradient <- tibble(
              time = c(0, 2, 4, 10),
              Contentration_B = c(0, 50, 100, 100),
              Flow_rate = 2.5)%>%list
            Mass_Spectrum <- tibble(
              "Instrument" = "SCIEX TripleTOF 6600",
              "MS_type" = "Q-TOF",
              "Data_aquisition" = "DDA TOP10",
              "Ion_mode" = "both",
            )%>%list
            Internal_Standard <- data.frame(
              "Compound_name" = "",
              "Exact_mass" = NA,
              "Retention time" = NA

            )%>%list()

            .Object@General <- General
            .Object@Pre_process <-Pre_process
            .Object@Moblie_phase_A <- Moblie_phase_A
            .Object@Moblie_phase_B <- Moblie_phase_B
            .Object@Chroma_column <- Chroma_column
            .Object@Chroma_gradient <- Chroma_gradient
            .Object@Mass_Spectrum <- Mass_Spectrum
            .Object@Internal_Standard <- Internal_Standard
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
          new.object@Internal_Standard <- x@Internal_Standard[i]
          new.object

          })
setMethod("c" , "MS_Exp" ,
          definition = function(x,y){
            new.object <- MS_Exp()
            new.object@General <- rbind(x@General,y@General)
            new.object@Pre_process <- c(x@Pre_process,y@Pre_process)
            new.object@Moblie_phase_A <- c(x@Moblie_phase_A,y@Moblie_phase_A)
            new.object@Moblie_phase_B <- c(x@Moblie_phase_B,y@Moblie_phase_B)
            new.object@Chroma_column <- c(x@Chroma_column,y@Chroma_column)
            new.object@Chroma_gradient <- c(x@Chroma_gradient,y@Chroma_gradient)
            new.object@Mass_Spectrum <-c(x@Mass_Spectrum,y@Mass_Spectrum)
            new.object@Internal_Standard <-c(x@Internal_Standard,y@Internal_Standard)
            new.object



          })

setMethod("length" ,
          "MS_Exp" ,
          definition = function(x){
            nrow(x@General)})





