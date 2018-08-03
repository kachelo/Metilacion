#' High level constructor for MethylationExperiment class
#'
#' @description This is a higher level function that allows to construct
#'    a MethylationExperimentLimma or MethylationExperimentBeta object with
#'    the colData and rawData.
#' @docType methods
MethylationExperiment <- function(ClassName){
    stopifnot(ClassName %in%
        c("MethylationExperimentLimma", "MethylationExperimentBeta"))
    object <- switch(ClassName,
        "MethylationExperimentLimma" = .MethylationExperimentLimma(),
        "MethylationExperimentBeta" = .MethylationExperimentBeta()
    )
}
object <- MethylationExperiment(ClassName = "MethylationExperimentLimma")
