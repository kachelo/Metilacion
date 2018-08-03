#' High level constructor for MethylationExperiment class
#'
#' @description This is a higher level function that allows to construct
#'    a MethylationExperimentLimma or MethylationExperimentBeta object with
#'    the colData and rawData.
#' @importFrom minfi getAnnotation
#' @docType methods
MethylationExperiment <- function(ClassName, ...){
    #Check for right ClassName
    stopifnot(ClassName %in%
        c("MethylationExperimentLimma", "MethylationExperimentBeta"))

    #Create an empty object
    object <- switch(ClassName,
        "MethylationExperimentLimma" = .MethylationExperimentLimma(),
        "MethylationExperimentBeta" = .MethylationExperimentBeta()
    )

    #Parse the args if present

    object <- MethylationExperiment(ClassName = "MethylationExperimentLimma")
    targets <- read.csv(
        paste0(
            system.file(package = "MethylCheckDistr"),
            "/extdata/SampleSheet.csv"
        )
    )
    # ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    # object
    # colData = DataFrame(targets)
    nrows <- 200
    ncols <- 5
    counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
    object <- SummarizedExperiment(
        colData = DataFrame(targets),
        assays = SimpleList(counts = counts)
    )
    object

}

