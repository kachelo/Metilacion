#' MethylationExperimentLimma-Class
#'
#' @description  This is a concrete class extends MethylationExperiment to impose a normal
#'    distribution using limma package to fit it. In this context, as it
#'    inherits from a SummarizedExperiment the fitted models will be included
#'    a "model" column, in the features (Rows), where each model is stored as a
#'    list.
#' @docType class
#' @exportClass MethylationExperimentLimma
setClass(
    Class = "MethylationExperimentLimma",
    contains = "MethylationExperiment"
)




