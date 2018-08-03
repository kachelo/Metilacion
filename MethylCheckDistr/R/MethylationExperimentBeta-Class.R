#' MethylationExperimentBeta-Class
#'
#' @description  This is a concrete class extends MethylationExperiment to
#'    impose a beta distribution using betareg package to fit it. In this
#'    context, as it inherits from a SummarizedExperiment the fitted models
#'    will be included a "model" column, in the features (rowData), where each
#'    model is stored as a list.
#' @docType class
#' @exportClass MethylationExperimentBeta
.MethylationExperimentBeta <- setClass(
    Class = "MethylationExperimentLimma",
    contains = "MethylationExperiment"
)




