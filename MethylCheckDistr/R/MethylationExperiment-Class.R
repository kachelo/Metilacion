#' MethylationExperiment-Class
#'
#' @description  This is a VIRTUAL class that extends SummarizedExperiment to
#'    to model a Methylation Experiment where the idea is to contain all the
#'    data in a single object instead of having multiple copies across the
#'    analysis pipeline, e. g., raw, normalized, filtered objects, etc. This
#'    class is then extended into concrete clases according to the distribution
#'    imposed in the analysis, i. e., normal (limma style), beta (betareg style)
#'    or any other distribution based implemented class. In this context,
#'    whatever thet fitted model is, it will be included in the features
#'    (rowData) as "model" column, where the fit itself is represented as a list.
#' @docType class
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods new
#' @exportClass MethylationExperiment
setClass(
    Class = "MethylationExperiment",
    contains = c("SummarizedExperiment", "VIRTUAL")
)
