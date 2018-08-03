#' MethylationExperiment-Class
#'
#' @description  This is a VIRTUAL class that extends SummarizedExperiment to
#'    to model a Methylation Experiment where the idea is to contain all the
#'    data in a single object instead of having multiple copies across the
#'    analysis pipeline, e. g., raw, normalized, filtered objects, etc. This
#'    class is then extended into concrete clases according to the distribution
#'    imposed in the analysis, i. e., normal (limma style), beta (betareg style)
#'    or any other distribution based implemented class.
#' @docType class
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @exportClass MethylationExperiment
setClass(
    Class = "MethylationExperiment",
    contains = c("SummarizedExperiment", "VIRTUAL")
)
