validate_ribo <- function(object) {
    errors <- character()
    
    num_experiments <- length(object@experiments)
    
    if (num_experiments <= 0) {
        msg <- paste("File appears to have no experiments.")
        errors <- c(errors, msg)
    }
    
    length.min <- object@length.min
    length.max <- object@length.max 
    
    if (length.min > length.max) {
        msg <- paste("The minimum read length should be less than or equal ",
                     "to the maximum read length.", sep = "")
        errors <- c(errors, msg)
    }
    
    if (length(errors) == 0) TRUE else errors
}

#' Class "ribo"
#' 
#' An S4 class to be used with the ribor package
#' 
#' @param object A 'ribo' object
#' @slot handle A handle to the ribo file of interest
#' @slot experiments A character vector of experiment names in the file 
#' @slot format.version The format version of the ribo file 
#' @slot reference The reference transcriptome used in the ribo file 
#' @slot length.min The minimum read length of the data in the file 
#' @slot length.max The maximum read length of the data in the file 
#' @slot left.span Left span of the junction regions
#' @slot right.span Right span of the junction regions
#' @slot length.offset Length offset of all transcripts 
#' @slot has.metadata Value denoting whether the root ribo file has metadata 
#' @slot experiment.info Data table of information on the experiments 
#' @slot transcript.info Hash of the lengths and offsets of each transcript
#' @slot transcript.alias Hash that goes from alias to original transcript name
#' @slot transcript.original Hash that goes from original to alias transcript 
#' @return A 'ribo' object
#' @seealso \code{\link{create_ribo}} to create a ribo file
ribo <- setClass(
    "ribo",
    
    representation(handle = "H5IdComponent",
                   experiments = "character",
                   format.version = "integer",
                   reference = "character",
                   length.min = "integer",
                   length.max = "integer",
                   left.span = "integer",
                   right.span = "integer",
                   length.offset = "numeric",
                   has.metadata = "logical",
                   experiment.info = "data.table",
                   transcript.info = "hash",
                   transcript.alias = "hash",
                   transcript.original = "hash"),
    
    validity = validate_ribo
)