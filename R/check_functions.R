check_ribo <- function(ribo.object, stop = TRUE) {
    # Helper method that checks the internal contents and class of a parameter
    # check_ribo takes in an object and checks for proper contents and class
    # Args:
    # ribo.object An S4 object of class "ribo"
    #
    # Return:
    # None

    is.ribo <- is.ribo(ribo.object)
    if (stop && !is.ribo) {
        stop("Param ribo.object should be of class ribo.", call. = FALSE)
    }
    return(is.ribo)
}

#' Checks for a valid 'ribo' object
#'
#' The function \code{\link{is.ribo}} checks whether or not a given ribo object
#' is valid.
#'
#' @param ribo.object Object in question
#' @importFrom methods is
#' @examples
#' file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")
#' sample <- create_ribo(file.path, rename = rename_default)
#' 
#' is.ribo(sample)
#' @return A boolean indicating whether or not a given object is a valid 'ribo' object
#' @export
is.ribo <- function(ribo.object) {
    return (all(is(ribo.object) == "ribo"))
}

check_alias <- function(ribo.object,
                        alias) {
    check_ribo(ribo.object)
    has.alias <- !is.empty(ribo.object@transcript.alias) &&
                 !is.empty(ribo.object@transcript.original)
    if (alias && !has.alias) {
        stop("Transcripts do not have any aliases.", call. = FALSE)
    }
}

check_lengths <- function(ribo.object, range.lower, range.upper) {
    # Helper method that checks for correct lengths
    #
    # check_lengths directly reads the .ribo file for its lowest and highest read
    # length and compares it to the corresponding parameters
    #
    # Args:
    # ribo.object S3 object of class "ribo"
    # range.lower lowest read length
    # range.upper highest read length
    #
    # Return:
    # none

    min.length <- get_attributes(ribo.object)$length_min
    max.length <- get_attributes(ribo.object)$length_max

    if ((range.lower < min.length | range.lower > range.upper)) {
        stop(
            "Param range.lower must be greater than or equal to the minimum
            length and less than range.upper.",
            call. = FALSE
        )
    } else if ((range.upper > max.length |
                range.upper < range.lower)) {
        stop(
            "Param range.upper must be less than or equal to the maximum length
            and greater than or equal to range.lower.",
            call. = FALSE
        )
    }
}

check_experiments <- function(ribo.object, experiments) {
    # Helper method that checks if the user-given experiments
    # are present in the current ribo file
    # Args:
    # ribo.object - S3 object of class ribo, contains the handle to the file
    # Return:
    # None
    ribo.experiments <- get_experiments(ribo.object)
    matched.experiments <- intersect(experiments, ribo.experiments)

    if (!length(matched.experiments)) {
        stop("Param 'experiments' contained no valid experiments.",
             call. = FALSE)
    }

    #deals with missing experiments
    check <- setdiff(experiments, matched.experiments)
    if (length(check)) {
        for (experiment in check) {
            warning("'", experiment, "'", " was not found.", call. = FALSE)
        }
        warning("Param 'experiments' contained experiments that were not ", 
                "found. The returned data frame ignores these experiments.",
                call. = FALSE)
    }
}
