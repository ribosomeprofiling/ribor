check_alias <- function(ribo.object,
                        alias) {
    has.alias <- !is.empty(alias_hash(ribo.object)) &&
                 !is.empty(original_hash(ribo.object))
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
    # ribo.object 
    # range.lower lowest read length
    # range.upper highest read length
    #
    # Return:
    # none

    min.length <- length_min(ribo.object)
    max.length <- length_max(ribo.object)

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

check_experiments <- function(ribo.object, experiment) {
    # Helper method that checks if the user-given experiments
    # are present in the current ribo file
    # Args:
    # ribo.object
    # Return:
    # None
    ribo.experiments <- experiments(ribo.object)
    matched.experiments <- intersect(experiment, ribo.experiments)

    if (!length(matched.experiments)) {
        stop("Param 'experiments' contained no valid experiments.",
             call. = FALSE)
    }

    #deals with missing experiments
    check <- setdiff(experiment, matched.experiments)
    if (length(check)) {
        for (exp in check) {
            warning("'", exp, "'", " was not found.", call. = FALSE)
        }
        warning("Param 'experiment' contained experiments that were not ", 
                "found. The returned data frame ignores these experiments.",
                call. = FALSE)
    }
}
