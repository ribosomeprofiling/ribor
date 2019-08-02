#' Retrieves the coverage data for a given transcript
#'
#' The function \code{\link{get_coverage}} generates a data.table of coverage
#' data over the length of a given transcript.
#'
#' The function \code{\link{get_coverage}} first checks the experiments in the
#' 'experiments' parameter to see if they are present in the .ribo file.
#' It will then check these experiments for coverage data which is an optional
#' dataset. As a result, this function safe guards against experiments that
#' do not have coverage data, but it also, by default, includes all of the
#' experiments in a file in the experiments' parameter.
#'
#' The function checks the coverage of one transcript at a time at
#' each read length from 'range.lower' to 'range.upper', inclusive. However,
#' the parameter 'length' allows the user to obtain the coverage
#' information of a transcript across the range of read lengths indicated by
#' 'range.lower' and 'range.upper'.
#'
#' If the ribo.object is generated with aliases, the 'alias' parameter, if set to TRUE,
#' allows the user to use the alias of the transcript as the 'name' parameter instead of
#' the original transcript name.
#'
#' @examples
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- create_ribo(file.path)
#'
#' #get the experiments of interest that also contain coverage data
#' experiments <- c("Hela_1", "Hela_2", "Hela_3", "WT_1")
#'
#' #the ribo file contains a transcript named 'MYC'
#' coverage.data <- get_coverage(ribo.object = sample,
#'                               name = "MYC",
#'                               range.lower = 2,
#'                               range.upper = 5,
#'                               length = TRUE,
#'                               experiments = experiments)
#'
#' @param ribo.object A 'ribo' object
#' @param name Name of the transcript
#' @param range.lower Lower bound of the read length
#' @param range.upper Upper bound of the read length
#' @param length Logical value that denotes if the coverage should be summed across read lengths
#' @param experiments List of experiments to obtain coverage information on
#' @param tidy Logical value denoting whether or not the user wants a tidy format
#' @param alias Option to report the transcripts as aliases/nicknames
#' @return A data table containing the coverage data
#' @seealso \code{\link{ribo}} to generate the necessary ribo.object parameter
#' @importFrom rhdf5 h5read
#' @importFrom tidyr gather
#' @importFrom hash has.key
#' @export
get_coverage <- function(ribo.object,
                         name,
                         range.lower = rangeLower(ribo.object),
                         range.upper = rangeUpper(ribo.object),
                         length = TRUE,
                         tidy = FALSE,
                         alias = FALSE,
                         experiments = get_experiments(ribo.object)) {
    matched.experiments <- initialize_coverage(ribo.object,
                                               alias,
                                               range.lower,
                                               range.upper,
                                               experiments)
    total.experiments   <- length(matched.experiments)

    info <- retrieve_transcript_info(ribo.object,
                                     name,
                                     alias)

    result <- fill_coverage(ribo.object,
                            total.experiments,
                            matched.experiments,
                            range.lower,
                            range.upper,
                            length,
                            alias,
                            info)

    return <- create_datatable(length,
                               range.lower,
                               range.upper,
                               matched.experiments,
                               result)

    return(check_tidy_coverage(return, tidy, length))
}


fill_coverage <- function(ribo.object,
                          total.experiments,
                          matched.experiments,
                          range.lower,
                          range.upper,
                          length,
                          alias,
                          info) {
    current.offset      <- info[1]
    transcript.length   <- info[2]
    length.offset       <- ribo.object@length.offset
    min.length          <- get_read_lengths(ribo.object)[1]
    read.range          <- range.upper - range.lower + 1

    result <- initialize_matrix_coverage(total.experiments,
                                         transcript.length,
                                         read.range,
                                         length)
    #fills in the entries of the matrix
    for (i in seq(total.experiments)) {
        #for the given experiment, get its coverage
        experiment <- matched.experiments[i]
        path <- get_dataset_path(experiment, "coverage")
        #get the coverage information for the specific transcript at each read length
        for (current.length in range.lower:range.upper) {
            correct.length <- 1 + (current.length - min.length) * length.offset
            coverage.start <- correct.length + current.offset
            coverage.stop  <- coverage.start + transcript.length - 1
            coverage <- t(h5read(ribo.object@handle,
                                 path,
                                 index = list(coverage.start:coverage.stop)))
            coverage <- as.integer(coverage)

            # compute the correct offset to store the coverage information in each case
            if (length) {
                result[i,] <- result[i,] + coverage
            } else {
                current.experiment     <- (i - 1) * read.range
                current.read           <- current.length - range.lower
                current.index          <- current.experiment + current.read + 1
                result[current.index,] <- coverage
            }
        }
    }
    return(result)
}

get_dataset_path <- function(experiment,
                             dataset) {
  return(paste("/experiments/",
               experiment,
               "/", dataset, "/", dataset, sep =""))
}

initialize_coverage <- function(ribo.object,
                                alias,
                                range.lower,
                                range.upper,
                                experiments) {
    #perform checks on the parameters
    check_alias(ribo.object, alias)
    check_lengths(ribo.object, range.lower, range.upper)
    check_experiments(ribo.object, experiments)
    ribo.experiments    <- get_experiments(ribo.object)

    #generate list of experiments also present in the ribo file that have coverage data
    filter.experiments  <- intersect(experiments, ribo.experiments)
    matched.experiments <- check_coverage(ribo.object, filter.experiments)
    matched.experiments <- intersect(matched.experiments, filter.experiments)

    return(matched.experiments)
}


retrieve_transcript_info <- function(ribo.object,
                                     name,
                                     alias) {
    #helper method that checks the alias parameter and retrieves the corresponding information
    if (alias) {
        if (!has.key(name, ribo.object@transcript.alias)) {
            stop("Alias name was not found.", call. = FALSE)
        }
        name <- ribo.object@transcript.alias[[name]]
    }

    #generate offsets
    info <- as.vector(unlist(ribo.object@transcript.info[[name]]))
    if (is.null(info)) {
        stop("Transcript name was not found. Check name and 'alias' parameter.",
             call. = FALSE)
    }
    return(info)
}

initialize_matrix_coverage <- function(total.experiments,
                                       transcript.length,
                                       read.range,
                                       length) {
    #create matrix of correct size based on param transcript
    result <- matrix()
    if (length) {
        result <- matrix(0L,
                         nrow = 1 * total.experiments,
                         ncol = transcript.length)
    } else {
        result <- matrix(nrow = read.range * total.experiments,
                         ncol = transcript.length)
    }
    colnames(result) <- seq_len(transcript.length)
    return(result)
}

check_tidy_coverage <- function(return, tidy, length) {
    #checks the tidy case for coverage
    if (tidy) {
        #determine which columns to remove
        tidy.columns <- "experiment"
        if (!length) {
            tidy.columns <- c("experiment", "read.length")
        }
        return <- setDT(gather(
            return,
            key = "position",
            value = "count",-tidy.columns
        ))
    }
    return(return)
}


create_datatable <- function (length,
                              range.lower,
                              range.upper,
                              matched.experiments,
                              matrix) {
# Given a matrix of coverage data, create_datatable generates the correct
# data.table based on a set of parameters
#
# Returns:
# Data table that wraps the matrix with the correct and appropriate labels

    matched.size <- length(matched.experiments)
    if (length) {
        return (data.table(experiment = matched.experiments,
                           matrix))
    }
    range <- range.upper - range.lower + 1
    return (data.table(
        experiment  = rep(matched.experiments, each = range),
        read.length = rep(c(range.lower:range.upper), matched.size),
        matrix
    ))
}


check_coverage <- function(ribo.object, experiments) {
    # helper function that both generates a list of experiments with coverage
    # data using get_info and produces a warning message for each experiment
    # (provided by the user) that does not have coverage data
    #
    # Args:
    # ribo.object: S3 object of class "ribo"
    # experiment.list: list of experiments inputted by the user
    #
    # Returns:
    # A list of experiments in the ribo.object that have coverage data

    handle <- ribo.object@handle
    #obtain the coverage data
    table <- get_content_info(handle)
    has.coverage <- table[table$coverage == TRUE,]
    has.coverage <- has.coverage$experiment

    #find the experiments in the experiment.list that do not have
    #coverage and print warnings
    check <- setdiff(experiments, has.coverage)
    if (length(check)) {
        for (experiment in check) {
            warning("'", experiment, "'",
                    " did not have coverage data.",
                    call. = FALSE)
        }
        warning(
            "Param 'experiments' contains experiment names that did not have coverage data.",
            "The return value ignores these experiments.",
            call. = FALSE
        )
    }

    #return a list of experiments with coverage
    return(has.coverage)
}
