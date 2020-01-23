#' Retrieves the coverage data for a given transcript
#'
#' The function \code{\link{get_coverage}} generates a DataFrame of coverage
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
#'                               experiment = experiments)
#'
#' @param ribo.object A 'ribo' object
#' @param name Name of the transcript
#' @param range.lower Lower bound of the read length
#' @param range.upper Upper bound of the read length
#' @param length Logical value that denotes if the coverage should be summed across read lengths
#' @param experiment List of experiments to obtain coverage information on
#' @param tidy Logical value denoting whether or not the user wants a tidy format
#' @param alias Option to report the transcripts as aliases/nicknames
#' @param compact Option to return a DataFrame with Rle and factor as opposed to a raw data.frame
#' @return A data table containing the coverage data
#' @seealso \code{\link{ribo}} to generate the necessary ribo.object parameter
#' @importFrom rhdf5 h5read
#' @importFrom tidyr gather
#' @importFrom S4Vectors Rle DataFrame
#' @importFrom methods as
#' @importFrom hash has.key
#' @export
get_coverage <- function(ribo.object,
                         name,
                         range.lower = length_min(ribo.object),
                         range.upper = length_max(ribo.object),
                         length = TRUE,
                         tidy = FALSE,
                         alias = FALSE,
                         compact = TRUE,
                         experiment = experiments(ribo.object)) {
    # check the parameters and prepare parameters for reading from ribo file
    if (missing(name)) stop("Please provide a transcript name.")
    validObject(ribo.object)
    matched.experiments <- initialize_coverage(ribo.object,
                                               alias,
                                               range.lower,
                                               range.upper,
                                               experiment)
    
    total.experiments   <- length(matched.experiments)
    info <- retrieve_transcript_info(ribo.object,
                                     name,
                                     alias)

    # read from the ribo file and generate a data.frame
    result <- fill_coverage(ribo.object,
                            total.experiments,
                            matched.experiments,
                            range.lower,
                            range.upper,
                            length,
                            alias,
                            info)
    
    # gather the data to make it tidy 
    if (tidy) {
      tidy.columns <- "experiment"
      if (!length) tidy.columns <- c("experiment", "length")
      
      result <- gather(result,
                       key = "position",
                       value = "count", 
                       -tidy.columns)
    }
    
    # convert to DataFrame with Rle encoding if necessary 
    if (compact) return(prepare_DataFrame(ribo.object, result))
    return(result)
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
    length.offset       <- length_offset(ribo.object)
    min.length          <- length_min(ribo.object)
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
            coverage <- t(h5read(path(ribo.object),
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
    matched.size <- length(matched.experiments)
    if (length) {
      return (data.frame(experiment = matched.experiments,
                         result, check.names = FALSE))
    }
    range <- range.upper - range.lower + 1
    return (data.frame(experiment  = rep(matched.experiments, each = range),
                       length      = rep(c(range.lower:range.upper), matched.size),
                       result, 
                       check.names = FALSE))
    
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
                                experiment) {
    #perform checks on the parameters
    check_alias(ribo.object, alias)
    check_lengths(ribo.object, range.lower, range.upper)
    # check_experiments(ribo.object, experiments)
    ribo.experiments    <- experiments(ribo.object)
    
    #generate list of experiments also present in the ribo file that have coverage data
    filter.experiments  <- intersect(experiment, ribo.experiments)
    matched.experiments <- check_coverage(ribo.object, filter.experiments)
    matched.experiments <- intersect(matched.experiments, filter.experiments)

    return(matched.experiments)
}


retrieve_transcript_info <- function(ribo.object,
                                     name,
                                     alias) {
    #helper method that checks the alias parameter and retrieves the corresponding information
    if (alias) {
        if (!has.key(name, alias_hash(ribo.object))) {
            stop("Alias name was not found.", call. = FALSE)
        }
        name <- alias_hash(ribo.object)[[name]]
    }

    #generate offsets
    info <- as.vector(unlist(transcript_info(ribo.object)[[name]]))
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
    colnames(result) <- as.character(seq_len(transcript.length))
    return(result)
}


create_dataframe <- function (length,
                              range.lower,
                              range.upper,
                              matched.experiments,
                              compact,
                              matrix) {
    # Given a matrix of coverage data, create_dataframe generates the correct
    # DataFrame based on a set of parameters
    #
    # Returns:
    # Data table that wraps the matrix with the correct and appropriate labels

    matched.size <- length(matched.experiments)
    if (length) {
        return (DataFrame(experiment = Rle(factor(matched.experiments)),
                          matrix, check.names = FALSE))
    }
    range <- range.upper - range.lower + 1
    return (DataFrame(
        experiment  = Rle(factor(rep(matched.experiments, each = range))),
        length = Rle(factor(rep(c(range.lower:range.upper), matched.size))),
        matrix, check.names = FALSE
    ))
}


check_coverage <- function(ribo.object, experiment) {
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

    path <- path(ribo.object)
    #obtain the coverage data
    table <- get_content_info(path)
    has.coverage <- table[table$coverage == TRUE,]
    has.coverage <- has.coverage$experiment

    #find the experiments in the experiment.list that do not have
    #coverage and print warnings
    check <- setdiff(experiment, has.coverage)
    if (length(check)) {
        for (exp in check) {
            warning("'", exp, "'",
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
