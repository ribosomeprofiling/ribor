#' Get information about the .ribo file
#'
#' The function \code{\link{get_info}} provides information on the attributes, metadata,
#' and datasets of the ribo file.
#'
#' The \code{\link{get_info}} first provides information on the format version, left_span, right_span,
#' longest read length, shortest read length, metagene_radius, and reference model. The last element of the
#' returned list contains the information about the presence of coverage and RNA-seq data which are
#' optional datasets to include in a .ribo file.
#'
#' @param ribo.object ribo.object is an S4 object of class "Ribo"
#' @examples
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- Ribo(file.path)
#'
#' #retrieve information
#' get_info(sample)
#'
#' @return Returns a list containing a nested list of file attributes, a logical
#' value denoting whether the root file has additional metadata, and a
#' data.frame of information on each experiment
#'
#' @seealso \code{\link{Ribo}} to generate the necessary ribo.object parameter
#'
#' @importFrom rhdf5 h5ls h5readAttributes
#' @export
get_info <- function(ribo.object) {
    validObject(ribo.object)
    path   <- path(ribo.object)
    
    #retrieve an experiment list
    exp.list <- experiments(ribo.object)
    result <- get_attributes(ribo.object)
    has.metadata <- ("metadata" %in% names(result))
    
    if (has.metadata) {
        result <- result[-which(names(result) == "metadata")]
    }
    
    experiment.info <- get_content_info(path)
    result <- list(
        "has.metadata"    = has.metadata,
        "attributes"      = result,
        "experiment.info" = experiment.info
    )
    return(result)
}

#' Retrieves the metadata of an experiment
#'
#' \code{\link{get_metadata}} provides information on all of the user-inputted
#' metadata of an experiment. If the experiment is not found, then the
#' attributes of the root .ribo file is returned instead.
#'
#' @param ribo.object object of class 'ribo'
#' @param name The name of the experiment
#' @param print Logical value indicating whether or not to neatly print the output
#' @examples
#' #ribo object use case
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- Ribo(file.path)
#'
#' #the ribo file contains an experiment named 'Hela_1'
#' get_metadata(sample, "Hela_1")
#'
#' @return
#' If a valid experiment name is provided, a list of elements providing all of the metadata of the
#' experiment is returned.
#'
#' If the name is not provided and the root file has metadata, then a list of elements
#' providing all of the metadata found in the root file is returnend.
#'
#' @seealso \code{\link{Ribo}} to generate the necessary ribo.object parameter
#' @importFrom rhdf5 h5readAttributes
#' @importFrom yaml yaml.load
#' @export
get_metadata <- function(ribo.object,
                         name = NULL,
                         print = TRUE) {
    path <- path(ribo.object)
    exp.list <- experiments(ribo.object)
    file_path = "/"
    
    if (!is.null(name)) {
        file_path = paste("experiments/", name, sep = "")
        if (!(name %in% exp.list)) {
            stop("'",
                 name,
                 "'",
                 " is not a valid experiment name.",
                 call. = FALSE)
        }
    }
    
    attribute <- h5readAttributes(path, file_path)
    if ("metadata" %in% names(attribute)) {
        result <- yaml.load(string = attribute[["metadata"]])
        if (print) {
            print_metadata(result, 0)
            invisible(result)
        } else {
            return(result)
        }
    } else {
        stop("File does not have metadata.")
    }
}

print_metadata <- function(metadata, index) {
    for (i in seq_along(metadata)) {
        if (length(metadata[[i]]) >= 2) {
            cat(rep(" ", index * 3), paste(names(metadata)[i], ": \n", sep = ""))
            info <- metadata[[i]]
            print_metadata(info, index + 1)
        } else {
            cat(rep(" ", index * 3),
                paste(names(metadata)[i], ": ", metadata[[i]], "\n", sep = ""))
        }
    }
}

#' Provides a list of experiments from a .ribo file
#'
#' The function \code{\link{get_experiments}} provides a list of experiment names in the .ribo file.
#'
#' \code{\link{get_experiments}} returns a list of strings denoting the experiments. It obtains this
#' by reading directly from the .ribo file through the path of the 'ribo.object' parameter. To generate
#' the param 'ribo.object', call the \code{\link{Ribo}} function and provide the path to the .ribo file of interest.
#'
#' The user can then choose to create a subset from this list for any specific experiments of interest
#' for later function calls. Many functions that have the param 'experiment.list'
#'  call \code{\link{get_experiments}} to generate a default list of all experiments in the
#' .ribo file.
#' @examples
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- Ribo(file.path)
#'
#' #get a list of the experiments
#' get_experiments(sample)
#'
#' @seealso \code{\link{Ribo}} to generate the necessary ribo.object parameter
#' @param ribo.object S4 object of class "Ribo"
#' @return A list of the experiment names
#' @importFrom rhdf5 h5ls
#' @export
get_experiments <- function(ribo.object) {
    validObject(ribo.object)
    result <- h5ls(ribo.object@path)
    result <- result[result$group == "/experiments",]
    return(result$name)
}