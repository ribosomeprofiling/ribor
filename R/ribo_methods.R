#'@include create_ribo.R
#'@include ribo_class.R
NULL

#' @describeIn ribo Displaying the file contents
#' @examples 
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- create_ribo(file.path)
#' 
#' show(sample)
#' @export
setMethod(f = "show",
          signature = "ribo",
          definition = function(object) {
              file.values <- list("format version"   = object@format.version,
                                  "reference"        = object@reference,
                                  "min read length"  = object@length.min,
                                  "max read length"  = object@length.max,
                                  "left span"        = object@left.span,
                                  "right span"       = object@right.span,
                                  "transcript count" = length(object@transcript.info),
                                  "has.metadata"     = object@has.metadata,
                                  "metagene radius"  = object@metagene.radius,
                                  "has.alias"        = !is.empty(object@transcript.alias))
              file.info <- data.frame(info = names(file.values),
                                      " " = unlist(unname(file.values)), check.names = FALSE,
                                                         stringsAsFactors = FALSE,
                                                         fix.empty.names = FALSE)
              
              #center by reformatting the output of the data tables
              file.info   <- format_file_info(file.info)
              experiment.info <- format_experiment_info(object@experiment.info)

              print_output(file.info, experiment.info)
              invisible(object)
          })

format_experiment_info <- function(experiment.info) {
    #helper method to format the experiment info
    experiment.info <- lapply(as.list(experiment.info), as.character)
    experiment.info <- as.data.frame(experiment.info, stringsAsFactors = FALSE)

    #format in the case that the experiment names are really long
    exp.val.name <- unlist(names(experiment.info)[-1])
    exp.val      <- unlist(experiment.info[, -1])
    name.size    <- max(nchar(exp.val))
    name.size    <- max(name.size,
                        max(nchar(exp.val.name)))


    names(experiment.info)[-1] <- format(exp.val.name,
                                         width = name.size,
                                         justify = "right")

    experiment.info <-  as.data.frame(format(experiment.info,
                                             width   = name.size,
                                             justify = "right"))

    exp.title <- unlist(names(experiment.info[1]))
    exp.name  <- unlist(experiment.info[, 1])
    name.size <- max(nchar(exp.name), nchar(exp.title))

    names(experiment.info)[1] <- format(names(experiment.info)[1],
                                        width = name.size,
                                        justify = "right")

    experiment.info[, 1]      <-  format(experiment.info[, 1],
                                         width = name.size,
                                         justify = "right")
    return(experiment.info)
}

format_file_info <- function(file.info) {
    #helper method to format the file info
    name.width <- max(nchar(unname(unlist(file.info))))
    file.info   <- format(file.info,
                          justify = "right",
                          width = name.width)

    return(file.info)
}

print_output <- function (file.info,
                          experiment.info) {
    # Final helper method that prints the output of the formatted attributes 
    cat("General File Information:\n")
    print(data.frame(file.info, check.names = FALSE),
          row.names = FALSE,
          quote = FALSE)
    cat("\n")
    cat("Dataset Information:\n")
    print(data.frame(experiment.info, check.names = FALSE),
          row.names = FALSE,
          quote = FALSE)
}

setGeneric("path", function(object) standardGeneric("path"))
setMethod("path", "ribo", function(object) object@path)

setGeneric("experiments", function(object) standardGeneric("experiments"))
setMethod("experiments", "ribo", function(object) object@experiments)

setGeneric("format_version", function(object) standardGeneric("format_version"))
setMethod("format_version", "ribo", function(object) object@format.version)

setGeneric("reference", function(object) standardGeneric("reference"))
setMethod("reference", "ribo", function(object) object@reference)

setGeneric("length_min", function(object) standardGeneric("length_min"))
setMethod("length_min", "ribo", function(object) object@length.min)

setGeneric("length_max", function(object) standardGeneric("length_max"))
setMethod("length_max", "ribo", function(object) object@length.max)

setGeneric("left_span", function(object) standardGeneric("left_span"))
setMethod("left_span", "ribo", function(object) object@left.span)

setGeneric("right_span", function(object) standardGeneric("right_span"))
setMethod("right_span", "ribo", function(object) object@right.span)

setGeneric("metagene_radius", function(object) standardGeneric("metagene_radius"))
setMethod("metagene_radius", "ribo", function(object) object@metagene.radius)

setGeneric("length_offset", function(object) standardGeneric("length_offset"))
setMethod("length_offset", "ribo", function(object) object@length.offset)

setGeneric("has_metadata", function(object) standardGeneric("has_metadata"))
setMethod("has_metadata", "ribo", function(object) object@has.metadata)

setGeneric("experiment_info", function(object) standardGeneric("experiment_info"))
setMethod("experiment_info", "ribo", function(object) object@experiment.info)

setGeneric("transcript_info", function(object) standardGeneric("transcript_info"))
setMethod("transcript_info", "ribo", function(object) object@transcript.info)

setGeneric("alias_hash", function(object) standardGeneric("alias_hash"))
setMethod("alias_hash", "ribo", function(object) object@transcript.alias)

setGeneric("transcript_original", function(object) standardGeneric("transcript_original"))
setMethod("transcript_original", "ribo", function(object) object@transcript.original)

setGeneric("rangeLower", function(object) standardGeneric("rangeLower"))
setMethod("rangeLower", "ribo", function(object) get_read_lengths(object)[1])

setGeneric("rangeUpper", function(object) standardGeneric("rangeUpper"))
setMethod("rangeUpper", "ribo", function(object) get_read_lengths(object)[2])

setGeneric("original_hash", function(object) standardGeneric("original_hash"))
setMethod("original_hash", "ribo", function(object) object@transcript.original)

# validity method 
setValidity("ribo", validate_ribo)
