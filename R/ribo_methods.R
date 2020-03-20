#'@include create_ribo.R
#'@include ribo_class.R
NULL

#' @rdname Ribo-class
#' @aliases Ribo-class
#' @examples 
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- Ribo(file.path)
#' 
#' show(sample)
#' 
#' @export
setMethod(f = "show",
          signature = "Ribo",
          definition = function(object) {
              file.values <- list("format version"   = format_version(object),
                                  "reference"        = reference(object),
                                  "min read length"  = length_min(object),
                                  "max read length"  = length_max(object),
                                  "left span"        = left_span(object),
                                  "right span"       = right_span(object),
                                  "transcript count" = length(transcript_info(object)),
                                  "has.metadata"     = has_metadata(object),
                                  "metagene radius"  = metagene_radius(object),
                                  "has.alias"        = !is.empty(alias_hash(object)))
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

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("path", "Ribo", function(object) object@path)

setGeneric("experiments", function(object) standardGeneric("experiments"))

#' @rdname Ribo-class
#' @aliases Ribo-class
#' @param object Ribo object
#' @docType methods
#' @export
setMethod("experiments", "Ribo", function(object) object@experiments)


setGeneric("format_version", function(object) standardGeneric("format_version"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("format_version", "Ribo", function(object) object@format.version)


setGeneric("reference", function(object) standardGeneric("reference"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("reference", "Ribo", function(object) object@reference)


setGeneric("length_min", function(object) standardGeneric("length_min"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("length_min", "Ribo", function(object) object@length.min)

setGeneric("length_max", function(object) standardGeneric("length_max"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("length_max", "Ribo", function(object) object@length.max)

setGeneric("left_span", function(object) standardGeneric("left_span"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("left_span", "Ribo", function(object) object@left.span)

setGeneric("right_span", function(object) standardGeneric("right_span"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("right_span", "Ribo", function(object) object@right.span)

setGeneric("metagene_radius", function(object) standardGeneric("metagene_radius"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("metagene_radius", "Ribo", function(object) object@metagene.radius)


setGeneric("length_offset", function(object) standardGeneric("length_offset"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("length_offset", "Ribo", function(object) object@length.offset)

setGeneric("has_metadata", function(object) standardGeneric("has_metadata"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("has_metadata", "Ribo", function(object) object@has.metadata)

setGeneric("experiment_info", function(object) standardGeneric("experiment_info"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("experiment_info", "Ribo", function(object) object@experiment.info)


setGeneric("transcript_info", function(object) standardGeneric("transcript_info"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("transcript_info", "Ribo", function(object) object@transcript.info)

setGeneric("alias_hash", function(object) standardGeneric("alias_hash"))
 
#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("alias_hash", "Ribo", function(object) object@alias.hash)


setGeneric("original_hash", function(object) standardGeneric("original_hash"))

#' @rdname Ribo-class
#' @param object Ribo object
#' @aliases Ribo-class
#' @docType methods
#' @export
setMethod("original_hash", "Ribo", function(object) object@original.hash)
