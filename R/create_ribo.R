#' @include ribo_class.R
#' @include ribo_methods.R
NULL

initialize_transcript_info <- function(transcript.names,
                                       transcript.lengths) {
    # helper method that generates a hash of the length offset and length
    # of each transcript
    
    num.transcripts <- length(transcript.names)
    if (num.transcripts <= 0) {
        stop("Please make sure that there is a positive amount of transcripts.")
    }
    
    hash.value <- rep(list(c("offset" = 0, "length" = 0)), 
                      length = num.transcripts)
    names(hash.value)         <- transcript.names
    hash.value[[1]]["length"] <-  transcript.lengths[[1]]
    
    i <- 2
    while (i <= num.transcripts) {
        hash.value[[i]][["offset"]] <- hash.value[[i - 1]][["length"]] + 
            hash.value[[i - 1]][["offset"]]
        hash.value[[i]][["length"]] <- transcript.lengths[[i]]
        i <- i + 1
    }
    
    transcript.info <- hash(hash.value)
    length.offset <- hash.value[[num.transcripts]][["offset"]] +
        hash.value[[num.transcripts]][["length"]]
    names(length.offset) <- NULL
    
    return(list(transcript.info = transcript.info,
                length.offset   = length.offset))
}

#' Set the aliases of a ribo object 
#'
#' The function \code{\link{set_aliases}} allows the user to add aliases to a valid
#' ribo object.
#'
#' If there is a different naming convention from the default appris transcriptome,
#' there may be no simple way to generate convenient aliases from the original reference names.
#' As a result, the user can first generate the ribo object and get the reference names, use custom
#' (and likely more intricate) functions to generate a list of aliases, and then pass in a character
#' vector of these aliases. The character vector should match the order of and correspond to the
#' list of reference names retrieved from \code{\link{get_reference_names}}
#'
#' @param ribo.object A 'ribo' object
#' @param rename A function that renames original transcript name into an alias
#' @export
#' @importFrom rhdf5 h5read
#' @importFrom hash hash is.empty
#' @importFrom methods validObject
#' @return A modified 'ribo' object that contains alias information
#' @examples
#' #generate a ribo object with transcript nicknames/aliases
#' file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")
#' sample <- Ribo(file.path)
#' sample <- set_aliases(ribo.object = sample,
#'                       rename = rename_default)
set_aliases <- function(ribo.object, rename) {
    if (!is.empty(alias_hash(ribo.object)) ||
        !is.empty(original_hash(ribo.object))) {
        warning("Ribo object already has a naming convention. Aliases will be overriden.",
                call. = FALSE)
    }
    path <- path(ribo.object)
    original   <- h5read(path,
                         name = "reference/reference_names")
    
    num.transcripts <- length(original)
    alias <- rename_transcripts(path, rename)
    ribo.object@alias.hash    <- hash(keys = alias, values = original)
    ribo.object@original.hash <- hash(keys = original, values = alias)
    return(ribo.object)
}

#' Creates an S4 object of class "Ribo"
#'
#' \code{\link{Ribo}} creates a "Ribo" object. It creates a path, extracts the root folder attributes,
#' and provides information about the reference transcript names and lengths
#'
#' An important option is the param 'rename' which allows the user to nickname long
#' transcript names. For the appris human transcriptome, a default function \code{\link{rename_default}}
#' has been provided. In subsequent calls of certain functions, the user can make use of these renamed
#' references with the 'alias' parameter.
#'
#' This object is required as an argument for almost all of the functions in this package, and all of the
#' functions in this package can accept the returned object of this function. This object is not meant to
#' be modified or changed by the user. It is meant to serve as an intermediary between the .ribo file adevnd
#' an R environment by creating an object that holds pertinent information.
#'
#' The information stored in this object include the .ribo file path, the list of experiments,
#' the format version, the reference model, the maximum read length, the minimum read length, the left span,
#' the right span, and other information about the transcript information.
#'
#'
#' @param name The path to the .ribo file
#' @param rename A function that renames the original transcript or an already generated
#' character vector of aliases
#' @return Returns an S4 object of class "Ribo" containing a path to the HDF5 file,
#'         various attributes in the root folder, and information about the transcripts 
#'         such as names and lengths
#' @importFrom rhdf5 H5Fopen h5readAttributes h5ls h5read
#' @importFrom hash hash is.empty
#' @importFrom methods new
#' @importFrom tools file_path_as_absolute
#' @examples
#' #generate a ribo object with transcript nicknames/aliases
#' file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")
#' sample <- Ribo(file.path, rename = rename_default )
#' @seealso
#' If a ribo object is already generated but aliases want to be added or updated, use the
#' \code{\link{set_aliases}} function.
#' @export
Ribo <- function(name, rename = NULL) {
    ribo.path   <- file_path_as_absolute(name)
    attributes <- h5readAttributes(ribo.path, name = "/")
    transcript.names   <- h5read(ribo.path,
                                 name = "reference/reference_names")
    
    transcript.lengths <- h5read(ribo.path, 
                                 name = "reference/reference_lengths")
    file_info          <- h5ls(ribo.path, recursive = TRUE, all = FALSE)
    transcript.info    <- initialize_transcript_info(transcript.names,
                                                     transcript.lengths)
    
    has.metadata   <- ("metadata" %in% names(attributes))
    
    ribo.object <- new("Ribo", 
                       path            = ribo.path,
                       experiments     = file_info[file_info$group == "/experiments",]$name,
                       format.version  = as.integer(attributes$format_version),
                       reference       = attributes$reference,
                       length.min      = attributes$length_min,
                       length.max      = attributes$length_max,
                       left.span       = attributes$left_span,
                       right.span      = attributes$right_span,
                       metagene.radius = attributes$metagene_radius,
                       length.offset   = transcript.info[['length.offset']],
                       has.metadata    = has.metadata,
                       experiment.info = get_content_info(ribo.path),
                       transcript.info = transcript.info[['transcript.info']])
    
    if (!is.null(rename)) {
        ribo.object <- set_aliases(ribo.object, rename)
    }
    return(ribo.object)
}
