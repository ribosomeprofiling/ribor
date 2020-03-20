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


#' Ribo Class
#' 
#' The Ribo object serves as the main utility vehicle for the ribor package.
#' Specifically, it allows the user to interface with a .ribo file in the R
#" environment through the use of the ribor package. Almost all functions in
#' ribor rely on the Ribo object to read, visualize, and inspect the
#' contents of the .ribo file. 
#' The information stored in this object include the .ribo file path, the 
#' list of experiments, the format version, the reference model, the minimum
#' read length, maximum read length, the left span, the right span, and 
#' other transcript information.
#' 
#' Note that the path parameter takes in a file path and stores it. While 
#' using the package, be sure to not to move or change the location of the 
#' .ribo file. The default names of the transcripts may be difficult to use depending on the 
#' settings used to generate the .ribo file. As a result, we have provided a 
#' rename parameter that integrates well with the Appris reference transcriptome.
#' Users may also define a simple function that processes a given default 
#' transcript name in a one-to-one manner to another custom alias. 
#' 
#' 
#' @param path The path to the .ribo file
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
#' @rdname Ribo-class
#' @export
#' @aliases experiments alias_hash experiment_info format_version has_metadata
#' @aliases left_span length_max length_min length_offset metagene_radius
#' @aliases length_offset metagene_radius original_hash path reference right_span transcript_info
Ribo <- function(path, rename = NULL) {
    ribo.path   <- file_path_as_absolute(path)
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
