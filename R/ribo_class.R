### =========================================================================
### Factor objects
### -------------------------------------------------------------------------
###
### The Factor class serves a similar role as factor in base R except that
### the levels of a Factor object can be any vector-like object.
###
#' @rdname Ribo-class
#' @aliases Ribo-class
#' @export
setClass(
    "Ribo",
    slots =      c(path            = "character",
                   experiments     = "character",
                   format.version  = "integer",
                   reference       = "character",
                   length.min      = "integer",
                   length.max      = "integer",
                   left.span       = "integer",
                   right.span      = "integer",
                   metagene.radius = "integer",
                   length.offset   = "numeric",
                   has.metadata    = "logical",
                   experiment.info = "data.frame",
                   transcript.info = "hash",
                   alias.hash      = "hash",
                   original.hash   = "hash"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

setValidity("Ribo", function(object) {
    # The validity method is to protect against function calls with ribo
    # objects that have been modified in a way that could generate incorrect
    # output
    
    # Generate all of the values found on the .ribo file itself 
    attributes         <- h5readAttributes(object@path, name = "/")
    file_info          <- h5ls(object@path, recursive = TRUE, all = FALSE)
    transcript.names   <- h5read(object@path,
                                 name = "reference/reference_names")
    transcript.lengths <- h5read(object@path, 
                                 name = "reference/reference_lengths")
    transcript.info    <- initialize_transcript_info(transcript.names,
                                                     transcript.lengths)
    
    errors <- character()
    expected_radius <- attributes$metagene_radius
    
    if (expected_radius != object@metagene.radius) {
        msg <- paste("Metagene radius is different from the .ribo file.")
        errors <- c(errors, msg)
    }
    
    num_experiments <- length(object@experiments)
    
    if (num_experiments <= 0) {
        msg <- paste("File appears to have no experiments.")
        errors <- c(errors, msg)
    }
    
    if (!all(file_info[file_info$group == "/experiments",]$name == object@experiments)) {
        msg <- paste("Experiments are different from the .ribo file.")
        errors <- c(errors, msg)
    }
    
    
    length.min <- object@length.min
    length.max <- object@length.max 
    
    if (length.min > length.max) {
        msg <- paste("The minimum read length should be less than or equal ",
                     "to the maximum read length.", sep = "")
        errors <- c(errors, msg)
    }
    
    if (length.min != attributes$length_min) {
        msg <- paste("The minimum read length differs from the .ribo file minimum read length.")
        errors < c(errors, msg)
    }
    
    if (length.max != attributes$length_max) {
        msg <- paste("The maximum read length differs from the .ribo file maximum read length.")
        errors <- c(errors, msg)
    }
    
    if (object@length.offset != transcript.info[['length.offset']]) {
        msg <- paste("The length offset differs from the .ribo file read length")
        errors <- c(errors, msg)
    }
    
    if (object@left.span != attributes$left_span) {
        msg <- paste("The left span differs from the .ribo file left span.")
        errors <- c(errors, msg)
    }
    
    if (object@right.span != attributes$right_span) {
        msg <- paste("The right span differs from the .ribo file right span.")
        errors <- c(errors, msg)
    }
    
    if (!all(object@experiment.info == get_content_info(object@path))) {
        msg <- paste("The experiment info differs from the .ribo file experiment info.")
        errors <- c(errors, msg)
    }
    
    if (length(errors) == 0) TRUE else errors
})