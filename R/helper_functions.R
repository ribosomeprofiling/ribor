#' Retrieves a list of reference names
#'
#' Gets a list of reference names by reading directly from the .ribo file
#'
#' @param ribo.object A 'ribo' object
#' @return a list of the reference names
#' @importFrom rhdf5 h5read
#' @export
#' @examples
#' #generate a ribo object with transcript nicknames/aliases
#' file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")
#' sample <- create_ribo(file.path)
#'
#' #get the reference names
#' names <- get_reference_names(sample)
get_reference_names <- function(ribo.object) {
    # Retrieves the reference transcript names
    check_ribo(ribo.object)
    return(h5read(ribo.object@path,
                  name = "reference/reference_names"))
}

get_reference_lengths <- function(ribo.object) {
    # Retrieves the reference transcript lengths
    check_ribo(ribo.object)
    row.names <- h5read(ribo.object@path,
                        name = "reference/reference_names")
    lengths   <- h5read(ribo.object@path,
                        name = "reference/reference_lengths")
    return(data.frame(transcript = row.names, length = lengths))
}


#' Renames the transcripts
#'
#' The function \code{\link{rename_transcripts}} strives to make the transcript names less
#' cumbersome to write and easier to use.
#'
#' Transcript names found in a .ribo file can often be long and inconvenient to use.
#' As a result, this function allows the user to rename the transcripts.
#'
#' Often times, a short function can be used on the ribo file reference names
#' to split and extract a more convenient name, and a function with a similar input and
#' output to {\code{\link{rename_default}}} can be passed in.
#'
#' However, if there is no simple function that takes the original name and renames it into
#' a unique alias, then the user can provide a character vector of the same length as the number of
#' transcripts in the ribo file. This character vector would provide aliases that match the order
#' of the original reference names returned by the {\code{\link{get_reference_names}}} function.
#'
#' @param ribo a path to the ribo file or a 'ribo' object
#' @param rename A function that renames the original transcript or an already generated
#' character vector of aliases
#' @importFrom rhdf5 h5read
#' @seealso
#' {\code{\link{rename_default}}} to view expected input and output of a 'rename' function
#' {\code{\link{create_ribo}}} to generate a ribo object
#' @examples
#' file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")
#' sample <- create_ribo(file.path, rename = rename_default)
#'
#' aliases <- rename_transcripts(sample, rename = rename_default)
#' @export
#' @return A character vector denoting the renamed transcript aliases
rename_transcripts <- function(ribo, rename) {
    #ensure that the ribo path is retrieved
    ribo.path <- ribo
    if (is.ribo(ribo)) {
        ribo.path <- ribo@path
    }

    #handle the function case and the vector case
    simplify <- NULL
    original <- h5read(ribo.path, 
                       name = "reference/reference_names")

    if (is.function(rename)) {
        simplify <- vapply(X = original,
                           FUN = rename,
                           FUN.VALUE = "character")
    } else if (is.vector(rename)) {
        #handle the case that rename is a list as well
        simplify <- unlist(rename)
    } else {
        stop("Param 'rename' must be a function or a vector.")
    }

    if (anyDuplicated(simplify)) {
        stop("Invalid param 'rename'. Redundant values found in the alias.",
             call. = FALSE)
    }

    # helper function to rename the long transcript names
    return(simplify)
}

#' Rename function for appris transcriptome naming convention
#'
#' The function {\code{\link{rename_default}}} is the default renaming function for the
#' appris human transcriptome. It takes one single transcript name and returns a simplified
#' alias.
#' @examples
#' original <- paste("ENST00000613283.2|ENSG00000136997.17|",
#'                   "OTTHUMG00000128475.8|-|MYC-206|MYC|1365|protein_coding|",
#'                   sep = "")
#' alias <- rename_default(original)
#' @param x Character denoting original name of the transcript
#' @return Character denoting simplified name of the object
#' @export
rename_default <- function(x) {
    return(unlist(strsplit(x, split = "|", fixed = TRUE))[5])
}

get_content_info <- function(ribo.path) {
    file_info     <- h5ls(ribo.path, recursive = TRUE, all = FALSE)
    experiments   <- file_info[file_info$group == "/experiments", ]$name
    length <- length(experiments)

    #creates the separate lists for reads, coverage, rna.seq, and metadata
    #to eventually put in a data frame
    reads.list    <- vector(mode = "integer", length = length)
    coverage.list <- vector(mode = "logical", length = length)
    rna.seq.list  <- vector(mode = "logical", length = length)
    metadata.list <- vector(mode = "logical", length = length)

    #ls function provides information about the contents of each experiment
    ls <- h5ls(ribo.path)

    #loop over all of the experiments
    for (i in seq(length)) {
        experiment <- experiments[i]
        #gathers information on the number of reads for each experiment by looking at
        #the attributes
        name           <-
            paste("/experiments/", experiment, sep = "")
        attribute      <- h5readAttributes(ribo.path, name)
        reads.list[i]     <- attribute[["total_reads"]]

        #creates separate logical lists to denote the presence of
        #reads, coverage, RNA-seq, metadata
        metadata.list[i]  <- ("metadata" %in% names(attribute))

        group.contents <- ls[ls$group == name,]
        group.names    <- group.contents$name

        coverage.list[i]  <- ("coverage" %in% group.names)
        rna.seq.list[i]   <- ("rnaseq" %in% group.names)
    }

    experiments.info       <- data.frame(
        experiment  = experiments,
        total.reads = reads.list,
        coverage    = coverage.list,
        rna.seq     = rna.seq.list,
        metadata    = metadata.list,
        stringsAsFactors = FALSE
    )
    return(experiments.info)
}


get_attributes <- function(ribo.object) {
    # Retrieves the attributes of the ribo.object
    path  <- ribo.object@path
    attribute <- h5readAttributes(path, "/")
    return(attribute[-which(names(attribute) == "time")])
}

get_read_lengths <- function(ribo.object) {
    # Retrieves the minimum and maximum read lengths
    #
    # get_read_lengths finds the minimum and maximum read lengths of the .ribo file
    attributes <- get_attributes(ribo.object)
    result <- c(attributes$length_min, attributes$length_max)
    return(result)
}



make_dataframe <- function(ribo.object,
                           matched.list,
                           range.info,
                           conditions,
                           matrix) {
    alias      <- conditions[["alias"]]
    ref.names <- get_reference_names(ribo.object)
    # helper function that creates a polished dataframe out of a filled in matrix
    if (alias) {
        original <- ref.names
        ref.names <- vector(mode = "character", length = length(original))
        for (i in seq(length(ref.names))) {
            ref.names[i] <- ribo.object@transcript.original[[original[[i]]]]
        }
    }
    return(help_make_dataframe(ref.names,
                               conditions,
                               matched.list,
                               range.info,
                               matrix))
}

help_make_dataframe <- function(ref.names,
                                conditions,
                                matched.list,
                                range.info,
                                matrix) {

    # helper that generates data frame of the correct size
    # the helper method is used for both the get_region_counts
    # and the get_metagene functions,
    range.lower <- range.info['range.lower']
    range.upper <- range.info['range.upper']

    transcript  <- conditions[["transcript"]]
    length      <- conditions[["length"]]

    ref.length  <- length(ref.names)
    num.reads   <- range.upper - range.lower + 1
    total.list  <- length(matched.list)
    
    # get the matrix data and pass create a data frame of the 
    # correct format
    if (transcript & length) {
        experiment.list   <- matched.list
        return (data.frame(experiment = matched.list,
                           matrix,
                           stringsAsFactors = FALSE,
                           check.names = FALSE))
    } else if (transcript) {
        #sum transcripts only
        experiment.column <- rep(matched.list, each = num.reads)
        read.column <- rep(c(range.lower:range.upper), total.list)
        return (data.frame(experiment = experiment.column,
                           length     = read.column,
                           matrix, stringsAsFactors = FALSE,
                           check.names = FALSE))
    } else if (length) {
        #length only
        experiment.column <- rep(matched.list, each = ref.length)
        transcript.column <- rep(ref.names, total.list)
        return (data.frame(experiment = experiment.column,
                           transcript = transcript.column,
                           matrix, stringsAsFactors = FALSE,
                           check.names = FALSE))
    }
    #!transcript and !length
    experiment.column <- rep(matched.list, each = num.reads * ref.length)
    transcripts       <- rep(ref.names, total.list * num.reads)
    ref.read          <- rep(c(range.lower:range.upper), each = ref.length)
    read.column <- rep(ref.read, total.list)
    return (data.frame(experiment  = experiment.column,
                       transcript  = transcripts,
                       length      = read.column,
                       matrix, stringsAsFactors = FALSE,
                       check.names = FALSE))
}


generate_matrix <- function(ribo.object,
                            transcript,
                            length,
                            normalize,
                            ncol,
                            file,
                            path,
                            index){
    #helper method that generates the matrix under different circumstances
    experiment <- strsplit(path, split="/")[[1]][3]
    output <- t(h5read(file=file, index=index, name=path))

    ref.length <- length(get_reference_names(ribo.object))
    single <- (ref.length == nrow(output))

    if (transcript & length) {
        output <- matrix(unlist(colSums(output)), ncol=ncol, byrow=TRUE)
    } else if (transcript) {
        # condense all transcripts together
        condense_transcripts <- seq(from = 1 ,
                                    to = nrow(output),
                                    by = ref.length)
        
        output <- lapply(condense_transcripts, sum_lengths, 
                                               ref.length = ref.length,
                                               mat        = output)

        output <- matrix(unlist(output), ncol=ncol, byrow=TRUE)
   } else if (length) {
        #only sum across lengths if there is more than one length 
        if (!single) {
            val <- seq(ref.length)
            condense_lengths <- lapply(val,
                                       seq,
                                       to = nrow(output),
                                       by = ref.length)
            output <- lapply(condense_lengths, sum_transcripts, mat=output)
            output <- matrix(unlist(output), ncol=ncol, byrow=TRUE)
      }
  }
  return(output)
}

sum_transcripts <- function (index, mat) {
  #in the case that there is only one region
  if (ncol(mat) == 1) {
    return(sum(mat[index, ]))
  }
  return (colSums(mat[index, ]))
}

sum_lengths <- function(index, ref.length, mat) {
  #in the case that there is only one region 
  if (ncol(mat) == 1) {
    return(sum(mat[index:(index+ref.length-1), ]))
  }
  return (colSums(mat[index:(index+ref.length-1), ]))
}


prepare_DataFrame <- function(ribo.object, DF) {
  #factor and Rle the columns of the metagene and region_count functions
  if (!is.null(DF$region)) DF$region <- Rle(factor(DF$region))
  if (!is.null(DF$transcript)) DF$transcript <- factor(DF$transcript)
  if (!is.null(DF$position)) DF$position <- Rle(DF$position)
  if (!is.null(DF$length)) DF$length <- Rle(DF$length)
  DF$experiment = Rle(factor(DF$experiment))
  
  # add the metadata
  DF@metadata[[1]] <- get_info(ribo.object)$experiment.info[, c("experiment", 
                                                                "total.reads")]
  return(DF)
}

strip_rlefactor <- function(DF) {
  if (!is.null(DF$region)) DF$region <- as.character(DF$region)
  if (!is.null(DF$transcript)) DF$transcript <- as.character(DF$transcript)
  if (!is.null(DF$position)) DF$position <- as.integer(DF$position)
  if (!is.null(DF$length)) DF$length <- as.integer(DF$length)
  DF$experiment <- as.character(DF$experiment)
  return(DF)
}
