#' Returns the overall length of each region with UTR Junctions 
#' 
#' The function {\code{\link{get_internal_region_coordinates}}} retrieves the 
#' lengths for the UTR5, UTR5 Junction, CDS, UTR3 Junction,
#' and UTR3 regions of every transcript. 
#' 
#' @examples 
#' # generate a ribo object 
#' file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")
#' sample <- Ribo(file.path, rename = rename_default)
#' 
#' # get the region lengths
#' region_lengths <- get_internal_region_lengths(sample, alias = TRUE)
#' 
#' @export
#' @importFrom rhdf5 h5read
#' @param ribo.object A 'Ribo' object 
#' @param alias Option to return the transcript names as aliases  
#' @return A data.frame of the region lengths 
get_internal_region_lengths <- function(ribo.object, alias = FALSE) {
    validObject(ribo.object)
    check_alias(ribo.object, alias)
    
    #generate the start and stop indices, get the matrix of positions
    start <- seq(1, 10, by = 2)
    stop  <- seq(2, 10, by = 2)
    result <- compute_boundaries(ribo.object)
    
    #in the case of NA, report the length as 0 
    result[, c(start)][is.na(result[, c(start)])] <- 0
    result[, c(stop)][is.na(result[, c(stop)])] <- -1
    
    #generate the lengths 
    diff <- lapply(start, function(c, result) result[, c+1] - result[, c] + 1, result)
    diff <- matrix(unlist(diff), byrow = FALSE, nrow = nrow(result))
    colnames(diff) <- c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3")
    
    #generate the data.frame 
    references <- change_reference_names(ribo.object, alias)
    return(data.frame(transcript=references, diff))
}


#' Returns the overall length of each region with UTR Junctions 
#' 
#' The function {\code{\link{get_region_lengths}}} retrieves the 
#' lengths for the UTR5, UTR5 Junction, CDS, UTR3 Junction,
#' and UTR3 regions of every transcript. 
#' 
#' This function is deprecated, and we recommend {\code{\link{get_internal_region_lengths}}}.
#' 
#' @export
#' @importFrom rhdf5 h5read
#' @param ribo.object A 'Ribo' object 
#' @param alias Option to return the transcript names as aliases  
#' @return A data.frame of the region lengths 
get_region_lengths <- function(ribo.object, alias = FALSE) {
    .Deprecated("get_internal_region_lengths")
    return(get_internal_region_lengths(ribo.object = ribo.object, alias = alias))
}



#' Returns the overall length of each region
#' 
#' The function {\code{\link{get_original_region_coordinates}}} retrieves the 
#' lengths for the UTR5, CDS, and UTR3 regions of every transcript. 
#' 
#' @examples 
#' # generate a ribo object 
#' file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")
#' sample <- Ribo(file.path, rename = rename_default)
#' 
#' # get the region coordinates
#' region_lengths <- get_original_region_lengths(sample, alias = TRUE)
#' 
#' @export
#' @importFrom rhdf5 h5read
#' @param ribo.object A 'Ribo' object 
#' @param alias Option to return the transcript names as aliases  
#' @return A data.frame of the region lengths 
get_original_region_lengths <- function(ribo.object, alias = FALSE) {
    validObject(ribo.object)
    check_alias(ribo.object, alias)
    
    #generate the start and stop indices, get the matrix of positions
    start <- seq(1, 6, by = 2)
    stop  <- seq(2, 6, by = 2)
    result <- compute_original_boundaries(ribo.object)
    
    #in the case of NA, report the length as 0 
    result[, c(start)][is.na(result[, c(start)])] <- 0
    result[, c(stop)][is.na(result[, c(stop)])] <- -1
    
    #generate the lengths based on the differences in coordinate position 
    diff <- lapply(start, function(c, result) result[, c+1] - result[, c] + 1, result)
    diff <- matrix(unlist(diff), byrow = FALSE, nrow = nrow(result))
    colnames(diff) <- c("UTR5", "CDS", "UTR3")
    
    #generate the data.frame 
    references <- get_reference_names(ribo.object)
    return(data.frame(transcript=references, diff))
}


#' Retrieves the region stop and start coordinates 
#' 
#' The function {\code{\link{get_internal_region_coordinates}}} retrieves the 
#' start and site positions for the UTR5, UTR5 Junction, CDS, UTR3 Junction,
#' and UTR3 regions of every transcript. 
#' 
#' To note, because of the R-specific 1-based indexing, the positions start at
#' 1 instead of 0 in other programming languages. The positions provided in 
#' the returned data.frame will correspond to the positions in the output of 
#' {\code{\link{get_coverage}}}.
#' 
#' Additionally, within the transcripts, there are edge cases. 
#' NA values found in the returned data.frame means that the region has no 
#' start and stop position and a length of zero after computing the boundaries
#' of the UTR5 and UTR3 junction.
#' 
#' @export 
#' @examples 
#' # generate a ribo object 
#' file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")
#' sample <- Ribo(file.path, rename = rename_default)
#' 
#' # get the region coordinates
#' coord <- get_internal_region_coordinates(sample, alias = TRUE)
#' 
#' @inheritParams get_region_lengths
#' @return A data.frame of start and stop coordinates for every region
get_internal_region_coordinates <- function(ribo.object, alias=FALSE) {
    validObject(ribo.object)
    check_alias(ribo.object, alias)
    references <- change_reference_names(ribo.object, alias)
    #generate the data.frame from the boundary positions
    return(data.frame(transcript=references, compute_boundaries(ribo.object)))
}


#' Retrieves the region stop and start coordinates 
#' 
#' The function {\code{\link{get_region_coordinates}}} retrieves the 
#' start and site positions for the UTR5, UTR5 Junction, CDS, UTR3 Junction,
#' and UTR3 regions of every transcript. 
#' 
#' To note, because of the R-specific 1-based indexing, the positions start at
#' 1 instead of 0 in other programming languages. The positions provided in 
#' the returned data.frame will correspond to the positions in the output of 
#' {\code{\link{get_coverage}}}.
#' 
#' Additionally, within the transcripts, there are edge cases. 
#' NA values found in the returned data.frame means that the region has no 
#' start and stop position and a length of zero after computing the boundaries
#' of the UTR5 and UTR3 junction.
#' 
#' @export 
#' 
#' @inheritParams get_region_lengths
#' @return A data.frame of start and stop coordinates for every region
get_region_coordinates <- function(ribo.object, alias=FALSE) {
    .Deprecated("get_internal_region_coordinates")
    return(get_internal_region_coordinates(ribo.object = ribo.object, alias = alias))
}




#' Retrieves the region stop and start coordinates 
#' 
#' The function {\code{\link{get_original_region_coordinates}}} retrieves the 
#' start and site positions for the UTR5, UTR5 Junction, CDS, UTR3 Junction,
#' and UTR3 regions of every transcript. 
#' 
#' To note, because of the R-specific 1-based indexing, the positions start at
#' 1 instead of 0 in other programming languages. The positions provided in 
#' the returned data.frame will correspond to the positions in the output of 
#' {\code{\link{get_coverage}}}.
#' 
#' Additionally, within the transcripts, there are edge cases. 
#' NA values found in the returned data.frame means that the region has no 
#' start and stop position and a length of zero after computing the boundaries
#' of the UTR5 and UTR3 junction.
#' 
#' @export 
#' @examples 
#' # generate a ribo object 
#' file.path <- system.file("extdata", "HEK293_ingolia.ribo", package = "ribor")
#' sample <- Ribo(file.path, rename = rename_default)
#' 
#' # get the region coordinates
#' coord <- get_original_region_coordinates(sample, alias = TRUE)
#' 
#' @inheritParams get_original_region_lengths
#' @return A data.frame of start and stop coordinates for every region
get_original_region_coordinates <- function(ribo.object, alias=FALSE) {
    validObject(ribo.object)
    check_alias(ribo.object, alias)
    references <- change_reference_names(ribo.object, alias)
    #generate the data.frame from the boundary positions
    return(data.frame(transcript=references, compute_original_boundaries(ribo.object)))
}




compute_boundaries <- function(ribo.object) {
    #helper method that computes the boundaries for each region
    annotation <- t(h5read(path(ribo.object),
                           name = "/reference/annotation"))
    left.span <- left_span(ribo.object)
    right.span <- right_span(ribo.object)
    UTR5_end <- annotation[, 1]
    CDS_end <-  annotation[, 2]
    UTR3_end <- annotation[, 3]
    
    #compute the junction start and stop
    UTR5J_start <- UTR5_end + 1 - left.span 
    UTR5J_start[UTR5J_start <= 0] <- 1 
    UTR5J_stop <- UTR5_end + right.span
    UTR5J_stop[UTR5J_stop > CDS_end] <- CDS_end[UTR5J_stop > CDS_end]
    UTR3J_start <- CDS_end + 1 - left.span
    UTR3J_start[UTR3J_start <= UTR5J_stop] <- 
        UTR5J_stop[UTR3J_start <= UTR5J_stop] + 1
    UTR3J_stop <- CDS_end + right.span
    UTR3J_stop[UTR3J_stop > UTR3_end] <- UTR3_end[UTR3J_stop > UTR3_end]
    UTR3J_start[UTR3J_start > UTR3J_stop] <- NA
    UTR3J_stop[UTR3J_start > UTR3J_stop] <- NA
    
    #compute the CDS start and stop
    CDS_start <- UTR5J_stop + 1 
    CDS_stop <- UTR3J_start - 1
    CDS_start[CDS_start > CDS_stop] <- NA
    CDS_stop[CDS_start > CDS_stop] <- NA
    
    #compute the UTR5 start and stop
    UTR5_start  <- rep(1, times = length(get_reference_names(ribo.object)))
    UTR5_stop   <- UTR5J_start - 1
    UTR5_start[UTR5_stop == 0] <- NA
    UTR5_stop[UTR5_stop == 0] <- NA
    
    #compute the UTR3 start and stop
    UTR3_start <- UTR3J_stop + 1 
    UTR3_stop  <- UTR3_end
    UTR3_start[UTR3_start > UTR3_stop] <- NA
    UTR3_stop[UTR3_start > UTR3_stop] <- NA
    info <- list(UTR5_start, UTR5_stop, UTR5J_start, UTR5J_stop,
                 CDS_start, CDS_stop, UTR3J_start, UTR3J_stop,
                 UTR3_start, UTR3_stop)
    result <- do.call(cbind, info)
    colnames(result) <- c("UTR5_start", "UTR5_stop", "UTR5J_start", 
                          "UTR5J_stop", "CDS_start", "CDS_stop", "UTR3J_start",
                          "UTR3J_stop", "UTR3_start", "UTR3_stop")
    return(result)
}



compute_original_boundaries <- function(ribo.object) {
    #helper method that computes the boundaries for each region
    annotation <- t(h5read(path(ribo.object),
                           name = "/reference/annotation"))
    UTR5_end <- annotation[, 1]
    CDS_end <-  annotation[, 2]
    UTR3_end <- annotation[, 3]
    
    #compute the CDS start and stop
    CDS_start <- UTR5_end + 1 
    CDS_stop <- CDS_end
    CDS_start[CDS_start > CDS_stop] <- NA
    CDS_stop[CDS_start > CDS_stop] <- NA
    
    #compute the UTR5 start and stop
    UTR5_start  <- rep(1, times = length(get_reference_names(ribo.object)))
    UTR5_stop   <- UTR5_end
    UTR5_start[UTR5_stop == 0] <- NA
    UTR5_stop[UTR5_stop == 0] <- NA
    
    #compute the UTR3 start and stop
    UTR3_start <- CDS_end + 1 
    UTR3_stop  <- UTR3_end
    UTR3_start[UTR3_start > UTR3_stop] <- NA
    UTR3_stop[UTR3_start > UTR3_stop] <- NA
    info <- list(UTR5_start, UTR5_stop,
                 CDS_start, CDS_stop, 
                 UTR3_start, UTR3_stop)
    result <- do.call(cbind, info)
    colnames(result) <- c("UTR5_start", "UTR5_stop",
                          "CDS_start", "CDS_stop",
                          "UTR3_start", "UTR3_stop")
    return(result)
}
