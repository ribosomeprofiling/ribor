#' ribor: A package for reading .ribo files 
#'
#' The 'ribor' package offers a suite of reading functions for the datasets
#' present in a .ribo file and also provides some rudimentary plotting 
#' functions.
#' 
#' @section Vignette:
#' To get started with the ribor package, please see the vignette page at
#' \url{https://ribosomeprofiling.github.io/ribor/walkthrough.html}.
#' 
#' @section Package Content:
#' \subsection{Generating a ribo object}{
#'  \code{\link{create_ribo}} to get started 
#' }
#' 
#' \subsection{Length Distribution}{
#'  \code{\link{get_length_distribution}} to get length distribution counts
#'  
#'  \code{\link{plot_length_distribution}} to plot the length distribution
#' }
#' 
#' \subsection{Region Counts}{
#'   \code{\link{get_region_counts}} to get region counts
#'    
#'   \code{\link{plot_region_counts}} to plot the region counts 
#' }
#' 
#' \subsection{Metagene Coverage}{
#'   \code{\link{get_metagene}} to get metagene site coverage
#'   
#'   \code{\link{get_tidy_metagene}} to get a tidy format of the metagene site coverage
#'   
#'   \code{\link{plot_metagene}} to plot the metagene site coverage
#' }
#' 
#' @docType package
#' @name ribor
#' @importFrom methods show setClass setGeneric setMethod is validObject new
NULL