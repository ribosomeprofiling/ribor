#' Retrieves the region counts from a .ribo file
#'
#' \code{\link{get_region_counts}} will return the particular region counts
#' of any subset of regions for a given set of experiments.
#'
#' This function will return a DataFrame of the counts at each specified region
#' for each specified experiment. The region options are "UTR5", "UTR5J", "CDS",
#' "UTR3J", and "UTR3". The user can specify any subset of regions in the form of a vector,
#' a list, or a single string if only one region is desired.
#'
#' The dimensions of the returned DataFrame depend on the parameters
#' range.lower, range.upper, length, and transcript.
#'
#' The param 'length' condenses the read lengths together.
#' When length is TRUE and transcript is FALSE, the
#' DataFrame presents information for each transcript across
#' all of the read lengths. That is, each transcript has a value
#' that is the sum of all of the counts across every read length.
#' As a result, information about the transcript at each specific
#' read length is lost.
#'
#' The param 'transcript' condenses the transcripts together.
#' When transcript is TRUE and length is FALSE data
#' table presents information at each read length between range.lower and
#' range.upper inclusive. That is, each separate read length denotes the
#' sum of counts from every transcript. As a result, information about the
#' counts of each individual transcript is lost.
#'
#' When 'transcript' is set to FALSE, the 'alias' parameter specifies whether
#' or not the returned DataFrame should present each transcript as an alias
#' instead of the original name. If 'alias' is set to TRUE, then the column
#' of the transcript names will contain the aliases rather than the original
#' reference names of the .ribo file.
#'
#' If both 'length' and 'transcript' are TRUE, then the resulting
#' DataFrame prints out one row for each experiment. This provides the metagene
#' information across all transcripts and all reads in a given experiment.
#'
#' If both length' and 'transcript' are FALSE, calculations are done to the data,
#' all information is preserved for both the read length and the transcript.
#' The DataFrame would just present the entire stored raw data
#' from the read length 'range.lower' to the read length 'range.upper' which in most
#' cases would result in a slow run time with a massive DataFrame returned.
#'
#' When 'transcript' is set to FALSE, the 'alias' parameter specifies whether
#' or not the returned DataFrame should present each transcript as an alias
#' instead of the original name. If 'alias' is set to TRUE, then the column
#' of the transcript names will contain the aliases rather than the original
#' reference names of the .ribo file.
#' @examples
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- create_ribo(file.path)
#'
#' #specify the regions and experiments of interest
#' regions <- c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3")
#' experiments <- c("Hela_1", "Hela_2", "WT_1")
#'
#' #obtains the region counts at each individual read length, summed across every transcript
#' region.counts <- get_region_counts(ribo.object = sample,
#'                                    region = regions,
#'                                    range.lower = 2,
#'                                    range.upper = 5,
#'                                    length = FALSE,
#'                                    transcript = TRUE,
#'                                    tidy = FALSE,
#'                                    alias = FALSE,
#'                                    experiments = experiments)
#'
#' @param ribo.object A 'ribo' object
#' @param range.lower Lower bound of the read length
#' @param range.upper Upper bound of the read length
#' @param length Option to condense the read lengths together, preserve the transcripts
#' @param transcript Option to condense the transcripts together, preserve the read lengths
#' @param tidy Option to return the DataFrame in a tidy format
#' @param region Specific region of interest
#' @param alias Option to report the transcripts as aliases/nicknames
#' @param experiments List of experiment names
#' @param normalize Option to normalize the counts as counts per million reads
#' @param compact Option to return a DataFrame with Rle and factor as opposed to a raw data.frame
#' @return A DataFrame of the region counts
#' @importFrom rhdf5 h5read
#' @importFrom methods as 
#' @importFrom S4Vectors DataFrame Rle
#' @importFrom tidyr gather
#' @importFrom dplyr left_join
#' @export
get_region_counts <- function(ribo.object,
                              range.lower = rangeLower(ribo.object),
                              range.upper = rangeUpper(ribo.object),
                              length = TRUE,
                              transcript = TRUE,
                              tidy = TRUE,
                              alias = FALSE,
                              normalize = FALSE,
                              region = c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3"),
                              compact = TRUE, 
                              experiments = get_experiments(ribo.object)) {
    range.info <- c(range.lower = range.lower, range.upper = range.upper)
    region     <- check_rc_input(ribo.object, region, range.info, experiments, alias)
    conditions <- c(transcript = transcript, length = length,
                    normalize  = normalize,   alias = alias)

    path              <- ribo.object@path
    matched.experiments <- intersect(experiments, get_experiments(ribo.object))
    total.experiments   <- length(matched.experiments)
    ref.names           <- get_reference_names(ribo.object)
    ref.length          <- length(ref.names)

    values    <- c("UTR5"  = 1, "UTR5J" = 2, "CDS" = 3, "UTR3J" = 4, "UTR3" = 5)
    columns   <- unname(values[region])
    ncol      <- length(region)
    range.min <- get_read_lengths(ribo.object)[1]
    row.start <- (range.lower - range.min) * ref.length + 1
    row.stop  <- row.start + ref.length* (range.upper - range.lower + 1) - 1
    rows      <- c(row.start:row.stop)
    matched.experiments <- intersect(experiments, get_experiments(ribo.object))
    
    # compute the file paths
    file_paths <- vapply(matched.experiments, get_rc_path, FUN.VALUE = "character")
    data  <- lapply(file_paths, generate_matrix, ribo.object = ribo.object,
                                                 transcript  = transcript,
                                                 length      = length,
                                                 normalize   = normalize,
                                                 file        = path,
                                                 index       = list(columns, rows),
                                                 ncol        = ncol)
    data           <- do.call(rbind, data)
    colnames(data) <- region
    
    # make_dataframe is a generic helper that returns a DataFrame so need to cast
    result         <- as.data.frame(make_dataframe(ribo.object, matched.experiments,
                                                   range.info, conditions, data))
    
    if (normalize) {
        # normalize the counts to reads/total reads * 1000000
        info <- get_info(ribo.object)$experiment.info[, c("experiment", "total.reads")]
        result <- left_join(result, info, by = "experiment")
        result[, region] <- result[, region] * 1000000/result$total.reads
        result[, "total.reads"] <- NULL
    }
    
    if (tidy) result <- gather(result, "region", "count", region)
    
    if (compact) {
      result <- as(result, "DataFrame")
      return(prepare_DataFrame(ribo.object, result))
    } 
    return (result)
}


get_rc_path <- function(experiment) {
    return(paste("/experiments/",
                 experiment,
                 "/region_counts/region_counts",
                 sep = ""))
}

check_rc_input <- function(ribo.object,
                           region,
                           range.info,
                           experiments,
                           alias) {
    #helper function that checks for valid parameters given by the user
    #calls error messages on any incorrect parameters
    region <- check_regions(ribo.object, region)

    range.lower <- range.info[["range.lower"]]
    range.upper <- range.info[["range.upper"]]

    check_alias(ribo.object, alias)
    check_lengths(ribo.object, range.lower, range.upper)
    check_experiments(ribo.object, experiments)
    return(region)
}

check_regions <- function(ribo.object,
                          region) {
    region <- toupper(region)
    region.options <- c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3")

    if (typeof(region) != 'list' & typeof(region) != "character") {
        stop("Please specify the regions as a single string, a vector, or a list.",
             call. = FALSE)
    }
    
    if (length(which(!region %in% region.options))) {
        stop("Please indicate the region(s) with ", 
             "the following: 'UTR5', 'UTR5J, 'CDS', 'UTR3J', 'UTR3'",
             call. = FALSE)
    }
    
    region <- factor(region, 
                     levels  = c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3"), 
                     ordered = TRUE)
    region <- as.vector(sort(region))
    return(region)
}


#' Retrieves the length distribution of a given region
#'
#' The function {\code{\link{get_length_distribution}}} retrieves the raw or normalized
#' counts at each read length from 'range.lower' to 'range.upper'.
#'
#' This function is a wrapper function of {\code{\link{get_region_counts}}}, and the
#' returned DataFrame is valid input for {\code{\link{plot_length_distribution}}}.
#'
#' @examples
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- create_ribo(file.path)
#'
#' #specify the experiments of interest
#' experiments <- c("Hela_1", "Hela_2", "WT_1")
#'
#' #gets the normalized length distribution from read length 2 to 5
#' length.dist <- get_length_distribution(ribo.object = sample,
#'                                        region = "CDS",
#'                                        range.lower = 2,
#'                                        range.upper = 5)
#'
#' @inheritParams get_region_counts
#' @importFrom methods as
#' @importFrom S4Vectors DataFrame Rle
#' @export
#' @seealso
#' {\code{\link{plot_length_distribution}}} to plot the output of this function
#' @return
#' A DataFrame of the counts at each read length
get_length_distribution <- function(ribo.object,
                                    region,
                                    range.lower = rangeLower(ribo.object),
                                    range.upper = rangeUpper(ribo.object),
                                    compact = TRUE,
                                    experiments = get_experiments(ribo.object)) {
  if (length(region) != 1) {
    stop("Please provide only one region.")
  }
  
  result <- get_region_counts(ribo.object,
                              range.lower = range.lower,
                              range.upper = range.upper,
                              region = region,
                              length = FALSE,
                              transcript = TRUE,
                              normalize = FALSE,
                              compact = compact,
                              experiments = experiments)[, -3]
  
  # by not using DataFrame, we will need to add total.reads back
  if (!compact) {
    result %>% 
      left_join(get_info(ribo.object)$experiment.info[, c("experiment", "total.reads")],
                by = "experiment") -> result
  }
  return (result)
}


#' Plots the length distribution
#'
#' The function \code{\link{plot_length_distribution}} can take either a DataFrame
#' or a "ribo" object to generate a line graph of the length distributions from
#' range.lower to range.upper.
#'
#' The param 'fraction' will plot the fractions of each length relative
#' to the total sum of the read length range provided by param 'range.lower'
#' and 'range.upper'. When fraction is set to FALSE, the total count of each
#' read length is plotted.
#'
#' When given a "ribo" object, \code{\link{plot_length_distribution}} calls
#' \code{\link{get_region_counts}} to retrieve the necessary information
#' for plotting.
#'
#' The user can instead provide a DataFrame with the same structure as the
#' output of the \code{\link{get_region_counts}} function where the 'transcript'
#' parameter is set to FALSE and 'length' parameters is the default value of
#' TRUE. This also means that the many of the remaining parameters of the
#' \code{\link{plot_length_distribution}} function are not necessary. The run
#' time becomes substantially faster when \code{\link{plot_region_counts}} is
#' given the direct DataFrame to plot. Note that there is no manipulation by
#' this function on the DataFrame. This responsibility is given to the user
#' and allows for more control.
#'
#' @examples
#' #ribo object use case
#'
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- create_ribo(file.path)
#'
#' #specify experiments of interest
#' experiments <- c("Hela_1", "Hela_2", "WT_1")
#'
#' plot_length_distribution(x = sample,
#'                          region = "CDS",
#'                          range.lower = 2,
#'                          range.upper = 5,
#'                          experiments = experiments,
#'                          fraction = TRUE)
#'
#'
#' #DataFrame use case
#' #obtains the region counts at each individual read length, summed across every transcript
#' region.counts <- get_length_distribution(ribo.object = sample,
#'                                          region      = "CDS",
#'                                          range.lower = 2,
#'                                          range.upper = 5,
#'                                          experiments = experiments)
#'
#' #the param 'length' must be set to FALSE and param 'transcript' must be set
#' #to TRUE to use a DataFrame
#' plot_length_distribution(region.counts)
#'
#'
#' @seealso \code{\link{get_region_counts}} to generate a DataFrame that can
#' be provided as input,
#' \code{\link{ribo}} to create a ribo.object that can be provided as input
#' @param x A 'ribo' object or a DataFrame generated from \code{\link{get_region_counts}}
#' @param region the region of interest
#' @param range.lower a lower bounds for a read length range
#' @param range.upper an upper bounds for a read length range
#' @param experiments a list of experiment names
#' @param fraction logical value that, if TRUE, presents the count as a fraction of the total reads in the given ranges
#' @param title a title for the generated plot
#' @importFrom tidyr gather
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot geom_line theme_bw theme labs aes_string
#' @importFrom stats aggregate
#' @importFrom rlang .data
#' @export
#' @return A 'ggplot' of the length distribution
plot_length_distribution <- function(x,
                                     region,
                                     experiments,
                                     range.lower,
                                     range.upper,
                                     fraction = FALSE,
                                     title = "Length Distribution") {
    x <- check_ld_input(x, region, range.lower, range.upper, experiments)
    y.axis <- "Count"
    y.value <- "count"
    
    if (fraction) {
        if (is(x, "DataFrame")) {
          info <- x@metadata[[1]]
          x %>% 
            strip_rlefactor() %>% 
            as.data.frame() %>% 
            left_join(info, by = "experiment") -> x
        } else {
          x %>% 
            left_join(get_info(x)$experiment.info[, c("experiment", "total.reads")],
                      by = "experiment") 
        }
        
        x %>% 
          mutate(fraction = .data$count/.data$total.reads) -> x
        y.axis <- "Fraction"
        y.value <- "fraction"
    } else if (is(x, "DataFrame")) {
        x <- as.data.frame(x)
    }
    ggplot(x, aes_string(x="length", y=y.value, color="experiment")) +
        geom_line() +
        theme_bw() +
        theme(plot.title=element_text(hjust=0.5)) +
        labs(title=title, x="Read Length", y=y.axis, color="Experiment")
}

check_ld_input <- function(x,
                           region,
                           range.lower,
                           range.upper,
                           experiments) {
      #check the plot_length_distribution output
    is.ribo <- check_ribo(x, stop=FALSE)
    if (is.ribo) {
        if (missing(region) || length(region) != 1){
            stop("Please indicate a single region.") 
        } 
        if (missing(experiments)) experiments = get_experiments(x)
        if (missing(range.lower)) range.lower <- rangeLower(x)
        if (missing(range.upper)) range.upper <- rangeUpper(x)
        
        x <- strip_rlefactor(get_length_distribution(ribo.object = x,
                                                     region      = region,
                                                     range.lower = range.lower,
                                                     range.upper = range.upper,
                                                     experiments = experiments))
    } else if (is(x, "DataFrame") || is(x, "DFrame")) {

        x <- strip_rlefactor(x)
        col.names <- c("experiment", "length", "count")
        types <- c("integer", "double")
        mismatch <- !all(names(x) == col.names, 
                         typeof(x[, "experiment"]) == "character",
                         typeof(x[, "length"]) %in% types, 
                         typeof(x[, "count"]) %in% types,
                         ncol(x) == 3,
                         length(x@metadata) > 0)
        if (mismatch) {
              stop("Please make sure that the DataFrame is of the correct format.
                    It requires a non-empty metadata field.",
                    call.=FALSE)
        }
    } else if (is.data.frame(x)) {
        col.names <- c("experiment", "length", "count", "total.reads")
        types <- c("integer", "double")
        mismatch <-  !all(names(x) == col.names,                    
                          typeof(x[, "experiment"]) == "character",
                          typeof(x[, "length"]) %in% types,
                          typeof(x[, "count"]) %in% types,    
                          typeof(x[, "total.reads"]) %in% types, 
                          ncol(x) == 4)
        if (mismatch) {
          stop("Please make sure that the data frame is of the correct format.",
               " It requires a 'total.reads' column.",
               call.=FALSE)
        }
    } else {
        stop("Please make sure that param 'x' is either", 
             "a DataFrame or a ribo object.",
             call.=FALSE)
    }
    return(x)
}


#' Plots the region counts of UTR5, CDS, and UTR3
#'
#' The function \code{\link{plot_region_counts}} can take either a DataFrame
#' or a "ribo" object to generate the a stacked bar plot of proportions that
#' correspond to the "UTR5", "CDS", and "UTR3" regions.
#'
#' When given a 'ribo' object, \code{\link{plot_region_counts}} calls
#' \code{\link{get_region_counts}} to retrieve the necessary information
#' for plotting. This option is in the case that a DataFrame of the
#' region count information is not required.
#'
#' The user can instead provide a DataFrame with the same structure as the
#' output of the \code{\link{get_region_counts}} function where the 'transcript'
#' and 'length' parameters are the default values of TRUE. This also means that
#' the remaining parameters of the \code{\link{plot_region_counts}} function are not necessary.
#' The run time becomes substantially faster when \code{\link{plot_region_counts}} is given
#' the direct DataFrame to plot. However, the DataFrame needs to follow the format and 
#' types in the output of the reading functions 
#'
#' @examples
#' #ribo object use case
#' #generate the ribo object
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- create_ribo(file.path)
#'
#' #specify the regions and experiments of interest
#' regions <- c("UTR5", "CDS", "UTR3")
#' experiments <- c("Hela_1", "Hela_2", "WT_1")
#'
#' plot_region_counts(sample,
#'                    range.lower = 2,
#'                    range.upper = 5,
#'                    experiments)
#'
#' #DataFrame use case
#' #obtains the region counts at each individual read length, summed across every transcript
#' region.counts <- get_region_counts(sample,
#'                                    region = regions,
#'                                    range.lower = 2,
#'                                    range.upper = 5,
#'                                    tidy = TRUE,
#'                                    length = TRUE,
#'                                    transcript = TRUE)
#'
#' #the params 'length' and 'transcript' must be set to true to use a DataFrame
#' plot_region_counts(region.counts)
#'
#' @seealso \code{\link{get_region_counts}} to generate a DataFrame that can be provided as input,
#' \code{\link{ribo}} to create a ribo.object that can be provided as input
#'
#' @param x A 'ribo' object or a DataFrame generated from \code{\link{get_region_counts}}
#' @param range.lower a lower bounds for a read length range
#' @param range.upper an upper bounds for a read length range
#' @param experiments a list of experiment names
#' @param title a title for the generated plot
#' @importFrom dplyr left_join mutate %>% group_by summarize arrange desc
#' @importFrom tidyr gather
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_col theme_bw theme ggtitle coord_flip theme
#' @importFrom ggplot2 labs scale_fill_discrete element_blank geom_text position_stack
#' @export
#' @return A 'ggplot' of the region counts for each of the experiments 
plot_region_counts <- function(x,
                               experiments,
                               range.lower,
                               range.upper,
                               title = "Region Counts") {
    rc <- check_plot_rc_input(x, range.lower, range.upper, experiments)
    rc <- as.data.frame(rc)
    
    #prepare data for visualization
    rc %>%
        group_by(.data$experiment) %>%
        summarize(sum=sum(.data$count)) -> rc.total

    rc.total %>%
        left_join(rc, by = "experiment") -> rc

    rc %>%
        mutate(percentage=round(100 * .data$count/sum, 1)) %>%
        mutate(region=factor(.data$region, 
               levels=c("UTR3", "CDS", "UTR5"))) -> rc

    #text label of percentages only in the "CDS" region
    percentages <- replace(rc$percentage, 
                           rc$region != "CDS", "")

    ggplot(rc, 
           aes(x=.data$experiment, y=.data$percentage, fill=.data$region)) +
        geom_col() +
        coord_flip() +
        geom_text(aes(x=.data$experiment, y=50, label=percentages), size=3) +
        theme_bw() +
        theme(plot.title   = element_text(hjust = 0.5),
              panel.border = element_blank(),
              panel.grid   = element_blank()) +
        scale_fill_discrete(breaks = c("UTR5", "CDS", "UTR3")) +
        labs(title=title, fill="Region", x="Experiment", y="Percentage")
}


check_plot_rc_input <- function(x,
                                range.lower,
                                range.upper,
                                experiments) {
    #helper method that checks the plot_region_counts input and returns data
    #for plotting
    is.ribo <- check_ribo(x, stop = FALSE)
    regions <- c("UTR5", "CDS", "UTR3")
    if (is.ribo) {
        if (missing(experiments) || is.null(experiments)) {
            experiments <- get_experiments(x)
        } 
        if (missing(range.lower)) range.lower <- rangeLower(x)
        if (missing(range.upper)) range.upper <- rangeUpper(x)
        check_lengths(x, range.lower, range.upper)
        x <- strip_rlefactor(get_region_counts(ribo.object = x,
                                                region      = regions,
                                                range.lower = range.lower,
                                                range.upper = range.upper,
                                                length      = TRUE,
                                                transcript  = TRUE,
                                                experiments = experiments))
    } else if (is(x, "DataFrame") || is.data.frame(x)) {
        if (is(x, "DataFrame")) x <- strip_rlefactor(x)

        col.names <- c("experiment", "region", "count")
        mismatch  <- !all(names(x) == col.names, 
                          typeof(x[, "experiment"]) == "character",
                          typeof(x[, "region"]) == "character",
                          typeof(x[, "count"]) %in% c("double", "integer"),
                          ncol(x) == 3)
        if (mismatch) {
            stop("Please make sure that the DataFrame is of ",  
                 "the correct format.", call. = FALSE)
        } else if (!identical(unique(x$region), regions)) {
            stop("Please make sure that the DataFrame only includes the ", 
                 "'UTR5','CDS', and 'UTR3' regions.", call. = FALSE)
        }
    } else {
        stop("Please provide a ribo object or a DataFrame of the",
             "correct format.", call. = FALSE)
    }
    return(x)
}