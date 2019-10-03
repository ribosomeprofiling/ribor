#' Retrieves the region counts from a .ribo file
#'
#' \code{\link{get_region_counts}} will return the particular region counts
#' of any subset of regions for a given set of experiments.
#'
#' This function will return a data.table of the counts at each specified region
#' for each specified experiment. The region options are "UTR5", "UTR5J", "CDS",
#' "UTR3J", and "UTR3". The user can specify any subset of regions in the form of a vector,
#' a list, or a single string if only one region is desired.
#'
#' The dimensions of the returned data table depend on the parameters
#' range.lower, range.upper, length, and transcript.
#'
#' The param 'length' condenses the read lengths together.
#' When length is TRUE and transcript is FALSE, the
#' data table presents information for each transcript across
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
#' or not the returned data.table should present each transcript as an alias
#' instead of the original name. If 'alias' is set to TRUE, then the column
#' of the transcript names will contain the aliases rather than the original
#' reference names of the .ribo file.
#'
#' If both 'length' and 'transcript' are TRUE, then the resulting
#' data table prints out one row for each experiment. This provides the metagene
#' information across all transcripts and all reads in a given experiment.
#'
#' If both length' and 'transcript' are FALSE, calculations are done to the data,
#' all information is preserved for both the read length and the transcript.
#' The data table would just present the entire stored raw data
#' from the read length 'range.lower' to the read length 'range.upper' which in most
#' cases would result in a slow run time with a massive data.table returned.
#'
#' When 'transcript' is set to FALSE, the 'alias' parameter specifies whether
#' or not the returned data.table should present each transcript as an alias
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
#' @param tidy Option to return the data table in a tidy format
#' @param region Specific region of interest
#' @param alias Option to report the transcripts as aliases/nicknames
#' @param experiments List of experiment names
#' @param normalize Option to normalize the counts as counts per million reads
#' @return A data table of the region counts
#' @importFrom rhdf5 h5read
#' @importFrom data.table data.table
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
                              experiments = get_experiments(ribo.object)) {
    range.info <- c(range.lower = range.lower, range.upper = range.upper)
    region     <- check_rc_input(ribo.object, region, range.info, experiments, alias)
    conditions <- c(transcript = transcript, length = length,
                    normalize  = normalize,   alias = alias)

    handle              <- ribo.object@handle
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
    paths <- vapply(matched.experiments, get_rc_path, FUN.VALUE = "character")
    data  <- lapply(paths, generate_matrix, ribo.object = ribo.object,
                                            transcript  = transcript,
                                            length      = length,
                                            normalize   = normalize,
                                            file        = handle,
                                            index       = list(columns, rows),
                                            ncol        = ncol)
    data           <- do.call(rbind, data)
    colnames(data) <- region
    result         <- make_datatable(ribo.object, matched.experiments,
                                     range.info, conditions, data)
    if (normalize) {
        info <- get_info(ribo.object)$experiment.info[, c("experiment", "total.reads")]
        result <- left_join(result, info, by = "experiment")
        result[, region] <- result[, region] * 1000000/result$total.reads
        result[, "total.reads"] <- NULL
    }
    if (tidy) result <- setDT(gather(result, "region", "count", region))
    return(result)
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
#' returned data table is valid input for {\code{\link{plot_length_distribution}}}.
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
#' @param total Include a column with the total reads, required for plotting
#' @importFrom dplyr left_join
#' @export
#' @seealso
#' {\code{\link{plot_length_distribution}}} to plot the output of this function
#' @return
#' A data table of the counts at each read length
get_length_distribution <- function(ribo.object,
                                    region,
                                    range.lower = rangeLower(ribo.object),
                                    range.upper = rangeUpper(ribo.object),
                                    experiments = get_experiments(ribo.object),
                                    total = TRUE) {
    if (length(region) != 1) {
      stop("Please provide only one region.")
    }

    rc <- get_region_counts(ribo.object,
                            range.lower = range.lower,
                            range.upper = range.upper,
                            region = region,
                            length = FALSE,
                            transcript = TRUE,
                            normalize = FALSE,
                            experiments = experiments)[, -"region"]

    if (total) {
        exp.info <- get_info(ribo.object)$experiment.info
        reads <- exp.info[, c("experiment", "total.reads")]
        rc %>% left_join(reads, by="experiment") %>% setDT() -> rc
  }
  return(rc)
}


#' Plots the length distribution
#'
#' The function \code{\link{plot_length_distribution}} can take either a data.table
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
#' The user can instead provide a data.table with the same structure as the
#' output of the \code{\link{get_region_counts}} function where the 'transcript'
#' parameter is set to FALSE and 'length' parameters is the default value of
#' TRUE. This also means that the many of the remaining parameters of the
#' \code{\link{plot_length_distribution}} function are not necessary. The run
#' time becomes substantially faster when \code{\link{plot_region_counts}} is
#' given the direct data.table to plot. Note that there is no manipulation by
#' this function on the data.table. This responsibility is given to the user
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
#' #data.table use case
#' #obtains the region counts at each individual read length, summed across every transcript
#' region.counts <- get_length_distribution(ribo.object = sample,
#'                                          region      = "CDS",
#'                                          range.lower = 2,
#'                                          range.upper = 5,
#'                                          experiments = experiments)
#'
#' #the param 'length' must be set to FALSE and param 'transcript' must be set
#' #to TRUE to use a data.table
#' plot_length_distribution(region.counts)
#'
#'
#' @seealso \code{\link{get_region_counts}} to generate a data.table that can
#' be provided as input,
#' \code{\link{ribo}} to create a ribo.object that can be provided as input
#' @param x A 'ribo' object or a data table generated from \code{\link{get_region_counts}}
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
        x %>% mutate(fraction = .data$count/.data$total.reads) -> x
        y.axis <- "Fraction"
        y.value <- "fraction"
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
        
        x <- get_length_distribution(ribo.object = x,
                                     region      = region,
                                     range.lower = range.lower,
                                     range.upper = range.upper,
                                     experiments = experiments)
    } else if (is.data.frame(x)){
        col.names <- c("experiment", "length", "count", "total.reads")
        types <- c("integer", "double")
        mismatch <- all(names(x) == col.names, 
                        typeof(x[[1]]) == "character",
                        typeof(x[[2]]) %in% types, 
                        typeof(x[[3]]) == "character",
                        typeof(x[[4]]) %in% types,
                        ncol(x) == 4)
        if (mismatch) {
              stop("Please make sure that the data table is of the correct format.",
                   call.=FALSE)
        }
    } else {
        stop("Please make sure that param 'x' is either", 
             "a data.table or a ribo object.",
             call.=FALSE)
    }
    return(x)
}


#' Plots the region counts of UTR5, CDS, and UTR3
#'
#' The function \code{\link{plot_region_counts}} can take either a data.table
#' or a "ribo" object to generate the a stacked bar plot of proportions that
#' correspond to the "UTR5", "CDS", and "UTR3" regions.
#'
#' When given a 'ribo' object, \code{\link{plot_region_counts}} calls
#' \code{\link{get_region_counts}} to retrieve the necessary information
#' for plotting. This option is in the case that a data.table of the
#' region count information is not required.
#'
#' The user can instead provide a data.table with the same structure as the
#' output of the \code{\link{get_region_counts}} function where the 'transcript'
#' and 'length' parameters are the default values of TRUE. This also means that
#' the remaining parameters of the \code{\link{plot_region_counts}} function are not necessary.
#' The run time becomes substantially faster when \code{\link{plot_region_counts}} is given
#' the direct data.table to plot. Note that there is no manipulation by this function on the
#' data.table, making this input option more error prone.
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
#' #data.table use case
#' #obtains the region counts at each individual read length, summed across every transcript
#' region.counts <- get_region_counts(sample,
#'                                    region = regions,
#'                                    range.lower = 2,
#'                                    range.upper = 5,
#'                                    tidy = TRUE,
#'                                    length = TRUE,
#'                                    transcript = TRUE)
#'
#' #the params 'length' and 'transcript' must be set to true to use a data.table
#' plot_region_counts(region.counts)
#'
#' @seealso \code{\link{get_region_counts}} to generate a data.table that can be provided as input,
#' \code{\link{ribo}} to create a ribo.object that can be provided as input
#'
#' @param x A 'ribo' object or a data table generated from \code{\link{get_region_counts}}
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
#' @return A 'ggplot' of the region counts
plot_region_counts <- function(x,
                               experiments,
                               range.lower,
                               range.upper,
                               title = "Region Counts") {
    rc <- check_plot_rc_input(x, range.lower, range.upper, experiments)

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
        x <- get_region_counts(ribo.object = x,
                               region      = regions,
                               range.lower = range.lower,
                               range.upper = range.upper,
                               length      = TRUE,
                               transcript  = TRUE,
                               experiments = experiments)
    } else if (is.data.frame(x)) {
        col.names <- c("experiment", "region", "count")
        mismatch  <- !all(names(x) == col.names, 
                     typeof(x[[1]]) == "character",
                     typeof(x[[2]]) == "character",
                     typeof(x[[3]]) == "double",
                     ncol(x) == 3)
        if (mismatch) {
            stop("Please make sure that the data table is of",  
                 "the correct format.", call. = FALSE)
        } else if (!identical(unique(x$region), regions)) {
            stop("Please make sure that the data table only includes the ", 
                 "'UTR5','CDS', and 'UTR3' regions.", call. = FALSE)
        }
    } else {
        stop("Please provide a ribo object or a data.table of the",
             "correct format.", call. = FALSE)
    }
    return(x)
}