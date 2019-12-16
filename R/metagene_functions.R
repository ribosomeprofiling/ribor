#' Retrieves the metagene data from a .ribo file
#'
#' The function \code{\link{get_metagene}} returns a data frame that provides
#' the coverage at the positions surrounding the metagene start or stop site.
#'
#' The dimensions of the returned data frame depend on the parameters
#' range.lower, range.upper, length, and transcript.
#'
#' The param 'length' condenses the read lengths together.
#' When length is TRUE and transcript is FALSE, the
#' data frame presents information for each transcript across
#' all of the read lengths. That is, each transcript has a value
#' that is the sum of all of the counts across every read length.
#' As a result, information about the transcript at each specific
#' read length is lost.
#'
#' The param 'transcripts' condenses the transcripts together.
#' When transcript is TRUE and length is FALSE, the data
#' table presents information at each read length between range.lower and
#' range.upper inclusive. That is, each separate read length denotes the
#' sum of counts from every transcript. As a result, information about the
#' counts of each individual transcript is lost.
#'
#' If both 'length' and 'transcript' are TRUE, then the resulting
#' data frame prints out one row for each experiment. This provides the metagene
#' information across all transcripts and all reads in a given experiment.
#'
#' If both length' and 'transcript' are FALSE, no calculations are done to the data,
#' all information is preserved for both the read length and the transcript.
#' The data frame would just present the entire stored raw data
#' from the read length 'range.lower' to the read length 'range.upper' which in most
#' cases would result in a slow run time with a massive DataFrame returned.
#'
#' When 'transcript' is set to FALSE, the 'alias' parameter specifies whether
#' or not the returned DataFrame should present each transcript as an alias
#' instead of the original name. If 'alias' is set to TRUE, then the returned
#' data frame will contain the aliases rather than the original
#' reference names of the .ribo file.
#'
#' @param ribo.object A 'ribo' object
#' @param site "start" or "stop" site coverage
#' @param range.lower Lower bound of the read length, inclusive
#' @param range.upper Upper bound of the read length, inclusive
#' @param transcript Option to condense the transcripts together, preserving information at each read length
#' @param length Option to condense the read lengths together, preserving information at each transcript
#' @param experiments List of experiment names
#' @param alias Option to report the transcripts as aliases/nicknames
#' @return A DataFrame of the metagene information
#' @examples
#'
#' #generate the ribo object by providing the file.path to the ribo file
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- create_ribo(file.path)
#'
#'
#' #extract the total metagene information for all experiments
#' #across the read lengths and transcripts of the start site
#' #from read length 2 to 5
#' metagene_info <- get_metagene(ribo.object = sample,
#'                               site = "start",
#'                               range.lower = 2,
#'                               range.upper = 5,
#'                               length = TRUE,
#'                               transcript = TRUE,
#'                               experiments = get_experiments(sample))
#'
#'
#' #Note that length, transcript, and experiments in this case are the
#' #default values and can be left out. The following generates the same output.
#'
#' metagene_info <- get_metagene(ribo.object = sample,
#'                               site = "start",
#'                               range.lower = 2,
#'                               range.upper = 5)
#'
#' @seealso
#' \code{\link{create_ribo}} to generate the necessary 'ribo' class object,
#' \code{\link{plot_metagene}} to visualize the metagene data,
#' \code{\link{get_tidy_metagene}} to obtain tidy metagene data under certain conditions
#' @importFrom rhdf5 h5read
#' @importFrom methods as 
#' @importFrom S4Vectors DataFrame Rle
#' @importFrom hash keys
#' @export
get_metagene <- function(ribo.object,
                         site,
                         range.lower = rangeLower(ribo.object),
                         range.upper = rangeUpper(ribo.object),
                         transcript = TRUE,
                         length = TRUE,
                         alias = FALSE,
                         experiments = get_experiments(ribo.object)) {
    range.info <- c(range.lower = range.lower, range.upper = range.upper)
    conditions <- c(transcript = transcript, length = length, alias = alias)
    site <- tolower(site)
    check_metagene_input(ribo.object, site, range.info, experiments, alias)

    path              <- ribo.object@path
    range.min           <- get_read_lengths(ribo.object)[1]
    metagene.radius     <- h5readAttributes(path, "/")[["metagene_radius"]]
    ncol                <- 2 * metagene.radius + 1
    columns             <- seq(ncol)
    ref.length          <- length(get_reference_names(ribo.object))

    row.start <- (range.lower - range.min) * ref.length + 1
    row.stop <- row.start + ref.length*(range.upper - range.lower + 1) - 1
    rows <- c(row.start:row.stop)

    #gather information to use in filling and labeling final data frame
    matched.experiments <- intersect(experiments, get_experiments(ribo.object))
    exp_paths <- vapply(matched.experiments, get_metagene_path, site = site,
                        FUN.VALUE = "character")
    data  <- lapply(exp_paths, generate_matrix,
                               ribo.object = ribo.object,
                               transcript = transcript,
                               length = length,
                               normalize = FALSE,
                               ncol = ncol,
                               file = path,
                               index = list(columns, rows))

    data <- as.data.frame(do.call(rbind, data))
    colnames(data) <- c(as.character(-metagene.radius:metagene.radius))
    
    result <- as(make_dataframe(ribo.object,
                                matched.experiments,
                                range.info,
                                conditions,
                                data), "DataFrame")
    
    return(prepare_DataFrame(ribo.object, result))
}


get_metagene_path <- function(experiment, site) {
    # helper method that generates a path within the ribo file
    dataset.name <- paste(site, "_site_coverage", sep = "")
    data.path <- paste("/experiments/", experiment, "/metagene/", sep = "")
    path <- paste(data.path, dataset.name, sep = "")
    return(path)
}
#' Retrieves the metagene data in a tidy format
#'
#' The function \code{\link{get_tidy_metagene}} provides the user with a tidy data format for easier
#' data cleaning and manipulation. In providing this functionality, the user must length the
#' transcripts together and is only provided the option to length the read lengths together.
#'
#' The dimensions of the returned data frame depend on the parameters
#' range.lower, range.upper, and length.
#'
#' The param 'length' condenses the read lengths together.
#' When length is TRUE, then the resulting data frame prints out one row
#' for each experiment. This provides a tidy format of the metagene information
#' across all transcripts and all read lengths in a given experiment. Each row
#' in the data frame represents the total metagene coverage count of a given experiment
#' at a given position.
#'
#' When the param  'length' is FALSE, then the resulting data frame prints out the
#' metagene coverage count at each position of the metagene radius for each read length.
#' This provides a tidy format of the metagene information across the transcripts, preserving
#' the metagene coverage count at each read length.
#'
#' @param ribo.object A 'ribo' object
#' @param site "start" or "stop" site coverage
#' @param range.lower Lower bound of the read length
#' @param range.upper Upper bound of the read length
#' @param length Option to condense the read lengths together
#' @param alias Option to report the transcripts as aliases/nicknames
#' @param experiments List of experiment names
#' @return
#' A tidy data frame of the metagene information
#' @examples
#' #generate the ribo object by loading in a ribo function and calling the \code{\link{ribo}} function
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- create_ribo(file.path)
#'
#' #extract the total metagene information in a tidy format
#' #for all experiments across the read lengths and transcripts
#' #of the start site from read length 2 to 5
#'
#' metagene_info <- get_tidy_metagene(ribo.object = sample,
#'                                    site = "start",
#'                                    range.lower = 2,
#'                                    range.upper = 5,
#'                                    length = TRUE,
#'                                    experiments = get_experiments(sample))
#'
#' #Note that length and experiments in this case are the
#' #default values and can be left out. The following generates the same output.
#' metagene_info <- get_tidy_metagene(ribo.object = sample,
#'                                    site = "start",
#'                                    range.lower = 2,
#'                                    range.upper = 5)
#'
#' @seealso
#' \code{\link{create_ribo}} to generate the necessary 'ribo' class object.
#' \code{\link{plot_metagene}} to visualize the metagene data,
#' \code{\link{get_metagene}} to obtain tidy metagene data under certain conditions
#' @importFrom rhdf5 h5read
#' @export
get_tidy_metagene <- function(ribo.object,
                              site,
                              range.lower = rangeLower(ribo.object),
                              range.upper = rangeUpper(ribo.object),
                              length = TRUE,
                              alias = FALSE,
                              experiments = get_experiments(ribo.object)) {
  site <- tolower(site)
  result <- get_metagene(ribo.object,
                         site,
                         range.lower,
                         range.upper,
                         length,
                         transcript = TRUE,
                         alias = alias,
                         experiments = experiments)
  
  result <- strip_rlefactor(result)
  metagene.radius <- as.integer((ncol(result) - 2) / 2)
  result <- data.frame(result, check.names = FALSE)
  tidy.data <- gather(result,
                      key = "position",
                      value = "count",
                      c(as.character(-metagene.radius:metagene.radius)))
  tidy.data$position <- as.integer(tidy.data$position)
  tidy.data <- as(tidy.data, "DataFrame")
  return(prepare_DataFrame(ribo.object, tidy.data))
}


check_metagene_input <- function(ribo.object,
                                 site,
                                 range.info,
                                 experiments,
                                 alias) {
  #check_metagene_input is a helper function that checks the validity of
  #the metagene function parameters
  #check param validity
  if (site != "start" & site != "stop") {
    stop("Please type 'start' or 'stop' to indicate the 'site' parameter value.")
  }

  range.lower <- range.info[["range.lower"]]
  range.upper <- range.info[["range.upper"]]

  check_alias(ribo.object, alias)
  check_ribo(ribo.object)
  check_lengths(ribo.object, range.lower, range.upper)
  check_experiments(ribo.object, experiments)
}

#' Plots the metagene coverage data
#'
#' The function \code{\link{plot_metagene}} plots the metagene site coverage,
#' separating by experiment.
#'
#' If a DataFrame is provided as param 'x', then the only additional parameter
#' is the optional title' parameter for the generated plot. If a ribo.object is
#' provided as param 'x', the rest of the parameters listed are necessary.
#'
#' When given a ribo class object, the \code{\link{plot_metagene}} function
#' generates a DataFrame by calling the \code{\link{get_tidy_metagene}}
#' function, so the run times in this case will be mostly comprised of a call
#' to the \code{\link{get_metagene}} function.
#'
#' This function uses ggplot in its underlying implementation.
#'
#' @param x A 'ribo' object or a data frame generated from \code{\link{get_metagene}}
#' @param site "start" or "stop" site
#' @param range.lower lower bound of the read length, inclusive
#' @param range.upper upper bound of the read length, inclusive
#' @param experiments list of experiments
#' @param normalize When TRUE, normalizes the data by the total reads.
#' @param title title of the generated plot
#' @param tick x-axis labeling increment
#' @examples
#' #a potential use case is to directly pass in the ribo object file as param 'x'
#'
#' #generate the ribo object to directly use
#' file.path <- system.file("extdata", "sample.ribo", package = "ribor")
#' sample <- create_ribo(file.path)
#'
#' #specify experiments of interest
#' experiments <- c("Hela_1", "Hela_2", "WT_1")
#'
#' #plot the metagene start site coverage for all experiments in 'sample.ribo'
#' #from read length 2 to 5
#' plot_metagene(x = sample,
#'               site = "start",
#'               range.lower = 2,
#'               range.upper = 5,
#'               experiments = experiments)
#'
#' #Note that the site, range.lower, range.upper, and experiments are only
#' #necessary if a ribo object is being passed in as param 'x'. If a ribo
#' #object is passed in, then the param 'experiments' will be set to all of
#' #the experiments by default.
#'
#' #If a DataFrame is passed in, then the plot_metagene function
#' #does not need any other information. All of the elements of the DataFrame
#' #will be used, assuming that it contains the same column names and number of
#' #columns as the output from get_tidy_metagene()
#'
#' #gets the metagene start site coverage from read length 2 to 5
#' #note that the data must be summed across transcripts and read lengths
#' #for the plot_metagene function
#' data <- get_tidy_metagene(sample,
#'                           site = "start",
#'                           range.lower = 2,
#'                           range.upper = 5)
#'
#' #plot the metagene data
#' plot_metagene(data)
#'
#' @importFrom dplyr left_join mutate
#' @importFrom ggplot2 ggplot geom_line theme_bw ggtitle aes expand_limits
#' @importFrom ggplot2 element_text theme labs scale_x_continuous
#' @importFrom rlang .data
#' @importFrom tidyr gather
#' @export
#' @return
#' A 'ggplot' of the metagene site coverage
plot_metagene <- function(x,
                          site,
                          experiments,
                          range.lower,
                          range.upper,
                          normalize = FALSE,
                          title = "Metagene Site Coverage",
                          tick = 10) {
    x <- check_plot_metagene(x,
                             site,
                             range.lower,
                             range.upper,
                             experiments)

    y.value <- "count"
    y.label <- "Count"
    
    if (normalize) {
        per.million <- 1000000
        info <- x@metadata[[1]]
        x %>% 
          strip_rlefactor() %>% 
          as.data.frame() %>% 
          left_join(info, by = "experiment") %>% 
          mutate(normalize = per.million * .data$count/.data$total.reads) -> x
        y.value <- "normalize"
        y.label <- "Counts per 1M Reads"
    } else {
      x <- as.data.frame(x)
    }

    metagene.radius <- max(x$position)
    ggplot(x,
           aes_string(x     = "position", 
                      y     = y.value, 
                      color = "experiment")) +
    scale_x_continuous(breaks=seq(-metagene.radius, metagene.radius, tick)) +
    expand_limits(y=0) +
    geom_line() +
    theme_bw() +
    theme(plot.title=element_text(hjust = 0.5)) +
    labs(title=title, x="Position", y=y.label, color="Experiment")
}

check_plot_metagene <- function(x,
                                site,
                                range.lower,
                                range.upper,
                                experiments) {
    is.ribo <- check_ribo(x, stop = FALSE)

    #x is a ribo object
    if (is.ribo) {
        if (missing(site)) {
          stop("Please indicate the 'site' parameter with either 'start' or 'stop'",
               call. = FALSE)
        } 
        if (missing(experiments)) experiments = get_experiments(x)
        if (missing(range.lower)) range.lower <- rangeLower(x)
        if (missing(range.upper)) range.upper <- rangeUpper(x)
        x <- strip_rlefactor(get_tidy_metagene(x,
                                               site,
                                               range.lower,
                                               range.upper,
                                               length = TRUE,
                                               experiments = experiments))
    } else if (class(x) == "DataFrame" || class(x) == "DFrame") {
        x <- strip_rlefactor(x)
        col.names <- c("experiment", "position", "count")
        types <- c("integer", "double")
        mismatch <- !all(names(x) == col.names,
                         typeof(x[[1]]) == "character",
                         typeof(x[[2]]) %in% types,
                         typeof(x[[3]]) %in% types,
                         ncol(x) == 3)
        if (mismatch) {
            stop("Please make sure that the data frame is of the correct format.",
                 call.=FALSE)
        }
    } else {
        #not a data frame
        stop("Please make sure that param 'x' is either ", 
             "a DataFrame or a ribo object.")
    }
    return(x)
}
