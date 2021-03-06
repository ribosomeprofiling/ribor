% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq_functions.R
\name{get_rnaseq}
\alias{get_rnaseq}
\title{Information on the RNA-Seq data of the experiments, if any}
\usage{
get_rnaseq(
  ribo.object,
  tidy = TRUE,
  region = c("UTR5", "UTR5J", "CDS", "UTR3J", "UTR3"),
  experiment = experiments(ribo.object),
  compact = TRUE,
  alias = FALSE
)
}
\arguments{
\item{ribo.object}{A 'Ribo' object}

\item{tidy}{Option to return the data frame in a tidy format}

\item{region}{Specific region(s) of interest}

\item{experiment}{List of experiment names}

\item{compact}{Option to return a DataFrame with Rle and factor as opposed to a raw data.frame}

\item{alias}{Option to report the transcripts as aliases/nicknames}
}
\value{
An annotated data frame containing the RNA-Seq counts for the regions in specified in the `region` parameter with the option of
presenting the data in a tidy format. Additionally, the function returns a DataFrame with Rle and factor applied if the `compact` parameter
is set to TRUE and a data.frame without any Rle or factor if the `compact` parameter is set to FALSE
}
\description{
\code{\link{get_rnaseq}} returns a data frame containing information on the transcript name, experiment, and
sequence abundance
}
\details{
As a default value, experiment.list is presumed to include all of the
experiments within a ribo file. RNA-Seq data is an optional dataset to
include in a .ribo file. The experiments in experiment.list are checked
for experiment existence in the ribo file and then checked for RNA-seq data.

The returned DataFrame can either be in the tidy format for easier data
cleaning or in a condensed non-tidy format. The data will present RNA-seq counts
for each transcript in each valid experiment in experiment.list.

The 'alias' parameter specifies whether or not the returned DataFrame
should present each transcript as an alias instead of the original name.
If 'alias' is set to TRUE, then the column of the transcript names will
contain the aliases rather than the original reference names of the .ribo
file.
}
\examples{
#generate the ribo object
file.path <- system.file("extdata", "sample.ribo", package = "ribor")
sample <- Ribo(file.path)

#list out the experiments of interest that have RNA-Seq data
experiments <- c("Hela_1", "Hela_2", "WT_1")
regions <- c("UTR5", "CDS", "UTR3")
rnaseq.data <- get_rnaseq(ribo.object = sample,
                          tidy = TRUE,
                          region = regions,
                          experiment = experiments)


}
\seealso{
\code{\link{Ribo}} to generate the necessary ribo.object parameter
}
