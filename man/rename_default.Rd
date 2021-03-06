% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{rename_default}
\alias{rename_default}
\title{Rename function for appris transcriptome naming convention}
\usage{
rename_default(x)
}
\arguments{
\item{x}{Character denoting original name of the transcript}
}
\value{
Character denoting simplified name of the object
}
\description{
The function {\code{\link{rename_default}}} is the default renaming function for the
appris human transcriptome. It takes one single transcript name and returns a simplified
alias.
}
\examples{
original <- paste("ENST00000613283.2|ENSG00000136997.17|",
                  "OTTHUMG00000128475.8|-|MYC-206|MYC|1365|protein_coding|",
                  sep = "")
alias <- rename_default(original)
}
