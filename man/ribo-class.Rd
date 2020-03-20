\name{Ribo-class}
\alias{Ribo-class}
\alias{Ribo}
\alias{class:Ribo}
\alias{show, Ribo-method}

\alias{path}
\alias{path, Ribo-method}

\alias{experiments}
\alias{experiments, Ribo-method}

\alias{format_version}
\alias{format_version, Ribo-method}

\alias{reference}
\alias{reference, Ribo-method}

\alias{length_min}
\alias{length_min, Ribo-method}

\alias{length_max}
\alias{length_max, Ribo-method}

\alias{left_span}
\alias{left_span, Ribo-method}

\alias{right_span}
\alias{right_span, Ribo-method}

\alias{metagene_radius}
\alias{metagene_radius, Ribo-method}

\alias{length_offset}
\alias{length_offset, Ribo-method}

\alias{has_metadata}
\alias{has_metadata, Ribo-method}

\alias{experiment_info}
\alias{experiment_info, Ribo-method}

\alias{transcript_info}
\alias{transcript_info, Ribo-method}

\alias{alias_hash}
\alias{alias_hash, Ribo-method}

\alias{original_hash}
\alias{original_hash, Ribo-method}

\title{Ribo Class}
\description{
  The Ribo object serves as the main utility vehicle for the ribor package.
  Specifically, it allows the user to interface with a .ribo file in the R
  environment through the use of the ribor package. Almost all functions in
  ribor rely on the Ribo object to read, visualize, and inspect the
  contents of the .ribo file. 
  
  The information stored in th is object include the .ribo file path, the 
  list of experiments, the format version, the reference model, the minimum
  read length, maximum read length, the left span, the right span, and 
  other transcript information.
}

\usage{
Ribo(file_path, rename)
}


\arguments{
  \item{name}{
    The path to the .ribo file
  }
  \item{rename}{
    A function that renames the original transcript or an already    
    generated character vector of aliases
  }
}
\details{
  Note that the name parameter takes in a file path and stores it. While 
  using the package, be sure to not to move or change the location of the 
  .ribo file.
  
  The default names of the transcripts may be difficult to use depending on the 
  settings used to generate the .ribo file. As a result, we have provided a 
  rename parameter that integrates well with the Appris reference transcriptome.
  Users may also define a simple function that processes a given default 
  transcript name in a one-to-one manner to another custom alias. 
}


\section{Accessors}{
  Note that no mutators are provided. The Ribo object stays valid if and only if 
  the original slot values declared in the constructor function. Subsequent 
  function calls with an illegally mutated Ribo object will raise errors.
  \itemize{
    \item \code{experiments(x)} to get a character vector of the experiments 

    \item \code{format_version(x)} to get the format version

    \item \code{reference(x)} to get the reference
    
    \item \code{length_min(x)} to get the lower bound of the read length, inclusive

    \item \code{length_max(x)} to get the upper bound of the read length, inclusive
    
    \item \code{left_span(x)} to get the left span
    
    \item \code{right_span(x)} to get the right span
    
    \item \code{metagene_radius(x)} to get the metagene radius
    
    \item \code{length_offset(x)} to get the length offset
    
    \item \code{has_metadata(x)} to get a logical value denoting the presence of 
    metadata on the .ribo file
    
    \item \code{experiment_info(x)} to get information on the experiments
    
    \item \code{transcript_info(x)} to get specific information on transcript 
    offsets and lengths
    
    \item \code{alias_hash(x)} to get the hash mapping from alias name to original 
    transcript name
    
    \item \code{original_hash(x)} to get the hash mapping from original name to 
    alias name
  }
}

\value{A Ribo Object}

\examples{
#generate the ribo object
file.path <- system.file("extdata", "sample.ribo", package = "ribor")
sample <- Ribo(file.path)
}