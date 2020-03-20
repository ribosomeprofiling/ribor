\name{riborNews}
\title{ribor News}
\encoding{UTF-8}

\section{Version 0.99.0}{\itemize{
\item
  Submitted to Bioconductor
}}

\section{Version 0.99.7}{\itemize{
\item {
  All parameter names are in singular version. The parameters `region` and  
  `experiment` are the parameter names for all applicable functions. The
  `experiment` parameter can still take in a list of multiple experiments 
  and the `region` parameter can as well except for the case of getting the   length distribution.
}

\item {
  Validator and accessor methods have been added. Use ?Ribo for more
  details on the available methods.
}

\item {
  The return types for the reader functions are currently either S4Vectors
  DataFrame or the R native data.frame, depending on the presence of the 
  added `compact` parameter.
}

\item {
  Documentation has been added to be more thorough.
}
}