\name{riborNews}
\title{ribor News}
\encoding{UTF-8}

\section{Version 0.99.0}{\itemize{
\item
Submitted to Bioconductor
}}

\section{Version 0.99.7}{\itemize{
\item
All parameter names are in singular version.
Specifically, the parameters `region` and `experiment` are the names
for all applicable functions. The `experiment` parameter can
still take in a list of multiple experiments, and the region
parameter can as well except for the case of getting the length
distribution.
Validator and accessor methods have been added. Use ?Ribo for more
details on the available methods.
The return types are currently either s4Vectors DataFrame or the 
R native data.frame, depending on the presence of the added `compact` parameter.
Documentation has been added to be more thorough.
}}