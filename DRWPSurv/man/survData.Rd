\name{survData}
\alias{survData}
\docType{data}
\title{
Survival data
}
\description{
The survival data of patients in \code{mRNA_matrix}.
}
\usage{data("survData")}
\format{
  The format is:
 num [1:100, 1:2] 0 0 0 0 0 1 0 0 1 1 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:100] "GSM79114" "GSM79115" "GSM79116" "GSM79118" ...
  ..$ : chr [1:2] "status" "time"
}
\examples{
data(survData)
}
\keyword{datasets}
