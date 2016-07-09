\name{predict.DRWPSurv}
\alias{predict.DRWPSurv}

\title{
Make predictions from a "DRWPSurv" object
}
\description{
This function predicts the risk of new samples from a fitted \code{DRWPSurv} object.
}
\usage{
predict.DRWPSurv(object, newx.mRNA, type = "link", s = "lambda.min")
}

\arguments{
  \item{object}{
Fitted \code{"DRWPSurv"} model object.
}
  \item{newx.mRNA}{
A matrix with new samples to predict.
}
  \item{type}{
Type of prediction required. Type "link" gives the linear predictors for "cox" models. See \code{predict.glmnet} for details.
}
  \item{s}{
Value(s) of the penalty parameter lambda at which predictions are required. See \code{predict.glmnet} for details.
}
}
\value{
the fitted relative-risk of new patients in \code{newx.mRNA}.
}
\references{
  Liu, W., et al., Topologically inferring pathway activity for precise survival outcome prediction of breast cancer. 
  Unpublished.
}
\author{
Wei Liu
}

\seealso{
\code{\link{fit.DRWPSurv}}
}
\examples{
data(dGMMirGraph)
data(pathSet)
data(mRNA_matrix)
data(survData)
trainSmpl.Idx <- sample(1:dim(mRNA_matrix)[2], floor(4/5*dim(mRNA_matrix)[2]))
testSmpl.Idx <- setdiff(1:dim(mRNA_matrix)[2], trainSmpl.Idx)
trainSmpl <- mRNA_matrix[ ,trainSmpl.Idx]
testSmpl <- mRNA_matrix[ ,testSmpl.Idx]
fit <- fit.DRWPSurv(x.mRNA = t(trainSmpl), y = survData[trainSmpl.Idx,], DEBUG=TRUE, standardize=TRUE, globalGraph = dGMMirGraph, pathSet = pathSet, Gamma=0.7, alpha= 1, nfolds = 5)
predict.DRWPSurv(object = fit, newx.mRNA = t(testSmpl), type="link",s="lambda.min")
}

