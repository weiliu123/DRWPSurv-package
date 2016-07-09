\name{DRWPSurv-package}
\alias{DRWPSurv-package}
\alias{DRWPSurv}
\docType{package}
\title{
\packageTitle{DRWPSurv}
}
\description{
\packageDescription{DRWPSurv}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{DRWPSurv}
\packageIndices{DRWPSurv}
Very simple to use. Accepts x,y data for Lasso-Cox models, and makes predictions for new samples. 
Only 2 functions:
\code{fit.DRWPSurv} topologically infers pathway activities and fits the Lasso-Cox model.
\code{predict.DRWPSurv} predicts the risk of new patients using the Lasso-Cox model fitted by \code{fit.DRWPSurv}.
}
\author{
\packageAuthor{DRWPSurv}

Maintainer: \packageMaintainer{DRWPSurv}
}
\references{
  Liu, W., et al., Topologically inferring pathway activity for precise survival outcome prediction of breast cancer. 
  Unpublished.
}

\keyword{ package }

\examples{
data(dGMMirGraph)
data(pathSet)
data(mRNA_matrix)
data(survData)
trainSmpl.Idx <- sample(1:dim(mRNA_matrix)[2], floor(4/5*dim(mRNA_matrix)[2]))
testSmpl.Idx <- setdiff(1:dim(mRNA_matrix)[2], trainSmpl.Idx)
trainSmpl <- mRNA_matrix[ ,trainSmpl.Idx]
testSmpl <- mRNA_matrix[ ,testSmpl.Idx]
fit <- fit.DRWPSurv(x.mRNA = t(trainSmpl), y = survData[trainSmpl.Idx,], DEBUG=TRUE, 
                    standardize=TRUE, globalGraph = dGMMirGraph, pathSet = pathSet, 
                    Gamma=0.7, alpha= 1, nfolds = 5)
predict.DRWPSurv(object = fit, newx.mRNA = t(testSmpl), type="link",s="lambda.min")
}
