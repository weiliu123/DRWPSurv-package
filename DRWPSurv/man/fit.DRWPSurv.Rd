\name{fit.DRWPSurv}
\alias{fit.DRWPSurv}
\title{
Fit a Lasso-Cox model using DRWPSurv
}
\description{
DRWPSurv is a pathway-based survival prediction method which topologically infers survival associated pathway activities 
and constructs Lasso-Cox models using the pathway activities as features. 
}
\usage{
fit.DRWPSurv(x.mRNA, y, DEBUG = FALSE, standardize = TRUE, globalGraph, pathSet, 
             Gamma = 0.7, alpha = 1, nfolds = 5)
}
\arguments{
  \item{x.mRNA}{
a n x p matrix of gene expression measurements with n genes and p samples.
}
  \item{y}{
a p x 2 matrix of survival data. The two columns represent survival time 'time' and disease status 'status' respectively. 
}
  \item{DEBUG}{
Logical. Should debugging information be plotted.
}
  \item{standardize}{
Logical flag for \code{x.mRNA} standardization, prior to fitting the model. Default is \code{TRUE}.
}
  \item{globalGraph}{
An \code{igraph} R object containing the global pathway graph.
}
  \item{pathSet}{
A list of pathways. Each pathway is represented as a vector of pathway member genes and metabolites.
}
  \item{Gamma}{
A numeric value. The restart probability in directed random walk. Default is \code{Gamma = 0.7}.
}
  \item{alpha}{
The elasticnet mixing parameter in \code{glmnet}.
}
  \item{nfolds}{
Number of folds - default is 10. Although \code{nfolds} can be as large as the sample size (leave-one-out CV), 
it is not recommended for large datasets. Smallest value allowable is \code{nfolds=3}.
}
}
\details{
DRWPSurv integrates gene expression profiles and prior gene interaction information to topologically infer survival 
associated pathway activities, and uses the pathway activities as features to construct Lasso-Cox model. 
It uses topological importance of genes evaluated by directed random walk to enhance the robustness of pathway activities 
and thereby improve the predictive performance.
}
\value{
Fitted \code{"DRWPSurv"} model object.
	\item{fit.cox}{An object of class "cv.glmnet"}
	\item{W}{The weights of nodes in \code{globalGraph}.}
	\item{geneCoxZP}{The z statistic and p-value (univariate Cox regression) of nodes in \code{globalGraph}.}
	\item{globalGraph}{An \code{igraph} R object containing the global pathway graph.}
	\item{pathSet}{A list of pathways. Each pathway is represented as a vector of pathway member genes and metabolites.}
	\item{features}{The selected pathway features to construct the Lasso-Cox model.}
	\item{sigGenes}{A vector of genes used to infer the pathway activities in \code{features}.}
	\item{sigPathGenes}{A list of pathways in \code{features}. Each pathway contains the member genes used to infer its pathway activity.}
}
\references{
  Liu, W., et al., Topologically inferring pathway activity for precise survival outcome prediction of breast cancer. 
  Unpublished.
}
\author{
Wei Liu
}

\seealso{
 \code{\link{predict.DRWPSurv}}
}
\examples{
data(dGMMirGraph)
data(pathSet)
data(mRNA_matrix)
data(survData)
trainSmpl.Idx <- sample(1:dim(mRNA_matrix)[2], floor(4/5*dim(mRNA_matrix)[2]))
# testSmpl.Idx <- setdiff(1:dim(mRNA_matrix)[2], trainSmpl.Idx)
trainSmpl <- mRNA_matrix[ ,trainSmpl.Idx]
# testSmpl <- mRNA_matrix[ ,testSmpl.Idx]
fit <- fit.DRWPSurv(x.mRNA = t(trainSmpl), y = survData[trainSmpl.Idx,], DEBUG=TRUE, 
                    standardize=TRUE, globalGraph = dGMMirGraph, pathSet = pathSet, 
                    Gamma=0.7, alpha= 1, nfolds = 5)
}

