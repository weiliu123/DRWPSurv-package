\name{getPathActivity}
\alias{getPathActivity}

\title{
Inferring pathway activity
}
\description{
This function infers pathway activities using the DRWPSurv method.
}
\usage{
getPathActivity(x, pathSet, w, vertexZP)
}

\arguments{
  \item{x}{
    A n x p matrix of gene expression measurements with n samples and p genes.
}
  \item{pathSet}{
	A list of pathways. Each pathway is represented as a vector of pathway member genes and metabolites.
}
  \item{w}{
	A numerical vector containing the topological weights of nodes in \code{globalGraph}.
}
  \item{vertexZP}{
	A p x 2 matrix which contains the z statistic and p-value of p genes in \code{x}.
}
}
\details{
For each pathway, we selected the genes with p-value < 0.05 (univariate cox regression) in the pathway to infer pathway activity. 
The expression values of genes are weighted by their topological weights obtained from directed random walk on the global pathway graph. 
Pathway activity inference transforms the gene expression profiles into pathway activity profiles, which are then used to fit the Lasso-Cox model.
}
\value{
	\item{pathwayActivity}{The pathway activities of pathways in \code{pathSet}.}
	\item{sigGenes}{The genes used to infer the pathway activity.}
}
\references{
Liu, W., et al., Topologically inferring risk-active pathways toward precise cancer classification by directed random walk. 
  Bioinformatics, 2013. 29(17): p. 2169-77.
}
\author{
Wei Liu
}

\examples{

}

