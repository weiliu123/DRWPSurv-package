getPathActivity <-
function(x, pathSet, w, vertexZP)
{
	# infer pathway expression profile
	# output:
	# the pathway expression profile
	
	# infer non metabolic pathway expression profile
	pathActivity <- c()
	sigGenes <- vector("list", length = length(pathSet))  # The differential genes to construct pathActivity
	names(sigGenes) <- names(pathSet)
	for (i in 1 : length(pathSet)){
		Vpathwayi <- pathSet[[i]]
		if (length(Vpathwayi) > 0){
			n <- 0    # the number of differential genes in ith pathway 			
			pathActivity_tmp <- matrix(nrow=dim(x)[1],ncol=1,data=0) 
			sigGenesi <- c()	
			Idx_pathwayi <- c()   
			for (j in 1 : length(Vpathwayi)){
				Idx <- which(colnames(x)==Vpathwayi[j])
				if (length(Idx) > 0){
					if ( colnames(x)[Idx] %in% names(w)){
						if(vertexZP[colnames(x)[Idx],"p-value"] < 0.05){ # p-value < 0.05
						# browser()
							pathActivity_tmp <- pathActivity_tmp + sign(vertexZP[colnames(x)[Idx],"ZScore"]) * w[colnames(x)[Idx]] * x[,Idx]
							n <- n + 1
							Idx_pathwayi <- rbind(Idx_pathwayi,Idx)
							sigGenesi <- c(sigGenesi, Vpathwayi[j])							
						}
					}
				}
			}
			if(n > 0){
				pathActivity_tmp <- pathActivity_tmp / sqrt(sum(w[colnames(x)[Idx_pathwayi]]^2))
				colnames(pathActivity_tmp) <- names(pathSet)[i] 
				pathActivity <- cbind(pathActivity, pathActivity_tmp)
				sigGenes[[i]] <- sigGenesi
			}#else
		}#else
	}
	rownames(pathActivity) <- rownames(x)
	return(list(pathActivity=pathActivity, sigGenes=sigGenes))		
}
