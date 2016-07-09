fit.DRWPSurv <-
function(x.mRNA, y,DEBUG=FALSE, standardize=TRUE, globalGraph, pathSet, Gamma=0.7, alpha= 1, nfolds = 5)
{
	if(standardize){
		x.mRNA <- scale(x.mRNA, center = TRUE, scale = TRUE)			
	}
	commonmRNA <- intersect(colnames(x.mRNA),V(globalGraph)$name)
	x.mRNA <- x.mRNA[, commonmRNA]
	x <- x.mRNA
		
	# calculate the p-value of each gene using Cox PH model
	if(DEBUG) cat('Calculating Cox score...')
	
	Survdata.mRNA <- data.frame(x.mRNA, y)
	geneCoxZP <- matrix(NA,nrow=ncol(x.mRNA),ncol=2)
	rownames(geneCoxZP) <- colnames(x.mRNA)
	colnames(geneCoxZP) <- c("ZScore","p-value")
	for(i in 1 : ncol(x.mRNA)){
		res.coxph <- coxph(as.formula(paste("Surv(time, status)~", colnames(Survdata.mRNA)[i])), Survdata.mRNA)
		geneCoxZP[i,] <- summary(res.coxph)$coefficients[c(4,5)]		
	}
	geneWeight <- -log(geneCoxZP[ ,2]+2.2e-16)
	geneWeight[which(is.na(geneWeight))] <- 0
	geneWeight <- (geneWeight - min(geneWeight)) / (max(geneWeight) - min(geneWeight))
			
	W <- getW(geneWeight=geneWeight, globalGraph=globalGraph)
		
	if(DEBUG) cat('Done\n')
	
	if(DEBUG) cat('Performing directed random walk...')
	vertexWeight <- DRW(igraphM = globalGraph, p0 = W, EdgeWeight=FALSE, gamma = Gamma)
	if(DEBUG) cat('Done\n')
	
	pA <- getPathActivity(x, pathSet, vertexWeight, geneCoxZP)
	x <- pA$pathActivity
			
	if(DEBUG) cat('Fitting DRWPSurv model...')
	
	fit.cox <- cv.glmnet(x,y,family="cox",alpha=alpha, nfolds=nfolds,maxit=10^6)
	if(DEBUG) cat('Done\n')
	pAcoef <- coef(fit.cox,s="lambda.min")
	sigPA <- intersect(colnames(pA$pathActivity),rownames(pAcoef)[which(pAcoef!=0)])
	# sigfactors <- setdiff(rownames(pAcoef)[which(pAcoef!=0)], sigPA)
	sigGenes <- unique(unlist(pA$sigGenes[sigPA]))
	features <- rownames(pAcoef)[which(pAcoef!=0)]

	fit <- list(fit.cox=fit.cox, W=W, geneCoxZP=geneCoxZP, globalGraph=globalGraph, pathSet=pathSet, features=features, sigGenes=sigGenes, sigPathGenes=pA$sigGenes[sigPA])
	class(fit) <- "DRWPSurv"
	fit
}
