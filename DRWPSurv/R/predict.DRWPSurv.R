predict.DRWPSurv <-
function(object,newx.mRNA, type="link",s="lambda.min")
{
	newx <- newx.mRNA
	newx <- scale(newx, center = TRUE, scale = TRUE)
	if(any(!(rownames(object$geneCoxZP) %in% colnames(newx)))){
		stop("The genes/miRNAs in the training set do not match the genes/miRNAs in the test set!")
	}	
	
	newx <- newx[,rownames(object$geneCoxZP)]
	pA <- getPathActivity(newx, object$pathSet, object$W, object$geneCoxZP)
	newx <- pA$pathActivity
		
	lpnew <- predict(object$fit.cox, newx, type=type, s=s)	
	lpnew
}
