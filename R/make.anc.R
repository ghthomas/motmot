make.anc <-
function(y, x, data=NULL, common.mean=FALSE) {
	
	if(is.factor(x) == "FALSE"){
		stop("The discrete trait must be a factor")	
	}
	
	m <- model.frame(y ~ x, data)
	x <- model.matrix(y ~ x, m)
	return(x)
}

