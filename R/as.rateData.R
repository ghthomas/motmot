as.rateData <-
function(y, x, rateMatrix=NULL, phy=NULL, data, meserr.col=NULL, meserr.propn=NULL, log.y=FALSE) {
	
		if (is.numeric(y)) { y <- colnames(data)[y] } else { y <- y }
		if (is.numeric(x)) { x <- colnames(data)[x] } else { x <- x }
		if (is.numeric(meserr.col) ) { meserr.col <- colnames(data)[meserr.col] } else { meserr.col <- meserr.col }

		if (!is.null(meserr.col) & !is.null(meserr.propn) ) {
				stop("Only one source of measurement error should be specified. Use messerr.col to specify a column of measurement error 
						OR use meserr.propn to specify a single numeric value as a proportion error (e.g. enter 0.05 for 5% measurement error.")
						}
	
		if (is.null(rateMatrix))	{ rateMatrix <- as.rateMatrix(phy=phy, x=x, data=data) } else { rateMatrix <- rateMatrix }
	
		if (is.null(meserr.col)) { dat <- data.frame(x=data[,x], y=data[,y], row.names = rownames(data)) }
					 
		if (!is.null(meserr.col)) { dat <- data.frame(x=data[,x], y=data[,y], meserr=data[,meserr.col], row.names = rownames(data)) }
		if (!is.null(meserr.propn)) { dat <- data.frame(x=data[,x], y=data[,y], meserr=rep(meserr.propn, dim(data)[1]), row.names = rownames(data)) }
					 
			 
	
	
		if(sum(sort(rownames(dat)) != sort(rownames(rateMatrix[[1]]))) > 0){warning("Length or names of phenotypic and of phylogenetic data do not match - non-matching taxa will be dropped")}	# check names

		
		# get names of species common to data frame and phylogeny
		sharedSpecies <- intersect(rownames(rateMatrix[[1]]), rownames(dat))
	
	
		# subset data to only those species in phylogeny and alphabetise
		dat <- dat[match(sharedSpecies ,rownames(dat)),] 
		dat<- dat[sort(rownames(dat), index.return = TRUE)$ix, ]
		
		# index missing data for y
		ccdat <- complete.cases(dat)
		idx <- which(ccdat)
		missing.dat <- rownames(dat)[which(!ccdat)]

	
		if(length(missing.dat > 0)) {cat("Dropping species due to missing data (x, y, or meserr):", "\n")
		
		cat("Dropped species: ", missing.dat, "\n")
		}
			
			
		# alphabetise and prune rate matrices to remove data with missing species
		Vmat <- vector(mode="list", length = length(rateMatrix))
		sumMat <- matrix(0, dim(rateMatrix[[1]])[1], dim(rateMatrix[[1]])[2])
		for(i in 1:length(rateMatrix)) {
			Vmatrix <- rateMatrix[[i]]
			sumMat <- sumMat + Vmatrix
			nms <- rownames(Vmatrix)
			snms <- sort(nms, index.return = TRUE)
			Vmatrix <- Vmatrix[snms$ix, snms$ix]
			Vmatrix <- Vmatrix[idx, idx]
			Vmat[[i]] <- Vmatrix
			
			}
		
	
			xdat <- as.matrix(dat$x)
			if (log.y==TRUE) {ydat <- log(dat[idx, "y"]) } else {ydat <- dat[idx, "y"]}
			xdat <- xdat[idx, ]
			
			if (is.null(meserr.col) & is.null(meserr.propn)) {traits <- list(y = ydat, x = xdat, Vmat = Vmat, meserr=rep(0,length(ydat)))} 
			if (!is.null(meserr.col)) {	meserr <- dat[idx,"meserr"]^2
				traits <- list(y = ydat, x = xdat, Vmat = Vmat, meserr=meserr)}
			
			if (!is.null(meserr.propn)) {	meserr <- (dat[idx,"meserr"]*ydat)^2
				traits <- list(y = ydat, x = xdat, Vmat = Vmat, meserr=meserr)}
	
	
			attr(traits, "class") <- "rateData"
			return(traits)
			}




