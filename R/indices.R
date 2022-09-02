#' Scalar indices of diversity 
#' 
#' This function includes some indices that characterize a species-abundance/occurrence distribution.
#' 
#' This set is not complete and does not intend to supercede additional R packages (e.g. vegan). However, some metrics are presented here as they are not
#' implemented elsewhere or because they are invoked more frequently. The following entries can be added to the \code{method} argument of the function, which are
#' also named accordingly in the output table/vector.
#' 
#' \code{"richness"}: The number of sampled species. 
#' 
#' \code{"shannon"}: The Shannon entropy.
#' 
#' \code{dom}: The Berger-Parker dominance index, the proportion of occurrences in the time bin that belong to the most frequent taxon.
#' 
#' \code{"hill2"}: The second order Hill number (Jost, 2006; q=2), which will be calculated by default. You can specify additional Hill numbers with adding \code{"hillXX"} to the \code{method}
#' argument, such as \code{"hill3"} for (q=3). The first Hill number is defined as the exponentiad version of Shannon entropy (Eq. 3 in Jost, 2006). 
#' 
#' \code{"squares"}: The 'squares' richness estimator of J. Alroy (2018).
#'
#' \code{"chao2"}: The Chao2 estimator for incidence-based data.
#' 
#' \code{"SCOR"}: The Sum Common Species Occurrence rate of Hannisdal et al. (2012). This method will only be calculated if the occurrence entries (vector)
#'  a collection vector is provided (see examples). 
#' 
#' @param x either a \code{character} vector of occurrences or a table of counts (\code{matrix}).
#' @param method (\code{character}): The type of metric that is to be calculated. The default value (\code{NULL}) will calculate all implemented metrics. 
#' @param samp (\code{vector}): Only applicable if \code{x} is vector of occurrence entries. The id of the sampling units. A necessary variable for SCOR. 
#' 
#' @return A named numeric vector.
#' @section References:
#' Alroy, J. 2018. Limits to species richness in terrestrial communities. Ecology Letters.
#' 
#' Hannisdal, B., Henderiks, J., & Liow, L. H. (2012). Long-term evolutionary and ecological responses of calcifying phytoplankton to changes in atmospheric CO2. Global Change Biology, 18(12), 3504–3516. https://doi.org/10.1111/gcb.12007
#'
#' Jost, L. (2006). Entropy and diversity. Oikos, 113, 363–375. https://doi.org/10.1111/j.2006.0030-1299.14714.x
#' @examples
#' # the coral data
#'   data(corals)
#' 
#' # Pleistocene subset
#'   plei <- corals[corals$stg==94,]
#' 
#' # calculate everything
#'   pleiIndex<-indices(plei$genus, plei$coll)
#' 
#' 
#' @export
indices <- function(x, samp=NULL, method=NULL){
	
	if("SCOR"%in%method & is.null(samp)) stop("If you want to calculate SCOR, you must provide a collection vector.")

	# which Hill numbers should be returned
	if(!is.null(method)){
		listMet<-strsplit(method, "hill")
		qs <- unlist(lapply(listMet, function(y) as.numeric(y[2])))
		qVec <- qs[is.finite(qs)]
	}else{
		qVec <- c(2)
	}

	# use occurrences - reformulate to 
	if(is.character(x) | is.factor(x)) ab <- matrix(tabulate(factor(x)),nrow=1)

	if(is.numeric(x) & !is.matrix(x)) ab<-matrix(x, nrow=1)

	# vectorized forms
	if(is.matrix(x)){
		if(!is.null(samp)) stop("You cannot specifiy a collection variable if the input data are already forming a table of sites.")
		ab <- x
	} 

	# EVERYTHING FROM HERE ON REQUIRES matrix collection/species abundance
	final <- matrix(ncol=0, nrow=nrow(ab))

	# calculate the proportions for further use
	props <- t(apply(ab,1, function(y){
		y/sum(y)
	}))

	# simple richness
	if("richness" %in% method | is.null(method)){
		richness<-apply(ab, 1, function(y){
			return(sum(y!=0))

		})
		richness <- matrix(richness, nrow=length(richness))
		colnames(richness) <- "richness"
		final <- cbind(final, richness)

	}

	# simple Shannon Entropy
	if("shannon"%in%method | is.null(method)){
		shannon <- matrix(apply(props, 1, function(y){
			-sum(y*log(y), na.rm=T)
		}), ncol=1)
		colnames(shannon) <- "shannon"
		final <- cbind(final, shannon)
	}


	# Hill numbers
	if(length(qVec)>0){
		
		hill <- matrix(NA, nrow=nrow(ab), ncol=length(qVec))
		
		for(i in 1:length(qVec)){
			
			if(qVec[i]==1){

				hill[,i] <- exp(-apply(props*log(props),1, sum, na.rm=T))

			}else{
				raised <- props^qVec[i]
				raised[props==0] <-NA
	
				hill[,i] <- apply(raised, 1,sum, na.rm=T)^(1/1-qVec[i])

			}
		}

		colnames(hill) <- paste("hill", qVec, sep="")
		final <- cbind(final, hill)
	}

	if("dominance"%in%method | is.null(method)){
		dominance <- matrix(apply(props, 1, function(y){
			max(y, na.rm=T)

		}), ncol=1)
		colnames(dominance) <- "dominance"
		final <- cbind(final, dominance)
	}

	# squares method
	if("squares"%in%method | is.null(method)){
		squares <- matrix(apply(ab, 1, squaresMethod), ncol=1)
		colnames(squares) <- "squares"
		final <- cbind(final, squares)

	}

	# Chao2 method
	if("chao2"%in%method | is.null(method)){
		chao2 <- matrix(apply(ab, 1, Chao2Method), ncol=1)
		colnames(chao2) <- "chao2"
		final <- cbind(final, chao2)
	}

	# calculate SCOR
	if(("SCOR"%in%method | is.null(method)) & !is.null(samp)){
		# create a cross tablulation
		crossTab<-table(samp, x)
		class(crossTab) <- "matrix"
		
		# for every species
		pI <- apply(crossTab, 2, function(y){
			# number of samples ("sites") in which species i is present/number of samples 
			sum(y>0)/nrow(crossTab)
		})
		lambdaI <- -log(1-pI)
		SCOR  <- matrix(sum(lambdaI), ncol=1)
		# add to the rest
		colnames(SCOR) <- "SCOR"
		final <- cbind(final, SCOR)

	}

	if(nrow(final)==1){
		final <-final[1,]
	}else{
		rownames(final) <- rownames(x)

	}

	return(final)
	
}

squaresMethod<-function(ab){
	ab<- ab[ab!=0]
    S <- length(ab)
    s1 <- length(which(ab == 1))
    if (s1 == S)
            return(NA)
    return(S + s1^2 * sum(ab^2) / (sum(ab)^2 - s1 * S))
}

Chao2Method<-function(ab){
	ab <- ab[ab!=0]
    s <- length(ab)
    return(s + (s - 1) / s * length(which(ab == 1)) * (length(which(ab == 1)) - 1) / (2 *  (length(which(ab == 2)) + 1)))
}
	