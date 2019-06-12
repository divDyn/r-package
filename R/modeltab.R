#' Origination/extinction response table for statistical modelling.
#' 
#' This function takes an occurrence dataset and reformats it to a table that can be used as input for logistic models.
#' 
#' Every entry in the output table corresponds to one cell in the \code{bin}/\code{tax} matrix. This function omits duplicates and concatenates two \code{logical} vectors (response variables) to the occurrence dataset:  
#' The \code{ori} vector is \code{TRUE} in the interval when the taxon first appeared, and \code{FALSE} in all others. The \code{ext} vector is \code{TRUE} in the interval the taxon appeared for the last time, and \code{FALSE} in the rest.
#' 
#' @param x \code{(data.frame)} Fossil occurrence table.
#' 
#' @param tax \code{(character)} Variable name of the occurring taxa (variable type: \code{factor} or \code{character} - such as \code{"genus"})
#' 
#' @param bin \code{(character)} Variable name of the bin numbers of the occurrences. This variable should be \code{numeric} and should increase as time passes by (use negative values for age estimates). The current version only supports discreet, non-negative integer interval numbers.
#' 
#' @param taxvars \code{(character)} Taxon-specific column names of the variables that should be included in the output table. Only one entry/taxon is used, make sure that the data are clean.
#' 
#' @param rt \code{(logical)} Should the range-through assumption be applied within the function? If set to \code{TRUE} then missing occurrences will be interpolated with \code{FALSE} values in both the \code{ext} and \code{ori} variables. .
#' 
#' @param singletons \code{(logical)} Should single-interval taxa be removed from the final table? This is recommended, as it is impossible to get a \code{FALSE} response for these taxa. 
#' @rdname modeltab
#' @examples
#' # load necessary data
#' data(corals)
#' # simple table
#' modTab<-modeltab(corals, bin="stg", tax="genus", taxvars=c("ecology", "family"))
#'
#' @export
modeltab<-function(x, tax, bin, taxvars=NULL,  rt=FALSE, singletons=FALSE){
#	x <- corals[corals$stg!=95,]
#	bin<- "stg"
#	tax <- "genus"
#	taxvars <- "ecology"

	# sub dataset
	# the cell-variable
	unVar<- paste(x[, tax, drop=TRUE], x[, bin, drop=TRUE], sep="")
	subDat <- x[,c(tax, bin, taxvars)]
	subDat <- subDat[!duplicated(unVar),]

	# omit NAs
	bNeed<- !(is.na(subDat[,tax, drop=TRUE]) | is.na(subDat[,bin, drop=TRUE]))
	# taxon vars
	subDat<- subDat[bNeed,]

	# the data are complied - reorder
	subDat <- subDat[order(subDat[,tax, drop=TRUE], subDat[, bin, drop=TRUE]),]
	

	# filter singleton taxa
	if(!singletons){
		taxInd <- as.numeric(factor(subDat[, tax, drop=TRUE]))

		#if two 1 follow each other that indicates that a singleton is present
		diffInd <- c(1,diff(taxInd))
		diffShift<-c(diffInd[2:length(diffInd)], 1)

		# subset
		subDat<- subDat[!(diffInd & diffShift),]

	}

	# for extinctions
	rows <- 1:nrow(subDat)

	# for every taxon
	respInd<-tapply(INDEX=subDat[, tax], X=rows, function(w){
		# sampled bins
		bins<- subDat[w, bin, drop=TRUE]
		# return the FAD,LAD row numbers
		c(w[bins==min(bins)], w[bins==max (bins)])
	})

	# FAD-LAD row numbers 
	numResp <- matrix(unlist(respInd), ncol=2, byrow=T)

	# construct response variables
	ori <- rep(FALSE, length(rows))
	ext <- rep(FALSE, length(rows))
	ori[numResp[,1]] <- TRUE
	ext[numResp[,2]] <- TRUE
		

	# the model table to be returned
	datMod <- cbind( ori,ext, subDat)

	# if range interpolation is allowed
	if(rt){
		occurrence <- rep(TRUE, length(ori))
		datMod <- cbind(occurrence, datMod)

		# the rows of original data, gaps are NAs
		indInNew <- tapply(subDat[, tax], X=rows, function(w){
			bins <- subDat[w, bin, drop=TRUE]
			ran <- (range(bins)[1]:range(bins)[2])

			one <- rep(NA, length(ran))

			one[ran%in%bins] <- w
			return(one)
		})

		# in a single vector
		indInNew <- unlist(indInNew)

		# use this to subset the sampled part
		datMod<-datMod[indInNew,]
		datMod$occurrence[is.na(datMod$occurrence)] <- FALSE
		# fill in the obvious gaps
		# taxon name
		datMod[, tax] <- fill(datMod[, tax, drop=TRUE])
		# bin - increment by one
		datMod[, bin] <- fill(datMod[, bin, drop=TRUE], inc=1)

		# extinction/origination response
		datMod[, "ext"] <- fill(datMod[, "ext", drop=TRUE])
		datMod[, "ori"] <- fill(datMod[, "ori", drop=TRUE], forward=FALSE)

		if(!is.null(taxvars)){
			# the taxon variables
			for(i in 1:length(taxvars)){
				datMod[, taxvars[i]]<- fill(datMod[, taxvars[i]])
			}
		}
	}

	# final object
	rownames(datMod) <- NULL
	return(datMod)

}

# It is very easy to append additional variables to this table. Such as the output of affinity()