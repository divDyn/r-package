#' Sampling statistics and diversity indices in every bin
#'
#' This function will return the basic sampling summaries of a dataset
#'
#' Secondary function of the package that calculates a number of sampling related variables and diversity estimators for each bin. 
#' In contrast to the (\code{\link{divDyn}}) function, the bins are treated independently in this function.
#' The function also returns the maximum subsampling quota for OxW subsampling 
#' (\code{\link{subtrialOXW}}) with a given \code{xexp} value.
#'
#' By setting \code{total} to \code{FALSE} (default), the following results are output:
#'
#'	\code{occs}: The number of occurrences in each time bin.
#'
#'	\code{colls}: The number of collections in each time bin.
#'
#'  \code{xQuota}: The maximum quota for OxW subsampling (\code{\link{subtrialOXW}}) with the given \code{xexp} value. 
#' 	The number of occurrences in each collection is tabulated, and is raised to the power of \code{xexp}. 
#' 	The \code{xQuota} value is the sum of these values across all collections in a time slice.
#'
#'	\code{refs}: The number of references in each time bin.
#'
#'	\code{SIBs}: The number of Sampled-In-Bin taxa in each time bin.
#'
#'	\code{occ1}: The number of taxa in each time bin, that occur in only 1 collection.
#'
#'	\code{ref1}: The number of taxa in each time bin, that occur in only 1 reference.
#'
#'	\code{occ2}: The number of taxa in each time bin, that occur in exactly 2 collections.
#'
#'	\code{ref2}: The number of taxa in each time bin, that occur in exactly 2 references.
#'
#'	\code{u}: Good's u, coverage estimator based on the number of single-collection taxa (occ1).
#'
#'	\code{uPrime}: Good's u, coverage estimator based on the number of single-reference taxa (ref1).
#'
#'	\code{chao1occ}: Chao1 extrapolation estimator, based on the the number of single-collection and two-collection taxa (occ1).
#'
#'	\code{chao1ref}: Chao1 extrapolation estimator, based on the the number of single-reference and two-reference taxa (occ2).
#'
#' @param x \code{(data.frame)}: The occurrence dataset.
#' @param tax \code{(character)}: The column name of taxon names.
#' @param bin \code{(character)}: The column name of bin names.
#' @param coll \code{(character)}: The column name of collection numbers. (optional)
#' @param ref \code{(character)}: The column name of reference numbers. (optional)
#' @param noNAStart (logical) Useful when the dataset does not start from bin no. 1, but positive integer bin numbers are provided. Then \code{noNAStart=TRUE} will cut the first part of the resulting table, so the first row will contain the estimates for the lowest bin number. In case of positive integer bin identifiers, and if \code{noNAStart=FALSE}, the index of the row will be the bin number. 
#' @param duplicates \code{(logical)}: The function will check whether there are duplicate occurrences (multiple species/genera). When set to \code{NULL}, nothing will happen, but the function will notify you if duplicates are present. If set to \code{TRUE}, the function will not do anything with these, if set to \code{FALSE}, the duplicates will be omitted. 
#' @param xexp (\code{numeric}): Argument of the OxW subsampling type (\code{\link{subtrialOXW}}).Setting this parameter to a valid numeric value will return the maximum quota for \code{xexp}.
#' @param indices (\code{logical}): Setting this value to \code{TRUE} will calculate all indices implemented in (\code{\link{indices}}). 
#' @examples
#'	data(corals)
#'	# slice-specific sampling
#'	basic <- binstat(corals, tax="genus", bin="stg")
#'
#' # subsampling diagnostic
#'  subStats <- subsample(corals, method="cr", tax="genus", FUN=binstat, 
#'    bin="stg", q=100,noNAStart=FALSE)
#'
#' # maximum quota with xexp
#'	more <- binstat(corals, tax="genus", bin="stg", coll="collection_no", xexp=1.4)
#'	
#' @export
binstat <- function(x, tax="genus", bin="stg", coll=NULL, ref=NULL, noNAStart=FALSE, duplicates=NULL, xexp=NULL, indices=FALSE){
	datUni<- unique(x[,c(tax, bin, coll, ref)])
	if(nrow(datUni)!=nrow(x)){
		if(is.null(duplicates)){
			if(!is.null(coll) | !is.null(ref))	message("The database contains duplicate occurrences (multiple species/genus).") 
		}else{
			if(!duplicates & (is.null(coll)& is.null(ref))) stop("duplicates=FALSE implies that a 'coll' or 'ref' value is present!")
			if(!duplicates){
				x <- datUni
			}
		}
	}

	# final output variables depend on what was entered
	needed <- "thebin"
	# factorize everything
	binVar<-x[,bin, drop=TRUE]
	
	# note if binVar is a positive integer
	uBin <- unique(binVar)
	uBin <- uBin[!is.na(uBin)]
	if(sum(uBin%%1)==0 & sum(uBin<1)==0){
		posint <- TRUE
	}else{
		posint <- FALSE
	}
	
	
	if(!is.null(coll)) collVar<-as.integer(factor(x[,coll, drop=TRUE]))
	taxVar<-as.integer(factor(x[,tax, drop=TRUE]))
	
	if(!is.null(ref)) refVar<-as.integer(factor(x[,ref, drop=TRUE]))
	
	# the number of occurrences
	occs<-table(binVar)
	thebin<-as.numeric(names(occs))
	needed<- c(needed, "occs")
	
	# the number of collections)
	if(!is.null(coll)){
		colls<-table((binVar[!duplicated(paste(binVar, collVar))]))
		needed<- c(needed, "colls")
		if(!is.null(xexp)){
			xQuota <- tapply(INDEX=x[,bin, drop=TRUE], X=x[, coll, drop=TRUE], function(y){
			#	y <- x[x[,bin]==85,coll]
				sum(table(y)^xexp)
			})
			needed<- c(needed, "xQuota")
			xQuota <- xQuota[names(occs)]
		}else{
			xQuota<-rep(NA, length(occs))
		}
	}else{
		colls<-rep(NA, length(occs))
		xQuota<-rep(NA, length(occs))
	}
	
#	table(unique(x[,c(bin, coll)])[,bin])
	
	
	# the number of references
	if(!is.null(ref)){ 
		refs<-table((binVar[!duplicated(paste(binVar, refVar))]))
		needed<- c(needed, "refs")
	}else{
		refs<-rep(NA, length(occs))
	}
	
	# number of sampled taxa
	SIBs<-table((binVar[!duplicated(paste(taxVar, binVar))]))
	needed<- c(needed, "SIBs")
	
	rows<-1:nrow(x)
	
	# calculate taxon counts
	cVar<-NULL
	tot<-tapply(INDEX=binVar, X=rows, FUN=function(w){
		# the basic stuff
		tVar<-taxVar[w]
		
		# no collection occ table
		occTable<-table(tVar)
		
		if(!is.null(ref)){
			rVar<-refVar[w]
			refTable<-table(unique(cbind(rVar, tVar))[, "tVar"])
			# single reference taxa
			p1<-sum(refTable==1)
			# doubletond
			p2<-sum(refTable==2)
		}else{
			p1 <- NA
			p2 <- NA
		}
		
		if(!is.null(coll)){
			cVar<-collVar[w]
		}
		
		
		# single collection (occurrence) taxa
		o1<-sum(occTable==1)
		
		o2<-sum(occTable==2)
	
		if(indices){
			
			# cVar defaults to NULL when 
			ind<-indices(x=as.character(tVar), cVar)
		}else{
			ind <- NULL
		}
		
		resThis <- as.numeric(c(o1, p1, o2, p2, ind))
		names(resThis) <- c("occ1", "ref1", "occ2", "ref2", names(ind))
		return(resThis)
	})
	# all returned var names
	outVars<-names(tot[[1]])
	# within indices
	indNames <- outVars[!outVars%in%c("occ1", "ref1", "occ2", "ref2")]
	# process the output
	tapRes<-matrix(nrow=length(tot), ncol=0)

	for(i in 1:length(outVars)){
		tapRes<- cbind(tapRes, unlist(lapply(tot, function(w)w[[i]])))
	}
	
	colnames(tapRes) <- outVars
	tapRes <- tapRes[names(occs),]
	needed <- c(needed, "occ1", "occ2")
	
	if(!is.null(ref)){ 
		needed <- c(needed, "ref1", "ref2")
	}
	
	# Good's u (occs)
	u<-1-(tapRes[names(occs), "occ1", drop=TRUE]/occs)
	# Chao' estime
	chao1occ<-SIBs+(tapRes[names(occs), "occ1", drop=TRUE]^2)/(2*tapRes[names(occs), "occ2", drop=TRUE])
	needed <- c(needed, "u","chao1occ")
	
	# Alroy's u' (references)
	if(!is.null(ref)){
		uPrime<-1-(tapRes[names(occs), "ref1", drop=TRUE]/occs)
		chao1ref<-SIBs+(tapRes[names(occs), "ref1", drop=TRUE]^2)/(2*tapRes[names(occs), "ref2", drop=TRUE])
		needed <- c(needed, "uPrime","chao1ref")
	}else{
		uPrime<-rep(NA, length(u))
		chao1ref<-rep(NA, length(u))
	}

	# concatenate everything (also emtpy variables)
	all<-cbind(
		thebin,
		occs, 
		colls, 
		refs, 
		SIBs, 
		tapRes, 
		u,
		uPrime, 
		chao1occ, 
		chao1ref, 
		xQuota)
	
	rownames(all)<-names(occs)
	
	all<- all[, c(needed, indNames)]
	
	# if the entries are positive integers, coerce indexing!
	if(posint){
		place<- matrix(NA, ncol=ncol(all), nrow=max(as.numeric(names(occs))))
		colnames(place)<-colnames(all)
		place[, "thebin"] <- 1:nrow(place)
		rownames(place)<- place[, "thebin"]
		place[as.numeric(names(occs)),]<-all
		all<-place
	}
	
	if(noNAStart){
		firstVal<-min(binVar, na.rm=T)
		all<-all[as.character(firstVal:max(as.numeric(names(occs)), na.rm=T)),]
		
	}
	all<-as.data.frame(all)
	colnames(all)[colnames(all)=="thebin"] <- bin



	return(all)
}

#' Occurrence database summary
#' 
#' The function calculates global statistics of the entire database
#' 
#' The function returns the following values.
#' 
#'  \code{bins}: The total number of bins sampled. 
#'
#'	\code{occs}: The total number of sampled occurrences.
#'
#'	\code{colls}: The total number of sampled collections.
#'
#'	\code{refs}:  The total number of sampled references.
#'
#'	\code{taxa}:  The total number of sampled taxa.
#'
#'  \code{gappiness}: The proportion of sampling gaps in the ranges of the taxa (without the range-endpoints).
#'
#' @param x \code{(data.frame)}: The occurrence dataset.
#' @param tax \code{(character)}: The column name of taxon names.
#' @param bin \code{(character)}: The column name of bin names.
#' @param coll \code{(character)}: The column name of collection numbers. (optional)
#' @param ref \code{(character)}: The column name of reference numbers. (optional)
#' @param duplicates \code{(logical)}: The function will check whether there are duplicate occurrences (multiple species/genera). When set to \code{NULL}, nothing will happen, but the function will notify you if duplicates are present. If set to \code{TRUE}, the function will not do anything with these, if set to \code{FALSE}, the duplicates will be omitted. 
#' @examples
#'	data(corals)
#'	  sumstat(corals, tax="genus", bin="stg", coll="collection_no", ref="reference_no")
#' @export
sumstat<- function(x, tax="genus", bin="stg", coll=NULL, ref=NULL, duplicates=NULL){
	# duplicate searching and omission part
	datUni<- unique(x[,c(tax, bin, coll, ref)])
	if(nrow(datUni)!=nrow(x)){
		if(is.null(duplicates)){
			if(!is.null(coll) | !is.null(ref))	message("The database contains duplicate occurrences (multiple species/genus).") 
		}else{
			if(!duplicates & (is.null(coll)& is.null(ref))) stop("duplicates=FALSE implies that a 'coll' or 'ref' value is present!")
			if(!duplicates){
				x <- datUni
			}
		}
	}
	
	nOcc<- nrow(x)
	nTaxa <- length(levels(factor(x[,tax, drop=TRUE])))
	nBin <- length(levels(factor(x[,bin, drop=TRUE])))
	if(!is.null(coll)){
		nColl <- length(levels(factor(x[,coll, drop=TRUE])))
	}else{
		nColl <- NA
	}
	if(!is.null(ref)){
		nRef <- length(levels(factor(x[,ref, drop=TRUE])))
	}else{
		nRef<-NA
	}
	
	# Paul's sampling statistic
	binVar<-x[,bin, drop=TRUE]
	taxVar<-x[,tax, drop=TRUE]
	
	pSamp<-taxPS(taxVar, binVar)
	
	all<-data.frame(bins=nBin, occs=nOcc, taxa=nTaxa,colls=nColl,refs=nRef, gappiness=pSamp)
	return(all)
}


taxPS<-function(taxVar, binVar)
{
	#taxa sampled
	cTax<-levels(factor(taxVar))
	
	#repeat procedure for all taxa
	sap<-sapply(cTax, function(w){
		#taxon specific dataset
		nSamp<-sort(unique(binVar[taxVar==w]))
		
		#exclude range endpoints from calculation
		if (length(nSamp)>2)
		{
			#sampled bins without endpoints
	    	nSamp2<-nSamp[2:(length(nSamp)-1)]
	    	observed<-length(nSamp2)
	    	
	    	#implied range without the endpoints
	    	expected<-length((nSamp[1]+1):(nSamp[length(nSamp)]-1))
	    	
			return(c(observed, expected))
	 
	    }else{
			return(c(0,0))
		}
			
		
	})
	
	#frequency
	ps<-sum(sap[1,])/sum(sap[2,])
	return(ps)
	
}

