#' Basic sampling summary for every bin
#'
#' This function will return the basic sampling summaries of a dataset
#'
#' A number of sampling related variables are calculated with the function.
#'
#' By setting \code{total} to \code{FALSE} (default), the following results are output:
#'
#'	\code{occs}: The number of occurrences in each time bin.
#'
#'	\code{colls}: The number of collections in each time bin.
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
#'	\code{dom}: Dominance, the proportion of occurrences in the time bin that belong to the most frequent taxon.
#'
#'	\code{u}: Good's u, coverage estimator based on the number of single-collection taxa (occ1).
#'
#'	\code{uPrime}: Good's u, coverage estimator based on the number of single-reference taxa (ref1).
#'
#'	\code{chao1occ}: Chao1 extrapolation estimator, based on the the number of single-collection and two-collection taxa (occ1).
#'
#'	\code{chao1ref}: Chao1 extrapolation estimator, based on the the number of single-reference and two-reference taxa (occ2).
#'
#' Setting the \code{total} to \code{TRUE} will lead to the following results
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
#'  \code{gappiness}: The proportion of sampling gaps in the ranges of the taxa §. 
#'
#' @param dat \code{(data.frame)}: The occurrence dataset.
#' @param tax \code{(character)}: The column name of taxon names.
#' @param bin \code{(character)}: The column name of bin names.
#' @param coll \code{(character)}: The column name of collection numbers. (optional)
#' @param ref \code{(character)}: The column name of reference numbers. (optional)
#' @param total \code{(logical)}:  If \code{FALSE}, then the function will provide sampling statistics for the individual bins. \code{TRUE} on the other hand, will provide them for the total dataset.
#' @param noNAStart (logical) Useful when the dataset does not start from bin no. 1, but positive integer bin numbers are provided. Then \code{noNAStart=TRUE} will cut the first part of the resulting table, so the first row will contain the estimates for the lowest bin number. In case of positive integer bin identifiers, and if \code{noNAStart=FALSE}, the index of the row will be the bin number. 
#' @param duplicates \code{(logical)}: The function will check whether there are duplicate occurrences (multiple species/genera). When set to \code{NULL}, nothing will happen, but the function will notify you if duplicates are present. If set to \code{TRUE}, the function will not do anything with these, if set to \code{FALSE}, the duplicates will be omitted. 
#' 
#' @examples
#'	data(corals)
#'	# slice-specific sampling
#'	basic<-sampstat(corals, tax="genus", bin="slc")
#'
#' # subsampling diagnostic
#'  subStats<-subsample(corals, method="cr", tax="genus", FUN=sampstat, 
#'    bin="slc", q=100,noNAStart=FALSE)
#'	
#' @export
sampstat <- function(dat, tax="genus", bin="slc", coll=NULL, ref=NULL, total=FALSE, noNAStart=FALSE, duplicates=NULL){
	datUni<- unique(dat[c(tax, bin, coll, ref)])
	if(nrow(datUni)!=nrow(dat)){
		if(is.null(duplicates)){
			if(!is.null(coll) | !is.null(ref))	message("The database contains duplicate occurrences (multiple species/genus).") 
		}else{
			if(!duplicates & (is.null(coll)& is.null(ref))) stop("duplicates=FALSE implies that a 'coll' or 'ref' value is present!")
			if(!duplicates){
				dat <- datUni
			}
		}
		
	}
	if(total){
		nOcc<- nrow(dat)
		nTaxa <- length(levels(factor(dat[,tax])))
		nBin <- length(levels(factor(dat[,bin])))
		if(!is.null(coll))	nColl <- length(levels(factor(dat[,coll])))
		if(!is.null(ref))  nRef <- length(levels(factor(dat[,ref])))
		
		# Paul's sampling statistic
		binVar<-dat[,bin]
		taxVar<-dat[,tax]
		
		pSamp<-taxPS(taxVar, binVar)
		
		all<-data.frame(bins=nBin, occs=nOcc, taxa=nTaxa,colls=nColl,refs=nRef, gappiness=pSamp)
	}
	
	if(!total){
		
		# final output variables depend on what was entered
		needed <- NULL
		# factorize everything
		binVar<-dat[,bin]
		
		# note if binVar is a positive integer
		uBin <- unique(binVar)
		uBin <- uBin[!is.na(uBin)]
		if(sum(uBin%%1)==0 & sum(uBin<1)==0){
			posint <- TRUE
		}else{
			posint <- FALSE
		}
		
		
		if(!is.null(coll)) collVar<-as.integer(factor(dat[,coll]))
		taxVar<-as.integer(factor(dat[,tax]))
		
		if(!is.null(ref)) refVar<-as.integer(factor(dat[,ref]))
		
		# the number of occurrences
		occs<-table(binVar)
		needed<- c(needed, "occs")
		
		# the number of collections)
		if(!is.null(coll)){
			colls<-table((binVar[!duplicated(paste(binVar, collVar))]))
			needed<- c(needed, "colls")
		}else{
			colls<-rep(NA, length(occs))
		}
		
	#	table(unique(dat[,c(bin, coll)])[,bin])
		
		
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
		
		rows<-1:nrow(dat)
		
		# calculate taxon counts
		tot<-tapply(INDEX=binVar, X=rows, FUN=function(x){
			# the basic stuff
			tVar<-taxVar[x]
			
			# no collection occ table
			occTable<-table(tVar)
			
			if(!is.null(ref)){
				rVar<-refVar[x]
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
				cVar<-collVar[x]
				occTable<-table(unique(cbind(cVar, tVar))[, "tVar"])
			}
			
			
			# single collection (occurrence) taxa
			o1<-sum(occTable==1)
			
			o2<-sum(occTable==2)
		
			
			if(length(occTable)>0){
				dom<-max(occTable)/length(x)
			}else{
				dom<-NA
			}
			
			return(as.numeric(c(o1, p1, o2, p2, dom)))
		})
		
		occ1 <- unlist(lapply(tot, function(x)x[[1]]))
		ref1 <- unlist(lapply(tot, function(x)x[[2]]))
		occ2 <- unlist(lapply(tot, function(x)x[[3]]))
		ref2 <- unlist(lapply(tot, function(x)x[[4]]))
		dom <- unlist(lapply(tot, function(x)x[[5]]))
		needed <- c(needed, "occ1", "occ2", "dom")
		
		if(!is.null(ref)){ 
			needed <- c(needed, "ref1", "ref2")
		}
		
		# Good's u (occs)
		u<-1-(occ1[names(occs)]/occs)
		# Chao' estime
		chao1occ<-SIBs+(occ1[names(occs)]^2)/(2*occ2[names(occs)])
		needed <- c(needed, "u","chao1occ")
		
		# Alroy's u' (references)
		if(!is.null(ref)){
			uPrime<-1-(ref1[names(occs)]/occs)
			chao1ref<-SIBs+(ref1[names(occs)]^2)/(2*ref2[names(occs)])
			needed <- c(needed, "uPrime","chao1ref")
		}else{
			uPrime<-rep(NA, length(u))
			chao1ref<-rep(NA, length(u))
		
		}
		occ1 <- occ1[names(occs)]
		ref1 <- ref1[names(occs)]
		occ2 <- occ2[names(occs)]
		ref2 <- ref2[names(occs)]
		dom <- dom[names(occs)]
		
		# concatenate everything (also emtpy variables)
		all<-cbind(
			occs, 
			colls, 
			refs, 
			SIBs, 
			occ1, 
			ref1, 
			occ2, 
			ref2, 
			dom,
			u,
			uPrime, 
			chao1occ, 
			chao1ref)
		rownames(all)<-names(occs)
		
		all<- all[, needed]
		
		# if the entries are positive integers, coerce indexing!
		if(posint){
			place<- matrix(NA, ncol=ncol(all), nrow=max(as.numeric(names(occs))))
			rownames(place)<-1:nrow(place)
			colnames(place)<-colnames(all)
			place[as.numeric(names(occs)),]<-all
			all<-place
		}
		
		if(noNAStart){
			firstVal<-min(binVar, na.rm=T)
			all<-all[firstVal:length(occs),]
			
		}
		all<-as.data.frame(all)
	}
	
	
	return(all)
}



taxPS<-function(taxVar, binVar)
{
	#taxa sampled
	cTax<-levels(factor(taxVar))
	
	#repeat procedure for all taxa
	sap<-sapply(cTax, function(x){
		#taxon specific dataset
		nSamp<-sort(unique(binVar[taxVar==x]))
		
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

