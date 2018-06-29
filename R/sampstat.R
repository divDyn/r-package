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
#' @param dat \code{(data.frame)}: the occurrence dataset.
#' @param tax \code{(character)}: the column name of taxon names.
#' @param bin \code{(character)}: the column name of bin names.
#' @param coll \code{(character)}: the column name of collection numbers.
#' @param ref \code{(character)}: the column name of reference numbers.
#' @param total \code{(logical)}:  If \code{FALSE}, then the function will provide sampling statistics for the individual bins. \code{TRUE} on the other hand, will provide them for the total dataset.
#' @param duplicates \code{(logical)}: The function will check whether there are duplicate occurrences (multiple species/genera). When set to \code{NULL}, nothing will happen, but the function will notify you if duplicates are present. If set to \code{TRUE}, the function will not do anything with these, if set to \code{FALSE}, the duplicates will be omitted. 
#'
#' @examples
#'	data(corals)
#'	# slice-specific sampling
#'	basic<-sampstat(corals, tax="genus", bin="slc")
#'
#' # subsampling diagnostic
#'  subStats<-subsample(corals, method="cr", tax="genus", FUN=sampstat, bin="slc", q=100,noNAstart=F)
#'	
#' @export
sampstat <- function(dat, tax="genus", bin="slc", coll="collection_no", ref="reference_no", total=FALSE, noNAstart=FALSE, duplicates=NULL){
	datUni<- unique(dat[c(tax, bin, coll, ref)])
	if(nrow(datUni)!=nrow(dat)){
		if(is.null(duplicates)){
			message("The database contains duplicate occurrences (multiple species/genus).") 
		}else{
			if(!duplicates){
				dat <- datUni
			}
		}
		
	}
	if(total){
		nOcc<- nrow(dat)
		nTaxa <- length(levels(factor(dat[,tax])))
		nBin <- length(levels(factor(dat[,bin])))
		nColl <- length(levels(factor(dat[,coll])))
		nRef <- length(levels(factor(dat[,ref])))
		
		# Paul's sampling statistic
		binVar<-dat[,bin]
		taxVar<-dat[,tax]
		
		pSamp<-taxPS(taxVar, binVar)
		
		all<-data.frame(bins=nBin, occs=nOcc, taxa=nTaxa,colls=nColl,refs=nRef, gappiness=pSamp)
	}
	
	if(!total){
		# factorize everything
		binVar<-as.integer(dat[,bin])
		collVar<-as.integer(factor(dat[,coll]))
		taxVar<-as.integer(factor(dat[,tax]))
		refVar<-as.integer(factor(dat[,ref]))
		
		# the number of collections)
		collections<-tabulate((binVar[!duplicated(paste(binVar, collVar))]))
		
	#	table(unique(dat[,c(bin, coll)])[,bin])
		
		# the number of occurrences
		occs<-tabulate(binVar)
		
		# the number of references
		refs<-tabulate((binVar[!duplicated(paste(binVar, refVar))]))
		
		# number of sampled taxa
		SIBs<-tabulate((binVar[!duplicated(paste(taxVar, binVar))]))
		
		rows<-1:nrow(dat)
		
		# calculate taxon counts
		tot<-tapply(INDEX=binVar, X=rows, FUN=function(x){
			# the basic stuff
			cVar<-collVar[x]
			tVar<-taxVar[x]
			rVar<-refVar[x]
			
			occTable<-table(unique(cbind(cVar, tVar))[, "tVar"])
			refTable<-table(unique(cbind(rVar, tVar))[, "tVar"])
			
			# single collection (occurrence) taxa
			o1<-sum(occTable==1)
			# single reference taxa
			p1<-sum(refTable==1)
			
			o2<-sum(occTable==2)
			# doubletond
			p2<-sum(refTable==2)
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
		
		# Good's u (occs)
		u<-1-(occ1[as.character(1:length(occs))]/occs)
		
	
		# Alroy's u' (references)
		uPrime<-1-(ref1[as.character(1:length(occs))]/occs)
		
		chao1occ<-SIBs+(occ1[as.character(1:length(occs))]^2)/(2*occ2[as.character(1:length(occs))])
		chao1ref<-SIBs+(ref1[as.character(1:length(occs))]^2)/(2*ref2[as.character(1:length(occs))])
		
		# concatenate everything
		all<-cbind(
			occs, 
			collections, 
			refs, 
			SIBs, 
			occ1[as.character(1:length(occs))], 
			ref1[as.character(1:length(occs))], 
			occ2[as.character(1:length(occs))], 
			ref2[as.character(1:length(occs))], 
			dom[as.character(1:length(occs))],
			u,
			uPrime, 
			chao1occ, 
			chao1ref)
		colnames(all)[c(2,3,5:9)]<-c("colls","refs", "occ1", "ref1", "occ2", "ref2", "dom")
		rownames(all)<-1:length(occs)
		
		if(noNAstart){
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

