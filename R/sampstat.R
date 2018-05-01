#' Basic sampling summary for every bin
#'
#' This function will return the basic sampling summaries of a dataset
#'
#' A number of sampling related variables are calculated with the function.
#'
#'	occs: The number of occurrences in each time bin.
#'
#'	colls: The number of collections in each time bin.
#'
#'	refs: The number of references in each time bin.
#'
#'	SIBs: The number of Sampled-In-Bin taxa in each time bin.
#'
#'	occ1: The number of taxa in each time bin, that occur in only 1 collection.
#'
#'	ref1: The number of taxa in each time bin, that occur in only 1 reference.
#'
#'	occ2: The number of taxa in each time bin, that occur in exactly 2 collections.
#'
#'	ref2: The number of taxa in each time bin, that occur in exactly 2 references.
#'
#'	dom: Dominance, the proportion of occurrences in the time bin that belong to the most frequent taxon.
#'
#'	u: Good's u, coverage estimator based on the number of single-collection taxa (occ1).
#'
#'	uPrime: Good's u, coverage estimator based on the number of single-reference taxa (ref1).
#'
#'	chao1occ: Chao1 extrapolation estimator, based on the the number of single-collection and two-collection taxa (occ1).
#'
#'	chao1ref: Chao1 extrapolation estimator, based on the the number of single-reference and two-reference taxa (occ2).
#'
#' @param dat (data.frame): the occurrence dataset.
#' @param tax (char): the column name of taxon names.
#' @param bin (char): the column name of bin names.
#' @param coll (char): the column name of collection numbers.
#' @param ref (char): the column name of reference numbers.
#'
#' @examples
#'	data(scleractinia)
#'	# slice-specific sampling
#'	basic<-sampstat(scleractinia, tax="genus", bin="slc")
#'
#' # subsampling diagnostic
#'  subStats<-subsample(scleractinia, method="cr", tax="genus", FUN=sampstat, bin="slc", q=100,noNAstart=F)
#'	
#' @export
sampstat <- function(dat, tax="occurrence.genus_name", bin="slc", coll="collection_no", ref="occurrence.reference_no", noNAstart=FALSE){

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
	return(all)
}



