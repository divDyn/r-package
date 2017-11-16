#' FAD - LAD matrix from occurrence data
#' 
#' Function to generate range data from an occurrence dataset.
#' 
#' The following variables are produced:
#'
#' taxon: the name of the taxon.
#'
#' FAD: first appearance date in the given time slices.
#'
#' LAD: last appearance date in the given time slices.
#'
#' range: the range of the taxo in the given time slices
#'
#' FAD_DATE: optionally, first appearance dates in numeric ages.
#'
#' LAD_DATE: optionally, last appearance dates in numeric ages.
#'
#' DUR: optionally, the duration of taxa in numeric ages. It is given for single-interval taxa as well assuming that the taxa's range span over the entire time slice.
#'
#' nThrough: the number of through-ranging taxa
#'
#' @param dat (data.frame): Occurrence data.
#' @param tax (char): the column name of taxon names.
#' @param bin (char): the column name of bin names.
#' @param maxdate (char): optional, the column names of maximum radiometric dates of occurrences.
#' @param mindate (char): optional, the column names of minimum radiometric dates of occurrences.
#' @examples
#'	data(scleractinia)
#'	fl <- fadLad(scleractinia, tax="genus", bin="slc", maxdate="ma_max", mindate="ma_min")
#' @export
fadLad<-function(dat, tax, bin, maxdate=F, mindate=F)
#dat<-scleractinia
#tax<-"genus"
#bin<-"slc"
#maxdate<-"ma_max"
#mindate<-"ma_min"
#
{
	#if there are ages
	if(maxdate!=F & mindate!=F)
	{
		cTaxa<-levels(factor(dat[,tax]))
		
		mFadLad<-matrix(NA, nrow=length(cTaxa), ncol=6)
		dFadLad<-as.data.frame(mFadLad)
		colnames(dFadLad)<-c("FAD", "LAD", "range", "FAD_DATE", "LAD_DATE", "DUR")
		#rownames(dFadLad)<-cTaxa
		
		
		for(i in 1:length(cTaxa))
		{
			nTax<-dat[dat[,tax]==cTaxa[i], bin]
			
			nTaxMaxDate<-dat[dat[,tax]==cTaxa[i], maxdate]
			nTaxMinDate<-dat[dat[,tax]==cTaxa[i], mindate]
			
			nTaxMaxDate<-nTaxMaxDate[!is.na(nTaxMaxDate)]
			nTaxMinDate<-nTaxMinDate[!is.na(nTaxMinDate)]
			
			dFadLad[i,1]<-min(nTax, na.rm=T)
			dFadLad[i,2]<-max(nTax, na.rm=T)
			dFadLad[i,3]<-max(nTax, na.rm=T)-min(nTax, na.rm=T)
			if(length(nTaxMaxDate)>0){
				dFadLad[i,4]<-max(nTaxMaxDate, na.rm=T)
			}
			if(length(nTaxMinDate)>0){
				dFadLad[i,5]<-min(nTaxMinDate, na.rm=T)
			}
			if(length(nTaxMaxDate)>0 & length(nTaxMinDate)>0){
				dFadLad[i,6]<-max(nTaxMaxDate, na.rm=T)-min(nTaxMinDate, na.rm=T)
			}
			
		}
		taxon<-cTaxa
		dFadLad<-cbind(taxon, dFadLad)
		return(dFadLad)
	}else{
		cTaxa<-levels(factor(dat[,tax]))
		
		mFadLad<-matrix(NA, nrow=length(cTaxa), ncol=3)
		dFadLad<-as.data.frame(mFadLad)
		colnames(dFadLad)<-c("FAD", "LAD", "range")
		
		for(i in 1:length(cTaxa))
		{
			nTax<-dat[dat[,tax]==cTaxa[i], bin]
			
			dFadLad[i,1]<-min(nTax, na.rm=T)
			dFadLad[i,2]<-max(nTax, na.rm=T)
			dFadLad[i,3]<-max(nTax, na.rm=T)-min(nTax, na.rm=T)
			
		}
		taxon<-cTaxa
		dFadLad<-cbind(taxon, dFadLad, stringsAsFactors=F)
		return(dFadLad)
	
	}
	
	
}
