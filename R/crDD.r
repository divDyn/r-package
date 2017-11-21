#' Diversity dynamics with classical rarefaction
#' 
#' This script calculates time series of diversity dynamics with iterated subsampling in every time slice.
#' 
#' 
#' The output variables are the same as in the divDyn function. Likely updates: c++ core, arguments: sd - list output. 
#' @param dat (data.frame): the data frame containing PBDB occurrences.
#' @param quota (num): the number of occurrences to rarefy (quota). Slices with lower number of occurrences will be dropped.
 
#' 
#' @param tax (char): variable  name of the occurring taxa (variable type: factor) - such as "occurrence.genus_name"
#' 
#' @param bin (char): variable name of the time slice numbers of the particular occurrences (variable type: int)- such as "slc" or whatever. Bin numbers should be in ascending order,can contain NA's, it can start from a number other than 1 and must not start with 0.
#' @param coll (char): variable name of the collection ids. 
#' @param iter (num): the number of resampling iterations.
#' @param noNAStart (bool): useful when the dataset does not start from bin No. 1. Then noNAStart=TRUE will cut the first part of the resulting table, 
#' 						so the first row will contain the estimates for the lowest bin number.
#' 
#' @param average (char): the way the counts are averaged after the trials are finished. Can be 'geom' or 'arit' meaning geometric or arithmetic averaging, respecitvely.
#' @param intactBins (num): Vector of those bin identifiers that must be included in the metric calculations, but should not be subsampled. Can be useful for including recent species list, without collection identifiers.
#' @examples
#'	# import data
#'	  data(scleractinia)
#'	  data(stages)
#'
#'	# calculate metrics of diversity dynamics
#'    dd <- divDyn(scleractinia, tax="genus", bin="slc")
#'	  rarefDD <- crDD(scleractinia, quota=100, iter=10, intactBins=95)
#'
#'	# plotting
#'	  library(geoscale)
#'	  geoscalePlot(stages$mid,dd$divRT, age.lim=c(260,0), 
#'	    data.lim=c(0,300), type="l", units="Period", vers="ICS2014", 
#'	    ts.width=0.1, tick.scale=50, boxes="Period",ts.col=F,
#'	    label="range-through diversity (genera)")
#'	  lines(stages$mid,rarefDD100$divRT, col="blue")
#'	  legend("topleft", legend=c("raw","rarefaction"),
#'	    col=c("black", "blue"), lwd=c(1,1))
#'
#' @export
crDD<-function (dat, quota, tax="genus", bin="slc", coll="collection_no",  iter=100, noNAStart, average="geom", intactBins=NULL){

#	#0. temporary variables
#		dat<-pb.orig
#		tax<-"clgen"
#		bin<-"slcFine"
#		coll<-"collection_no"
#		quota<-100
#		iter<-10
#		noNAStart<-F
	
	#1. preparation
	
		# the time slice value of the first bin in the analysed interval
			nStart<-min(dat[,bin])
			
		# the time slice value of the last bin in the analysed interval
			nEnd <- max(dat[,bin])
			
		# length of the new vectors
			nVectorLength <- nEnd-nStart+1	
			
		# displacement vector for the new variables: start filling from vector[1]
			nSubtract<- nStart-1
		
		#the vector of time slice numbers
			nTimeSlice<-c(rep(NA, nStart-1), nStart:nEnd)
			
			# simplify analysis matrix, omit as much data as possible
				dSubDat<-subset(dat, select=c(coll, tax, bin))
		
			
			#count genus occurrences of different species as one
				dSubDat<-unique(dSubDat)
		
		# Determine the number of unique occurrences per stage
			nOcs <- numeric(nVectorLength)
			for (i in nStart:nEnd) {
				tTemp <- table(dSubDat[,coll][dSubDat[,bin]==i])
				nOcs[i-nSubtract] <- sum(tTemp)
			}
							
		#the bins where the quorum is larger than nQuota
			bBinSQ<-nOcs>=quota
			
			bBinSQ<-c(rep(FALSE, nStart-1), bBinSQ)
			
		# bins which will be added in every subsampling trials
			dIntBin<-dSubDat[dSubDat[,bin]%in%intactBins,]
			
			
	#2. subsampling trials
		#count vector declaration
		
			#the taxon counts	
			mSIBTrials<-matrix(ncol=iter, nrow=nEnd)
			m1tTaxaTrials<-matrix(ncol=iter, nrow=nEnd)
			m3tTaxaTrials<-matrix(ncol=iter, nrow=nEnd)
			mPtTaxaTrials<-matrix(ncol=iter, nrow=nEnd)
			mU2tTaxaTrials<-matrix(ncol=iter, nrow=nEnd)
			mD2tTaxaTrials<-matrix(ncol=iter, nrow=nEnd)
			mUGfTaxaTrials<-matrix(ncol=iter, nrow=nEnd)
			mDGfTaxaTrials<-matrix(ncol=iter, nrow=nEnd)
			
			mThroughTrials<--matrix(ncol=iter, nrow=nEnd)
			mExtTrials<-matrix(ncol=iter, nrow=nEnd)
			mOriTrials<-matrix(ncol=iter, nrow=nEnd)
			mBoundaryCrosserTrials<-matrix(ncol=iter, nrow=nEnd)
			mRTTrials<-matrix(ncol=iter, nrow=nEnd)
		#the trials
		for(k in 1:iter)
		{
			#the empty database
				dTrialDat<-data.frame()
				
			#the subsample dataset (dTrialDat)
			for (i in nStart:nEnd)
			{
				
				#subsample only if the i-th bin has enough occurrences
				if (bBinSQ[i])
				{
					#dataset of the timeslice
					dSliceDat<-subset(dSubDat, dSubDat[,bin]==i)
					
					#a random subset corresponding to quota number of occurrences
					dTrialDatPlus<-dSliceDat[sample(1:length(dSliceDat[,tax]), quota, replace=F),]
					
				
					#concatenate the subset of data that corresponds to the time slice to the general database
					dTrialDat<-rbind(dTrialDat, dTrialDatPlus)
		
				
				}
			
			}
			
			# the intact bins
				if(nrow(dIntBin)>0){
					dTrialDat<-dTrialDat[!dTrialDat[bin,]%in%intactBins,]
					dTrialDat<-rbind(dTrialDat,dIntBin)
				}
				
			#calculate the counts of the subsample dataset
				dDivDyn<-divDyn(dTrialDat, tax, bin, noNAStart=F, inf=FALSE)
				
				#augment with NAs the end of the table				
				dDivDyn[((length(dDivDyn[,1])):nEnd), ]<-NA
				
			#the counts
				
				#export the counts	
				mSIBTrials[,k] <- dDivDyn$divSIB
				mRTTrials[,k] <- dDivDyn$divRT
				m1tTaxaTrials[,k]<-dDivDyn$t1
				m3tTaxaTrials[,k]<-dDivDyn$t3
				mPtTaxaTrials[,k]<-dDivDyn$tPart
				mU2tTaxaTrials[,k]<-dDivDyn$t2u
				mD2tTaxaTrials[,k]<-dDivDyn$t2d
				mDGfTaxaTrials[,k]<-dDivDyn$tGFd
				mUGfTaxaTrials[,k]<-dDivDyn$tGFu				
				mThroughTrials[,k]<-dDivDyn$tThrough
				mExtTrials[,k]<-dDivDyn$tExt
				mOriTrials[,k]<-dDivDyn$tOri
				mBoundaryCrosserTrials[,k]<-dDivDyn$divBC
			
			
			#loop counting
				
				cat(k, paste("of ", iter, " iterations, ", "\r",sep="") )
				flush.console()
		}
	
	#3. average the counts (using the geometric mean)
	if(average=="geom"){
		#the geometric meaning
		logm1t<-log(m1tTaxaTrials)
		logm1t[is.infinite(logm1t)]<-NA
		n1tTaxaGeo<- apply(logm1t, 1, mean, na.rm=T) 
		n1tTaxaGeo<- exp(n1tTaxaGeo)
		
		logm3t<-log(m3tTaxaTrials)
		logm3t[is.infinite(logm3t)]<-NA
		n3tTaxaGeo<- apply(logm3t, 1, mean, na.rm=T) 
		n3tTaxaGeo<- exp(n3tTaxaGeo)
		
		logmPt<-log(mPtTaxaTrials)
		logmPt[is.infinite(logmPt)]<-NA
		nPtTaxaGeo<- apply(logmPt, 1, mean, na.rm=T) 
		nPtTaxaGeo<- exp(nPtTaxaGeo)
		
		logmD2t<-log(mD2tTaxaTrials)
		logmD2t[is.infinite(logmD2t)]<-NA
		nD2tTaxaGeo<- apply(logmD2t, 1, mean, na.rm=T) 
		nD2tTaxaGeo<- exp(nD2tTaxaGeo)
		
		logmU2t<-log(mU2tTaxaTrials)
		logmU2t[is.infinite(logmU2t)]<-NA
		nU2tTaxaGeo<- apply(logmU2t, 1, mean, na.rm=T) 
		nU2tTaxaGeo<- exp(nU2tTaxaGeo)
		
		logmUGf<-log(mUGfTaxaTrials)
		logmUGf[is.infinite(logmUGf)]<-NA
		nUGfTaxaGeo<- apply(logmUGf, 1, mean, na.rm=T) 
		nUGfTaxaGeo<- exp(nUGfTaxaGeo)
		
		logmDGf<-log(mDGfTaxaTrials)
		logmDGf[is.infinite(logmDGf)]<-NA
		nDGfTaxaGeo<- apply(logmDGf, 1, mean, na.rm=T) 
		nDGfTaxaGeo<- exp(nDGfTaxaGeo)
		
		#for the per capita rates
		logmThrough<-log(mThroughTrials)
		logmThrough[is.infinite(logmThrough)]<-NA
		nThroughGeo<- apply(logmThrough, 1, mean, na.rm=T) 
		nThroughGeo<- exp(nThroughGeo)
		
		logmExt<-log(mExtTrials)
		logmExt[is.infinite(logmExt)]<-NA
		nExtGeo<- apply(logmExt, 1, mean, na.rm=T) 
		nExtGeo<- exp(nExtGeo)
		
		logmOri<-log(mOriTrials)
		logmOri[is.infinite(logmOri)]<-NA
		nOriGeo<- apply(logmOri, 1, mean, na.rm=T) 
		nOriGeo<- exp(nOriGeo)
	
	
		logmSIB<-log(mSIBTrials)
		logmSIB[is.infinite(logmSIB)]<-NA
		nSIBGeo<- apply(logmSIB, 1, mean, na.rm=T) 
		nSIB<- exp(nSIBGeo)
		
		logmBC<-log(mBoundaryCrosserTrials)
		logmBC[is.infinite(logmBC)]<-NA
		nBoundaryCrosserGeo<- apply(logmBC, 1, mean, na.rm=T) 
		nBoundaryCrosserGeo<- exp(nBoundaryCrosserGeo)
		
		logmRT<-log(mRTTrials)
		logmRT[is.infinite(logmRT)]<-NA
		nRTGeo<- apply(logmRT, 1, mean, na.rm=T) 
		nRTGeo<- exp(nRTGeo)
	}
	if(average=="arit"){
	#the geometric meaning
		n1tTaxaGeo<- apply(m1tTaxaTrials, 1, mean, na.rm=T) 
		
		n3tTaxaGeo<- apply(m3tTaxaTrials, 1, mean, na.rm=T) 
		
		nPtTaxaGeo<- apply(mPtTaxaTrials, 1, mean, na.rm=T) 
		
		nD2tTaxaGeo<- apply(mD2tTaxaTrials, 1, mean, na.rm=T) 
		
		nU2tTaxaGeo<- apply(mU2tTaxaTrials, 1, mean, na.rm=T) 
		
		nUGfTaxaGeo<- apply(logmUGf, 1, mean, na.rm=T) 
		
		nDGfTaxaGeo<- apply(mDGfTaxaTrials, 1, mean, na.rm=T) 
		
		#for the per capita rates
		nThroughGeo<- apply(mThroughTrials, 1, mean, na.rm=T) 
		
		nExtGeo<- apply(mExtTrials, 1, mean, na.rm=T) 
		
		nOriGeo<- apply(mOriTrials, 1, mean, na.rm=T) 
		
		nSIBGeo<- apply(mSIBTrials, 1, mean, na.rm=T) 
	
		nBoundaryCrosserGeo<- apply(mBoundaryCrosserTrials, 1, mean, na.rm=T) 
		
		nRTGeo<- apply(mRTTrials, 1, mean, na.rm=T) 
		
	}
	#4. calculate the rates and the metrics
		#Per Capita rates
			nFooteExt<- -log(nThroughGeo/(nThroughGeo+nExtGeo))
			nFooteOri<- -log(nThroughGeo/(nThroughGeo+nOriGeo))
		
		
		#Three-timer sampling completeness (Alroy, 2008)	
			n3tSampComp <- n3tTaxaGeo/(n3tTaxaGeo+nPtTaxaGeo)
		
			#Sampling completeness of the entire time series
				nTot3tSampCompGeo<- sum(n3tTaxaGeo, na.rm=T)/(sum(n3tTaxaGeo, na.rm=T)+sum(nPtTaxaGeo, na.rm=T))
						
		
		#Three-timer rates by Alroy (2008)
			#uncorrected:
			n3tExt <- log(nD2tTaxaGeo/n3tTaxaGeo) # two-timer/three-timer ratio (bottom)
			n3tOri <- log(nU2tTaxaGeo/n3tTaxaGeo) # two-timer/three-timer ration (top)
			
			#corrected:
				#extinction rates:
					thtSampCompNext <- c(n3tSampComp[2:nEnd],NA) # Sampling probability in subsequent bin (BIN5)
					nC3tExt <- n3tExt + log(thtSampCompNext)
					#nC3tExtGeo[nC3tExtGeo<0] <- 0 # omit negative values
				
				#origination rates:
					thtSampCompPrev <- c(NA, n3tSampComp[1:nEnd-1]) # Sampling probability in previous bin (BIN5)
					nC3tOri <- n3tOri + log(thtSampCompPrev)
					#nC3tOriGeo[nC3tOriGeo<0] <- 0 #omit negative values
	 
		#corrected sampled-in-bin diversity
			nCorrSIB<-nSIB*nTot3tSampCompGeo/n3tSampComp
		
		#Gap filler estimates (Alroy, 2013)
			nGfExt<-log((nD2tTaxaGeo+nPtTaxaGeo)/(n3tTaxaGeo+nPtTaxaGeo+nUGfTaxaGeo))
			nGfOri<-log((nU2tTaxaGeo+nPtTaxaGeo)/(n3tTaxaGeo+nPtTaxaGeo+nDGfTaxaGeo))
	

	 
		
	
	
	#5. return the results
		#augment this variable
		nOcs<-c(rep(NA, nStart-1), nOcs[1:nVectorLength])
		
		
		dResults<-cbind(bin=nTimeSlice, occs=nOcs, t1=n1tTaxaGeo, t2d=nD2tTaxaGeo,
						t2u=nU2tTaxaGeo, t3=n3tTaxaGeo,  tPart=nPtTaxaGeo, tExt=nExtGeo, tOri=nOriGeo, divBC=nBoundaryCrosserGeo,
						tThrough=nThroughGeo,divRT=nRTGeo,
						extPC=nFooteExt, oriPC=nFooteOri, samp3t=n3tSampComp,  
						ext3t=n3tExt,ori3t=n3tOri, extC3t=nC3tExt,  oriC3t=nC3tOri,  divSIB=nSIB, 
						divCSIB=nCorrSIB, extGF=nGfExt, oriGF=nGfOri)
		
		dResults<-as.data.frame(dResults, stringsAsFactors=F)
		#want to see the NA's at the beginnning? (when the time series does not start with bin 1)
			if (missing(noNAStart)) {}
			else
			{
				if (noNAStart==TRUE)
				{
					dResults<-dResults[nStart:nEnd,]
				}
				
				if (noNAStart!=TRUE & noNAStart!=FALSE)
				{
					print("You have entered an invalid argument for noNAStart, no cropping will occurr in the final table.")
				}
			}	
		
		
		#return
		return(dResults)

}

