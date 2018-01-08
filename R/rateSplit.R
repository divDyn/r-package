#' Test of rate split (selectivity testing)
#'	
#' This function will determine whether there are meaningful differences between the taxonomic rates in the individual time slices of two subsets of an occurrence database.
#'	
#' Splitting an occurrence database to its subsets secreases the amount of information passed to the rate calculations and therefore the precision of the individual estimates. In order to assess the subsets individually and compare them, it is advised to test whether they the split is meaningful, given the total data. Examples of this use can be found in Kiessling and Simpson (2011), Global Change Biology 17, 56-67 and Kiessling and Kocsis (2015), Paleobiology 41, 402-414.
#'	
#' The meaningfulness of the split is dependent on the estimate accurracy and the magnitude of the difference.
#' @param sel (character): Variable name to do the splitting of the dataset. Can have only two levels.
#'	
#' @param output (character): Either "simple" or "full". Simple returns the indices of the series where selectivity can be suggested. "full" returns a matrix of Akaike weights, or binomial probabilities.
#'	
#' @param method (character): Either "AIC", "binom" or "combine". The "AIC" method calculates the Akaike weights of the single and dual rate models. The "binom" method assumes a binomial error distribution of the counts that are necessary for the rate calculations. The "combine" method shows slices that pass both tests, the "AIC" being usually the stronger.
#'	
#' @param alpha (num): Threshold for the between meaningful and meaningless split. If method=="AIC", the value corresponds to the minimum weight value the dual model should have , by default it is 0.89, which corresponds to the likelihood ratio of 8. If method=="binom", the value corresponds to the alpha value of the binomial test.
#'	
#' @param AICc (logical): Only applicable for the "AIC" method. Toggles whether the small sample corrected AIC (AICc) should be used instead of the regular one.
#'	
#' @param rate (character): The rate metric. Currently only the per capita rates of Foote (2000) are available.
#' @param dat (data.frame): the data frame containing PBDB occurrences.
#'	 
#' @param tax (char): variable  name of the occurring taxa (variable type: factor) - such as "occurrence.genus_name"
#'	 
#' @param bin (char): variable name of the time slice numbers of the particular occurrences (variable type: int)- such as "slc" or whatever. Bin numbers should be in ascending order,can contain NA's, it can start from a number other than 1 and must not start with 0.
#' @examples
#'	# example with the coral dataset of Kiessling and Kocsis (2015)
#'	data(scleractinia)
#'	data(stages)
#'	
#'	# split by ecology
#'	  z<-scleractinia[scleractinia$ecology=="z",]
#'	  az<-scleractinia[scleractinia$ecology=="az",]
#'	
#'	# calculate diversity dynamics
#'	ddZ<-divDyn(z, tax="genus", bin="slc")
#'	ddAZ<-divDyn(az, tax="genus", bin="slc")
#'	
#'	# origination rate plot
#'	plotTS(stages, boxes="per", shading="series", xlim=54:95, 
#'	  ylab="raw per capita originations")
#'	lines(stages$mid, ddZ$oriPC, lwd=2, lty=1, col="blue")
#'	lines(stages$mid, ddAZ$oriPC, lwd=2, lty=2, col="red")
#'	legend("topright", inset=c(0.1,0.1), legend=c("z", "az"), 
#'	  lwd=2, lty=c(1,2), col=c("blue", "red"), bg="white")
#'	
#'	# The ratesplit function
#'	rs<-ratesplit(rbind(z, az), sel="ecology", tax="genus", bin="slc")
#'	rs
#'	
#'	# display selectivity with points
#'	# select the higher rates
#'	selIntervals<-cbind(ddZ$oriPC[rs$ori], ddAZ$oriPC[rs$ori])
#'	groupSelector<-apply(selIntervals, 1, function(x) x[1]<x[2])
#'	# draw the points
#'	points(stages$mid[rs$ori[groupSelector]], ddAZ$oriPC[rs$ori[groupSelector]],
#'	  pch=16, col="red", cex=2)
#'	points(stages$mid[rs$ori[!groupSelector]], ddZ$oriPC[rs$ori[!groupSelector]],
#'	  pch=16, col="blue", cex=2)
#'	
#'	
#' @export
ratesplit<-function(dat,  sel, tax="genus", bin="slc", rate="pc", method="AIC",AICc=T, na.rm=T, alpha=NULL, output="simple"){
#	 #dummy data:
# 	
#	 dat<-scleractinia
#	 tax<-"genus"
#	 bin<-"slc"
#	 cRate<-"nFooteOri"
#	 
#	
#	na.rm=T
#	
	# the types of entries
	entries<-unique(dat[,sel])
	
	# omit NAs automatically
	if(!na.rm){
		if(sum(is.na(entries))>0){
			return(NA)
		}
	}else{
		entries<-entries[!is.na(entries)]
		dat<-dat[!is.na(dat[,sel]),]
	}
	
	if(length(entries)>2){
		stop("More than two types of entries are present in the selector variable. Function was designed to operate with 2.")
	}
	 
	if(is.null(alpha)){
		if(method=="AIC") alpha <- 0.89
		if(method=="binom") alpha <- 0.05
	}
	 
	#data selection:
	dFirst<-subset(dat, dat[,sel]==entries[1])
	dSecond<-subset(dat, dat[,sel]==entries[2])
	dComb<-dat
	
	 
	# the diversity dynamics
	dd1<-divDyn(dFirst, tax=tax, bin=bin)
	dd2<-divDyn(dSecond, tax=tax, bin=bin)
	ddBoth<-divDyn(dComb, tax=tax, bin=bin)
	
	# combine them to matrices
	newOri <- cbind(ddBoth$tOri,dd1$tOri, dd2$tOri)
	newExt <- cbind(ddBoth$tExt,dd1$tExt, dd2$tExt)
	newBC<- cbind(ddBoth$tThrough+ddBoth$tExt,dd1$tThrough+dd1$tExt, dd2$tThrough+dd2$tExt)
	newThr<- cbind(ddBoth$tThrough, dd1$tThrough, dd2$tThrough)
	
	# method for the per capita rates
	if(rate=="pc"){
		fext<- cbind(ddBoth$extPC, dd1$extPC, dd2$extPC)     #Foote extinction rates based on total ranges
		fori<- cbind(ddBoth$oriPC,dd1$oriPC, dd2$oriPC)           #Foote origination rates based on total ranges
		
		if(method=="AIC"){
		
		
			# Log-likelihoods and AIC
			loglikExt <- newThr*log(newThr/newBC) + newExt*log(newExt/newBC) # log-likelihood extinction
			
			
		
			# Log-likelihoods and AIC
			loglikOri <- newThr*log(newThr/newBC) + newOri*log(newOri/newBC) # log-likelihood extinction
			
			
			# number of parameters
			K1 <- 1
			K2 <- 2
			
			# correction parameter for AIC
			cp <- newBC[,1] 
			
			if (AICc==F){
				# AIC test simple and corrected (specific for three columns)
				# simple AIC
				# extinction
				extAIC12 <- 2*K2-2*(loglikExt[,2]+loglikExt[,3])
				extAIC1 <- 2*1-2*(loglikExt[,1])
				
				# origination
				oriAIC12 <- 2*K2-2*(loglikOri[,2]+loglikOri[,3])
				oriAIC1 <- 2*1-2*(loglikOri[,1])
			}
			
			if(AICc==T){
				# corrected AIC (small samples)
				# extinction
				extAIC12 <- 2*K2-2*(loglikExt[,2]+loglikExt[,3])+(2*K2*(K2+1))/(cp-K2-1)
				extAIC1 <- 2*K1-2*(loglikExt[,1])+ (2*K1*(K1+1))/(cp-K1-1)
				
				# origination
				oriAIC12 <- 2*K2-2*(loglikOri[,2]+loglikOri[,3])+(2*K2*(K2+1))/(cp-K2-1)
				oriAIC1 <- 2*K1-2*(loglikOri[,1])+ (2*K1*(K1+1))/(cp-K1-1)
			}
			
			# difference
				zOri <- oriAIC1-oriAIC12
				zOri1 <- replace(zOri, zOri < 0, 0)
				zOri2 <- abs(replace(zOri, zOri > 0, 0))
				
				zExt <- extAIC1-extAIC12
				zExt1 <- replace(zExt, zExt < 0, 0)
				zExt2 <- abs(replace(zExt, zExt > 0, 0))
				
				# change to weights
				# extinctions
				wExt1 <- exp(-0.5*zExt1)/(exp(-0.5*zExt1)+exp(-0.5*zExt2)) # single rate preferred
				wExt2 <- exp(-0.5*zExt2)/(exp(-0.5*zExt1)+exp(-0.5*zExt2)) # dual rate preferred
			
				# extinctions
				wOri1 <- exp(-0.5*zOri1)/(exp(-0.5*zOri1)+exp(-0.5*zOri2)) # single rate preferred
				wOri2 <- exp(-0.5*zOri2)/(exp(-0.5*zOri1)+exp(-0.5*zOri2)) # dual rate preferred
			
			# likelihood ratio larger than alpha
				AkaikeSelExt <- which(wExt2 > alpha)
				AkaikeSelOri <- which(wOri2 > alpha)
				
				if(output=="simple"){
					res<-list(ext=AkaikeSelExt, ori=AkaikeSelOri)
					
					if(method=="combine"){
						aicRes<-res
					}
					
				}
				if(output=="full"){
					if(method=="combine"){
						stop("The full output option does not work with the combine method")
					}
					
					res<-data.frame(singleOri=wOri1, dualOri=wOri2,singleExt=wExt1, dualExt=wExt2)
					return(res)
				}
		}
		
		if(method=="binom"){
			
			# binomials 
			dualBC <- apply(newBC[,2:3], 1,  sum) # for boundary crossers
			dualExt <- apply(newExt[,2:3], 1, sum) 
			dualOri <- apply(newOri[,2:3], 1, sum) 
			
			# successes in vulnerable category
			kExt <- newExt[,2] 
			kOri <- newOri[,2] 
			
			nExt <- dualExt
			nOri <- dualOri
			
			# diversity proportion of vulnerable category
				pr <- newBC[,2]/dualBC # boundary crossers 
		
			# Binomial test (specific for three columns) - extinctions
				sig <- pbinom(kExt,nExt,pr,lower.tail=F)+dbinom(kExt,nExt,pr)
				sig2 <- pbinom(nExt-kExt,nExt,1-pr,lower.tail=F)+dbinom(nExt-kExt,nExt,1-pr)
				
				# which are significant
				binExt <- which(sig <= alpha | sig2 <= alpha)
			
			# Binomial test (specific for three columns) - originations
				sig <- pbinom(kOri,nOri,pr,lower.tail=F)+dbinom(kOri,nOri,pr)
				sig2 <- pbinom(nOri-kOri,nOri,1-pr,lower.tail=F)+dbinom(nOri-kOri,nOri,1-pr)
				
				# which are significant
				binOri <- which(sig <= alpha | sig2 <= alpha)	
				
			# output
				res<-list(ext=binExt, ori=binOri)
		
		}
	}
	if(method=="combine"){
		res$ext<-res$ext[res$ext%in%aicRes]
		res$ori<-res$ori[res$ori%in%aicRes]
	}
	
	
	
	return(res)
}

