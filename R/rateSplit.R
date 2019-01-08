#' Test of rate split (selectivity)
#'	
#' This function will determine whether there are meaningful differences between the taxonomic rates in the individual time bins of two subsets of an occurrence database.
#'	
#' Splitting an occurrence database to its subsets secreases the amount of information passed to the rate calculations and therefore the precision of the individual estimates. Therefore, our ability to tell apart two similar values decreases with the number of sampled taxa. In order to assess the subsets individually and compare them, it is advised to test whether the split into two subsets is meaningful, given the total data. Examples of this use can be found in Kiessling and Simpson (2011) and Kiessling and Kocsis (2015).
#' The meaningfulness of the split is dependent on the estimate accurracy and the magnitude of the difference. Two different methods are implemented: \code{binom} and \code{combine}.
#'
#' \strong{References}
#'
#' Foote, M. (1999) Morphological Diversity In The Evolutionary Radiation Of Paleozoic and Post-Paleozoic Crinoids. Paleobiology 25, 1–115. doi:10.1017/S0094837300020236.
#'
#' Kiessling, W., & Simpson, C. (2011). On the potential for ocean acidification to be a general cause of ancient reef crises. Global Change Biology, 17(1), 56-67.
#'
#' Kiessling, W., & Kocsis, A. T. (2015). Biodiversity dynamics and environmental occupancy of fossil azooxanthellate and zooxanthellate scleractinian corals. Paleobiology, 41(3), 402-414.
#'
#' @param sel \code{(character)}: Variable name to do the splitting of the dataset. Can have only two levels.
#'	
#' @param output \code{(character)}: Either \code{"simple"} or \code{"full"}. \code{"simple"} returns the indices of the series where selectivity can be suggested. \code{"full"} returns a \code{matrix} of Akaike weights, or binomial probabilities.
#'	
#' @param method \code{(character)}: Either \code{"AIC"}, \code{"binom"} or \code{"combine"}. The \code{"AIC"} method calculates the Akaike weights of the single and dual rate models. The \code{"binom"} method assumes a binomial error distribution of the counts that are necessary for the rate calculations. The \code{"combine"} method shows slices that pass both tests, the \code{"AIC"} being usually the stronger.
#'	
#' @param alpha \code{(numeric)}: Threshold to discriminate between meaningful and meaningless split. If \code{method="AIC"}, the value corresponds to the minimum weight value the dual model should have. By default it is \code{0.89}, which corresponds to the likelihood ratio of 8. If \code{method="binom"}, the value corresponds to the alpha value of the binomial test (default: 0.05). If \code{method="combine"} than two alpha values are required (1st for the AIC, 2nd for the binomial test). If alpha is \code{NULL}, than the default values will be used.
#'	
#' @param AICc \code{(logical)}: Only applicable for the \code{"AIC"} method. Toggles whether the small sample corrected AIC (AICc) should be used instead of the regular one.
#'	
#' @param rate \code{(character)}: The rate metric. Currently only the per capita rates of Foote (1999) are available (\code{rate="pc"}).
#' @param dat \code{(data.frame)}: The fossil occurrence data.
#'	 
#' @param na.rm \code{(logical)}: Argument indicating whether the function should proceede when \code{NA}s are found in the \code{sel} column. Setting this argument to \code{TRUE} will proceede with the omission of these entries, while \code{FALSE} will coerce the function to output a single \code{NA} value.
#' @param tax \code{(character)}: Variable name of the occurring taxa (variable type: \code{factor} or \code{character}).
#'	 
#' @param bin \code{(character)}: Variable name of the bin numbers of the particular occurrences (\code{numeric}). Bin numbers should be in ascending order,can contain \code{NA}s, it can start from a number other than 1 and must not start with 0.
#'
#' @examples
#'	# example with the coral dataset of Kiessling and Kocsis (2015)
#'	data(corals)
#'	data(stages)
#'	
#'	# split by ecology
#'	  z<-corals[corals$ecology=="z",]
#'	  az<-corals[corals$ecology=="az",]
#'	
#'	# calculate diversity dynamics
#'	ddZ<-divDyn(z, tax="genus", bin="stg")
#'	ddAZ<-divDyn(az, tax="genus", bin="stg")
#'	
#'	# origination rate plot
#'	tsplot(stages, boxes="sys", shading="series", xlim=54:95, 
#'	  ylab="raw per capita originations")
#'	lines(stages$mid, ddZ$oriPC, lwd=2, lty=1, col="blue")
#'	lines(stages$mid, ddAZ$oriPC, lwd=2, lty=2, col="red")
#'	legend("topright", inset=c(0.1,0.1), legend=c("z", "az"), 
#'	  lwd=2, lty=c(1,2), col=c("blue", "red"), bg="white")
#'	
#'	# The ratesplit function
#'	rs<-ratesplit(rbind(z, az), sel="ecology", tax="genus", bin="stg")
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
ratesplit<-function(dat,  sel, tax="genus", bin="stg", rate="pc", method="AIC",AICc=TRUE, na.rm=TRUE, alpha=NULL, output="simple"){
#	 #dummy data:
# 	
#	 dat<-corals
#	 tax<-"genus"
#	 bin<-"stg"
#	 cRate<-"nFooteOri"
#	 
#	
#	na.rm=T
#	
	# defense
	if(length(output)!=1 | length(method)!=1) stop("Enter only one value for the output and method arguments.")
	if(!method%in%c("AIC", "binom", "combine")) stop("Invalid method.")
	if(!output%in%c("simple", "full")) stop("Invalid output.")


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
		if(method=="AIC") alphaAIC <- 0.89
		if(method=="binom") alphaBinom <- 0.05
		if(method=="combine") {
			alphaAIC <- 0.89
			alphaBinom <- 0.05
		}
		
	}else{
		if(method=="AIC") alphaAIC<-alpha[1]
		if(method=="binom") alphaBinom<-alpha[1]
		if(method=="combine"){
			if(length(alpha)!=2) stop("With the combine method you need at least two alpha values.")
			alphaAIC<-alpha[1]
			alphaBinom<-alpha[2]
		} 
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
		
		if(method=="AIC" | method=="combine"){
		
		
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
				AkaikeSelExt <- which(wExt2 > alphaAIC)
				AkaikeSelOri <- which(wOri2 > alphaAIC)
				
			if(output=="simple"){
				res<-list(ext=AkaikeSelExt, ori=AkaikeSelOri)
				
				if(method=="combine"){
					aicRes<-res
				}
				
			}
			if(output=="full"){
				
				res<-data.frame(singleOri=wOri1, dualOri=wOri2,singleExt=wExt1, dualExt=wExt2)
				if(method=="combine"){
					aicRes<-res
				}else{
					return(res)
				}
			}
		}
		
		if(method=="binom" | method=="combine"){
			
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
				sigExt <- stats::pbinom(kExt,nExt,pr,lower.tail=F)+stats::dbinom(kExt,nExt,pr)
				sigExt2 <- stats::pbinom(nExt-kExt,nExt,1-pr,lower.tail=F)+stats::dbinom(nExt-kExt,nExt,1-pr)
				
				# which are significant
				binExt <- which(sigExt <= alphaBinom | sigExt2 <= alphaBinom)
			
			# Binomial test (specific for three columns) - originations
				sigOri <- stats::pbinom(kOri,nOri,pr,lower.tail=F)+stats::dbinom(kOri,nOri,pr)
				sigOri2 <- stats::pbinom(nOri-kOri,nOri,1-pr,lower.tail=F)+stats::dbinom(nOri-kOri,nOri,1-pr)
				
				# which are significant
				binOri <- which(sigOri <= alphaBinom | sigOri2 <= alphaBinom)	
				
			# output
			if(output=="simple"){
				res<-list(ext=binExt, ori=binOri)
			}
			if(output=="full"){
				res<-cbind(sigExt, sigExt2, sigOri, sigOri2)
			
			}
		}
	}
	if(method=="combine"){
		if(output=="simple"){
			res$ext<-binExt[binExt%in%AkaikeSelExt]
			res$ori<-binOri[binOri%in%AkaikeSelOri]
		}
		if(output=="full"){
			res<-cbind(res, aicRes)
			names(res)<- c("binSignExt1", "binSignExt2", "binSignOri1", "binSignOri2", 
			"aicSingleOri", "aicDualOri", "aicSingleExt", "aicDualExt")
		}
	}
	
	
	
	return(res)
}

