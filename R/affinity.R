#' Calculate environmental affinities of taxa
#'
#' This function will return the environmental affinities of taxa, given the sampling conditions implied by the supplied dataset.
#'
#' Sampling patterns have an overprinting effect on the frequency of taxon occurrences in different environments. The environmental affinity (Foote, 2006; Kiessling and Aberhan, 2007; Kiessling and Kocsis, 2015) expresses whether the taxa are more likely to occur in an environment, given the sampling patterns of the dataset at hand. The function returns the likely preferred environment for each taxon as a vector. \code{NA} outputs indicate that the environmental affinity is equivocal based on the selected method.
#'
#' \strong{The following methods are implemented:}
#'
#' \code{'majority'}: Environmental affinity will be assigned based on the number of occurrences of the taxon in the different environments, without taking sampling of the entire dataset into account. If the taxon has more occurrences in \emph{environment 1}, the function will return \emph{environment 1} as the preferred habitat. 
#'
#' \code{'binom'}: The proportion of occurrences of a taxon in \emph{environment 1} and \emph{environment 2} will be compared to a null model, which is based on the distribution of all occurrences from the stratigraphic range of the taxon. The \code{alpha} value indicates the significance of the binomial tests, setting \code{alpha} to \code{1} will effectively switch the testing off: if the ratio of occurrences for the taxon is different from the ratio observed in the dataset, an affinity will be assigned. This is the default method.
#' 
#' \strong{References}
#'
#' Foote, M. (2006). Substrate affinity and diversity dynamics of Paleozoic marine animals. Paleobiology, 32(3), 345-366.
#'
#' Kiessling, W., & Aberhan, M. (2007). Environmental determinants of marine benthic biodiversity dynamics through Triassic–Jurassic time. Paleobiology, 33(3), 414-434.
#'
#' Kiessling, W., & Kocsis, Á. T. (2015). Biodiversity dynamics and environmental occupancy of fossil azooxanthellate and zooxanthellate scleractinian corals. Paleobiology, 41(3), 402-414.
#' 
#' @param dat \code{(data.frame)} The occurrence dataset containing the taxa with unknown environmental affinities.
#' @param env \code{(character)} The environmental variable of the occurrences.
#' @param method \code{(character)} The method used for affinity calculations. Can be either \code{"binom"} or \code{"majority"}.
#' @param tax \code{(character)} The column name of taxon names.
#' @param bin \code{(character)} The column name of bin names.
#' @param coll \code{(character)} The column name of collection numbers.
#' @param alpha \code{(numeric)} The alpha value of the binomial tests. By default binomial testing is off (\code{alpha=1}).
#' @param reldat \code{(data.frame)} Database with the same structure as \code{dat}. Typically \code{dat} is the subset of \code{reldat}. If given, the occurrence distribution of \code{reldat} is used 
#' as the null model of sampling.
#'
#' @examples
#'	data(corals)
#'	# omit values where no occurrence environment entry is present, or where unknown
#'	  fossils<-subset(corals, stg!=95)
#'	  fossilEnv<-subset(fossils, bath!="uk")
#'	# calculate affinities
#'	  aff<-affinity(fossilEnv, env="bath", tax="genus", bin="stg", alpha=1)
#'	
#' @export
affinity<-function(dat, env, tax="genus",  bin="stg", coll="collection_no", method="binom", alpha=1,reldat=NULL){

	if(is.null(alpha)& !is.null(reldat)) warning("Majority rule selected, reldat will be ignored.")
	
	match.arg(method, c("binom", "majority"))
	
	# omit everything from dat that is not necessary
		dat<-unique(dat[,c(coll, tax,bin, env)])
	
	# the affinity variable
		affDat<-dat[, env]
		
		if(0!=sum(is.na(affDat))) stop("The 'env' variable contains NAs, please omit these entries first.")
		affLevels<-unique(affDat)
		if(length(affLevels)!=2) stop("The 'env' variable contains more or less than 2 levels of entries. ")

	# create an FAD-LAD matrix first
		dFL<-fadlad(dat, tax, bin)
	
	# add the names to the matrix so that apply can process it
		dFL$taxon<-rownames(dFL)
	
	#the bins
		allBins<-sort(unique(dat[,bin]))
		allBins<-allBins[!is.na(allBins)]
	
	# tables of occurences
	# occs of the first env.
		firstDat<-dat[dat[,env]==affLevels[1],]
		fTab<-table(firstDat[,bin], firstDat[,tax])
		class(fTab)<-"matrix"
		
		# redo
		firstTab<-matrix(0, nrow=length(allBins), ncol=nrow(dFL))
		rownames(firstTab)<-allBins
		colnames(firstTab)<-rownames(dFL)
		firstTab[rownames(fTab), colnames(fTab)]<-fTab
		
	# occs of the second env.
		secondDat<-dat[dat[,env]==affLevels[2],]
		sTab<-table(secondDat[,bin], secondDat[,tax])
		class(sTab)<-"matrix"

		# redo
		secondTab<-matrix(0, nrow=length(allBins), ncol=nrow(dFL))
		rownames(secondTab)<-allBins
		colnames(secondTab)<-rownames(dFL)
		secondTab[rownames(sTab), colnames(sTab)]<-sTab
		
	# in case a relative dataset is added
		if(!is.null(reldat)){
			relLev<-reldat[,env]
			if(sum(!affLevels%in%relLev)==2) stop("The 'env' variable in 'reldat' does not contain the 'env' entries of 'dat'. ")
			
			# not all bins are represented!
			if(sum(allBins%in%reldat[,bin])!=length(allBins)) stop("The 'bin' variable of 'reldat' does cover the range of 'bin' in 'dat'.")
			dRel<-unique(reldat[,c(coll, tax,bin, env)])
			
		# occs of the first env.
			firstRelDat<-dRel[dRel[,env]==affLevels[1],]
			firstRel<-table(firstRelDat[,bin])
			class(firstRel)<-"numeric"
	
		# occs of the second env.
			secondRelDat<-dRel[dRel[,env]==affLevels[2],]
			secondRel<-table(secondRelDat[,bin])
			class(secondRel)<-"matrix"

		}else{
			firstRel<-apply(firstTab, 1, sum)
			secondRel<-apply(secondTab, 1, sum)
		}
			
	# prepare for bayesian
	if(method=="bayesian"){
		# prior probability of belonging to First
		pH1<- 0.5
		
		# prior probability of affinity to second
		pH2<- 0.5
	}
	
	# calculate the affinity of every taxon
	affVarTaxon<-apply(dFL, 1, FUN=function(x){

#	affVarTaxon<-character(length(dFL$taxon))
#	for (i in 1:length(affVarTaxon)){
	#	x<- dFL[i,,drop=F]
	
	
	#subset of the taxons ranges
		# bins where the taxon is present
		vectRange<-as.numeric(x[1]):as.numeric(x[2])
		
		# total number of occurrences from the first environment
		firstRelSum<-sum(firstRel[as.character(vectRange)], na.rm=T)
		
		# total number of occurrences from the second environment
		secondRelSum<-sum(secondRel[as.character(vectRange)], na.rm=T)
	
		# expected probabilities of occurrence based on the total dataset (reldat)
		pFirst<-firstRelSum/(firstRelSum+secondRelSum)
		pSecond<-secondRelSum/(firstRelSum+secondRelSum)

		# occurrences of taxon in the first environment
		first<-sum(firstTab[as.character(vectRange), as.character(x[length(x)])], na.rm=T)
		
		# occurrences of taxon in the second environment
		second<-sum(secondTab[as.character(vectRange), as.character(x[length(x)])], na.rm=T)
		
		# total
		all<-first+second
		
		#no occurrences from known habitats
		if (first==0 & second==0)
		{
		#	affVarTaxon[i]<-"non"
			return(NA)
		}
		
		# otherwise, continue
		
		# majority rule
		if(method=="majority"){
			if(first>second){
				return(affLevels[1])
			}else{
				if(second>first){
					return(affLevels[2])
				}else{
					return(NA)
				}
			}
		} # end majority
		
		# binomial basic
		if(method=="binom"){
		
			#first preference
			if (first/all>pFirst)
			{
				if(stats::binom.test(first, all, pFirst, "greater")$p.val<=alpha)
				{
				#	affVarTaxon[i]<-affLevels[1]
					return(affLevels[1])
				}else{
				#	affVarTaxon[i]<-"non"
					return(NA)
				}
			
			}
			
			#second preference
			if (second/all>pSecond)
			{
				if(stats::binom.test(second, all, pSecond, "greater")$p.val<=alpha)
				{
				#	affVarTaxon[i]<-affLevels[2]
					return(affLevels[2])
				}else{
				#	affVarTaxon[i]<-"non"
					return(NA)
				}
			
			}
			
			if((second/all>pSecond)==(first/all>pFirst)){
			#	affVarTaxon[i]<-"non"
				return(NA)
			}
		} # end binom
		
	#	if(method=="bayesian"){
	#		
	#		pEgivH1 <- dbinom(first, all, pFirst)
	#	#	pEgivH2 <- dbinom(second, all, pSecond)
	#		pEgivH2 <- dbinom(second, all, pFirst)
	#		
	#
	#		# Bayes' theorem
	#		pPost<- (pEgivH1*pH1)/(pEgivH1*pH1+pEgivH2*pH2)
	#		
	#		if(round(pPost,10)<0.5) return(affLevels[1])
	#		if(round(pPost,10)>0.5) return(affLevels[2])
	#		if(round(pPost,10)==0.5) return(NA)
	#	
	#	}		
	}
	)
	#table(affVarTaxon)
	names(affVarTaxon)<-rownames(dFL)
	
	return(affVarTaxon)
}


