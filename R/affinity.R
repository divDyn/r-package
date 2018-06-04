#' Calculate environmental affinities of taxa
#'
#' This function will return the environmental affinities of taxa, given the sampling conditions implied by the supplied dataset.
#'
#' Sampling patterns have an overprinting effect on the frequency of taxon occurrences in different environments. The environmental affinity (sensu Kiessling and Aberhan, 2007 - Paleobiology 33, 414-434 and Kiessling and Kocsis, 2015 - Paleobiology 41, 402-414) expresses whether the taxa are more likely to occur in an environment, given the sampling patterns of the dataset at hand. NA output indicates that the environmental affinity is equivocal based on the selected method.
#'
#' The following methods are implemented: 
#'
#' 'majority': environmental affinity will be assigned based on the sampled proportions, without assuming a equal sampling probabilities in the two environments (more occurrence in env. 1 means that taxon prefers env. 1)
#'
#' 'binom': The proportion of occurrences of a taxon in environment 1 and environment 2 will be compared to a null model which is based on the distribution of all occurrences from the stratigraphic range of the taxon. The alpha value indicates the significance of the binomial tests. 
#' 
#' 'bayesian': The same as 'binom', except that instead of a binomial test, Bayesian inference is used (Simpson and Harnik, 2009). Note: using choose(n, k)*p^k*(1-p)^(n-k) to calculate the conditional probability P(E|H1) (as it was published in Simpson and Harnik, 2009) did not produce a meaningful output. This is likely the case as this formula calculates the probability of the exact event of sampling k success out of n trials, which decreases systematically as n increases. However, using the CDF (pbinom(n, k, p)) to estimate P(E|H1) produces plausible final output. This corresponds to the event that given H1, out of n trials, k is the maximum number of env1 occurrences given its sampling conditions characterized by p (total dataset). The function is implemented thusly.
#'
#' @param dat (data.frame): the occurrence dataset containing the taxa with unknown environmental affinities.
#' @param env (char): The name of the column with the occurrences' environmental values.
#' @param method (character): The method used for affinity calculation. Can be either "binom", "bayesian" or "majority".
#' @param tax (char): the column name of taxon names.
#' @param bin (char): the column name of bin names.
#' @param coll (char): the column name of collection numbers.
#' @param alpha (num): the alpha value of the binomial tests. By default (alpha=1) binomial testing is off.
#' @param reldat (data.frame): database with the same structure as 'dat'. Typically 'dat' is the subset of 'reldat'. If given, the occurrence distribution of 'reldat' is used 
#' as the null model of sampling.
#'
#' @examples
#'	data(corals)
#'	# omit values where no occurrence environment entry is present, or where unknown
#'	fossils<-subset(corals, slc!=95)
#'	fossilEnv<-subset(fossils, bath!="uk")
#'	# calculate affinities
#'	  aff1<-affinity(fossilEnv, env="bath", tax="genus", bin="slc", alpha=1)
#'    aff2<-affinity(fossilEnv, env="bath", tax="genus", bin="slc", method="bayesian")
#'	
#' @export
affinity<-function(dat, env, tax="genus",  bin="slc", coll="collection_no", method="binom", alpha=1,reldat=NULL){

	if(is.null(alpha)& !is.null(reldat)) warning("Majority rule selected, reldat will be ignored.")
	
	# omit everything from dat that is not necessary
		dat<-unique(dat[,c(coll, tax,bin, env)])
	
	# the affinity variable
		affDat<-dat[, env]
		
		if(0!=sum(is.na(affDat))) stop("The 'env' variable contains NAs, please omit these entries first.")
		affLevels<-unique(affDat)
		if(length(affLevels)!=2) stop("The 'env' variable contains more or less than 2 levels of entries. ")

	# create an FAD-LAD matrix first
		dFL<-fadLad(dat, tax, bin)
	
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
				if(binom.test(first, all, pFirst, "greater")$p.val<=alpha)
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
				if(binom.test(second, all, pSecond, "greater")$p.val<=alpha)
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
		
		if(method=="bayesian"){
	#		pEgivH1<- choose(all, first) * (pFirst^first)*(1-pFirst)^(all-first)
			# le merném fogadni, bassza meg, hogy ezt használták...
			pEgivH1 <- pbinom(first, all, pFirst)
			pEgivH2 <- pbinom(second, all, pSecond)
			
		#	pEgivH2<- choose(all, second) * (pSecond^second)*(1-pSecond)^(all-second)
			# Bayes' theorem
			pPost<- (pEgivH1*pH1)/(pEgivH1*pH1+pEgivH2*pH2)
			
			if(round(pPost,10)>0.5) return(affLevels[1])
			if(round(pPost,10)<0.5) return(affLevels[2])
			if(round(pPost,10)==0.5) return(NA)
		
		}		
	}
	)
	#table(affVarTaxon)
	names(affVarTaxon)<-rownames(dFL)
	
	return(affVarTaxon)
}


