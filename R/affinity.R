#' Calculate environmental affinities of taxa
#'
#' This function will return the environmental affinities of taxa, given the sampling conditions implied by the supplied dataset.
#'
#' Sampling patterns have an overprinting effect on the frequency of taxon occurrences in different environments. The environmental affinity (sensu Kiessling and Aberhan, 2007 - Paleobiology 33, 414-434 and Kiessling and Kocsis, 2015 - Paleobiology 41, 402-414) expresses whether the taxa are more likely to occur in an environment, given the sampling patterns of the dataset at hand. NA output indicates that the environmental variable 
#'
#' @param dat (data.frame): the occurrence dataset containing the taxa with unknown environmental affinities.
#' @param env (char): The name of the column with the occurrences' environmental values.
#' @param tax (char): the column name of taxon names.
#' @param bin (char): the column name of bin names.
#' @param coll (char): the column name of collection numbers.
#' @param alpha (num): the alpha value of the binomial tests. By default (alpha=1) binomial testing is off.
#' @param reldat (data.frame): database with the same structure as 'dat'. Typically 'dat' is the subset of 'reldat'. If given, the occurrence distribution of 'reldat' is used 
#' as the null model of sampling.
#'
#' @examples
#'	data(scleractinia)
#'	# omit values where no occurrence environment entry is present, or where unknown
#'	fossils<-subset(scleractinia, slc!=95)
#'	fossilEnv<-subset(fossils, bath!="uk")
#'	# calculate affinities
#'	  aff<-affinity(fossilEnv, env="bath", tax="genus", bin="slc", alpha=1)
#'	
#' @export
affinity<-function(dat, env, tax="occurrence.genus_name",  bin="SLC", coll="collection_no", alpha=1,reldat=NULL){
#	dat<-scleractinia[scleractinia$slc!=95,]
#	dat<-dat[dat$envnow!="uk",]
#	env<- "envnow"
#	tax<-"genus"
#	bin<- "slc"
#	coll<-"collection_no"
#	alpha<-"0.1"
#	output<-"occs"
	
	# omit everything from dat that is not necessary
		dat<-dat[,c(coll, tax,bin, env)]
	
	# the affinity variable
		affDat<-dat[, env]
		
		if(0!=sum(is.na(affDat))) stop("The 'env' variable contains NAs, please omit these entries first.")
		affLevels<-unique(affDat)
		if(length(affLevels)!=2) stop("The 'env' variable contains more or less than 2 levels of entries. ")

		# in case a relative dataset is added
		if(!is.null(reldat)){
			relLev<-reldat[,env]
			if(sum(!affLevels%in%relLev)==2) stop("The 'env' variable in 'reldat' does not contain the 'env' entries of 'dat'. ")
			dRel<-reldat[,c(coll, tax,bin, env)]
		}else{
			dRel<-dat
		}
		
	# create an FAD-LAD matrix first
		dFL<-fadLad(dat, tax, bin)
	
	# add the names to the matrix so that apply can process it
		dFL$taxon<-rownames(dFL)
	

	# calculate the affinity of every taxon
	affVarTaxon<-apply(dFL, 1, FUN=function(x){

#	affVarTaxon<-character(length(dFL$taxon))
#	for (i in 1:length(affVarTaxon)){

	#subset of the taxons ranges
	#	dSub<-subset(dat, dat[,bin]>=dFL$FAD[i] & dat[,bin]<=dFL$LAD[i])
		dSub<-subset(dRel, dRel[,bin]>=as.numeric(x[1]) & dRel[,bin]<=as.numeric(x[2]))
		dSub<-dSub[, c(tax, coll, env)]
		dSub<-unique(dSub)
		
		dSub2<-dSub[,c(coll, env)]
		
		#throw away unknown occurrences
		#dSub2<-unique(dSub2)
		
		#ratio - total ratio of sampling in the corresponding interval
		pFirst<-length(dSub2[dSub2[, env]==affLevels[1],1])/length(dSub2[, env])
		pSecond<-length(dSub2[dSub2[, env]==affLevels[2],1])/length(dSub2[, env])
		
		#dataset corresponding to the taxon at hand
	#	dTax<-subset(dat, dat[,tax]==as.character(dFL$taxon[i]))
		dTax<-subset(dat, dat[,tax]==as.character(x[length(x)]))
		dTax<-unique(dTax[, c(tax, coll, env)])
		
		dTax2<-dTax[, c(coll, env)]

		#ratio of reefal to non-reefal collections
		#dTax2<-unique(dTax2)
		
		#reefal counts
		first<-length(dTax2[dTax2[, env]==affLevels[1],env])
		#nonreefal counts
		second<-length(dTax2[dTax2[, env]==affLevels[2],env])
		all<-length(dTax2[, env])
		
		#no occurrences from known habitats
		if (first==0 & second==0)
		{
		#	affVarTaxon[i]<-"non"
			return(NA)
		}else{
				
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
			
		}
	}
	)
	#table(affVarTaxon)
	names(affVarTaxon)<-rownames(dFL)
	
	return(affVarTaxon)
}


