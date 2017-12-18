# Subsampling scripts

# Wrapper function to perform a specified task with Classical Rarefaction of bins
subseries <- function(dat, FUN=NULL, quota,iter=10,  bin="SLC", tax="genus", intact=NULL, output="arithmetic", implementation="for", method="cr"){

#	bin <- "slc"
#	tax<- "genus"
#	dat <- scleractinia
#	quota<-40
#	output<- "arithmetic"
#	intact<-c(94,95)

	# check the presences of vectors

	#a. given that the final output is not a list
	if (output!="list"){
		#a. perform the function on the total dataset to assess final structure
		if(is.null(FUN)){
			wholeRes<-divDyn(dat, tax=tax, bin=bin)
		}else{
			if(is.function(FUN)){
				wholeRes<-FUN(dat, ...)
			}else{
				stop("Invalid FUN argument.")
			}
		}
		
		# default
		out<-"list"
		
		# if the output is a vector
		if(is.vector(wholeRes)){
			holder <- matrix(NA, ncol=iter, nrow=length(wholeRes))
			rownames(holder)<-names(wholeRes)
			out<-"vector"
		}
		# if the output is matrix
		if(is.matrix(wholeRes)){
			holder<-array(NA, dim=c(dim(wholeRes),iter))
			colnames(holder) <- colnames(wholeRes)
			rownames(holder) <- rownames(wholeRes)
			out<-"matrix"
		
		}
		#if the output is a data.frame
		if(is.data.frame(wholeRes)){
			nVar<-ncol(wholeRes)
			# repeat original structure
			holder <- matrix(NA, ncol=nVar*iter, nrow=nrow(wholeRes))
			holder<-as.data.frame(holder)
			colnames(holder)<-rep(1:iter, each=nVar)
			out<-"data.frame"
		}

	}
	
	# if it is a list, the averaging has to be done by the urser
	if(output=="list" | out=="list"){
		message("U have to do the averaging yourself")
		# force list output
		out <-"list"
		
		# create a final container
		containList<-list()
	}
	
	#b. iteration
	if(implementation=="for"){
	
		for(k in 1:iter){
			#1. produce one subsample
			if(method=="cr"){
				
				trialDat<-dat[subsampleCR(dat[,bin], quota, intact),]
			}
			
			#2. run the function on subsample
				if(is.null(FUN)){
					trialRes<- divDyn(trialDat, tax=tax, bin=bin)
				}else{
					trialRes <- FUN(trialDat, ...)
				}
			
			#3. save the results depending on the output type
			if(out=="vector"){
				holder[,k]<-trialRes
			}
			if(out=="matrix"){
				holder[,,k] <- trialRes
			}
			if(out=="data.frame"){
				holder[,(1:nVar)+((k-1)*nVar)] <- trialRes
			}
			if(out=="list"){
				containList<-c(containList, list(trialRes))
			}
				
		}
		
	}
	
	
	#c. averaging and return
	if(out=="vector"){
		if(output=="arit"){
			totalResult<-apply(holder,2,mean, na.rm=T)
		}
		if(output=="geom"){
			loggedVars<-log(holder)
			loggedVars[is.infinite(loggedVars)]<-NA
			totalResult[i]<- apply(loggedVars, 1, mean, na.rm=T) 
			totalResult[i]<- exp(totalResult[i])
		}
	}
	
	if(out=="matrix"){
	
		if(output=="arit"){
			totalResult<-apply(holder,c(1,2),mean, na.rm=T)
		}
		
		if(output=="geom"){
			loggedVars<-log(holder)
			loggedVars[is.infinite(loggedVars)]<-NA
			totalResult[i]<- apply(loggedVars, c(1,2), mean, na.rm=T) 
			totalResult[i]<- exp(totalResult[i])
		}
	}
	
	if(out=="data.frame"){
		totalResult <- matrix(NA, ncol=ncol(wholeRes), nrow=nrow(wholeRes))
		totalResult<-as.data.frame(totalResult)
		colnames(totalResult) <- colnames(wholeRes)
		rownames(totalResult) <- rownames(wholeRes)
		for(i in 1:nVar){
			# if it is numeric at all!!
			if(!is.character(wholeRes[,i])){
				# variable specific columns
				varCol<-seq(i,ncol(holder), by=nVar)
				varDat<-holder[,varCol]
				
				# which meaning should be used
				if(output=="arit"){
					totalResult[,i]<-apply(varDat, 1, mean, na.rm=T)
				}
				
				if(output=="geom"){
					loggedVars<-log(varDat)
					loggedVars[is.infinite(loggedVars)]<-NA
					totalResult[,i]<- apply(loggedVars, 1, mean, na.rm=T) 
					totalResult[,i]<- exp(totalResult[,i])
				}
				
				# some hacking required for the divDyn for the 'outer rates'
				
			}else{
				# some functionality can be added
				totalResult[,i] <- wholeRes[,i]
			}
			

		}
	
	}
	if(out=="list"){
		totalResult<- containList
	}
	
	return(totalResult)
}

# generate one subsample with classical rarefaction
subsampleCR <- function(binVar,quota,intact){
	binVar <- dat[,bin]
	rows<- 1:length(binVar)
	
	tap<-tapply(INDEX=binVar, X=rows, FUN=function(x, quota){
		if(length(x)>=quota){
			return(sample(x, quota, replace=F))
		}else{
			return(NULL)
		}
	
	}, quota=quota)
	
	# the rows that should be passed
	trialRows <- rows%in%unlist(tap)
	
	# these rows should be present regardless of the subsampling
	intactRows<-binVar%in%intact
	
	# combine the two
	subsampleRows<-trialRows | intactRows

	# return
	return(subsampleRows)
}