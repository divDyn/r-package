#' Cleanse Species Name Vector
#' 
#' This function will take a vector of binomial names with various qualifiers of open nomenclatures, and removes them form the vector entries. Only the the genus and species names will remain.
#'
#' This version will not keep subgenera, and will assign species to the base genus. The following qualifiers will be omitted:
#' \emph{"n."}, \emph{"sp."}, \emph{"?"}, \emph{"gen."}, \emph{"aff."}, \emph{"cf."}, \emph{"ex gr."}, \emph{"subgen."}, \emph{"spp"} and informal species designated with letters. Entries with \emph{"informal"} and \emph{"indet."} in them will also be invalidated. 
#' 
#' Functions called by the \code{misspells} and \code{stems} arguments were written by Gwen Antell. 
#' @param vec \code{(character)}: the vector containing species names with qualifiers of open taxonomy.
#'
#' @param debug \code{(logical)}: \code{FALSE} will return the cleaned species name vector, \code{TRUE} returns a data table that allows one by one checking.
#' @param collapse \code{(character)}: This argument will be passed to the paste function's argument of the same name. The character value to be inserted between the genus and species names.
#' @param subgenera \code{(logical)}: \code{FALSE} omits subgenus information (in parentheses) and will construct a unique binomen based on the genus and species names alone. \code{TRUE} will promote the subgenus names and it will create a new binomen based on the subgenus rather than the genus name.
#' @param stems \code{(logical)}: Setting this to \code{TRUE} will omit the adjective declination suffices from the species names. 
#' @param misspells \code{logical}: Resolution of common spelling mistakes, such as diphtongs and alternate spellings: 'ue' is replaced with 'u', 'ae' is replaced with 'e', 'll' with 'l', 'ss' with 's'and 'y' with 'i'. 
#' @examples
#' examp <- c("Genus cf. species", "Genus spp.", "Family indet.", 
#'   "Mygenus yourspecies", "Okgenus ? questionsp", 
#'   "Genus (cf. Subgenus) aff. species")
#' cleansp(examp) 
#' @export
# function to cleanse a noisy species name vector
cleansp <- function(vec, debug=FALSE, collapse="_", subgenera=FALSE, misspells=TRUE, stems=TRUE){
	
	# keep the original
	vecOrig<-vec
	# presplit for parenthesis error
	vec<-gsub("\\( ", "\\(", vec)
	
	# presplit the the double space error
	vec<-gsub("  ", " ", vec)
	
	# split the string
	split<-strsplit(vec, " ")
	
	# excluded parts
		# these entries will be omitted without further issues	
		exclude <- c("n.", "sp.", "?", "gen.", "aff.", "cf.", "ex", "gr.", "subgen.", paste(LETTERS, ".", sep=""), "spp.", "(n.", "subgen.", "?)", "(?", "(cf.")
		
		# if elements of this list are found in pairs, those will be excluded (but only if both are found, so species name "lato" can be left)
		jointExclude<-list(c("sensu", "lato"), c("sensu", "stricto"))
	
		# if these entries are around, the species name will be invalid
		special <- c("sp.1","sp.2", "informal", "indet.", letters)
	
	
	dual<-lapply(split, function(x){
	# missing entries
		if(sum(is.na(x))==length(x)){
			return(NA, NA)
		}
	
	#is a name starting with quotes - remove quotes
		
		quotes<-sapply(x, function(y){
			substr(y, 1,1)%in%c("\"")
		})
		if(sum(quotes)>0){
			tem<-unlist(strsplit(x[quotes], "\""))[2]
			x[quotes]<-tem
		}
		
	# omit the prefixes and suffixes
		x<-x[!x%in%exclude]

	# omit the jointly occurring notes (e.g. 'sensu lato')
		jointOcc<-unlist(lapply(jointExclude, function(y){
			sum(y%in%x)==length(y)
		
		}))
		if(sum(jointOcc)>0){
			je<-jointExclude[jointOcc]
			for(i in 1:length(je)){
				x<-x[!x%in%je[[i]]]
			}
			
		}
		
	# if there is a non-valid species name indicator - remove the entry
		if(sum(x%in%special)>0){
			return(c(NA, NA))
		}
		numConvert<-suppressWarnings(as.numeric(x))
		if(sum(!is.na(numConvert))>0){
			return(c(NA, NA))
		}
		
		if(length(x)==1){
			return(c(NA, NA))
		}
	
	# is there a subgenus name - omit it
		# first character is parenthesis- omit
		parenth1<-sapply(x, function(y){
			substr(y, 1,1)=="("
		})
		
		if(sum(parenth1)>0){
			# exclude
			if(!subgenera){
				x<-x[!parenth1]
			}else{
				subgen <- x[parenth1]
				subgen <-substr(subgen, 2,nchar(subgen))
				x<-x[!parenth1]
				x[1] <- subgen

				# in case it is a single entry, omit the unneded 
				if(substr(subgen, nchar(subgen),nchar(subgen))==")"){
					x[1]<-substr(subgen, 1,nchar(subgen)-1)
				}
			}
		}
		
		#last character is parenthesis- because of other words/qualifiers before
		parenth2<-sapply(x, function(y){
			substr(y, nchar(y),nchar(y))==")"
		})
		
		if(sum(parenth2)>0){
			if(!subgenera){
				x<-x[!parenth2]
			}else{
				subgen <- x[parenth2]
				subgen <-substr(subgen, 1,nchar(subgen)-1)
				x<-x[!parenth2]
				x[1] <- subgen
			}
		}
		

	# roots
		
		
	# return genus and species name - potentially subspecies and crap at the end
		return(x)
	
	})
	
#	# not two
#	len<-unlist(lapply(dual, length))!=2
#	
#	
#	View(cbind(gen,sp)[len,])
#	
#	prob<-vec[len]
#	View(prob)
	
	# merge the genus and species in one column
	singular<-unlist(lapply(dual, function(x){
		if(is.na(x[1])){
			return(NA)
		}else{
			paste(x[1:2], collapse=collapse)
		}
	}))
	if(stems){
		singular <- sapply(singular, get_root)
	}

	if(misspells){
		singular <- sapply(singular, misspell)
	}

	# if names start with " omit those as well
	# omit parentheses
	if(debug){
		
		gen<-unlist(lapply(dual, function(x){
			x[1]
		}))
		
		sp<-unlist(lapply(dual, function(x){
			x[2]
		
		}))
	
		dat<-data.frame(original=vecOrig,genus=gen,species=sp, omitted=rep(FALSE, length(gen)))
		dat$omitted[is.na(singular)] <- TRUE
		dat$binomen <- singular
		return(dat)
	}else{
		return(singular)
	}
}


# good, but this has to be rewritten in Rcpp
get_root <- function(s){
	last2 <- substr(s, nchar(s)-1, nchar(s))
	last1 <- substr(s, nchar(s), nchar(s))
	root <- s
	if (last2 %in% c('us', 'um', 'ae', 'is', 'es', 'ii')){ # check for 2-letter suffix
		root <- substr(s, 1, nchar(s)-2)
	} else { # check for 1-letter suffix
		if (last1 %in% c('a', 'i', 'e')){ root <- substr(s, 1, nchar(s)-1) }
	}
	return(root)
}


misspell <- function(term){
  term <- gsub('ue','u', term, fixed=TRUE)
  term <- gsub('ae','e', term, fixed=TRUE)
  term <- gsub('ll','l', term, fixed=TRUE)
  term <- gsub('ss','s', term, fixed=TRUE)
  term <- gsub('y', 'i', term, fixed=TRUE)
  term
}
