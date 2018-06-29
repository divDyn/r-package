#' Time series plotting using a custom time scale 
#'
#' This function allows the user to quickly plot a time scale data table
#'
#' As most analysis use an individually compiled time scale object, in order to ensure compatibility between the analyzed and plotted values, the time scale table used for the analysis could be plotted rather than a standardized table. Two example tables have been included in the package (\code{\link{stages}} and \code{\link{bins}}) that can serve as templates.
#' @param tsdat \code{(data frame)}: The time scale data frame.
#' @param boxes \code{(character)}: Column name indicating the names that should be plotted as boxes of the timescale.
#' @param ylim \code{(numeric)}: The vertical extent of the plot, analogous to the same argument of \code{\link[graphics]{plot}}. By default it is set to the \code{[0,1]} interval.
#' @param xlim \code{(numeric)}: The horizontal extent of the plot, analogous to the same argument of \code{\link[graphics]{plot}}. By default it is set to plot the entire table. If a numeric vector with two values is supplied, it will be interpreted as the standard \code{xlim} argument and the plot will be displayed based on numerically constrained ages. If it is an integer vector with more than two values are plotted, the interval corresponding to the row indices of the table will be plotted.
#' @param prop \code{(numeric)}: Proportion of the vertical extent of the plot to display the the time scale at the bottom.
#' @param gap \code{(numeric)}: Proportion of the vertical extent of the plot that should be a gap betwen the time scale and the plot.
#' @param bottom \code{(character)}: Column name of the table for the variable that contains the older ages of intervals.
#' @param top \code{(character)}: Column name of the table for the variable that contains the earliest ages of intervals.
#' @param shading \code{(character)}: Column name used for the shading. By default, no shading will be drawn (\code{shading = NULL}).
#' @param shading.col \code{(character)}: Name of colors that will be used for the shading, if shading is set. The provided colors will be repeated as many times as necessary.
#' @param plot.args \code{(list)}: Arguments that will be passed to the main \code{\link[graphics]{plot}} function.  Can be useful for the suppression of axes, font change etc.
#' @param labels.args \code{(list)}: Arguments that will be passed to the \code{\link[graphics]{text}} function that draws the labels. 
#' @param boxes.args \code{(list)}: Arguments that will be passed to the \code{\link[graphics]{rect}} function that draws the rectangles of time intervals.
#' @examples
#'	data(stages)
#'	  plotTS(stages, boxes="per", shading="series")
#'
#'	# only the Mesozoic, custom axes
#'	  plotTS(stages, boxes="period", shading="stage", xlim=52:81, 
#'	    plot.args=list(axes=F, main="Mesozoic"))
#'	  axis(1, at=seq(250, 75, -25), labels=seq(250, 75, -25))
#'	  axis(2)
#'	
#'	# only the Triassic, use the supplied abbreviations
#'	  plotTS(stages, boxes="short", shading="stage", xlim=c(250,199), 
#'	    ylab="variable", labels.args=list(cex=1.5, col="blue"), 
#'	    boxes.args=list(col="gray95"))
#' @export
plotTS<-function(tsdat,  boxes, ylim=c(0,1), xlim=NULL, prop=0.05, gap=0,
	bottom="bottom", top="top",
	xlab="age (Ma)", ylab="",
	shading=NULL,shading.col=c("white", "gray80"),
	plot.args=NULL,
	boxes.args=NULL,
	labels.args=NULL){
	
#	tsdat<-stages
#	boxes<-"per"
#	bottom="bottom"
#	top="top"
#	ylim=c(0,1)
#	xlab="age (Ma)"
#	ylab=""
#	shading="series"
#	shading.col=c("white", "gray80")
#	
#	boxes.args<-NULL
#	plot.args<-NULL
#	labels.args<-NULL
	tsCols<-colnames(tsdat)
	if(!bottom%in%tsCols) stop("The 'bottom' column is not found.")
	if(!top%in%tsCols) stop("The 'top' column is not found.")
	
	
	# length(ylim) should be 2, and numeric, default
	if(!is.numeric(ylim)){
		stop("ylim should be a numeric vector of length 2")
	}
	
	if(length(boxes)>1){
		stop("Only one column can be used to plot the boxes.")
		
	}else{
		if(!boxes%in%colnames(tsdat)){
			stop("The referenced 'boxes' column does not exist.")
		}
	}
	if(!is.null(shading)){
		if(length(shading)>1){
			stop("Only one column can be used to plot the shades.")
			
		}else{
			if(!shading%in%colnames(tsdat)){
				stop("The referenced 'shading' column does not exist.")
			}
		}
	}
	
	# if the log axis is planned
	yLog<-FALSE
	if(!is.null(plot.args)){
		if(!is.null(plot.args$log)){
			if(is.na(plot.args$log)) stop("Invalid log argument in plot.args.")
			if(plot.args$log=="y"){
				yLog<-TRUE
			
			}
		}
	
	}
	
	# xlim defense
	if(is.null(xlim)){
		xlim <- c(max(tsdat[,bottom],na.rm=T),min(tsdat[,top],na.rm=T))
	} else{
		if(is.numeric(xlim)){
			if(length(xlim)>2){
				if(sum(xlim%%1)==0 & sum(xlim<=0)==0 & sum(xlim>nrow(tsdat))==0){
					xlim <- c(max(tsdat[xlim,bottom],na.rm=T),min(tsdat[xlim,top],na.rm=T))
				}
			}
			if(length(xlim)==1){
				stop("Invalid xlim values.")
			}
		}
	}
	
	# the ylim adjustment
	if(!yLog){
		boxesTop<-ylim[1]-diff(range(ylim))*gap
		plotBottom<-ylim[1]
	
		ylim[1]<-ylim[1]-diff(range(ylim))*prop-diff(range(ylim))*gap
	}else{
		boxesTop<-exp(log(ylim[1])-diff(range(log(ylim)))*gap)
		plotBottom<-ylim[1]
	
		ylim[1]<-exp(log(ylim[1])-diff(log(range(ylim)))*prop-diff(log(range(ylim)))*gap)

	
	}
	
	# check whether the additional arguments work ok or not
	if(!is.null(boxes.args)){
		if(!is.list(boxes.args)){
			stop("Invalid additional arguments supplied for the boxes.args/rect().")
		}
	}
	
	if(!is.null(plot.args)){
		if(!is.list(plot.args)){
			stop("Invalid additional arguments supplied for plot.args/plot().")
		}
	}
	
	if(!is.null(labels.args)){
		if(!is.list(labels.args)){
			stop("Invalid additional arguments supplied for labels.args/text() .")
		}
	}
	
		

	# the empty plot
		# invoke plot
		plotArgs<-list(
			x=1,
			y=1,
			xlim=xlim, 
			ylim=ylim,
			xlab=xlab, 
			ylab=ylab,
			xaxs="i",
			yaxs="i",
			type="n"
		)
		plotArgs<-c(plotArgs, plot.args)
		do.call(plot, plotArgs)
		
	# the box drawing
	boxLev<-unique(tsdat[,boxes])
	xLeft<-rep(NA, length(boxLev))
	xRight<-rep(NA, length(boxLev))
	yTop<-rep(NA, length(boxLev))
	yBottom<-rep(NA, length(boxLev))
	labMid<-rep(NA, length(boxLev))
	for(i in 1:length(boxLev)){
		# the rectangles
			boolSel<-boxLev[i]==tsdat[,boxes]
			curBottom<-max(tsdat[boolSel,bottom], na.rm=T)
			curTop<-min(tsdat[boolSel,top], na.rm=T)
			
			# export the positions
			xLeft[i]<-curBottom
			xRight[i]<-curTop
			yTop[i]<-boxesTop
			yBottom[i]<-ylim[1]
			labMid[i]<-mean(c(curBottom, curTop))
		
	}
	
	#boxes
	# inner arguments
	boxArgs<-list(
		xleft=xLeft, 
		xright=xRight, 
		ytop= yTop,
		ybottom=yBottom
	)
	
	# process the supplied arguments, in case they are in the table
	boolUnique<-boxes.args%in%colnames(tsdat)
	if(sum(boolUnique)>0){
		ind<-which(boolUnique)
		for(i in ind){
			boxes.args[[i]]<-unique(tsdat[,boxes.args[[i]]])
		}
	}
	
	# combine with outer arguments
	boxArgs<-c(boxArgs, boxes.args)
	do.call(rect, boxArgs)
				
		
		
	# the labels
	# inner arguments
	labArgs<-list(
		label=boxLev, 
		y=(boxArgs$ytop+boxArgs$ybottom)/2, 
		x= labMid
	)
		
	# combine inner and outer arguments
	# process the supplied arguments, in case they are in the table
	boolUnique<-labels.args%in%colnames(tsdat)
	if(sum(boolUnique)>0){
		ind<-which(boolUnique)
		for(i in ind){
			labels.args[[i]]<-unique(tsdat[,labels.args[[i]]])
		}
	}
	
	labArgs<-c(labArgs, labels.args)
	do.call(text,labArgs)
	
	# shading
	if(!is.null(shading)){
		shadeLev<-unique(tsdat[,shading])
		for(i in 1:length(shadeLev)){
			
			# the rectangles
			boolSel<-shadeLev[i]==tsdat[,shading]
			curBottom<-max(tsdat[boolSel,bottom], na.rm=T)
			curTop<-min(tsdat[boolSel,top], na.rm=T)
			if(i%%length(shading.col)){
				rect(xleft=curBottom, xright=curTop, ybottom=plotBottom, ytop=ylim[2], border=NA, col=shading.col[i%%length(shading.col)+1])
			}else{
				rect(xleft=curBottom, xright=curTop, ybottom=plotBottom, ytop=ylim[2], border=NA, col=shading.col[i%%length(shading.col)+1])
			}
		
		}
		rect(xright=xlim[2], xleft=xlim[1], ytop=ylim[2], ybottom=plotBottom)
	}
	# force box()	
	if(gap>0){
		rect(xleft=xlim[1], xright=xlim[2], ytop=plotBottom, ybottom=boxesTop,border="white")
		abline(h=c(boxesTop, plotBottom))
	}
	
}


#' Quantile plot of time series 
#'
#' This intermediate-level function will plot a time series with the quantiles shown with transparency values.
#'
#' @param interpolate \code{(logical)}: In case the symmetric method is chosen, the series of quantile values can be interpolated with a LOESS function. 
#' @param x \code{(numeric)}: The x coordinates.
#' @param y \code{(numeric matrix)}: The series of distributions to be plotted. Every row represents a distribution of values. The number of rows must equal to the length of \code{x}. 
#' @param res \code{(numeric)}: If a single value is entered, than this argument represents the number of quantiles to be shown (coerced to 150, if higher is entered). If it is vector of values, it will be interpreted as the vector of quantiles to be shown. If \code{method="symmetric"}, only an odd number of quantiles are plotted. 
#' @param border \code{(character)}: The color of the quantile lines. A single value, by default, no lines are drawn (\code{border=NA}).
#' @param col \code{(character)}: The color of the quantiles, currently just a single color is allowed.
#' @param method \code{(character)}: The default \code{"symmetric"} method will plot the mid quantile range with highest opacity and the shades will be more translucent at the tails of the distributions. The \code{"decrease"} method will decrease the opacity with higher quantiles, which can make more sense for bottom-bounded distributions (e.g. exponential).
#' @param na.rm \code{(logical)}: If set to \code{FALSE}, than rows that are missing from the dataset will be plotted as gaps in the shading. If set to \code{TRUE}, than these gaps will be skipped. 
#' @examples
#' # some random values accross the Phanerozoic
#'	data(stages)
#'	plotTS(stages, boxes="per", shading="series", ylim=c(-5,5), ylab=c("normal distributions"))
#'	  randVar <- t(sapply(1:95, FUN=function(x){rnorm(150, 0,1)}))
#'	  shades(stages$mid, randVar, col="blue", res=10,method="symmetric")
#'	  
#' # a bottom-bounded distribution (log normal)
#'	plotTS(stages, boxes="per", shading="series", ylim=c(0,30), ylab="log-normal distributions")
#'	  randVar <- t(sapply(1:95, FUN=function(x){rlnorm(150, 0,1)}))
#'	  shades(stages$mid, randVar, col="blue", res=c(0,0.33, 0.66, 1),method="decrease")	 
#' @export
shades <- function(x, y, col, res=10, border=NA,interpolate=F, method="symmetric",na.rm=FALSE){
	if(nrow(y)!=length(x)) stop("length of x and y don't match")
	
	# omit missing?
	if(na.rm){
		logVector<- apply(y, 1, function(a){
			if(sum(is.na(a))==length(a)) FALSE else TRUE
		})
		
		x<-x[logVector]
		y<-y[logVector,]
	}
	
	if(length(res)==1){
		if(res>150) res<-150
	}else{
		# coerce interpolation
		interpolate <- FALSE
	}
	
	if(length(col)>1 | length(border)>1) stop("Please enter only one color as the col and border arguments.")
	
	
	if(interpolate){
		
		colRes<-100*2
		colQ<-seq(0,1,length.out=colRes)
		
		
		resolved<-apply(y, 1, function(x){
			x<-x[!is.na(x)]
			if(length(x)>3){
				quant<-quantile(x,colQ)
				lmod<-loess(quant ~ colQ, span=0.4)
				tap<-predict(lmod, colQ)
			return(tap)
			}else{
				return(rep(NA, colRes))
			}
		
		})
		
		
		
		for(i in 1:ncol(resolved)){
			if(!is.na(sum(resolved[,i]))){
				ran<-range(y[i,], na.rm=T)
				resolved[resolved[,i]<ran[1],i]<-ran[1]
				resolved[resolved[,i]>ran[2],i]<-ran[2]
			}
		}
			
		
	
		hex<-t(col2rgb(col))/255
		hex<-rgb(red=hex[1],green=hex[2],blue=hex[3],alpha=1/(colRes/2))

	}else{
		
		
		if(length(res)>1){
			qS<-res
			res<-length(qS)
		}else{
			if(res%%2==1){
				res<-res-1
			}
			qS<-seq(0,1,length.out=res)
			
			
		}
		
		
		
		resolved<-apply(y, 1, function(x){
			x<-x[!is.na(x)]
			if(length(x)>3){
				quant<-quantile(x,qS)
				return(quant)
			}else{
				return(rep(NA, res))
			}
		
		})
		
		hex<-t(col2rgb(col))/255
		if(method=="symmetric"){
	
			hex<-rgb(red=hex[1],green=hex[2],blue=hex[3],alpha=1/(res/2))
		}
		
		if(method=="decrease"){
			hex<-rgb(red=hex[1],green=hex[2],blue=hex[3],alpha=1/(res))
		
		}
	
	}
	
	
	colSums<-apply(resolved, 2, sum)
	streaks<-!is.na(colSums)
	sL<-streaklog(streaks)
	
	for(j in 1:sL$runs){
		bSelect<-rep(FALSE, length(x))
		bSelect[sL$starts[j]:(sL$starts[j]+sL$streaks[j]-1)] <- TRUE
		
		xSub<-x[bSelect]
		resSub<-resolved[,bSelect, drop=FALSE]
		
		# in case a single row is subsetted
		if(ncol(resSub)==1){
			# the x coordinate
			actInd<-which(bSelect)
			xPrev<-x[actInd-1]
			xNext<-x[actInd+1]
			
			xPrev<-xSub-0.2*(xSub-xPrev)
			xNext<-xSub-0.2*(xSub-xNext)
			
			xSub<-c(xPrev, xNext)
			
			# double
			resSub<-cbind(resSub, resSub)
		
		
		}
		
		if(method=="symmetric"){
			for(i in 1:(nrow(resSub)/2)){
				xNew<-c(xSub,rev(xSub))
				yNew<-c(resSub[i,],rev(resSub[nrow(resSub)-i+1,])	)
			#	polygon(xNew,yNew, col=allCols[i+1], border=F)
				polygon(xNew,yNew, col=hex, border=border)
			}
		}
		
		if(method=="decrease"){
		for(i in nrow(resSub):1){
				xNew<-c(xSub,rev(xSub))
				yNew<-c(resSub[1,],rev(resSub[nrow(resSub)-i+1,])	)
				polygon(xNew,yNew, col=hex, border=border)
			}
		
		}
	}
	
}
	


#' Plot time series counts or proportions as polygons
#' 
#' This function plots the changing shares of categories in association with an independent variable. 
#' 
#' This function is useful for displaying the changing proportions of a category as time progresses. Check out the examples for the most frequent implementations.
#' 
#' To be added: missing portions are omitted in this version, but should be represented as gaps in the polygons. 
#' 
#' @param x \code{(numeric)}: The independent variable through which the proportion is tracked. Identical entries are used to assess which values belong together to a set. Their values represent the x coordinate over the plot.
#' 
#' @param b (\code{character} or \code{factor}): A single vector with the category designations. This vector will be segmented using the entries of \code{x}.
#' 
#' @param ord \code{(character)}: The parameter of the variable order. Either \code{"up"} (increasing alphabetical order), \code{"down"} (decreasing alphabetical order) or the vector of categories in the desired order.
#' 
#' @param col \code{(character)}: The color of polygons, has to a be a vector with as many entries as there are categories in \code{b}. By default \code{(col=NULL)} this is grayscale.
#' 
#' @param border \code{(character)}: The a single color of the polygon borders. By default (\code{border=NA}), no borders are drawn. 
#' 
#' @param prop \code{(logical)}: Should the diagram show proportions (\code{TRUE}) or counts (\code{FALSE})?
#' 
#' @param xlim \code{(numeric)}: Two values, analogous to the \code{xlim} argument of \code{\link[graphics]{plot}}, and has to exceed the range of \code{x}. The polygons that represent non-zero values with the lowest and highest values of \code{x} will be extended to these \code{x} coordinates. 
#' 
#' @param ylim \code{(numeric)}: If \code{prop=TRUE}, then the argument controls the position of the proportions in the plotting area (useful to show proportions as a sub plot in a plot). If \code{prop=FALSE}, then the entire plotting area will be shifted by a single \code{ylim} value.
#' 
#' @param labs \code{(logical)}: Should the category names be plotted?
#' 
#' @param na.valid \code{(logical)}: If \code{TRUE}, than the missing values will be treated as an independent category. Entries where \code{x} is \code{NA} will be omitted either way.
#' 
#' @param labs.args \code{(list)}: Arguments for the \code{\link[graphics]{text}} function. If one entry for each argument is provided, then it will be applied to all labels. If the number of elements in an argument equals the number of categories to be plotted, then one to one assignment will be used. For example, for 4 categories in total, if the \code{labs.args} \code{list} contains a \code{col} vector element with 4 values, see examples).
#' 
#' @param plot \code{(logical)}: If set to \code{TRUE}, then the function will plot the output. If set to \code{FALSE}, then a matrix with the relevant values will be returned. This output is similar to the output of \code{\link[base]{table}}, but handles proportions instantly.
#' 
#' @param vertical \code{(logical)}: Horizontal or vertical plotting? If \code{FALSE}, the independent variable will be horizontal, if \code{TRUE}, the count/proportion variable will be horizontal. In the latter case \code{xlim} and \code{ylim} has reversed roles.
#' @examples
#' 
#' # dummy examples 
#'   # independent variable
#'   slc<-c(rep(1, 5), rep(2,7), rep(3,6))
#' 
#'   # the categories as they change
#'   v1<-c("a", "a", "b", "c", "c") # 1
#'   v2<-c("a", "b", "b", "b", "c", "d", "d") # 2
#'   v3<-c("a", "a", "a", "c", "c", "d") #3
#'   va<-c(v1, v2,v3)
#' 
#'   # basic function
#'     plot(NULL, NULL, ylim=c(0,1), xlim=c(0.5, 3.5))
#'     parts(slc, va, prop=T)
#'  
#'   # vertical plot
#'     plot(NULL, NULL, xlim=c(0,1), ylim=c(0.5, 3.5))
#'     parts(slc, va, col=c("red" ,"blue", "green", "orange"), xlim=c(0.5,3.5), labs=T, prop=T, vertical=T)
#' 
#'   # intensive argumentation
#'     plot(NULL, NULL, ylim=c(0,10), xlim=c(0.5, 3.5))
#'     parts(slc, va, ord=c("b", "c", "d", "a"), col=c("red" ,"blue", "green", "orange"), 
#' 	  xlim=c(0.5,3.5), labs=T, prop=F, 
#' 	  labs.args=list(cex=1.3, col=c("black", "orange", "red", "blue")))
#' 
#'   # just the values
#'     parts(slc, va, prop=T,plot=F)
#' 	
#' # real example
#'   # the proportion of coral occurrences through time in terms of bathymetry
#'   data(corals)
#'   data(stages)
#' 
#'   # time scale plot
#'   plotTS(stages, shading="series", boxes="per", xlim=c(250,0), 
#'     ylab="proportion of occurrences", ylim=c(0,1))
#'   
#'   # plot of proportions	
#'   cols <- c("#55555588","#88888888", "#BBBBBB88")
#'   types <- c("uk", "shal", "deep")
#'   
#'   parts(x=stages$mid[corals$slc], b=corals$bath, 
#'    ord=types, col=cols, prop=T,border=NA, labs=F)
#'    
#'   # legend
#'   legend("left", inset=c(0.1,0), legend=c("unknown", "shallow", "deep"), fill=cols, bg="white", cex=1.4) 
#' 
#' @export
parts<-function(x, b=NULL, ord="up", prop=F, plot=TRUE,  col=NULL, xlim=NULL, border=NULL, ylim=c(0,1), na.valid=FALSE, labs=T, labs.args=NULL, vertical=F){
	
	# starting arguments
	if(is.matrix(x)){
		if(ncol(x==2)){
			b<-x[,2]
			x<-x[,1]
		}
	}
	if(length(x)!=length(b)) stop("x and b have to have the same length.")
	
	#filter NAs
	if(!na.valid){
		bothNA<-is.na(b) | is.na(x)
		b<-b[!bothNA]
		x<-x[!bothNA]
	}else{
		# remove NAs in x
		b<-b[!is.na(x)]
		x<-x[!is.na(x)]
		# an assign another
		b[is.na(b)]<-"N/A"
		
	}
		
	# the unique entries
	bLevs<-unique(b)
	
	# the x values
	# increasing order
	if(!is.null(xlim)){
		if(xlim[1]<xlim[2]){
			x2<-sort(x)
		}
		if(xlim[1]>xlim[2]){
			x2<-sort(x, decreasing=T)
		}
	}else{
		x2<-sort(x)
		xlim<-range(x2)
		xR<-diff(xlim)
		xlim[1]<-xlim[1]-0.1*xR
		xlim[2]<-xlim[2]+0.1*xR
	
	}
	# the appropriate order of the bins
	xLev<-unique(x2)
	
	
	# the order of plotting of categoreies
	if(sum(bLevs%in%ord)==length(bLevs)){
		bLevs<-ord
	}else{
		if(length(ord)==1){
			if(ord=="up"){
				bLevs<-sort(bLevs)
			}
			if(ord=="down"){
				bLevs<-sort(bLevs, decreasing=T)
			}
		}else{
			stop("Invalid ord argument.")
		}
	}
	
	# arguments of the labels
	if(labs){
		if(!is.null(labs.args)){
			labs.args<-lapply(labs.args, function(x){
				if(length(x)!=length(bLevs) & length(x)!=1) stop("Provide each category lab an individual entry in arguments or one.")
				if(length(x)==1){
					return(rep(x, length(bLevs)))
				}
				if(length(x)==length(bLevs)){
					return(x)
				}
			
			})
		}
	}
	
	# default colors
	if(is.null(col)){
		dif<-floor(1/(length(bLevs)-1)*100)
		grayVal<-seq(0, 100, by=dif)
		col<-paste("gray", grayVal, sep="")
		names(col)<-bLevs
	}else{
		if(!length(col)==length(bLevs)) stop("The number of colour in 'col' doesn't match the number of categories.")
		names(col)<-bLevs
	
	}
	if(!is.null(border)){
		if(length(border)==1){
			border<-rep(border, length(bLevs))
			names(border)<-bLevs
		}else{
			if(!length(border)==length(ord)) stop("The number of colours in 'border' doesn't match the number of categories.")
		}
	}
	
	allProps<-tapply(INDEX=x, X=b, function(y){
		oneres<-sapply(bLevs, function(z){
			sum(z==y)
		})
		
		#if plotted, return the cumulative
		if(plot){
			if(prop){
				cumsum(oneres/length(y))
			}else{
				cumsum(oneres)
			}
		# if not plotted, then return the values
		}else{
			if(prop){
				oneres/length(y)
			}else{
				oneres
			}
		}
	})
	
	# reorder the by x
	allProps<-allProps[as.character(xLev)]
	
	# the labels
	if(labs){
		theProps<-tapply(INDEX=x, X=b, function(y){
			oneres<-sapply(bLevs, function(z){
				sum(z==y)
			})
			
			if(prop){
				oneres/length(y)
			}else{
				oneres			
			}
		})
		
		theProps<-theProps[as.character(xLev)]
		
		labPos <-sapply(bLevs, function(y){
			lap<-lapply(theProps, function(z){
				z[y]
			})
			res<-unlist(lap)
			xLev[which(res==max(res))[1]]
		})
	}
	
	if(prop){
		yRange<-range(ylim)
		if(length(diff(yRange))>1) stop('ylim argument must have two entries')
		if(diff(yRange)<0) stop("ylim values have to be ascending")
	}
	
	# get the values out
	allCum <-sapply(bLevs, function(y){
		lap<-lapply(allProps, function(z){
			z[y]
		})
		res<-unlist(lap)
		names(res)<-names(lap)
		res
	})
	
	# return values if plot is FALSE
	if(!plot){
		return(allCum)
	}
	
	# rescale if necessary	
	if(prop){
		allCum<-allCum*diff(yRange)
	}
	
	# and shift
	allCum <- allCum+ylim[1]
		
	for(i in length(bLevs):1){
		if(i>1){
			plotY<-c(
				allCum[1,bLevs[i]], 
				allCum[,bLevs[i]], 
				allCum[nrow(allCum),bLevs[i]],
				allCum[nrow(allCum),bLevs[i-1]],
				rev(allCum[,bLevs[i-1]]), 
				allCum[1,bLevs[i-1]]
			)
		
		}else{
			plotY<-c(
				allCum[1,bLevs[i]], 
				allCum[,bLevs[i]], 
				allCum[nrow(allCum),bLevs[i]],
				rep(ylim[1], length(xLev)+2))
		}
		# the xlims
		plotX<-c(xlim[1], xLev, xlim[2], xlim[2], rev(xLev), xlim[1])
	
		if(vertical){
			polygon(plotY, plotX, col=col[bLevs[i]], border=border[bLevs[i]])
		}else{
			polygon(plotX, plotY, col=col[bLevs[i]], border=border[bLevs[i]])
		}
	}
	
	if(labs){
		
		for(i in 1:length(labPos)){
			xCoord<-labPos[i]
			
			yVect<-allCum[as.character(labPos[i]), ]
			
			yPos<-which(names(yVect)==names(labPos)[i])
			if(yPos==1){
				yCoord<-mean(c(yVect[yPos],ylim[1]))
			}else{
				yCoord<-mean(c(yVect[yPos],yVect[yPos-1]))
			}
			if(vertical){
				newArgs<-list(
					label=names(labPos)[i],
					y=xCoord,
					x=yCoord
				)	
			}else{
				newArgs<-list(
					label=names(labPos)[i],
					x=xCoord,
					y=yCoord
				)	
			}
			
			plusArgs<-lapply(labs.args, function(u){
				u[i]
			})
			
			# plot the labels
			do.call(text, c(newArgs, plusArgs))
		}
	
	}
}

	
	
	