#' Time series plotting using a custom timescale object
#'
#' This function allows the user to quickly plot a time scale data table
#'
#' As most analysis use an individually compiled time scale object, in order to ensure compatibility between the analyzed and plotted values, the time scale table used for the analysis should be plotted rather than a standardized table.
#' @param tsdat (data frame): time scale data frame
#' @param boxes (character): column name indicating the names that should be plotted as a timescale
#' @param ylim (numeric): the vertical extent of the plot, analogous to the same argument of plot(). By default it is set to the [0,1] interval.
#' @param xlim (numeric): the horizontal extent of the plot, analogous to the same argument of plot(). By default it is set to  plot the entire table. If a numeric vector of length=2 is supplied, it will be interpreted as the standard xlim argument and the plot will be displayed based on numerically constrained ages. If it an integer vector with length>2, the interval corresponding to the row indices of the table will be plotted.
#' @param prop (numeric value): proportion of the vertical extent of the plot for the timescale
#' @param gap (numeric value): proportion of the vertical extent of the plot that should be a gap betwen the timescale and the plot.
#' @param bottom (character): column name of the table for the variable that contains the older ages of intervals.
#' @param mid (character): column name of the table for the variable that contains the mean ages of intervals.
#' @param top (character): column name of the table for the variable that contains the earliest ages of intervals.
#' @param shading (character): column name used for the shading.
#' @param shading.col (character): name of colors that will be used for the shading, if shading is set.
#' @param plot.args (list): list of arguments that will be passed to the main plot() function.  Can be useful for the suppression of axes, font change etc.
#' @param labels.args (list): list of arguments that will be passed to the text() function that draws the labels. 
#' @param boxes.args (list): list of arguments that will be passed to the rect() function that draws the rectangles of time intervals.
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
	bottom="bottom", mid="mid",top="top",
	xlab="age (Ma)", ylab="",
	shading=NULL,shading.col=c("white", "gray80"),
	plot.args=NULL,
	boxes.args=NULL,
	labels.args=NULL){
	
#	tsdat<-stages
#	boxes<-"per"
#	bottom="bottom"
#	mid="mid"
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
	
	if(length(shading)>1){
		stop("Only one column can be used to plot the shades.")
		
	}else{
		if(!shading%in%colnames(tsdat)){
			stop("The referenced 'shading' column does not exist.")
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
	boxesTop<-ylim[1]-diff(range(ylim))*gap
	plotBottom<-ylim[1]
	
	ylim[1]<-ylim[1]-diff(range(ylim))*prop-diff(range(ylim))*gap
	
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
			x=NULL,
			y=NULL,
			xlim=xlim, 
			ylim=ylim,
			xlab=xlab, 
			ylab=ylab,
			xaxs="i",
			yaxs="i"
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
#' This intermediate-level function will plot a time series with the quantiles shown as shades (using alpha values) around the central tendency.
#'
#' @param interpolate (logical value): In case the symmetric method is chosen, the series of quantile values can be interpolated with a LOESS function. 
#' @param boxes (character): column name indicating the names that should be plotted as a timescale
#' @param x (numeric): The x coordinates.
#' @param y (numeric matrix): The series of distributions to be plotted. Every row represents a distribution of values. The number of rows must equal to the length of x. 
#' @param res (numeric value): If a single value is entered, than it represents the number of quantiles to be shown (coerced to 150, if higher is entered). If it is vector of values, it will be interpreted as the vector of quantiles to be shown. If method="symmetric", only an odd number of quantiles are plotted. 
#' @param border (character value): the color of the quantile lines
#' @param col (character value): the color of the quantiles, currently just a single color is allowed.
#' @param method (character value): The default "symmetric" method will plot the mid quantile range with highest opacity and the shades will be more translucent at the tails of the distributions. The "decrease" method will decrease the opacity with higher quantiles, which makes the plotts of low-bounded distributions easier to interpret.
#' @examples
#'	data(stages)
#'	plotTS(stages, boxes="per", shading="series", ylim=c(-5,5), ylim=c(normal distributions))
#'	  randVar <- t(sapply(1:95, FUN=function(x){rnorm(150, 0,1)}))
#'	  shades(stages$mid, randVar, col="blue", res=20,method="symmetric")	  
#'	  
#'	plotTS(stages, boxes="per", shading="series", ylim=c(0,30), ylab="log-normal distributions")
#	  randVar <- t(sapply(1:95, FUN=function(x){rlnorm(150, 0,1)}))
#'	  shades(stages$mid, randVar, col="blue", res=c(0,0.33, 0.66, 1),method="decrease")	 
#' @export
shades <- function(x, y, col, res=100, border=NA,interpolate=F, method="symmetric"){
	if(nrow(y)!=length(x)) stop("length of x and y don't match")
	if(length(res)==1){
		if(res>150) res<-150
	}else{
		# coerce interpolation
		interpolate <- FALSE
	}
	
	if(length(col)>1) stop("Please enter only one color.")
	
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
		resSub<-resolved[,bSelect]
		
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
	