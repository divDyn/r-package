## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(divDyn)

## ----package, echo= TRUE, eval=FALSE------------------------------------------
#  library(divDyn)

## ----dat, echo= TRUE----------------------------------------------------------
# attach the time scale object
data(stages)

## ----tsplot1, echo= TRUE, plot=TRUE, fig.height=5.5---------------------------
tsplot(stages, boxes="sys", shading= "series")

## ----tsplot2, echo= TRUE, plot=TRUE, fig.height=5.5---------------------------
tsplot(stages, boxes="system", shading="stage", xlim=59:81)

## ----tsplot3, echo= TRUE, plot=TRUE, fig.height=5.5---------------------------
tsplot(stages, boxes="sys", shading="series", 
  labels.args=list(col="red", font=3), shading.col=c("white", "wheat"))

## ----tsplot4, echo= TRUE, plot=TRUE, fig.height=5.5---------------------------
tsplot(stages, boxes=c("sys"), shading="sys", boxes.col="systemCol")

## ----tsplot5, echo= TRUE, plot=TRUE, fig.height=5.5---------------------------
tsplot(stages, boxes=c("short","system"), shading="short", 
  xlim=59:69, boxes.col=c("col","systemCol"), labels.args=list(cex=0.5))

## ----tsplot6, echo= TRUE, plot=TRUE, fig.height=5.5---------------------------
tsplot(stages, boxes=c("short","system"), shading="short", 
  xlim=59:69, boxes.col=c("col","systemCol"),
  labels.args=list(list(cex=0.5),list(cex=1)),
  boxes.args=list(list(),list(density=80)))

## ----examp, echo= FALSE, results=TRUE-----------------------------------------
structure <- data.frame(
  tax= c("Sp1", "Sp1", "Sp1", "Sp2","Sp3", "Sp2"),
  bin= c(1,1,2, 2,2,3),
  locality = c("first", "second", "second", "second", "first", "second"))
  
  structure

## ----corDat, echo= TRUE-------------------------------------------------------
data(corals)

## ----tibbleOrNot, echo= FALSE, eval=FALSE-------------------------------------
#  library(tibble)
#  corals <- as_tibble(corals)

## ----fossils, echo= TRUE, results=TRUE----------------------------------------
fossils <- corals[corals$stg!=95,]
# the number of occurrences
nrow(fossils)

## ----fadlad, echo= TRUE, results=FALSE----------------------------------------
fl <- fadlad(fossils, bin="stg", tax="genus")

## ----ranplot, echo= TRUE, plot=TRUE, fig.height=5.5---------------------------
fossils$mid <- stages$mid[fossils$stg]

tsplot(stages, shading="series", boxes="sys",xlim=c(260,0))

ranges(fossils, tax="genus", bin="mid")

## ----ranplot2, echo= TRUE, plot=TRUE, fig.height=5.5--------------------------
tsplot(stages, shading="series", boxes="series",xlim=c(100,81))

ranges(fossils, tax="genus", bin="mid", labs=T, labels.args=list(cex=0.2))

## ----ranplot3, echo= TRUE, plot=TRUE, fig.height=5.5--------------------------
tsplot(stages, shading="stage", boxes=c("stage","series"),xlim=c(100,80))

ranges(fossils, tax="genus", bin="mid", labs=T, 
  labels.args=list(cex=0.6), filt="orig", occs=T)

## ----ranplot4, echo= TRUE, plot=TRUE, fig.height=5.5--------------------------
tsplot(stages, shading="stage", boxes=c("stage","series"),xlim=c(100,80))

ranges(fossils, tax="genus", bin="mid", labs=T, 
  labels.args=list(cex=0.6), filt="orig", occs=T, group="ecology")

## ----surv, echo= TRUE, results=FALSE------------------------------------------
surv <- survivors(corals, bin="stg", tax="genus")

## ----survTsplotblank, echo=TRUE, plot=TRUE, fig.height=5.5--------------------
# time scale plot
tsplot(stages, shading="series", boxes="sys", 
  xlim=c(260,0), ylab="proportion of survivors present",
  ylim=c(0.01,1),plot.args=list(log="y"))

## ----survTsplotShow, echo=TRUE, eval=FALSE------------------------------------
#  # lines for every cohort
#  for(i in 1:ncol(surv)) lines(stages$mid, surv[,i])

## ----survTsplot1, echo=FALSE, plot=TRUE, fig.height=5.5-----------------------
# time scale plot
tsplot(stages, shading="series", boxes="sys", xlim=c(260,0),
  ylab="proportion of survivors present", 
  ylim=c(0.01,1),plot.args=list(log="y"))
# lines for every cohort
for(i in 1:ncol(surv)) lines(stages$mid, surv[,i])

## ----survTsplot2, echo=FALSE, plot=TRUE, fig.height=5.5-----------------------
survFor<-survivors(corals, tax="genus", bin="stg", method="backward")
# plot
tsplot(stages, shading="series", boxes="sys", 
  xlim=c(260,0), ylab="proportion of constituents present",
  ylim=c(0.001,1),plot.args=list(log="y"))
# for every cohort
for(i in 1:ncol(surv)) lines(stages$mid, survFor [,i])

## ----sampTOT, echo=TRUE, result=TRUE------------------------------------------
samp <-sumstat(fossils, tax="genus", bin="stg", 
  coll="collection_no", ref="reference_no", duplicates=FALSE)
samp

## ----binstat, echo=TRUE, result=FALSE-----------------------------------------
samp <-binstat(fossils, tax="genus", bin="stg", 
  coll="collection_no", ref="reference_no")

## ----binstatInd, echo=TRUE, result=FALSE--------------------------------------
samp <-binstat(fossils, tax="genus", bin="stg", 
  coll="collection_no", ref="reference_no", indices=TRUE)
colnames(samp)

## ----occsAndColls, echo=TRUE, result=FALSE, plot=TRUE, fig.height=5.5---------
oldPar <- par(mar=c(4,4,2,4))
# basic plot
tsplot(stages, shading="series", boxes="sys", xlim=52:95,
    ylab="number of occurrences", ylim=c(0,3000))
lines(stages$mid[1:94], samp$occs)
# the collections (rescaled, other axis)
  lines(stages$mid[1:94], samp$colls*5, col="blue")
  axis(4, col="blue",col.ticks="blue",col.axis="blue",
    at=seq(0,3000,500), labels=seq(0,600,100))
  mtext(4, text="number of collections", col="blue", line=2)
par(oldPar)

## ----parts, echo=TRUE, result=FALSE, plot=TRUE, fig.height=5.5----------------
# numerical ages, as bins
fossils$stgMid <- stages$mid[fossils$stg]
#plotting
tsplot(stages, shading="series", boxes="sys", xlim=52:95,
  ylab="number of occurences", ylim=c(0,3000))
parts(fossils$stgMid, fossils$bath)

## ----parts2, echo=TRUE, result=FALSE, plot=TRUE, fig.height=5.5---------------
cols <- c("#FF0000AA", "#00FF00AA", "#0000FFAA")
# reorder too
reord <- c("shal","deep","uk")
plotnames <-c("shallow", "deep", "unknown")
tsplot(stages, shading="series", boxes="sys", xlim=52:95,
    ylab="the number of occurrences", ylim=c(0,3000))
parts(fossils$stgMid, fossils$bath, col=cols, ord=reord, labs=F)
legend("topleft", inset=c(0.01, 0.01), 
  legend= plotnames, fill=cols, bg="white")

## ----parts3, echo=TRUE, result=FALSE, plot=TRUE, fig.height=5.5---------------
tsplot(stages, shading="series", boxes="sys", xlim=52:95,
    ylab="proportion of occurrences", ylim=c(0,1))
parts(fossils$stgMid, fossils$bath, prop=T, col=cols, ord=reord, labs=F)
legend("bottomleft", inset=c(0.01, 0.1), 
  legend= plotnames, fill=cols, bg="white")

## ----omit, echo=TRUE, result=TRUE---------------------------------------------
omitColl <- omit(corals, tax="genus", om="coll", coll="collection_no")
omitRef <- omit(corals, tax="genus", om="ref", ref="reference_no")
omitBinref <- omit(corals, bin="stg", tax="genus", om="binref", ref="reference_no")
# the conserved number of occurrences will be
sum(!omitColl)
sum(!omitRef)
sum(!omitBinref)

## ----ddFirst, echo=TRUE, result=FALSE-----------------------------------------
ddFirst<-divDyn(corals, bin="stg", tax="genus", noNAStart=TRUE)

## ----ddRec, echo=TRUE, result=FALSE, plot=TRUE, fig.height=5.5----------------
# metrics
ddRec <-divDyn(corals, bin="stg", tax="genus")
# basic plot
  tsplot(stages, shading="series", boxes="sys", xlim=52:95,
    ylab="range-through richness (diversity)", ylim=c(0,250))
# lines
  lines(stages$mid, ddRec$divRT, col="black", lwd=2)

## ----ddFossils, echo=TRUE, result=FALSE, plot=TRUE, fig.height=5.5------------
# metrics
dd <-divDyn(fossils, bin="stg", tax="genus")
# basic plot
  tsplot(stages, shading="series", boxes="sys", xlim=52:95,
    ylab="richness (diversity)", ylim=c(0,250))
# lines
  lines(stages$mid, ddRec$divRT, col="black", lwd=2)
  lines(stages$mid[1:94], dd$divRT, col="blue", lwd=2)
# legend
  legend("topleft", legend=c("with recent", "without recent"), col=c("black", "blue"), lwd=c(2,2), bg="white")

## ----ddDiv, echo=TRUE, result=FALSE, plot=TRUE, fig.height=5.5----------------
# metrics
dd <-divDyn(fossils, bin="stg", tax="genus")
# basic plot
  tsplot(stages, shading="series", boxes="sys", xlim=52:95,
    ylab="richness (diversity)", ylim=c(0,250))
# lines
  lines(stages$mid[1:94], dd$divRT, col="red", lwd=2)
  lines(stages$mid[1:94], dd$divBC, col="blue", lwd=2)
  lines(stages$mid[1:94], dd$divSIB, col="green", lwd=2)
# legend
  legend("topleft", legend=c("RT", "BC", "SIB"), 
    col=c("red", "blue",   "green"), lwd=c(2,2,2), bg="white")

## ----samp3t, echo=TRUE, result=FALSE, plot=TRUE, fig.height=5.5---------------
# basic plot
  tsplot(stages, shading="series", boxes="sys", xlim=52:95,
    ylab="three-timer sampling completeness")
  # lines
  lines(stages$mid[1:94], dd$samp3t, col="black", lwd=2)

## ----ddCorr, echo=TRUE, result=FALSE, plot=TRUE, fig.height=5.5---------------
# basic plot
  tsplot(stages, shading="series", boxes="sys", xlim=52:95,
    ylab="three-timer sampling completeness", ylim=c(0,250))
  lines(stages$mid[1:94], dd$divSIB, col="black", lwd=2)
  lines(stages$mid[1:94], dd$divCSIB, col="blue", lwd=2)
# legend
  legend("topleft", legend=c("SIB", "corrected SIB" ), 
    col=c("black", "blue"), lwd=c(2,2), bg="white")

## ----ddExt, echo=TRUE, result=FALSE, plot=TRUE, fig.height=5.5----------------
# basic plot
  tsplot(stages, shading="series", boxes="sys", xlim=52:95,
    ylab="extinction rates", ylim=c(0,2))
  lines(stages$mid[1:94], dd$extPC, col="black", lwd=2)
  lines(stages$mid[1:94], dd$extGF, col="blue", lwd=2)
  lines(stages$mid[1:94], dd$ext2f3, col="green", lwd=2)
  # legend
  legend("topright", legend=c("per capita", "gap-filler", "second-for-third"),
    col=c("black", "blue", "green"), lwd=c(2,2,2), bg="white")

