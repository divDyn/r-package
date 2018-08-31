# divDyn
R package for analyzing diversity dynamics from fossil occurrence data

If you are interested in what this package does, or have questions about use, please check out the first vignette: 'Handout to the R package 'divDyn' v0.5.0 for diversity dynamics from fossil occurrence data' here:
https://github.com/adamkocsis/divDyn/blob/master/_archive/vignettes/Handout_0.5.1.pdf

As the package is still getting developed, if you have any requirements or recommendations about what to add (or if you find a mistake), please contact me at adam.kocsis@fau.de.


## Installing and update
As the package now uses compiled code to make things run faster, there are two ways to install.
1. The first option is to use the binaries I have built. 

I have updated the windows binaries so that they work with the latest internals (3.5.1). If you want to use the package, please update your R.

If you have a windows computer, you can install the package with the following R command:
`install.packages("https://github.com/adamkocsis/divDyn/raw/master/_bin/Win_x64_x86/divDyn_0.5.1.zip", repos=NULL)`

If you have a Mac running OS X, then use the following link (still older version, for R <3.5):
`install.packages("https://github.com/adamkocsis/divDyn/raw/master/_bin/Mac_OSX/divDyn_0.5.0.tgz", repos=NULL)`

If you have Linux computer you probably know how to solve these problems. The sources of different versions are in the /_archive folder. Otherwise you can try option no. 2. 

You can access earlier versions by changing the version number appropriately.

2. The second option is to compile the code for yourself. To do this:
- Install a compiler. For Windows, this would be included in Rtools (https://cran.r-project.org/bin/windows/Rtools/). 
- Make sure that the 'devtools' package is installed
- Run `devtools::install_github("adamkocsis/divDyn")`

 '*' if you use either a Mac or Linux, and point no. 2 doesn't work for you for some reason, then contact me  and I will compile a binary for you.

# Change log

## [0.5.1] - 2018.08.31
### Changed
- example bin identifier changes (from 'slc' to 'stg')
- the subsample(), subtrialCR(), subtrialOXW(), subtrialSQS() functions now automatically omit rows that have NA entries in the 'bin' variables.
- the 'stages' object now contains dates from Gradstein et al. 2016
- bin durations were added to both the 'stages' and 'bins' objects
- the 'stratkeys' and 'keys' objects were updated to version 0.9.2

## [0.5.0-1] - 2018.08.29
### Changed
- stratigraphic assignment objects (in 'keys') now have a -1 entries for empty strings


## [0.5.0] - 2018.08.20
### Added
- bug fix that made R crash randomly
- subsampling trial functions (subtrialCR(), subtrialOXW(), subtrialSQS()) and cleaned help files for 'subsample()'
- multiple bin entry support for 'ranges()'
- proper Description in description file

### Changed
- the functions 'plotTS()' and 'fadLad()' were renamed to 'tsplot()' and 'fadlad()' for easier typing. 
- 'stratkeys' and the corresponding 'keys' updated to v0.9.1
- the 'sampstat()' function no longer outputs empty variables when the argument columns are set to NULL

### Deleted
- the 'exact' method implementation for SQS as it was slow and and incorrect. The 'inexact' method was confirmed with the iNEXT package.


## [0.4.2] - 2018.07.24
### Added
- the ranges() function to plot occurrence data and stratigraphic ranges over time (see examples)
- fadLad(): function defenses 
- fadLad(): two bin columns are now allowed as input (for age uncertainty estimation)
- fadLad(): the ages argument for inverting the time axis
- fixes to conform to the CRAN requirements

### Changed
- plotTS() now allows "no-boxes" plotting, and multiple layers of boxes


## [0.4.1] - 2018.07.07
### Added
- function defenses
- fixes to conform to the CRAN requirements

### Changed
- default settings for the sampstat() and subsample() functions to increase flexibility and user-friendliness
- the intact argument of subsample() was separated to the rem and keep arguments, to enable negative time-slice identifiers.

## [0.4.0] - 2018.06.29
### Added
- the first vignette is now included
- bug fix for the parts() function. Due to the overplotting of polygons, the RGBA colours were not appropriately visualized
- bug fix to the omit() function's binref method.
- the "total" method for the sampstat() function
- documentation to the built in datasets
- the map() function to resolve information downloaded from the PaleoDB
- 'keys' to the stratigraphic and environmental resolution, based on the FossilWorks dynamic timescale ?keys
- the 10 million year timescale of the Paleobiology Database, ?bin

### Changed
- The sampstat() function now outputs a data.frame, instead of a matrix.
- the 'method' argument of subsample() was renamed to 'type' to avoid problems with argument distribution with the frequent 'method' argument name

### Deleted 
- The previously added bayesian method from the affinity() function. The method had problems and is under investigation.


## [0.3.1] - 2018.06.04
### Added
- method argument for the affinity() function, allows implementation of a simple majority rule (method="majority") or the Bayesian approach of Simpson and Harnik (2009, method="bayesian")

### Changed
- the affinity() function was also completely reprogrammed

## [0.3.0] - 2018.06.01
### Added
- continuous time entry is now possible for the divDyn() function. (does not yet work with subsampling!!!!)
- the breaks argument to the divDyn() function

### Changed
- the examples and defaults match the current download from the Paleobiology Database
- the previous 'scleractinia' example dataset was renamed to 'corals' for the sake of simplicity

## [0.2.12] - 2018.05.17
### Added
- na.rm argument for the shades() function to omit gaps in plotting of the time series
- bug fix for the shades() function, single rows bounded by rows of NAs (gaps) from above and below do not crash the function anymore


## [0.2.11] - 2018.05.17
### Added
- bug fix for the plotTS() function, enables log scaling of y axis with plot.args=list(log="y")
- the parts() function for efficient plotting of proportion/count data over time 


## [0.2.10] - 2018.05.13
### Added
- proportional extinction and origination rates to the divDyn() function
- the omit() function for culling occurrence datasets (kind of slow yet)

### Changed
- argumentation of the divDyn function was extended to accomodate the omission of occurrences within the function (optional)


## [0.2.9] - 2018.05.09
### Changed
- the 'cr' method of the subsample() function was rewritten in C++ and with Rcpp classes for faster performance
- fixed bug in the subsample() function, when custom functions produce vector output


## [0.2.8] - 2018.05.01
### Added
- multiple versions of the shareholder quorum subsampling (SQS) routine to subsample()
- the sampstat() function, that outputs basic sampling stats for the time slices
- useFailed argument to subsample(). Changes output depending on whether the subsampling quota/quorum is reached.
- bug fix for the CR method in subsampling 

### Changes
- argument distribution in the subsample() function, FUN argument simplified to a single function

### Planned
- the subsampling routines will be remodularized so that their argument becomes the dataset, rather than the individual columns. This will allow the support for custom subsampling functions.
- 

## [0.2.7] - 2018.02.01
### Changes
- bug (due to output changes in fadLad) and speed fixes for the affinity() function. The interface also changed, it now outputs a single vector instead of a data.frame.
- bug fix in the ratesplit() function argument distribution. The "combine" version works properly now.

# Change log
## [0.2.6] - 2018.01.26
### Added
- John Alroy's (2015) second-for-third extinction/origination proportions/rates to the divDyn() function

## [0.2.5] - 2018.01.22
### Added
- survivors() function to calculate survivorship probabilities

### Changes
- the 'subseries()' function was renamed to 'subsample()' for easier application
- the 'fadLad()' function was redrafted for speed

## [0.2.4] - 2018.01.15
### Added
- the the "oxw" method for the subseries() function. This implements occurrence-weighted by list subsampling.

## [0.2.3] - 2018.01.08
### Added
- the subseries() function to execute certain functions with subsampling applied.


## [0.2.2] - 2017.12.18
### Added
- the shades() function for plotting distributions of values arranged in a time series

## [0.2.1] - 2017.12.12
### Added
- the ratesplit() function for selectivity tests
- period abbreviations for the 'stages' object

### Changes
- plotTS() function default changes


## [0.2.0] - 2017.12.11
### Changes
- added Rccp import
- divDyn() updated to v2, ca. 150x faster
- installation changes

## [0.1.6] - 2017-12-04
### Added
- the plotTS() time scale plotting function
- examples were updated

## [0.1.5] - 2017-11-27
### Added
- the streaklog() and whichmaxstreak() functions

## [0.1.4] - 2017-11-22
### Changes
- documentation

## [0.1.3] - 2017-11-21
### Changes
- output variable names in the divDyn() function

### Added
- the first version of the crDD() function is added.

### Notes
Pre-open versions were not registered.