# divDyn
R package for analysing diversity dynamics with Paleobiology Database occurrences

## installing and update
As the package now uses compiled code to make things run faster, there are two ways to install.
1. I have built binaries for Windows x86 and x64. If you have one of these architectures, you can install the package with:
`install.packages("https://github.com/adamkocsis/divDyn/raw/master/_bin/divDyn_0.2.2.zip", repos=NULL)`

2. You have to compile the code for yourself. For this:
- Install a compiler. For Windows, this would be included in Rtools.
- Make sure that the devtools package is installed
- Run `devtools::install_github("adamkocsis/divDyn")`




# Change log
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
Pre-alpha versions were not registered.