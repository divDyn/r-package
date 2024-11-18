# divDyn

R package for quantifying diversity dynamics using fossil sampling data

![](https://github.com/divDyn/assets/raw/master/logo/divDyn_logo_min.png)

## News

- The next update of the package (V0.8.2) is now availble from the CRAN servers! As this is a major update, I suggest everyone to update their copy. See the change log below for the changes. 

- The paper describing the package is available from the Wiley Online Library. You can download it from this link:
https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13161 
(or contact me, if you have problems). If you decide to use the methods implemented in this package, please refer to this document. The examples implemented can be reproduced with the files deposited at
https://github.com/divDyn/ddPhanero
using this (updated) vignette:
https://github.com/divDyn/ddPhanero/blob/master/doc/dd_phanero_web.html


## About the package

If you are interested in what this package does, or have questions about its use, please check out the first vignette: 'Handout to the R package 'divDyn' v0.8.2 for diversity dynamics from fossil occurrence data' here:
https://github.com/divDyn/r-package/blob/master/vignettes/handout.pdf

As the package is still getting developed, please note that some interface changes might occurr based on the feedback of people and my own experience on what is easier to use. If you have any requirements or recommendations about what to add (or more importantly, if you find a mistake), do not hesitate to contact me at adam.kocsis@fau.de.


## Installing 

### 1. From CRAN (v0.8.2)

You can install the appropriate binaries normally, by running 
`install.packages("divDyn")`

### 2. Using the source tarball to install (v0.8.2)
This you can do with running
```
install.packages(
  "https://github.com/divDyn/assets/raw/master/r_archive/source/divDyn_0.8.2.tar.gz", 
  repos=NULL, type="source")
```

from the R console. Note that there is some code in the package that requires compilation. To do this, the most straightforward way is to install/use Rtools (https://cran.r-project.org/bin/windows/Rtools/) on Windows, and XCode https://developer.apple.com/xcode/ on Mac.

The sources of the older versions are also in the _archive/source folder. You can access earlier versions by changing the version number in the command above.

### 3. Using the repository files and 'devtools' to install (v0.8.2)

To do this:
- You need a compiler, as for method 2
- Make sure that the 'devtools' package is installed
- Run `devtools::install_github("divDyn/r-package")`

If the frist method is not working for you, and you do not want to mess around with a compiler and such, then contact me and I will find a way to compile binaries for you.

# Change log

See the change log here: [Change Log](https://github.com/divDyn/r-package/blob/master/inst/NEWS)
