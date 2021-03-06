riot
=======

**Development of this package moved to
[Gitlab](https://gitlab.com/zauster/riot), this repository is not
updated regularly (if at all).**

Installation
------------

As this package is still in development, it can not yet be installed
through CRAN.

This development version can, however, be installed using the
`devtools` package:

```r
if (!require('devtools')) install.packages('devtools')
install_github("zauster/riot")
```


Usage
-----

After installation, load the package with:

```r
library(riot)
```

See the help-files for some examples of the usage of the package.

```r
help("riot")
```

To Do
-----

1. [x] fix DOCUMENTATION!
1. [ ] class system for SUTs (own class for supply/use tables?) and/or IO, preferably using the R6 class system
1. [ ] getter/setter methods for data points
1. [ ] compute IO table from supply and use tables
1. [ ] replace missings in sut tables
1. [ ] adjust rows/cols in SUT before SUTRASing, eg, adjust exports/imports according to some external data
1. [ ] doGRAS.long tests (in general more tests)
1. [ ] Decomposing bilateral flows accordings to WWZ in RcppArmadillo
