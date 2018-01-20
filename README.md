
riot
=======
[![Travis-CI Build Status](https://travis-ci.org/zauster/riot.png?branch=master)](https://travis-ci.org/zauster/riot)
[![codecov](https://codecov.io/gh/zauster/riot/branch/master/graph/badge.svg)](https://codecov.io/gh/zauster/riot)

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
1. [ ] class system for SUTs and/or IO, preferably using the R6 class system
1. [ ] compute IO table from supply and use tables
1. [ ] replace missings in sut tables
