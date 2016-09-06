#' Denmark SUTs
#'
#' Danish supply and use tables. The use table is split into domestic
#' and import use.
#' @name dk
#' @docType data
#' @format A list with a supply table (V, a matrix) and a separate
#'     import vector (m0), domestic use table (Ud, matrix) and import
#'     use table (Um, matrix). Furthermore, the country and time of
#'     the SUT system.
#' @author Oliver Reiter \email{oliver.reiter@snapdragon.cc}
#' @references \url{http://ec.europa.eu/eurostat/data/database}
#' @keywords data
"dk"

#' Danish national accounts data 
#'
#' Danish national accounts data, aka projection data that is used to
#' update the SUT.
#' @name dk.2011
#' @docType data
#' @format A list with data on sector-level: gross output (xbar,
#'     vector), intermediate consumption (ubar, vector). Additionally
#'     a vector "c0" with the absolute values of the negative entries
#'     of the trade and transport margins (see SUTRAS paper), as well
#'     as total imports.
#' @author Oliver Reiter \email{oliver.reiter@snapdragon.cc}
#' @references \url{http://ec.europa.eu/eurostat/data/database}
#' @keywords data
"dk.2011"


## setwd("/home/reitero/Arbeit/Rprogramming/riot/data")
## load("DK_SUTs_2010t2011.RData")
## load("DK_ProjData.RData")

## save(dk, file = "DK_SUT_2010.RData")
## save(dk.2011, file = "DK_ProjData_2011.RData")
