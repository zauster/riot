#' Function to check the consistency of a SUT system
#'
#' Given a SUT system containing a supply table and a domestic and a
#' import use table, this function checks if the rows and columns are
#' consistent.
#' @param sut the SUT system to be tested
#' @param tol the tolerance (epsilon), default = 1e-4
#' @param VA.included is a value-added included in the SUT system? If
#'     yes, then the function also checks for consistent
#'     columns. Default = FALSE.
#' @return If the SUT is consistent, returns only TRUE. If not,
#'     returns FALSE and a list of inconsistent rows and columns
isConsistent <- function(sut, tol = 1e-4,
                         VA.included = FALSE) {
    ##
    ## consistent row sums ?
    ##
    V.rowsums <- rowSums(sut$V) + sut$m0
    U.rowsums <- rowSums(sut$Ud) + rowSums(sut$Um)

    row.diff <- abs(V.rowsums - U.rowsums)

    consistent.rows <- all(row.diff < tol)


    ## 
    ## consistent column sums ?
    ##
    consistent.cols <- TRUE
    col.diff <- 0
    if(VA.included) {
        V.colsums <- colSums(sut$V)
        U.colsums <- colSums(sut$Ud) + colSums(sut$Um) + colSums(sut$VA)

        col.diff <- abs(V.colsums - U.colsums)

        consistent.cols <- all(col.diff < tol)
    }

    consistent <- consistent.rows & consistent.cols

    ## if not consistent, give more verbose output
    if(consistent == FALSE) {
        consistent <- list(consistent = consistent,
                           row.diff = row.diff[row.diff > tol],
                           col.diff = col.diff[col.diff > tol])
    }

    return(consistent)
}
## load("../data/DK_SUTs_mat.RData")

## dk.2010 <- list(V = V,
##                 Ud = Ud,
##                 Um = Um,
##                 year = 2010,
##                 country = "DK")
