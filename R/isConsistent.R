#' Function to check the consistency of a SUT system
#'
#' Given a SUT system containing a supply table and a domestic and a
#' import use table, this function checks if the rows and columns are
#' consistent.
#' @param SUT the SUT system to be tested
#' @param tol the tolerance (epsilon), default = 1e-4
#' @param VA.included is a value-added included in the SUT system? If
#'     yes, then the function also checks for consistent
#'     columns. Default = FALSE.
#' @return If the SUT is consistent, returns only TRUE. If not,
#'     returns FALSE and a list of inconsistent rows and columns
isConsistent <- function(SUT, tol = 1e-4,
                         VA.included = FALSE) {
    ##
    ## consistent row sums ?
    ##
    V.rowsums <- rowSums(SUT$V) + SUT$m0
    U.rowsums <- rowSums(SUT$Ud) + rowSums(SUT$Um)

    cbind(V.rowsums, U.rowsums)

    row.diff.abs <- abs(V.rowsums - U.rowsums)
    row.diff.percent <- row.diff.abs / pmax(V.rowsums, U.rowsums)
    row.diff.percent[!is.finite(row.diff.percent)] <- 0
    consistent.rows <- all(row.diff.percent < tol)

    res <- list(consistent = consistent.rows)

    if(consistent.rows == FALSE) {
        res$row.diff.percent <- row.diff.percent[row.diff.percent > tol] * 100
        res$row.diff.abs <- row.diff.abs[row.diff.percent > tol]
    }
    
    ## 
    ## consistent column sums ?
    ##
    if(VA.included) {
        V.colsums <- colSums(SUT$V)
        U.colsums <- colSums(SUT$Ud) + colSums(SUT$Um) + colSums(SUT$VA)

        col.diff.abs <- abs(V.colsums - U.colsums)
        col.diff.percent <- col.diff.abs / pmax(V.colsums, U.colsums)
        col.diff.percent[!is.finite(col.diff.percent)] <- 0
        consistent.cols <- all(col.diff.percent < tol)

        res$consistent <- res$consistent & consistent.cols

        if(consistent.cols == FALSE) {
            res$col.diff.percent <- col.diff.percent[col.diff.percent > tol] * 100
            res$col.diff.abs <- col.diff.abs[col.diff.percent > tol]
        }
    }

    return(res)
}
