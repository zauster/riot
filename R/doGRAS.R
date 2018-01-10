#' Do the GRAS algorithm for a "long" data.table
#'
#' This function calculates an updated "long" data.table, based on a data.table dt,
#' which meets the given row and columns totals. Note: dt can contain
#' negative elements.
#' @param dt the "long" data.table
#' @param rowcol the column which holds the row-dimension of the reshaped matrix
#' @param colcol the column which holds the column-dimension of the reshaped matrix
#' @param rowsum vector with the row totals
#' @param colsum vector with the column totals
#' @param ... parameter to be passed to doGRAS
#' @return the updated data.table in the same format as before
#' @author Oliver Reiter
#' @references Junius T. and J. Oosterhaven (2003), The solution of
#'     updating or regionalizing a matrix with both positive and
#'     negative entries, Economic Systems Research, 15, pp. 87-96.
#'
#'     Lenzen M., R. Wood and B. Gallego (2007), Some comments on the
#'     GRAS method, Economic Systems Research, 19, pp. 461-465.
#'
#'     Temurshoev, U., R.E. Miller and M.C. Bouwmeester (2013), A note
#'     on the GRAS method, Economic Systems Research, 25, pp. 361-367.
#' @keywords gras, matrix updating
#' @export
#' @import data.table
doGRAS.long <- function(dt, rowcol, colcol,
                        rowsum, colsum, ...) {

  mat.wide <- dcast.data.table(dt, get(rowcol) ~ get(colcol), fill = 0)

  ## save the rowvalues for later on
  ## rowvalues <- mat.wide[, 1, with = FALSE] # can't get it to work with the name of the column...
  rowvalues <- mat.wide[, rowcol]

  ## extract the matrix
  mat <- as.matrix(mat.wide[, 2:ncol(mat.wide), with = FALSE])

  ## do the RASing
  res <- doGRAS(mat, rowsum, colsum, ...)

  ## feed the matrix back into a data.table
  res <- data.table(prod.na = rowvalues,
                    res)

  ## and convert it back to a long data.table
  dt <- melt.data.table(res, id.vars = "prod.na",
                        variable.name = "induse",
                        variable.factor = FALSE)

  dt
}
