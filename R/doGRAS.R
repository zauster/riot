#' Do the GRAS algorithm for updating matrizes
#'
#' This function calculates an updated matrix X, based on a matrix A,
#' which meets the given row and columns totals. Note: A can contain
#' negative elements.
#' @param A the "base" matrix
#' @param u vector with the row totals
#' @param v vector with the column totals
#' @param epsilon the error tolerance level, default is 1e-10.
#' @param max.iter maximum number of iterations, default is 10000.
#' @param verbose should some information of the iterations be
#'     displayed? Default is FALSE
#' @return the updated matrix X
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
#' @examples
#' ## example from the papers
#' A <- matrix(c(7, 3, 5, -3, 2, 9, 8, 0, -2, 0, 2, 0),
#'             ncol = 4, nrow = 3, byrow = TRUE)
#' u <- c(15, 26, -1)
#' v <- c(9, 16, 17, -2)
#' doGRAS(A, u, v)
#' @export
doGRAS <- function(A, u, v,
                   epsilon = 1e-10,
                   max.iter = 10000,
                   verbose = FALSE) {

    ## TODO:
    ## check for correct dimensions
    ## conversion to vectors/matrices
    ## what about NAs?

    ## get col/row numbers
    m <- nrow(A)
    n <- ncol(A)

    if(n != length(v)) {
        stop(paste0("Unequal number of columns! ncol(A) = ", n, " and length(v) = ", length(v),
                    " should be equal!"))
    }
    if(m != length(u)) {
        stop(paste0("Unequal number of rows! nrow(A) = ", m, " and length(u) = ", length(u),
                    " should be equal!"))
    }

    ## separate A into positive (P) and negative (N) matrizes
    P <- A
    N <- abs(A)
    N[A >= 0] <- 0 ## set positive numbers to 0
    P[A < 0] <- 0  ## set negative numbers to 0

    ## initialize r
    r <- rep(1, m)
    s <- rep(1, n)
    error <- 1
    iter <- 1

    while((error > epsilon) & (iter <= max.iter)) {

        s.old <- s

        ##
        ## calculate multiplicator s
        ## 
        pj.r <- t(P) %*% r
        nj.r <- t(N) %*% sinv(r)

        s <- sinv(2 * pj.r) * (v + sqrt(v^2 + 4 * pj.r * nj.r))
        s.alt <- -sinv(v) * nj.r
        s[s == 0] <- s.alt[s == 0]

        ##
        ## calculate multiplicator r
        ##
        pi.s <- P %*% s
        ni.s <- N %*% sinv(s)

        r <- sinv(2 * pi.s) * (u + sqrt(u^2 + 4 * pi.s * ni.s))
        r.alt <- -sinv(u) * ni.s
        r[r == 0] <- r.alt[r == 0]

        ##
        ## calculate the remaining error
        diff <- abs(s.old - s)
        error <- max(diff)
        error.pos <- which.max(diff)

        if(verbose) {
            message("iter: ", iter, "  -> error: ", error, " at ", error.pos)
            }
        
        iter <- iter + 1
    }

    if(iter >= max.iter) {
        warning(paste0(" => Stopped after ", max.iter, " iterations. Optimal solution is NOT guaranteed!"))
    }

    if(verbose) {
        message(" => Iterations needed: ", iter - 1)
        message(" => Resulting error: ", error, " at ", error.pos)
        }

    r <- as.vector(r)
    s <- as.vector(s)

    ## X.1 <- diag(r) %*% P %*% diag(s) - diag(sinv(r)) %*% N %*% diag(sinv(s))
    X <- sweep(r * P, 2, s, "*") - sweep(sinv(r) * N, 2, sinv(s), "*")
    
    return(X)
}

   
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
doGRAS.long <- function(dt, rowcol, colcol,
                        rowsum, colsum, ...) {

    mat.wide <- dcast.data.table(dt, get(rowcol) ~ get(colcol), fill = 0)

    ## save the rowvalues for later on
    rowvalues <- mat.wide[, 1] # can't get it to work with the name of the column...

    ## extract the matrix
    ## mat <- as.matrix(mat.wide[, columnstoKeep, with = FALSE])
    mat <- as.matrix(mat.wide[, -1])
    
    ## do the RASing
    res <- doGRAS(mat, rowsum, colsum, ...)

    ## feed the matrix back into a data.table
    res <- data.table(prod.na = rowvalues, 
                      res)
    
    ## and convert it back to a long data.table
    dt <- melt(res, id.vars = "prod.na",
               variable.name = "induse", variable.factor = FALSE)
    dt
}
