#' Simple vector inverse function
#'
#' Calculate the inverse of a vector (1 / x) and set
#' non-finite values back to 1
#' @param x the vector
sinv <- function(x) {
    res <- 1 / as.vector(x)
    res[!is.finite(res)] <- 1
    return(res)
}

#' Vector inverse and diagonal function
#'
#' Calculate the inverse of a vector (1 / x) and set
#' non-finite values back to 1 and diagonalize it
#' @param x the vector
dinv <- function(x) {
    res <- 1 / as.vector(x)
    res[!is.finite(res)] <- 1
    return(diag(res))
}
