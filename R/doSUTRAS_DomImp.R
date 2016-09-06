#' do a SUTRAS update 
#'
#' do a SUTRAS update when the use table is separated in domestic and
#' imported products. Row and column totals are met.
#' @param SUT a SUT object containing supply and use matrizes
#' @param ProjData Projection (forecast) data containing intermediate
#'     consumption and gross output
#' @param epsilon the convergence tolerance, default = 1e-6.
#' @param max.iterations the maximal number of iterations that would
#'     be carried out
#' @param verbose should diagnostic output of the iterations be
#'     displayed? Default = FALSE.
#' @import data.table
doSUTRAS.DomImp <- function(SUT, #V, Ud, Um, m0,
                            ProjData, #ubar, xbar, M, c0,
                            epsilon = 1e-6,
                            max.iterations = 10000, 
                            verbose = FALSE) {

    ##
    ## use own variables for better readability
    V <- SUT$V
    Ud <- SUT$Ud
    Um <- SUT$Um
    m0 <- SUT$m0

    ubar <- ProjData$ubar
    xbar <- ProjData$xbar
    M <- ProjData$M
    c0 <- ProjData$c0

    ## 
    ## separate positive from negative entries
    ## V (make matrix)
    positiveV <- V >= 0
    Pv0 <- V * positiveV
    Nv0 <- abs(V * !positiveV)

    ## U matrices, Use table
    positiveUd <- Ud >= 0
    Pd0 <- Ud * positiveUd
    Nd0 <- abs(Ud * !positiveUd)

    positiveUm <- Um >= 0
    Pm0 <- Um * positiveUm
    Nm0 <- abs(Um * !positiveUm)
    
    iterations <- 0L

    rm <- rep(1, dim(Pm0)[1])
    rd <- rep(1, dim(Pd0)[1])
    su <- rep(1, dim(Pd0)[2])
    rv <- rep(1, dim(Pv0)[2])
    r <- 1

    repeat {
        rdm1 <- rd
        rmm1 <- rm

        ## rd
        pd <- Pd0 %*% su + rowSums(Nv0 %*% dinv(rv))
        nd <- rowSums(Nd0 %*% dinv(su)) + Pv0 %*% rv
        rd <- 0.5 * dinv(pd) %*% (-c0 + sqrt(c0 * c0 + 4 * (pd * nd)))
        
        ## rm
        rm <- sqrt(dinv(Pm0 %*% su) %*% (rowSums(Nm0 %*% dinv(su)) + r * m0))

        ## rv
        Y <- rowSums(t(Pv0) %*% dinv(rd))
        X <- xbar + sqrt(xbar * xbar + 4 * (Y * (t(Nv0) %*% rd)))
        rv <- 0.5 * dinv(Y) %*% X

        ## su
        ps <- t(Pd0) %*% rd + t(Pm0) %*% rm
        ns <- rowSums(t(Nd0) %*% dinv(rd)) + rowSums(t(Nm0) %*% dinv(rm))
        X <- ubar + sqrt(ubar * ubar + 4 * ps * ns)
        su <- 0.5 * dinv(ps) %*% X

        ## r
        r <- M / sum(m0 * sinv(rm))

        ## test for convergence 
        rd.test <- all(abs(rdm1 - rd) < epsilon)
        rm.test <- all(abs(rmm1 - rm) < epsilon)

        iterations <- iterations + 1L

        if(verbose == TRUE) {
            domdiff <- abs(rdm1 - rd)
            impdiff <- abs(rmm1 - rm)
            cat("\n----- Iteration:", iterations, " ------")
            cat("\nSum of abs diff:\n\tDomestic:", sum(domdiff),
                "| Imports:", sum(impdiff))
            cat("\nNumber of diff > epsilon:\n\tDomestic:",
                sum(domdiff > epsilon),
                "| Imports:", sum(impdiff > epsilon), "\n")

            if(iterations > 1) {
                if(sum(domdiff > epsilon) <= 10 & sum(domdiff > epsilon) > 0) {
                    dom.culprits <- as.vector(domdiff > epsilon)
                    dom.index <- 1:length(dom.culprits)
                    domvalues <- data.table(index = dom.index[dom.culprits],
                                            old = rdm1[dom.culprits],
                                            new = rd[dom.culprits],
                                            diff = rdm1[dom.culprits] - rd[dom.culprits])
                    cat("\nDomestic table:\n")
                    print.data.frame(domvalues, digits = 3, row.names = FALSE)
                }
                if(sum(impdiff > epsilon) <= 10 & sum(impdiff > epsilon) > 0) {
                    imp.culprits <- as.vector(impdiff > epsilon)
                    imp.index <- 1:length(imp.culprits)
                    impvalues <- data.table(index = imp.index[imp.culprits],
                                            old = rdm1[imp.culprits],
                                            new = rd[imp.culprits],
                                            diff = rdm1[imp.culprits] - rd[imp.culprits])
                    cat("\nImport table:\n")
                    print.data.frame(impvalues, digits = 3, row.names = FALSE)
                }
            }
        }
        else if(iterations %% 50 == 0) {
            cat(".")
        }

        if((rd.test & rm.test) == TRUE) {
            message("\n ==> SUTRAS Algorithm finished with SUCCESS \n")
            break
        } else if((iterations >= max.iterations) == TRUE) {
            stop("\n ==> Algorithm did not converge! <==\n")
            break
        }
    }

    ## convert them to vectors, so that multiplications below are easier
    rm <- as.vector(rm)
    rd <- as.vector(rd)
    su <- as.vector(su)
    rv <- as.vector(rv)

    ## 
    ## calculate the resulting matrizes
    ##

    ## V
    V <- t(rv * t(Pv0 * sinv(rd)) - sinv(rv) * t(Nv0 * rd))
    V[, 66] <- V[, 66] - c0

    ## Ud and Um
    Ud <- rd * t(t(Pd0) * su) - sinv(rd) * t(t(Nd0) * sinv(su))
    Um <- rm * t(t(Pm0) * su) - sinv(rm) * t(t(Nm0) * sinv(su))

    ## m vector
    m <- r * (m0 * sinv(rm))

    return(structure(list(country = SUT$country, year = SUT$year,
                          V = V, Ud = Ud, Um = Um, m0 = m),
                     class = "SUT")) 
}
