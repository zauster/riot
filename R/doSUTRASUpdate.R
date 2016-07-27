## 
## function to generate updated SUTs using the SUTRAS algorithm
## 
doSUTRASUpdate <- function(lst, epsilon = NULL, ## 0.0000005,
                           max.iterations = 50000,
                           verbose = FALSE,
                           TLSindex = 65, TTMindex = 66) {

    cat("\n>>---------------------------<<\n")
    cat("   [", lst[["country"]], "] [",
        lst[["SUTyear"]], " -> ", lst[["Projyear"]], "]\n")

    if(is.null(epsilon)) {
        epsilon <- lst[["epsilon"]]
    }
    
    ## verbose <- TRUE
    ## epsilon <- 0.00005
    ## max.iterations <- 10000L

    start <- Sys.time()

    tryResult <- try({
        c0 <- lst[["SUTs"]][["c0"]]
        m0 <- lst[["SUTs"]][["m0"]]
        Pv0 <- lst[["SUTs"]][["Pv0"]]
        Nv0 <- lst[["SUTs"]][["Nv0"]]
        Pd0 <- lst[["SUTs"]][["Pd0"]]
        Nd0 <- lst[["SUTs"]][["Nd0"]]
        Pm0 <- lst[["SUTs"]][["Pm0"]]
        Nm0 <- lst[["SUTs"]][["Nm0"]]

        M <- lst[["ProjData"]][["M"]]
        ubar <- lst[["ProjData"]][["ubar"]]
        xbar <- lst[["ProjData"]][["xbar"]]


        ## doSUTRAS begin here ====>
        
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
            r <- M / sum(m0 * inv(rm))

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
        V <- t(rv * t(Pv0 * inv(rd)) - inv(rv) * t(Nv0 * rd))
        V[, 66] <- V[, 66] - c0

        ## Ud and Um
        Ud <- rd * t(t(Pd0) * su) - inv(rd) * t(t(Nd0) * inv(su))
        Um <- rm * t(t(Pm0) * su) - inv(rm) * t(t(Nm0) * inv(su))

        ## doSUTRAS until here <====

        ## and Udb and Umb
        TLSvector <- V[, TLSindex]
        TTMvector <- V[, TTMindex]

        U <- Ud + Um

        TLSmatrix <- TLSvector * dinv(rowSums(U)) %*% U
        TTMmatrix <- TTMvector * dinv(rowSums(U)) %*% U

        Ushare <- abs(Ud) / (abs(Ud) + abs(Um))
        Ushare[!is.finite(Ushare)] <- 0

        ##
        ## replace following with doGRAS()
        ## 
        TLSmatrix.domestic <- TLSmatrix * Ushare
        TTMmatrix.domestic <- TTMmatrix * Ushare
        TLSmatrix.import <- TLSmatrix * (1 - Ushare)
        TTMmatrix.import <- TTMmatrix * (1 - Ushare)

        Udb <- Ud - TTMmatrix.domestic - TLSmatrix.domestic
        Umb <- Um - TTMmatrix.import - TLSmatrix.import

        ## m vector
        m <- r * (m0 * inv(rm))

        updatedSUTs <- list(V = V,
                            Ud = Ud, Um = Um,
                            Udb = Udb, Umb = Umb)

        
        ##
        ## Diagnostic output
        ##

        cat("Iterations until convergence:\t", iterations)

        ## test if all condition are met
        ## V Colsums == xbar
        Vdiff <- sum(abs(colSums(V[, 1:65]) - xbar[1:65]))
        cat("\nSum of abs diff (V - xbar):\t",
            round(Vdiff, digits = 5))
        if(Vdiff > 0.05) {
            Vcolsums <- colSums(V[, 1:65])
            cat("\n")
            print(cbind(Vcolsums = Vcolsums, xbar = xbar[1:65], diff = Vcolsums - xbar[1:65]))
            stop("\nSuppy table: Large differences in the column sums!!")
        }

        ## U Colsums == ubar
        Udiff <- sum(abs(colSums(Ud) + colSums(Um) - ubar))
        cat("\nSum of abs diff (U - ubar):\t",
            round(Udiff, digits = 5))
        if(Udiff > 0.05) {
            Ucolsums <- colSums(Ud) + colSums(Um)
            cat("\n")
            print(cbind(Ucolsums = Ucolsums, ubar = ubar, diff = Ucolsums - ubar))
            stop("\nUse table: Large differences in the column sums!!")
        }


        ## m sum = M
        mdiff <- abs(sum(m) - M)
        cat("\nSum of abs diff (M - m):\t", mdiff)
        if(mdiff > 0.05) {
            stop("\nImport vector: Large differences in the sum!!")
        }


        ##
        ## Purchaser prices
        ## 
        ## V Rowsums = U Rowsums
        Vsum <- rowSums(V) + m
        Usum <- rowSums(Ud) + rowSums(Um)

        Productdiff <- Vsum - Usum
        Productmax <- apply(cbind(Vsum, Usum), 1, max)
        relProductdiff <- Productdiff / Productmax
        relProductdiffmean <- mean(relProductdiff[is.finite(relProductdiff)], na.rm = TRUE)
        Productdiffsum <- abs(sum(Vsum) - sum(Usum))

        cat("\nPurchaser prices: ")
        cat("\nSum of abs diff product sums:\t", Productdiffsum)
        cat("\nSum of rel diff product mean:\t", relProductdiffmean, "\n")

        if(abs(relProductdiffmean) > 0.016) {
            stop("\n\n ==> Large differences in the total (purchaser prices) product sums! <== ")
        }

        ##
        ## Basic prices
        ## 
        ## V Rowsums = U Rowsums
        Vbsum <- rowSums(V[, 1:64]) + m
        Ubsum <- rowSums(Udb) + rowSums(Umb)

        Productdiff <- Vbsum - Ubsum
        Productmax <- apply(cbind(Vbsum, Ubsum), 1, max)
        relProductdiff <- Productdiff / Productmax
        relProductdiffmean <- mean(relProductdiff[is.finite(relProductdiff)], na.rm = TRUE)
        Productdiffsum <- abs(sum(Vbsum) - sum(Ubsum))

        cat("\nBasic prices: ")
        cat("\nSum of abs diff product sums:\t", Productdiffsum)
        cat("\nSum of rel diff product mean:\t", relProductdiffmean, "\n")

        if(abs(relProductdiffmean) > 0.016) {
            stop("\n\n ==> Large differences in the total (basic prices) product sums! <== ")
        }
        
        ## save some parameters
        parameters <- list(iterations.taken = iterations,
                           epsilon = epsilon,
                           duration = Sys.time() - start,
                           Vdiff = Vdiff,
                           Udiff = Udiff,
                           Proddiff = Productdiff)

        res <- list(updatedSUTs = updatedSUTs,
                    parameters = parameters)
        return(res)
    }, silent = TRUE)

    if(inherits(tryResult, "try-error")) {
        return(NULL)
    }
}
