
data(TestSUTs_AUT2010)
data(TestNA_AUT2011)
c <- rep(0, length(m0))

res <- doSUTRAS(V = V, Ud = Ud, Um = Um, m0 = m0, u_bar = ubar,
                x_bar = xbar, M = M, c = c, epsilon = 1e-10,
                maxiter = 100, verbose = FALSE)

names(xbar) <- NULL
names(ubar) <- NULL
names(M) <- NULL
expect_equal(colSums(res$V), xbar, tolerance = 0.001)
expect_equal(colSums(res$Ud) + colSums(res$Um), ubar, tolerance = 0.001)
expect_equal(sum(res$m), M, tolerance = 0.001)
