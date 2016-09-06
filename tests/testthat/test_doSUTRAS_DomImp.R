
data(DK_SUT_2010)
data(DK_ProjData_2011)

dk.2011$M <- dk.2011$xbar[65]
dk.2011$xbar <-  dk.2011$xbar[-65]
dk.2011$c0 <- c(dk.2011$c0, 0, 0, 0)

sut <- doSUTRAS.DomImp(dk.sut, dk.2011, epsilon = 0.0001)

isConsistent(sut)
