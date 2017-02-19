test_fisherMethod <- function() {
    checkEquals(fisherMethod(rep(0,10)),0)
    checkEquals(fisherMethod(rep(1,10)),1)
    checkEquals(fisherMethod(c(0,1)),0)
}

test_stoufferMethod <- function() {
    checkEquals(stoufferMethod(rep(0,10)),0)
    checkEquals(stoufferMethod(rep(1,10)),1)
}

test_addCLT <- function() {
    checkEquals(addCLT(rep(0,10)),0)
    checkEquals(addCLT(rep(1,10)),1)
    checkEquals(addCLT(c(rep(0,10),rep(1,10))),0.5)
}
