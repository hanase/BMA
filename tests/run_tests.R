library(BMA)
source('test_functions.R')

cran <- TRUE

## disable the following tests when submitting to CRAN
## to speed-up the checking procedure
if(!cran) {
    test.bic.glm.gamma()
    test.bic.glm.poisson()
    test.bic.glm.gamma.crime()
    test.bic.glm.logistic()
    test.bic.surv.veteran()
    test.bic.surv.pbc()
    test.bicreg()
    test.glib.vaso()
    test.glib.yates()
    test.glib.logreg()
    test.iBMA.glm()
    test.iBMA.surv()
    test.iBMA.bicreg()
    test.iBMA.bicreg.large()
    test.predict.gaussian1()
    test.predict.binomial2()
    test.predict.gaussian3()
    test.predict.gamma()
    test.predict.poisson()
    #test.predict.binomial6() # need to fix
    test.predict.bicreg()
    test.mc3.reg()
}

