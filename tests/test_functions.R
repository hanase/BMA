library(survival)
library(MASS)

start.test <- function(name) cat('\n<=== Starting test of', name, '====\n')
test.ok <- function(name) cat('\n==== Test of', name, 'OK.===>\n')

# These tests are mostly taken from examples that are set to dontrun.
# When making changes here, make sure these changes are also made in the examples.
test.bic.glm.gamma <- function(){
    test.name <- 'bic.glm Gamma family'
    start.test(test.name)
    data(cancer)
    surv.t<- veteran$time
    x<- veteran[,-c(3,4)]
    x$celltype<- factor(as.character(x$celltype))
    sel<- veteran$status == 0
    x<- x[!sel,]
    surv.t<- surv.t[!sel]
    
    glm.out.va <- bic.glm(x, y=surv.t, glm.family=Gamma(link="inverse"),
                          factor.type=FALSE)
    summary(glm.out.va)
    imageplot.bma(glm.out.va)
    plot(glm.out.va)
    test.ok(test.name)
}

test.bic.glm.poisson <- function(){
    test.name <- 'bic.glm Poisson family'
    x<- rbind(
        c(0, 0, 0),
        c(0, 1, 0),
        c(1, 0, 0),
        c(1, 1, 1))
    
    y<-c(4, 16, 1, 21)
    n<-c(1,1,1,1)
    
    models<- rbind(
        c(1, 1, 0),
        c(1, 1, 1))
    
    glm.out.yates <- bic.glm( x, y, n, glm.family = poisson(), 
                              factor.type=FALSE) 
    summary(glm.out.yates)
    test.ok(test.name)
}

test.bic.glm.gamma.crime <- function(){
    test.name <- 'bic.glm Gamma family with crime data'
    start.test(test.name)
    data(UScrime)
    f <- formula(log(y) ~  log(M)+So+log(Ed)+log(Po1)+log(Po2)+log(LF)+
                     log(M.F)+ log(Pop)+log(NW)+log(U1)+log(U2)+
                     log(GDP)+log(Ineq)+log(Prob)+log(Time))
    glm.out.crime <- bic.glm(f, data = UScrime, glm.family = gaussian()) 
    summary(glm.out.crime)
    test.ok(test.name)
}

test.bic.glm.logistic <- function(){
    test.name <- 'bic.glm logistic regression'
    start.test(test.name)
    
    data(birthwt)
    y<- birthwt$lo
    x<- data.frame(birthwt[,-1])
    x$race<- as.factor(x$race)
    x$ht<- (x$ht>=1)+0
    x<- x[,-9]
    x$smoke <- as.factor(x$smoke)
    x$ptl<- as.factor(x$ptl)
    x$ht <- as.factor(x$ht)
    x$ui <- as.factor(x$ui)
    
    glm.out.FT <- bic.glm(x, y, strict = FALSE, OR = 20, 
                          glm.family="binomial", factor.type=TRUE)
    summary(glm.out.FT)
    imageplot.bma(glm.out.FT)
    
    glm.out.FF <- bic.glm(x, y, strict = FALSE, OR = 20, 
                          glm.family="binomial", factor.type=FALSE)
    summary(glm.out.FF)
    imageplot.bma(glm.out.FF)
    
    glm.out.TT <- bic.glm(x, y, strict = TRUE, OR = 20, 
                          glm.family="binomial", factor.type=TRUE)
    summary(glm.out.TT)
    imageplot.bma(glm.out.TT)
    
    glm.out.TF <- bic.glm(x, y, strict = TRUE, OR = 20, 
                          glm.family="binomial", factor.type=FALSE)
    summary(glm.out.TF)
    imageplot.bma(glm.out.TF)
    test.ok(test.name)
}

test.bic.surv.veteran <- function(){
    test.name <- 'bic.surv with veteran data'
    start.test(test.name)

    data(cancer)

    test.bic.surv<- bic.surv(Surv(time,status) ~ ., data = veteran, 
                         factor.type = TRUE)
    summary(test.bic.surv, conditional=FALSE, digits=2)
    plot(test.bic.surv)

    imageplot.bma(test.bic.surv)
    test.ok(test.name)
}

test.bic.surv.pbc <- function(){
    test.name <- 'bic.surv with pbc data'
    start.test(test.name)
    
    data(cancer)
    x<- pbc[1:312,]
    surv.t<- x$time
    cens<- as.numeric((x$status == 2))

    x<- x[,c("age", "albumin", "alk.phos", "ascites", "bili", "edema", 
         "hepato", "platelet", "protime", "sex", "ast", "spiders", 
         "stage", "trt", "copper")]

    x$bili<- log(x$bili)
    x$alb<- log(x$alb)
    x$protime<- log(x$protime)
    x$copper<- log(x$copper)
    x$ast<- log(x$ast)
    
    test.bic.surv<- bic.surv(x, surv.t, cens, 
                             factor.type=FALSE, strict=FALSE)
    summary(test.bic.surv)
    test.ok(test.name)
}

test.bicreg <- function(){
    test.name <- 'bicreg'
    start.test(test.name)
    
    data(UScrime)
    x<- UScrime[,-16]
    y<- log(UScrime[,16])
    x[,-2]<- log(x[,-2])
    lma<- bicreg(x, y, strict = FALSE, OR = 20) 
    summary(lma)
    plot(lma)
    
    imageplot.bma(lma)
    test.ok(test.name)
}

test.glib.vaso <- function(){
    test.name <- 'glib with vaso data'
    start.test(test.name)
    ### Finney data
    data(vaso)
    x<- vaso[,1:2]
    y<- vaso[,3]
    n<- rep(1,times=length(y))
    
    finney.models<- rbind(
        c(1, 0),
        c(0, 1),
        c(1, 1))
    
    finney.glib <- glib (x,y,n, error="binomial", link="logit", 
                         models=finney.models, glimvar=TRUE, 
                         output.priorvar=TRUE, output.postvar=TRUE)
    summary(finney.glib)
    
    finney.bic.glm<- as.bic.glm(finney.glib)
    plot(finney.bic.glm,mfrow=c(2,1))
    test.ok(test.name)
}

test.glib.yates <- function(){
    test.name <- 'glib with Yates data'
    start.test(test.name)
    x<- rbind(
        c(0, 0, 0),
        c(0, 1, 0),
        c(1, 0, 0),
        c(1, 1, 1))
    
    y<-c(4, 16, 1, 21)
    n<-c(1,1,1,1)
    
    models<- rbind(
        c(1, 1, 0),
        c(1, 1, 1))
    
    glib.yates <- glib ( x, y, n, models=models, glimvar=TRUE,
                         output.priorvar=TRUE, output.postvar=TRUE) 
    summary(glib.yates)
    test.ok(test.name)
}

test.glib.logreg <- function(){
    test.name <- 'glib for logistic regression'
    start.test(test.name)
    data(birthwt)
    y<- birthwt$lo
    x<- data.frame(birthwt[,-1])
    x$race<- as.factor(x$race)
    x$ht<- (x$ht>=1)+0
    x<- x[,-9]
    x$smoke <- as.factor(x$smoke)
    x$ptl<- as.factor(x$ptl)
    x$ht <- as.factor(x$ht)
    x$ui <- as.factor(x$ui)
    
    glib.birthwt<- glib(x,y, error="binomial", link = "logit")
    summary(glib.birthwt)
    
    glm.birthwt<- as.bic.glm(glib.birthwt)
    
    imageplot.bma(glm.birthwt)
    plot(glm.birthwt)
    test.ok(test.name)
}

test.iBMA.glm <- function(){
    test.name <- 'iBMA for glm'
    start.test(test.name)
    data(birthwt)
    y<- birthwt$lo
    x<- data.frame(birthwt[,-1])
    x$race<- as.factor(x$race)
    x$ht<- (x$ht>=1)+0
    x<- x[,-9]
    x$smoke <- as.factor(x$smoke)
    x$ptl<- as.factor(x$ptl)
    x$ht <- as.factor(x$ht)
    x$ui <- as.factor(x$ui)
    
    ### add 41 columns of noise
    noise<- matrix(rnorm(41*nrow(x)), ncol=41)
    colnames(noise)<- paste('noise', 1:41, sep='')
    x<- cbind(x, noise)
    
    iBMA.glm.out<- iBMA.glm( x, y, glm.family="binomial", 
                             factor.type=FALSE, verbose = TRUE, 
                             thresProbne0 = 5 )
    summary(iBMA.glm.out)
    test.ok(test.name)
}

test.iBMA.surv <- function(){
    test.name <- 'iBMA for surv'
    start.test(test.name)
    
    data(cancer)
    
    surv.t<- veteran$time
    cens<- veteran$status
    veteran$time<- NULL
    veteran$status<- NULL
    lvet<- nrow(veteran)
    invlogit<- function(x) exp(x)/(1+exp(x))
    # generate random noise, 34 uniform variables 
    # and 10 factors each with 4 levels
    X <- data.frame(matrix(runif(lvet*34), ncol=34), 
                    matrix(letters[1:6][(rbinom(10*lvet, 3, .5))+1], 
                           ncol = 10))
    colnames(X) <- c(paste("u",1:34, sep=""),paste("C",1:10, sep=""))
    for(i in 35:44) X[,i] <- as.factor(X[,i])
    veteran_plus_noise<- cbind(veteran, X)
    
    
    test.iBMA.surv <- iBMA.surv(x = veteran_plus_noise, 
                                surv.t = surv.t, cens = cens, 
                                thresProbne0 = 5, maxNvar = 30, 
                                factor.type = TRUE, verbose = TRUE, 
                                nIter = 100)
    
    test.iBMA.surv
    summary(test.iBMA.surv)
    test.ok(test.name)
}

test.iBMA.bicreg <- function(){
    test.name <- 'iBMA for bicreg'
    start.test(test.name)
    data(UScrime)
    UScrime$M<- log(UScrime$M); UScrime$Ed<- log(UScrime$Ed); 
    UScrime$Po1<- log(UScrime$Po1); UScrime$Po2<- log(UScrime$Po2); 
    UScrime$LF<- log(UScrime$LF); UScrime$M.F<- log(UScrime$M.F)
    UScrime$Pop<- log(UScrime$Pop); UScrime$NW<- log(UScrime$NW); 
    UScrime$U1<- log(UScrime$U1); UScrime$U2<- log(UScrime$U2); 
    UScrime$GDP<- log(UScrime$GDP); UScrime$Ineq<- log(UScrime$Ineq)
    UScrime$Prob<- log(UScrime$Prob); UScrime$Time<- log(UScrime$Time) 
    noise<- matrix(rnorm(35*nrow(UScrime)), ncol=35)
    colnames(noise)<- paste('noise', 1:35, sep='')
    UScrime_plus_noise<- cbind(UScrime, noise)
    
    y<- UScrime_plus_noise$y
    UScrime_plus_noise$y <- NULL
    
    # run 2 iterations and examine results
    iBMA.bicreg.crime <- iBMA.bicreg( x = UScrime_plus_noise, 
                                      Y = y, thresProbne0 = 5, verbose = TRUE, maxNvar = 30, nIter = 2)
    summary(iBMA.bicreg.crime)
    orderplot(iBMA.bicreg.crime)
    
    # run from current state until completion
    iBMA.bicreg.crime <- iBMA.bicreg( iBMA.bicreg.crime, nIter = 200)
    summary(iBMA.bicreg.crime)
    orderplot(iBMA.bicreg.crime)
    
    test.ok(test.name)
}

test.iBMA.bicreg.large <- function(){
    test.name <- 'iBMA for bicreg with many predictors '
    start.test(test.name)
    set.seed(0)
    x <- matrix( rnorm(50*30), 50, 30)
    lp <- apply( x[,1:5], 1, sum)
    iBMA.bicreg.ex <- iBMA.bicreg( x = x,  Y = lp, thresProbne0 = 5, maxNvar = 20)
    
    explp <- exp(lp)
    prob <- explp/(1+explp)
    y=rbinom(n=length(prob),prob=prob,size=1)
    iBMA.glm.ex <- iBMA.glm( x = x, Y = y, glm.family = "binomial",
                             factor.type = FALSE, thresProbne0 = 5, maxNvar = 20)
    
    test.ok(test.name)
}

test.predict.gaussian1 <- function(){
    test.name <- 'predict example 1 (Gaussian)'
    start.test(test.name)
    data(UScrime)
    
    f <- formula(log(y) ~  log(M)+So+log(Ed)+log(Po1)+log(Po2)+
                     log(LF)+log(M.F)+log(Pop)+log(NW)+log(U1)+log(U2)+
                     log(GDP)+log(Ineq)+log(Prob)+log(Time))
    
    bic.glm.crimeT <- bic.glm(f, data = UScrime, 
                              glm.family = gaussian())
    predict(bic.glm.crimeT, newdata = UScrime)
    
    bic.glm.crimeF <- bic.glm(f, data = UScrime, 
                              glm.family = gaussian(),
                              factor.type = FALSE)
    predict(bic.glm.crimeF, newdata = UScrime)
    test.ok(test.name)
}

test.predict.binomial2 <- function(){
    test.name <- 'predict example 2 (binomial)'
    start.test(test.name)
    data(birthwt)
    
    y <- birthwt$lo
    x <- data.frame(birthwt[,-1])
    x$race <- as.factor(x$race)
    x$ht <- (x$ht>=1)+0
    x <- x[,-9]
    x$smoke <- as.factor(x$smoke)
    x$ptl <- as.factor(x$ptl)
    x$ht  <- as.factor(x$ht)
    
    x$ui <- as.factor(x$ui)
    
    bic.glm.bwT <- bic.glm(x, y, strict = FALSE, OR = 20,
                           glm.family="binomial",  
                           factor.type=TRUE)
    predict( bic.glm.bwT, newdata = x)
    
    bic.glm.bwF <- bic.glm(x, y, strict = FALSE, OR = 20,
                           glm.family="binomial",  
                           factor.type=FALSE)
    predict( bic.glm.bwF, newdata = x)
    test.ok(test.name)
}

test.predict.gaussian3 <- function(){
    test.name <- 'predict example 3 (Gaussian)'
    start.test(test.name)
    
    data(anorexia)
    
    anorexia.formula <- formula(Postwt ~ Prewt+Treat+offset(Prewt))
    
    bic.glm.anorexiaF <- bic.glm( anorexia.formula, data=anorexia,
                                  glm.family="gaussian", factor.type=FALSE)
    predict( bic.glm.anorexiaF, newdata=anorexia)
    
    bic.glm.anorexiaT <- bic.glm( anorexia.formula, data=anorexia,
                                  glm.family="gaussian", factor.type=TRUE)
    predict( bic.glm.anorexiaT, newdata=anorexia)
    test.ok(test.name)
}

test.predict.gamma <- function(){
    test.name <- 'predict example 4 (Gamma)'
    start.test(test.name)
    data(cancer)
    
    surv.t <- veteran$time
    x <- veteran[,-c(3,4)]
    x$celltype <- factor(as.character(x$celltype))
    sel<- veteran$status == 0
    x <- x[!sel,]
    surv.t <- surv.t[!sel]
    
    bic.glm.vaT <- bic.glm(x, y=surv.t, 
                           glm.family=Gamma(link="inverse"),
                           factor.type=TRUE)
    predict( bic.glm.vaT, x)
    
    bic.glm.vaF <- bic.glm(x, y=surv.t, 
                           glm.family=Gamma(link="inverse"),
                           factor.type=FALSE)
    predict( bic.glm.vaF, x)
    
    test.ok(test.name)
}

test.predict.poisson <- function(){
    test.name <- 'predict example 5 (Poisson)'
    start.test(test.name)
    
    x <- rbind.data.frame(c(0, 0, 0),
                          c(0, 1, 0),
                          c(1, 0, 0),
                          c(1, 1, 1))
    
    y <- c(4, 16, 1, 21)
    n <- c(1,1,1,1)
    
    bic.glm.yatesF <- bic.glm( x, y, glm.family=poisson(),
                               weights=n, factor.type=FALSE)
    
    predict( bic.glm.yatesF, x)
    test.ok(test.name)
}

test.predict.binomial6 <- function(){
    test.name <- 'predict example 6 (binomial)'
    start.test(test.name)
    ldose <- rep(0:5, 2)
    numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
    sex <- factor(rep(c("M", "F"), c(6, 6)))
    SF <- cbind(numdead, numalive=20-numdead) 
    
    budworm <- cbind.data.frame(ldose = ldose, numdead = numdead,
                                sex = sex, SF = SF)
    budworm.formula <- formula(SF ~ sex*ldose)
    
    bic.glm.budwormF <- bic.glm( budworm.formula, data=budworm,
                                 glm.family="binomial", factor.type=FALSE)
    predict(bic.glm.budwormF, newdata=budworm)
    
    bic.glm.budwormT <- bic.glm( budworm.formula, data=budworm,
                                 glm.family="binomial", factor.type=TRUE)
    predict(bic.glm.budwormT, newdata=budworm)
    
    test.ok(test.name)
}

test.predict.bicreg <- function(){
    test.name <- 'predict bicreg'
    start.test(test.name)
    
    # Example 1
    data(UScrime)
    
    x <- UScrime[,-16]
    y <- log(UScrime[,16])
    x[,-2]<- log(x[,-2])
    
    crimeBMA <- bicreg(x, y, strict = FALSE, OR = 20)
    predict( crimeBMA, x)
    
    # Example 2 (Venables and Ripley)
    
    npkBMA <- bicreg( x = npk[, c("block","N","K")], y=npk$yield)
    predict( npkBMA, newdata = npk)
    
    # Example 3 (Venables and Ripley)
    
    gasPRbma <- bicreg( x = whiteside[,c("Insul", "Temp")], 
                        y = whiteside$Gas)
    predict( gasPRbma, newdata = whiteside)
    
    test.ok(test.name)
}

test.mc3.reg <- function(){
    test.name <- 'MC3.REG'
    start.test(test.name)
    
    data(race)
    b<- out.ltsreg(race[,-1], race[,1], 2)
    races.run1<-MC3.REG(race[,1], race[,-1], num.its=20000, c(FALSE,TRUE), 
                        rep(TRUE,length(b)), b, PI = .1, K = 7, nu = .2, 
                        lambda = .1684, phi = 9.2)
    races.run1
    summary(races.run1)
    test.ok(test.name)
}
