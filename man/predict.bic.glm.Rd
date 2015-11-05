\name{predict.bic.glm}
\alias{predict.bic.glm}
\title{Predict function for Bayesian Model Averaging for generalized
  linear models.}
\description{
Bayesian Model Averaging (BMA) accounts for the model uncertainty
inherent in the variable selection problem by averaging over the best
models in the model class according to approximate posterior model
probability. This function predicts the response resulting from a BMA
generalized linear model from given data.
}
\usage{
\method{predict}{bic.glm}( object, newdata, \dots)
}
\arguments{
  \item{object}{a fitted object inheriting from class \code{bic.glm}.}
  \item{newdata}{a data frame containing observations on variables from 
                 which the predictor variables are to be selected or 
                 constructed from a formula.}
  \item{\dots}{ignored (for compatibility with generic function).}
}
\value{
  The predicted values from the BMA model for each observation in newdata.
}
\seealso{ \code{\link{bic.glm}}  }
\examples{

\dontrun{
# Example 1 (Gaussian)

     library(MASS)
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
}

\dontrun{
# Example 2 (binomial)

     library(MASS)
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
}

\dontrun{
# Example 3 (Gaussian)

     library(MASS)
     data(anorexia)

     anorexia.formula <- formula(Postwt ~ Prewt+Treat+offset(Prewt))

     bic.glm.anorexiaF <- bic.glm( anorexia.formula, data=anorexia,
                       glm.family="gaussian", factor.type=FALSE)
     predict( bic.glm.anorexiaF, newdata=anorexia)

     bic.glm.anorexiaT <- bic.glm( anorexia.formula, data=anorexia,
                           glm.family="gaussian", factor.type=TRUE)
     predict( bic.glm.anorexiaT, newdata=anorexia)
}

\dontrun{
# Example 4 (Gamma)

     library(survival)
     data(veteran)

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
}

# Example 5 (poisson - Yates teeth data)

     x <- rbind.data.frame(c(0, 0, 0),
                           c(0, 1, 0),
                           c(1, 0, 0),
                           c(1, 1, 1))

     y <- c(4, 16, 1, 21)
     n <- c(1,1,1,1)

     bic.glm.yatesF <- bic.glm( x, y, glm.family=poisson(),
                               weights=n, factor.type=FALSE)

     predict( bic.glm.yatesF, x)

\dontrun{
# Example 6 (binomial - Venables and Ripley)

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
}
}
\keyword{regression}
\keyword{models}
