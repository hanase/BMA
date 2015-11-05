\name{iBMA}
\alias{iBMA}
\alias{iBMA.glm}
\alias{iBMA.bicreg}
\alias{iBMA.surv}
\alias{iBMA.glm.data.frame}
\alias{iBMA.glm.matrix}
\alias{iBMA.glm.iBMA.intermediate.glm}
\alias{iBMA.bicreg.data.frame}
\alias{iBMA.bicreg.matrix}
\alias{iBMA.bicreg.iBMA.intermediate.bicreg}
\alias{iBMA.surv.data.frame}
\alias{iBMA.surv.matrix}
\alias{iBMA.surv.iBMA.intermediate.surv}
\title{ Iterated Bayesian Model Averaging variable selection for generalized linear models, linear models or survival models.}
\description{
  This function implements the iterated Bayesian Model Averaging method for variable selection. 
  This method works by making repeated calls to a Bayesian model averaging procedure, iterating through the variables in a fixed order. After each call to the Bayesian model averaging procedure only those variables 
  which have posterior probability greater than a specified threshold are retained, those variables whose 
  posterior probabilities do not meet the threshold are replaced with the next set of variables. The order in which the variables are to be considered
  is usually determined on the basis of the some measure of goodness of fit calculated univariately for each variable.
}

\usage{
iBMA.glm(x, ...)
iBMA.bicreg(x, ...)
iBMA.surv(x, ...)

\method{iBMA.glm}{matrix}(x, Y, wt = rep(1, nrow(X)), 
       thresProbne0 = 5, glm.family, maxNvar = 30, 
       nIter = 100, verbose = FALSE, sorted = FALSE, 
       factor.type = TRUE, ...)

\method{iBMA.glm}{matrix}(x, Y, wt = rep(1, nrow(X)), 
       thresProbne0 = 5, glm.family, maxNvar = 30, 
       nIter = 100, verbose = FALSE, sorted = FALSE, 
       factor.type = TRUE, ...)

\method{iBMA.glm}{iBMA.intermediate.glm}(x, nIter = NULL, 
        verbose = NULL, ...) 

\method{iBMA.bicreg}{matrix}(x, Y, wt = rep(1, nrow(X)), 
        thresProbne0 = 5, maxNvar = 30, nIter = 100, 
        verbose = FALSE, sorted = FALSE, ...)

\method{iBMA.bicreg}{data.frame}(x, Y, wt = rep(1, nrow(X)), 
        thresProbne0 = 5, maxNvar = 30, nIter = 100, 
        verbose = FALSE, sorted = FALSE, ...)

\method{iBMA.bicreg}{iBMA.intermediate.bicreg}(x, 
        nIter = NULL, verbose = NULL, ...) 

\method{iBMA.surv}{matrix}(x, surv.t, cens, 
        wt = rep(1, nrow(X)), thresProbne0 = 5, 
        maxNvar = 30, nIter = 100, verbose = FALSE, 
        sorted = FALSE, factor.type = TRUE, ...)

\method{iBMA.surv}{data.frame}(x, surv.t, cens, 
        wt = rep(1, nrow(X)), thresProbne0 = 5, 
        maxNvar = 30, nIter = 100, verbose = FALSE, 
        sorted = FALSE, factor.type = TRUE, ...)

\method{iBMA.surv}{iBMA.intermediate.surv}(x, nIter = NULL,verbose = NULL, ...) 
}
\arguments{
  \item{x}{a matrix or data.frame of independent variables,
    or else an object of class \code{iBMA.glm.intermediate},
    \code{iBMA.bicreg.intermediate} or \code{iBMA.surv.intermediate}
     that contains the current state of an incomplete selection.}
  \item{Y}{a vector of values for the dependent variable.}
  \item{surv.t}{a vector of survival times.}
  \item{cens}{a vector of indicators of censoring (0=censored 1=uncensored)}
  \item{wt}{an optional vector of weights to be used.}
  \item{thresProbne0}{ a number giving the probability threshold for including variables as a percent.}
  \item{glm.family}{glm family.}
  \item{maxNvar}{ a number giving the maximum number of variables to be
considered in a model. }  
\item{nIter}{ a number giving the maximum number of iterations that should be run. }
  \item{verbose}{ a logical value specifying if verbose output should be produced or not }
  \item{sorted}{ a logical value specifying if the variables have been sorted or not. If \code{FALSE} then \code{iBMA.glm} will sort the variables prior to running any iterations.}
  \item{factor.type}{  a logical value specifying how variables of class "factor"
          are handled.  A factor variable with d levels is turned into
          (d-1) dummy variables using a treatment contrast.   If
          'factor.type = TRUE', models will contain either all or none
          of these dummy variables.   If 'factor.type = FALSE', models
          are free  to select the dummy variables independently.   In
          this case, factor.prior.adjust determines the prior on these
          variables.}
  \item{\dots}{ other parameters to be passed to \code{bic.glm}, \code{bicreg} or \code{bic.surv} }
}
\details{
These methods can be run in a 'batch' mode by setting \code{nIter} to be larger than the number of variables. 
Alternatively, if \code{nIter} is set to be small, the procedure may return before all of the variables have been examined. 
In this case the returned result of the call will be of class 'iBMA.X.intermediate', and if iBMA.X is called with this result as the input, \code{nIter} more iterations will be run.

If on any iteration there are no variables that have posterior probability less than the threshold, the variable with the lowest posterior probability is dropped.
}
\value{
  An object of either type iBMA.X, or of type iBMA.X.intermediate, where 'X' is either 'glm', 'bicreg' or 'surv'. 
  Objects of type 'iBMA.X.intermediate' consist of a list with components for each parameter passed into iBMA.X as well as the following components:
  \item{sortedX }{a matrix or data.frame containing the sorted variables.}
  \item{call }{the matched call.}
  \item{initial.order }{the inital ordering of the variables.}
  \item{nVar }{the number of variables.}
  \item{currentSet}{a vector specifying the set of variables currently selected.} 
  \item{nextVar }{the next variable to be examined}
  \item{current.probne0}{the posterior probabilities for inclusion for each of the variables in the current set of variables.} 
  \item{maxProbne0 }{the maximum posterior probability calculated for each variable}
  \item{nTimes }{the number of times each variable has been included in the set of selected variables}
  \item{currIter }{the current iteration number}
  \item{new.vars }{the set of variables that will be added to the current set during the next iteration}
  \item{first.in.model}{a vector of numbers giving the iteration number that each variable was first examined in. A value of NA indicates that a variable has not yet been examined.} 
  \item{iter.dropped }{a vector giving the iteration number in which each variable was dropped from the current set. A value of NA indicates that a variable has not yet been dropeed.}
  
  Objects of the type iBMA.glm contain in addition to all of these elements the following components:
    \item{nIterations}{the total number of iterations that were run}
    \item{selected}{the set of variables that were selected (in terms of the initial ordering of the variables)}
    \item{bma}{an object of type 'bic.X' containing the results of the Bayesian model averaging run on the selected set of variables.}

}
\references{Yeung, K.Y., Bumgarner, R.E. and Raftery, A.E. (2005). ` Bayesian Model Averaging: Development of an improved multi-class, gene selection and classification tool for microarray data.' Bioinformatics, 21(10), 2394-2402 }

\author{
Ka Yee Yeung, \email{kayee@AT@u.washington.edu},
Adrian Raftery \email{raftery@AT@stat.washington.edu},
 Ian Painter \email{ian.painter@AT@gmail.com}
}

\note{ 
 The parameters \code{verbose} and \code{nIter} can be changed between sets of iterations. 

  The parameter \code{sorted} specifies if the variables should be sorted prior to iteration, if \code{sorted} is set to \code{FALSE} then the variables are sorted according to the decreasing single variable model R2 values for iBMA.bicreg or the single variable model increasing Chi-sq P-values for iBMA.glm and iBMA.surv.
  Subsequent reference to variables is in terms of this ordered set of variables.
  
  It is possible to obtain degenerate results when using a large number of predictor variables in linear regression. This problem is much less common with logistic regression and survival analysis.
}

\seealso{ 
\code{\link{bic.glm}}, 
\code{\link{bicreg}}, 
\code{\link{bic.surv}},
\code{\link{summary.iBMA.bicreg}}, 
\code{\link{print.iBMA.bicreg}}, 
\code{\link{orderplot.iBMA.bicreg}}
}
 
\examples{

\dontrun{
############ iBMA.glm
library("MASS")
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
}

\dontrun{
################## iBMA.surv
library(survival)
data(veteran)

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
veteran_plus_noise<- cbind(veteran, X)


test.iBMA.surv <- iBMA.surv(x = veteran_plus_noise, 
                            surv.t = surv.t, cens = cens, 
                            thresProbne0 = 5, maxNvar = 30, 
                            factor.type = TRUE, verbose = TRUE, 
                            nIter = 100)

test.iBMA.surv
summary(test.iBMA.surv)
}

\dontrun{
############ iBMA.bicreg ... degenerate example
library(MASS)
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
}

\dontrun{
# run from current state until completion
iBMA.bicreg.crime <- iBMA.bicreg( iBMA.bicreg.crime, nIter = 200)
summary(iBMA.bicreg.crime)
orderplot(iBMA.bicreg.crime)
}

set.seed(0)
x <- matrix( rnorm(50*30), 50, 30)
lp <- apply( x[,1:5], 1, sum)
iBMA.bicreg.ex <- iBMA.bicreg( x = x,  Y = lp, thresProbne0 = 5, maxNvar = 20)

explp <- exp(lp)
prob <- explp/(1+explp)
y=rbinom(n=length(prob),prob=prob,size=1)
iBMA.glm.ex <- iBMA.glm( x = x, Y = y, glm.family = "binomial",
                         factor.type = FALSE, thresProbne0 = 5, maxNvar = 20)

cat("\n\n CAUTION: iBMA.bicreg can give degenerate results when")
cat(" the number of predictor variables is large\n\n")


}
\keyword{ regression }
\keyword{ survival } 
