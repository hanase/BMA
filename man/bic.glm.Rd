\name{bic.glm}
\alias{bic.glm}
\alias{bic.glm.data.frame}
\alias{bic.glm.matrix}
\alias{bic.glm.formula}
\title{Bayesian Model Averaging for generalized linear models.}
\description{
Bayesian Model Averaging accounts for the model uncertainty inherent in the variable selection problem by averaging over the best models in the model class according to approximate posterior model probability.
}

\usage{
bic.glm(x, ...)

\method{bic.glm}{matrix}(x, y, glm.family, wt = rep(1, nrow(x)),
    strict = FALSE, prior.param = c(rep(0.5, ncol(x))), OR = 20, 
    maxCol = 30, OR.fix = 2, nbest = 150, dispersion = NULL, 
    factor.type = TRUE, factor.prior.adjust = FALSE, 
    occam.window = TRUE, call = NULL, ...)

\method{bic.glm}{data.frame}(x, y, glm.family, wt = rep(1, nrow(x)),
    strict = FALSE, prior.param = c(rep(0.5, ncol(x))), OR = 20, 
    maxCol = 30, OR.fix = 2, nbest = 150, dispersion = NULL, 
    factor.type = TRUE, factor.prior.adjust = FALSE, 
    occam.window = TRUE, call = NULL, ...)

\method{bic.glm}{formula}(f, data, glm.family, wt = rep(1, nrow(data)),
    strict = FALSE, prior.param = c(rep(0.5, ncol(x))), OR = 20, 
    maxCol = 30, OR.fix = 2, nbest = 150, dispersion = NULL, 
    factor.type = TRUE, factor.prior.adjust = FALSE, 
    occam.window = TRUE, ...)
}
\arguments{
  \item{x}{a matrix or data.frame of independent variables.}
  \item{y}{a vector of values for the dependent variable.}
  \item{f}{a formula}
  \item{data}{a data frame containing the variables in the model. }
  \item{glm.family}{a description of the error distribution and link function to
            be used in the model.  This can be a character string naming a
            family function, a family function or the result of a call to
            a family function.  (See 'family' for details of family
            functions.)}
  \item{wt}{an optional vector of weights to be used.}
  \item{strict}{a logical indicating whether models with more likely submodels are 
    eliminated. \code{FALSE} returns all models whose posterior model probability is
     within a factor of \code{1/OR} of that of the best model.}
  \item{prior.param}{a vector of values specifying the prior weights for each variable.}
  \item{OR}{ a number specifying the maximum ratio for excluding models in Occam's window }
  \item{maxCol}{a number specifying the maximum number of columns in design matrix 
                (excluding intercept) to be kept. }
  \item{OR.fix}{width of the window which keeps models after the leaps approximation
                is done.  
                Because the leaps and bounds gives only an approximation to BIC, 
                there is a need to increase the window at this first "cut" so as to 
                assure that no good models are deleted. 
                The level of this cut is at \code{1/(OR^OR.fix)}; the default value 
                for \code{OR.fix} is 2.}
  \item{nbest}{a number specifying the number of models of each size returned to 
       \code{bic.glm} by the modified leaps algorithm.}
  \item{dispersion}{a logical value specifying whether dispersion should be 
        estimated or not. Default is \code{TRUE} unless glm.family is poisson 
        or binomial}
  \item{factor.type}{a logical value specifying how variables of class "factor" are 
        handled. 
        A factor variable with d levels is turned into (d-1) dummy variables using a
         treatment contrast.  
        If \code{factor.type = TRUE}, models will contain either all or none of 
        these dummy variables.  
        If \code{factor.type = FALSE}, models are free  to select the dummy 
        variables independently.  
        In this case, factor.prior.adjust determines the prior on these variables.}
  \item{factor.prior.adjust}{a logical value specifying whether 
        the prior distribution on dummy variables for factors 
        should be adjusted when \code{factor.type=FALSE}.
        When \code{factor.prior.adjust=FALSE}, all dummy variables 
        for variable \code{i} have prior equal to \code{prior.param[i]}.
        Note that this makes the prior probability of the union of these variables 
        much higher than \code{prior.param[i]}.  
        Setting \code{factor.prior.adjust=T} corrects for this so that the union of 
        the dummies equals \code{prior.param[i]}
        (and hence the deletion of the factor has a prior of 
        \code{1-prior.param[i]}).  
        This adjustment changes the individual priors on each dummy variable to '
        \code{1-(1-pp[i])^(1/(k+1))}.}
  \item{occam.window}{a logical value specifying if Occam's window should be used. 
        If set to \code{FALSE} then all models selected by the modified leaps 
        algorithm are returned.}
      \item{call}{used internally}
\item{...}{unused}
}
\details{
Bayesian Model Averaging accounts for the model uncertainty inherent in the variable selection problem by averaging over the best models in the model class according to approximate posterior model probability.  
}
\value{
  \code{bic.glm} returns an object of class \code{bic.glm}
  
The function \code{summary} is used to print a summary of the results. 
The function \code{plot} is used to plot posterior distributions for the coefficients. 
The function \code{imageplot} generates an image of the models which were averaged over.

An object of class \code{bic.glm} is a list containing at least the following components:

  \item{postprob}{the posterior probabilities of the models selected}
  \item{deviance}{the estimated model deviances}
  \item{label}{labels identifying the models selected}
  \item{bic}{values of BIC for the models}
  \item{size}{the number of independent variables in each of the models}
  \item{which}{a logical matrix with one row per model and one column per variable indicating whether that variable is in the model}
  \item{probne0}{the posterior probability that each variable is non-zero (in percent)}
  \item{postmean}{the posterior mean of each coefficient (from model averaging)}
  \item{postsd}{the posterior standard deviation of each coefficient (from model averaging) }
  \item{condpostmean}{the posterior mean of each coefficient conditional on the variable being included in the model}
  \item{condpostsd}{the posterior standard deviation of each coefficient conditional on the variable being included in the model}
  \item{mle}{matrix with one row per model and one column per variable giving the maximum likelihood estimate of each coefficient for each model}
  \item{se}{matrix with one row per model and one column per variable giving the standard error of each coefficient for each model}
  \item{reduced}{a logical indicating whether any variables were dropped before model averaging}
  \item{dropped}{a vector containing the names of those variables dropped before model averaging}
  \item{call}{the matched call that created the bma.lm object}
}

\references{ Raftery, Adrian E. (1995). Bayesian model selection in social research (with Discussion). Sociological Methodology 1995 (Peter V. Marsden, ed.), pp. 111-196, Cambridge, Mass.: Blackwells. 

 An earlier version, issued as Working Paper 94-12, Center for Studies in Demography 
and  Ecology, University of Washington (1994) is available as a technical report
from the Department of Statistics, University of Washington.}

\author{
Chris Volinsky \email{volinsky@research.att.com},
Adrian Raftery \email{raftery@stat.washington.edu}, and
Ian Painter \email{ian.painter@gmail.com}
}

\note{If more than \code{maxcol} variables are supplied, then bic.glm does stepwise 
elimination of variables until \code{maxcol} variables are reached.
\code{bic.glm} handles factor variables according to the \code{factor.type} 
parameter. If this is true then factor variables are kept in the model or dropped in 
entirety. If false, then each dummy variable can be kept or dropped independently. 
If \code{bic.glm} is used with a formula that includes interactions between factor 
variables, then \code{bic.glm} will create a new factor variable to represent that 
interaction, and this factor variable will be kept or dropped in entirety if 
\code{factor.type} is true. 
This can create interpretation problems if any of the corresponding main effects are
dropped.
Many thanks to Sanford Weisberg for making source code for leaps available.
}

\seealso{\code{\link{summary.bic.glm}}, 
         \code{\link{print.bic.glm}}, 
         \code{\link{plot.bic.glm}}}

\examples{

\dontrun{
### logistic regression
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
}

\dontrun{
### Gamma family 
library(survival)
data(veteran)
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
}

### Poisson family
### Yates (teeth) data. 

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

\dontrun{
### Gaussian
library(MASS)
data(UScrime)
f <- formula(log(y) ~  log(M)+So+log(Ed)+log(Po1)+log(Po2)+log(LF)+
                       log(M.F)+ log(Pop)+log(NW)+log(U1)+log(U2)+
                       log(GDP)+log(Ineq)+log(Prob)+log(Time))
glm.out.crime <- bic.glm(f, data = UScrime, glm.family = gaussian()) 
summary(glm.out.crime)
# note the problems with the estimation of the posterior standard 
# deviation (compare with bicreg example)
}

}
\keyword{regression}
\keyword{models}
