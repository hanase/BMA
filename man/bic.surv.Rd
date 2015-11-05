\name{bic.surv}
\alias{bic.surv}
\alias{bic.surv.data.frame}
\alias{bic.surv.matrix}
\alias{bic.surv.formula}
\title{Bayesian Model Averaging for Survival models.}
\description{
Bayesian Model Averaging for Cox proportional hazards models for
censored survival data. This accounts for the model uncertainty inherent
in the variable selection problem by averaging over the best models in
the model class according to approximate posterior model probability.
}

\usage{
bic.surv(x, ...)

\method{bic.surv}{matrix}(x, surv.t, cens, strict = FALSE, 
      OR = 20, maxCol = 30, prior.param = c(rep(0.5, ncol(x))), 
      OR.fix = 2, nbest = 150, factor.type = TRUE, 
      factor.prior.adjust = FALSE, call = NULL, ...)

\method{bic.surv}{data.frame}(x, surv.t, cens, 
      strict = FALSE, OR = 20, maxCol = 30, 
      prior.param = c(rep(0.5, ncol(x))), OR.fix = 2, 
      nbest = 150, factor.type = TRUE, 
      factor.prior.adjust = FALSE, call = NULL, ...)

\method{bic.surv}{formula}(f, data, strict = FALSE, 
     OR = 20, maxCol = 30, prior.param = c(rep(0.5, ncol(x))), 
     OR.fix = 2, nbest = 150, factor.type = TRUE, 
     factor.prior.adjust = FALSE, call = NULL, ...)
}

\arguments{
  \item{x}{a matrix or data frame of independent variables.}
  \item{surv.t}{a vector of values for the dependent variable.}
  \item{cens}{a vector of indicators of censoring (0=censored 1=uncensored)}
  \item{f}{a survival model formula}
  \item{data}{a data frame containing the variables in the model.}
  \item{strict}{ logical indicating whether models with more likely submodels are eliminated. 
        \code{FALSE} returns all models whose posterior model probability is within a factor of \code{1/OR} of that of the best model.}
  \item{OR}{ a number specifying the maximum ratio for excluding models in Occam's window }
  \item{maxCol}{ a number specifying the maximum number of columns in design matrix (including intercept) to be kept. }
 
  \item{prior.param}{a vector of prior probabilities that parameters are non-zero. Default puts a prior of .5 on all parameters.  Setting to 1 forces the variable into the model.   }
   \item{OR.fix}{width of the window which keeps models after the leaps approximation is done.  
        Because the leaps and bounds gives only an approximation to BIC, there is a need to increase the window at 
        this first "cut" so as to ensure that no good models are deleted. 
        The level of this cut is at \code{1/(OR^OR.fix)}; the default value for \code{OR.fix} is 2.}
  \item{nbest}{a value specifying the number of models of each size returned to bic.glm by the modified leaps algorithm.}
  \item{factor.type}{a logical value specifying how variables of class "factor" are handled. 
        A factor variable with d levels is turned into (d-1) dummy variables using a treatment contrast.  
        If \code{factor.type = TRUE}, models will contain either all or none of these dummy variables.  
        If \code{factor.type = FALSE}, models are free  to select the dummy variables independently.  
        In this case, factor.prior.adjust determines the prior on these variables.}
  \item{factor.prior.adjust}{a logical value specifying if the prior distribution on
        dummy variables for factors should be adjusted when \code{factor.type=FALSE}.
        When \code{factor.prior.adjust=FALSE}, all dummy variables for variable 
 \code{i} have prior equal to \code{prior.param[i]}.
        Note that this makes the prior probability of the union of these variables much higher than \code{prior.param[i]}.  
        Setting \code{factor.prior.adjust=T} corrects for this so that the union of the dummies equals \code{prior.param[i]}
        (and hence the deletion of the factor has a prior of \code{1-prior.param[i]}).  
        This adjustment changes the individual priors on each dummy variable to \code{1-(1-pp[i])^(1/(k+1))}.}
  \item{call}{used internally}
\item{...}{unused}
}
\details{
Bayesian Model Averaging accounts for the model uncertainty inherent in the variable selection problem by averaging over the best models in the model class according to approximate posterior model probability.  
\code{bic.surv} averages of Cox regression models.
}
\value{
 \code{bic.surv} returns an object of class \code{bic.surv}
  
The function \code{summary} is used to print a summary of the results. 
The function \code{plot} is used to plot posterior distributions for the 
coefficients. 
The function \code{imageplot} generates an image of the models which were 
averaged over.

An object of class \code{bic.glm} is a list containing at least the following 
components:

  \item{postprob}{the posterior probabilities of the models selected}
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
\note{If more than \code{maxcol} variables are supplied, then bic.surv does 
stepwise elimination of variables until \code{maxcol} variables are reached.
Many thanks to Sanford Weisberg for making source code for leaps available.
}

\references{ 
Volinsky, C.T., Madigan, D., Raftery, A.E. and Kronmal, R.A. (1997). 
"Bayesian Model Averaging in Proportional Hazard Models:  
Assessing the Risk of a Stroke."  Applied Statistics 46: 433-448
}

\author{
Chris Volinsky 
\email{volinsky@AT@research.att.com};
Adrian Raftery 
\email{raftery@AT@stat.washington.edu};
Ian Painter 
\email{ian.painter@AT@gmail.com}
}

\seealso{ \code{\link{summary.bic.surv}}, 
          \code{\link{print.bic.surv}}, 
          \code{\link{plot.bic.surv}} }

\examples{

\dontrun{
## veteran data
library(survival)
data(veteran)

test.bic.surv<- bic.surv(Surv(time,status) ~ ., data = veteran, 
                         factor.type = TRUE)
summary(test.bic.surv, conditional=FALSE, digits=2)
plot(test.bic.surv)

imageplot.bma(test.bic.surv)
}


## pbc data
data(pbc)
x<- pbc[1:312,]
surv.t<- x$time
cens<- as.numeric((x$status == 2))

x<- x[,c("age", "albumin", "alk.phos", "ascites", "bili", "edema", 
         "hepato", "platelet", "protime", "sex", "ast", "spiders", 
         "stage", "trt", "copper")]

\dontrun{
x$bili<- log(x$bili)
x$alb<- log(x$alb)
x$protime<- log(x$protime)
x$copper<- log(x$copper)
x$ast<- log(x$ast)

test.bic.surv<- bic.surv(x, surv.t, cens, 
                         factor.type=FALSE, strict=FALSE)
summary(test.bic.surv)
}

}
\keyword{regression}
\keyword{survival}
