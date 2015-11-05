\name{bicreg}
\alias{bicreg}
\title{Bayesian Model Averaging for linear regression models.}
\description{Bayesian Model Averaging accounts for the model uncertainty inherent in the variable selection problem by averaging over the best models in the model class according to approximate posterior model probability.
}
\usage{
bicreg(x, y, wt = rep(1, length(y)), strict = FALSE, OR = 20, 
       maxCol = 31, drop.factor.levels = TRUE, nbest = 150)
}
\arguments{
  \item{x}{ a matrix of independent variables }
  \item{y}{ a vector of values for the dependent variable }
  \item{wt}{ a vector of weights for regression }
  \item{strict}{ logical. FALSE returns all models whose posterior model probability is within a factor of 1/OR of that of the best model. TRUE returns a more parsimonious set of models, where any model with a more likely submodel is eliminated.  }
  \item{OR}{ a number specifying the maximum ratio for excluding models in Occam's window }
  \item{maxCol}{ a number specifying the maximum number of columns in the design matrix (including the intercept) to be kept. }
  \item{drop.factor.levels}{ logical. Indicates whether factor levels can be individually dropped in the stepwise procedure to reduce the number of columns in the design matrix, or if a factor can be dropped only in its entirety. }
  \item{nbest}{ a value specifying the number of models of each size
    returned to bic.glm by the leaps algorithm. The default is 150
    (replacing the original default of 10).}
}
\details{
Bayesian Model Averaging accounts for the model uncertainty inherent in the variable selection problem by averaging over the best models in the model class according to the approximate posterior model probabilities.  
}
\value{
  \code{bicreg} returns an object of class \code{bicreg}
  
 The function 'summary' is used to print a summary of the results. The function 'plot' is used to plot posterior distributions for the coefficients.

An object of class \code{bicreg} is a list containing at least the following components:

  \item{postprob}{the posterior probabilities of the models selected}
  \item{namesx}{the names of the variables}
  \item{label}{labels identifying the models selected}
  \item{r2}{R2 values for the models}
  \item{bic}{values of BIC for the models}
  \item{size}{the number of independent variables in each of the models}
  \item{which}{a logical matrix with one row per model and one column per variable indicating whether that variable is in the model}
  \item{probne0}{the posterior probability that each variable is non-zero (in percent)}
  \item{postmean}{the posterior mean of each coefficient (from model averaging)}
  \item{postsd}{the posterior standard deviation of each coefficient (from model averaging) }
  \item{condpostmean}{the posterior mean of each coefficient conditional on the variable being included in the model}
  \item{condpostsd}{the posterior standard deviation of each coefficient conditional on the variable being included in the model}
  \item{ols}{matrix with one row per model and one column per variable giving the OLS estimate of each coefficient for each model}
  \item{se}{matrix with one row per model and one column per variable giving the standard error of each coefficient for each model}
  \item{reduced}{a logical indicating whether any variables were dropped before model averaging}
  \item{dropped}{a vector containing the names of those variables dropped before model averaging}
  \item{residvar}{residual variance for each model}
  \item{call}{the matched call that created the bicreg object}
}
\references{
Raftery, Adrian E. (1995). Bayesian model selection in social research (with Discussion). Sociological Methodology 1995 (Peter V. Marsden, ed.), pp. 111-196, Cambridge, Mass.: Blackwells.
 }
\author{ Original Splus code developed by Adrian Raftery (\email{raftery@AT@stat.washington.edu}) and revised by Chris T. Volinsky. Translation to R by Ian Painter. }
\seealso{ \code{\link{summary.bicreg}}, \code{\link{print.bicreg}}, \code{\link{plot.bicreg}} 
}
\examples{
library(MASS)
data(UScrime)
x<- UScrime[,-16]
y<- log(UScrime[,16])
x[,-2]<- log(x[,-2])
lma<- bicreg(x, y, strict = FALSE, OR = 20) 
summary(lma)
plot(lma)

imageplot.bma(lma)

}
\keyword{regression}
\keyword{models}
