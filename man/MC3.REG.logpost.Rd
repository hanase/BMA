\name{MC3.REG.logpost}
\alias{MC3.REG.logpost}
\title{ Helper function to MC3.REG  }
\description{
  Helper function to MC3.REG that calculates the posterior model probability (up to a constant).
}
\usage{
MC3.REG.logpost(Y, X, model.vect, p, i, K, nu, lambda, phi)
}
\arguments{
  \item{Y}{the vector of scaled responses.}
  \item{X}{the matrix of scaled covariates.}
  \item{model.vect}{ logical vector indicating which variables are to be included in the model }
  \item{p}{number of variables in model.vect}
  \item{i}{ vector of possible outliers}
  \item{K}{ a hyperparameter indicating the outlier inflation factor }
  \item{nu}{ regression hyperparameter. Default value is 2.58 if r2 for the full model is less than 0.9 
             or 0.2 if r2 for the full model is greater than 0.9. }
  \item{lambda}{ regression hyperparameter. Default value is 0.28 if r2 for the full model is less than 0.9 
             or 0.1684 if r2 for the full model is greater than 0.9. }
  \item{phi}{ regression hyperparameter. Default value is 2.85 if r2 for the full model is less than 0.9 
             or 9.2 if r2 for the full model is greater than 0.9.  }

}
\value{
  The log-posterior distribution for the model (up to a constant).
}
\references{    


Bayesian Model Averaging for Linear Regression Models 
Adrian E. Raftery, David Madigan, and Jennifer A. Hoeting (1997). 
Journal of the American Statistical Association, 92, 179-191. 

A Method for Simultaneous Variable and Transformation Selection in Linear Regression 
Jennifer Hoeting, Adrian E. Raftery and David Madigan (2002). 
Journal of Computational and Graphical Statistics 11 (485-507) 

A Method for Simultaneous Variable Selection and Outlier Identification in Linear Regression 
Jennifer Hoeting, Adrian E. Raftery and David Madigan (1996). 
Computational Statistics and Data Analysis, 22, 251-270 

Earlier versions of these papers are available via the World Wide Web using the url:
\url{http://www.stat.colostate.edu/~jah/papers/}

}
\author{Jennifer Hoeting \email{jah@AT@stat.colostate.edu} with the assistance of Gary Gadbury. Translation from Splus to R by Ian Painter \email{ian.painter@AT@gmail.com}.}
\note{ The implementation here differs from the Splus implentation. The Splus implementation uses global variables to contain the state of the current model and the history of the Markov-Chain. This implentation passes the current state and history to the function and then returns the updated state. }
\seealso{ \code{\link{MC3.REG}}, \code{\link{For.MC3.REG}}, \code{\link{MC3.REG.choose}} }
\keyword{regression}
