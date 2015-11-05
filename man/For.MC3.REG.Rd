\name{For.MC3.REG}
\alias{For.MC3.REG}
\title{ Helper function for MC3.REG }
\description{
 Helper function for MC3.REG which implements each step of the Metropolis-Hastings algorithm.
}
\usage{
For.MC3.REG(i, g, Ys, Xs, PI, K, nu, lambda, phi, outs.list)
}
\arguments{
  \item{i}{the current iteration number.}
  \item{g}{a list containing the current state and the history of the
    Markov-Chain. This list is in the same form as the return value (see
    the 'value' section below):
    \describe{
      \item{M0.var}{a logical vector specifying the variables in the
	current model.}
      \item{M0.out}{a logical vector specifying the outliers in the
	current model.}
      \item{M0.1}{a number representing the variables in the current
	model in binary form.}
    \item{M0.2}{a number represnting the outliers in the current model
      in binary form.}
    \item{outcnt}{the number of potential outliers} 
    }}
  \item{Ys}{the vector of scaled responses.}
  \item{Xs}{the matrix of scaled covariates.}
  \item{PI}{a hyperparameter indicating the prior probability of an outlier.  
        The default values are 0.1 if the data set has less than 50 observations, 0.02 otherwise. }
  \item{K}{ a hyperparameter indicating the outlier inflation factor }
  \item{nu}{ regression hyperparameter.
    Default value is 2.58 if r2 for the full model is less than 0.9 
             or 0.2 if r2 for the full model is greater than 0.9. }
	   \item{lambda}{ regression hyperparameter.
	     Default value is 0.28 if r2 for the full model is less than 0.9 
             or 0.1684 if r2 for the full model is greater than 0.9. }
	   \item{phi}{ regression hyperparameter.
	     Default value is 2.85 if r2 for the full model is less than 0.9 
             or 9.2 if r2 for the full model is greater than 0.9.  }
	   \item{outs.list}{ a vector of all potential outlier locations
	     (e.g. \code{c(10,12)} means the 10th and 12th points are
	     potential outliers).
	     If \code{NULL} and if \code{outliers} is \code{TRUE},
	     then potential outliers are estimated using the
	     \code{\link{out.ltsreg}} function.}
} 
\details{ This function implements a single Metropolis-Hastings step, choosing a proposal model, calculating the Bayes Factor between the current model and proposal model, and updating the current model to the proposal model if the step results in an update. 
}
\value{
  a list containing the current state and the history of the Markov-Chain, with components
    \item{flag}{a 0/1 number specifying whether the previous Metropolis-Hastings step resulted in a changed state or not.}
    \item{big.list}{a matrix containing the history of the Markov-Chain. Each row represents a unique model (combination of variables and outliers). The first column is the set of variables in the model (in binary form), the second column is the set of outliers in the model (in binary form), the third column is the log-posterior for the model (up to a constant) and the fourth column is the number of times that model has been visited.}
    \item{M0.var}{a logical vector specifying the variables in the current model.}
        \item{M0.out}{a logical vector specifying the outliers in the current model.}
    \item{M0.1}{a number representing the variables in the current model in binary form.}
    \item{M0.2}{a number represnting the outliers in the current model in binary form.}
    \item{outcnt}{the number of potential outliers}
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
\seealso{ \code{\link{MC3.REG}}, \code{\link{MC3.REG.choose}}, \code{\link{MC3.REG.logpost}} }
\keyword{regression}
