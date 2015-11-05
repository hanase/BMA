\name{MC3.REG}
\alias{MC3.REG}
\alias{as.data.frame.mc3}
\alias{[.mc3}
\title{Bayesian simultaneous variable selection and outlier identification}
\description{
  Performs Bayesian simultaneous variable selection and outlier identification (SVO) via Markov chain Monte Carlo model
   composition (MC3).  
}
\usage{
MC3.REG(all.y, all.x, num.its, M0.var= , M0.out= , outs.list= , 
        outliers = TRUE, PI=.1*(length(all.y) <50) + 
        .02*(length(all.y) >= 50),  K=7, nu= , lambda= , phi= )

}
\arguments{
  \item{all.y}{a vector of responses}
  \item{all.x}{ a matrix of covariates}
  \item{num.its}{the number of iterations of the Markov chain sampler}
  \item{M0.var}{a logical vector specifying the starting model. For example, if you  have 3 predictors and the starting model is X1 and X3, then \code{M0.var} would be \code{c(TRUE,FALSE,TRUE)}. The default is a logical vector of \code{TRUE}s. NOTE: the starting predictor model cannot be the null model. }
  \item{M0.out}{a logical vector specifying the starting model outlier set. The default value is a logical vector of \code{TRUE}'s the same length as \code{outs.list}. This can be \code{NULL} only if outs.list is \code{NULL}, otherwise it must be the same length as \code{outs.list} (but can be a vector of all \code{FALSE}) }
  \item{outs.list}{ a vector of all potential outlier locations  (e.g. \code{c(10,12)} means the 10th and 12th points are potential outliers). If \code{NULL} and if \code{outliers} is \code{TRUE}, then potential outliers are estimated using the \code{\link{out.ltsreg}} function.}
  \item{outliers}{a logical parameter indicating whether outliers are to be included. If \code{outs.list} is non null then this \code{outliers} is ignored. If \code{outs.list} is \code{NULL} and outliers is \code{TRUE}, potential outliers are estimated as described above.}
  \item{PI}{a hyperparameter indicating the prior probability of an outlier. The default values are 0.1 if the data set has less than 50 observations, 0.02 otherwise. }
  \item{K}{ a hyperparameter indicating the outlier inflation factor }
  \item{nu}{ regression hyperparameter. Default value is 2.58 if r2 for the full model is less than 0.9 or 0.2 if r2 for the full model is greater than 0.9. }
  \item{lambda}{ regression hyperparameter. Default value is 0.28 if r2 for the full model is less than 0.9 or 0.1684 if r2 for the full model is greater than 0.9. }
  \item{phi}{ regression hyperparameter. Default value is 2.85 if r2 for the full model is less than 0.9 or 9.2 if r2 for the full model is greater than 0.9.  }
}
\details{
Performs Bayesian simultaneous variable and outlier selection using Monte Carlo Markov Chain Model Choice (MC3). 
Potential models are visited using a Metropolis-Hastings algorithm on the integrated likelihood. 
At the end of the chain exact posterior probabilities are calculated for each model visited.
}
\value{

 An object of class \code{mc3}. Print and summary methods exist for this class.
 Objects of class \code{mc3} are a list consisting of at least
\item{post.prob}{The posterior probabilities of each model visited.}
\item{variables}{An indicator matrix of the variables in each model.}
\item{outliers}{An indicator matrix of the outliers in each model, if outliers were selected.}
\item{visit.count}{The number of times each model was visited.}
\item{outlier.numbers}{An index showing which outliers were eligable for selection.}
\item{var.names}{The names of the variables.}
\item{n.models}{The number of models visited.}
\item{PI}{The value of PI used.}
\item{K}{The value of K used.}
\item{nu}{The value of nu used.}
\item{lambda}{The value of lambda used.}
\item{phi}{The value of phi used.}
\item{call}{The function call.}

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
\author{Jennifer Hoeting \email{jah@AT@stat.colostate.edu} with the assistance of Gary Gadbury. Translation from Splus to R by Ian S. Painter.}
\note{

The default values for \code{nu}, \code{lambda} and \code{phi} are recommended when the R2 value for the full model with all outliers is less than 0.9.

If \code{PI} is set too high it is possible to generate sub models which are singular, at which point the function will crash.

The implementation of this function is different from that used in the Splus function. In particular, variables which were global are now passed between functions.}

\seealso{\code{\link{out.ltsreg}} \code{\link{as.data.frame.mc3}}}
\examples{

\dontrun{
# Example 1:   Scottish hill racing data.

data(race)
b<- out.ltsreg(race[,-1], race[,1], 2)
races.run1<-MC3.REG(race[,1], race[,-1], num.its=20000, c(FALSE,TRUE), 
                    rep(TRUE,length(b)), b, PI = .1, K = 7, nu = .2, 
                    lambda = .1684, phi = 9.2)
races.run1
summary(races.run1)
}

# Example 2: Crime data
library(MASS)
data(UScrime)

y.crime.log<- log(UScrime$y)
x.crime.log<- UScrime[,-ncol(UScrime)]
x.crime.log[,-2]<- log(x.crime.log[,-2])
crime.run1<-MC3.REG(y.crime.log, x.crime.log, num.its=2000, 
                     rep(TRUE,15), outliers = FALSE)
crime.run1[1:25,]
summary(crime.run1)




}
\keyword{regression}
\keyword{models}
