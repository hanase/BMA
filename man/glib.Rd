\name{glib}
\alias{glib}
\alias{glib.matrix}
\alias{glib.data.frame}
\alias{glib.bic.glm}
\alias{as.bic.glm}
\alias{as.bic.glm.glib}
\title{Model uncertainty in generalized linear models using Bayes factors}
\description{Function to evaluate Bayes factors and account for model uncertainty in generalized linear models. 
}

\usage{
glib(x, ...)

\method{glib}{matrix}(x, y, n = rep(1, nrow(x)), 
     error = "poisson", link = "log", scale = 1, 
     models = NULL, phi = c(1, 1.65, 5), psi = 1, 
     nu = 0, pmw = rep(1, nrow(models)), glimest = TRUE, 
     glimvar = FALSE, output.priorvar = FALSE, 
     post.bymodel = TRUE, output.postvar = FALSE, 
     priormean = NULL, priorvar = NULL, 
     nbest = 150, call = NULL, ...)

\method{glib}{data.frame}(x, y, n = rep(1, nrow(x)), 
     error = "poisson", link = "log", scale = 1, 
     models = NULL,  phi = c(1, 1.65, 5), 
     psi = 1, nu = 0, pmw = rep(1, nrow(models)), 
     glimest = TRUE, glimvar = FALSE, output.priorvar = FALSE, 
     post.bymodel = TRUE, output.postvar = FALSE, 
     priormean = NULL, priorvar = NULL, 
     nbest = 150, call = NULL, ...)

\method{glib}{bic.glm}(x, scale = 1, phi = 1, psi = 1, nu = 0, 
     glimest = TRUE, glimvar = FALSE, output.priorvar = FALSE, 
     post.bymodel = TRUE, output.postvar = FALSE, 
     priormean = NULL, priorvar = NULL, call = NULL, ...)

as.bic.glm(g, ...)

\method{as.bic.glm}{glib}( g, index.phi=1, ...)
}
\arguments{
  \item{x}{ an \code{n} x \code{p} matrix of independent variables  }
  \item{g}{ an object of type \code{bic.glm}}
  \item{y}{ a vector of values for the dependent variable }
  \item{n}{ an optional vector of weights to be used. }
  \item{error}{a string indicating the error family to use. Currently "gaussian", "gamma", "inverse gaussian", "binomial" and "poisson" are implemented. }
  \item{link}{a string indicating the link to use. Currently "identity", "log", "logit", "probit", "sqrt", "inverse" and "loglog" are implemented.}
  \item{scale}{the scale factor for the model. May be either a numeric constant or a string specifying the estimation, either "deviance" or "pearson". The default value is 1 for "binomial" and "poisson" error structures, and "pearson" for the others.}
  \item{models}{ an optional matrix representing the models to be averaged over. 
        \code{models} is a \code{n} x \code{p} matrix in which each row represents a model. 
        The corresponding entry in the row is 1 if that variable is included in the model; 0 if not. 
        The default value is \code{NULL} which will cause \code{glib} to call \code{bic.glm} with the parameter \code{occam.window} set to \code{FALSE} to obtain the models to average over.}
  \item{phi}{a vector of phi values. Default: \code{1}.}
  \item{psi}{a scalar prior parameter. Default: \code{1}. }
  \item{nu}{a scalar prior parameter. Default: 0 }
  \item{pmw}{a vector of prior model weights. These must be positive, but do not have to sum to one. 
        The prior model probabilities are given by  \code{pmw/sum(pmw)}. 
        The default is a vector of 1's of length \code{nrow(models)}}
  \item{glimest}{ a logical value specifying whether to output estimates and standard errors for each model. }
  \item{glimvar}{ a logical value specifying whether glim variance matrices are output for each model.}
  \item{output.priorvar}{ a logical value specifying whether the prior variance is output for each model and value of phi combination.}
  \item{post.bymodel}{a logical value specifying whether to output the posterior mean and sd for each model and value of phi combination.}
  \item{output.postvar}{a logical value specifying whether to output the posterior variance matrix for each model and value of phi combination.}
  \item{priormean}{an optional vector of length p+1 containing a user specified prior mean on the variables (including the intercept),
        where p=number of independent variables. }
  \item{priorvar}{an optional matrix containing a user specified prior variance matrix, a (p+1) x (p+1) matrix. 
        Default has the prior variance estimated as in Raftery(1996). }
  \item{nbest}{ an integer giving the number of best models of each size to be returned by bic.glm if \code{models == NULL}}
  \item{call}{the call to the function}
  \item{index.phi}{an index to the value of phi to use when converting a \code{glib} object to a \code{bic.glm} object}
\item{...}{unused}
}
\details{Function to evaluate Bayes factors and account for model
  uncertainty in generalized linear models. 
  This also calculates posterior distributions from a set of reference
  proper priors. 
\code{as.bic.glm} creates a 'bic.glm' object from a 'glib' object.
}
\value{
  \code{glib} returns an object of type \code{glib}, which is a list
  containing the following items:
    \item{inputs }{a list echoing the inputs
           (x,y,n,error,link,models,phi,psi,nu)}
  \item{bf }{a list containing the model comparison results:
    \describe{
    \item{twologB10}{an \code{nmodel} x \code{nphi} matrix whose
      \code{[i,j]} element is 2logB10 for model i against the null model
      with \code{phi=phi[j]}. A Laplace approximation (one-step Newton)
      is used. }
    \item{postprob}{ a matrix containing the posterior probabilities
                    of the models for each value of phi. }
    \item{deviance}{ a vector containing the deviances for the models. }
    \item{chi2}{ a vector containing the (DV0-DV1)/scale for the models }
    \item{npar}{ a vector containing the number of parameters estimated
      for each model. }
    \item{scale}{ the estimated or assigned scale used }
}}
  \item{posterior }{a list containing the Bayesian model mixing results:
    \describe{
      \item{prob0}{ an \code{ncol(x)} x \code{nphi} matrix whose [k,j]
        element is the posterior probability that the
        parameter corresponding to the k-th column of x is zero,
        for the j-th value of phi. }
      \item{mean}{ a \code{ncol(x)} x \code{nphi} matrix whose [k,j]
        element is the posterior mean of the parameter corresponding to
        the k-th column of x, for the j-th value of phi.}
    \item{sd}{ as for mean, but for the posterior standard deviation.
      NOTE: Both mean and sd are CONDITIONAL on the parameter being
      non-zero. They do not include the intercept. }
  }}
\item{glim.est }{a list containing the GLIM estimates for the different
   models:
   \describe{
     \item{coef}{An nmodel-list, each of whose elements is the coef
       value from "glim" for one of the models.}
    \item{se}{ as coef, but contains standard errors.}
    \item{var}{ as coef, but contains variance matrices of the
      estimates.}
    }}
  \item{posterior.bymodel }{a list containing model-specific posterior
    means and sds:
     \describe{
    \item{mean}{a list with \code{nmodel} elements, whose \code{i}th element
      is a \code{npar[i]xnphi} matrix, containing the posterior means of
      the npar[i] parameters of model \code{i},  for each value of phi.}
    \item{sd}{ as for mean, but for posterior standard deviations.}
    \item{var}{ a list with \code{nmodel} elements, whose \code{i}th element
      is a \code{npar[i]} by \code{npar[i]} by \code{nphi} array,
      containing the posterior variance matrix of the parameters of
      model i for each value of phi.}
}}  \item{prior }{a list containing the prior distributions:
 \describe{
   \item{mean}{ prior mean for the biggest model
     (this doesn't depend on phi)}
    \item{var}{ similar to corresponding member of posterior.bymodel.}
  }}
\item{models }{ an array containing the models used.}
\item{glm.out }{an object of type 'bic.glm' containing the results of
  any call to \code{bic.glm} }
\item{call }{the call to the function}
}
\references{ Raftery, A.E. (1988). Approximate Bayes factors for generalized linear  models. Technical Report no. 121, Department of Statistics, University of Washington. 
 
Raftery, Adrian E. (1995). Bayesian model selection in social research (with Discussion). Sociological Methodology 1995 (Peter V. Marsden, ed.), pp. 111-196, Cambridge, Mass.: Blackwells. 
 
Raftery, A.E. (1996). Approximate Bayes factors and accounting for model uncertainty in generalized linear models. Biometrika (83: 251-266). 
} 
\author{ Original Splus code developed by Adrian Raftery \email{raftery@AT@stat.washington.edu} and revised by Chris T. Volinsky. Translation to R by Ian S. Painter. }
\note{The outputs controlled by glimvar, output.priorvar and output.postvar can take up a lot of space, 
which is why these control parameters are F by default. }
\seealso{ \code{\link{bic.glm}}, \code{\link{summary.glib}} }
\examples{

\dontrun{
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
}

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

glib.yates <- glib ( x, y, n, models=models, glimvar=TRUE,
                     output.priorvar=TRUE, output.postvar=TRUE) 
summary(glib.yates)

\dontrun{
### logistic regression with no models specified
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

glib.birthwt<- glib(x,y, error="binomial", link = "logit")
summary(glib.birthwt)

glm.birthwt<- as.bic.glm(glib.birthwt)

imageplot.bma(glm.birthwt)
plot(glm.birthwt)
}


}
\keyword{regression}
\keyword{models}
