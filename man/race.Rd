\name{race}
\alias{race}
\docType{data}
\title{ Scottish Hill Racing data}
\description{
The record-winning times for 35 hill races in Scotland, as reported by Atkinson (1986). 
}
\usage{data(race)}
\format{
data.frame


}
\details{
  The distance travelled and the height climbed in each race is also given. The data contains a known error - Atkinson (1986) reports that the record for Knock Hill (observation 18) should actually be 18 minutes rather than 78 minutes.
\describe{
\item{Variable}{Description}
\item{Race}{Name of race} 
\item{Distance}{Distance covered in miles} 
\item{Climb}{Elevation climbed during race in feet}
\item{Time}{Record time for race in minutes}
}}
\references{
Atkison, A.C., Comments on "Influential Observations, High Leverage Points, and Outliers in Linear Regression", Statistical Science, 1 (1986) 397-402

}
\source{
 \url{http://www.statsci.org/data/general/hills.html}
}
\keyword{datasets}
