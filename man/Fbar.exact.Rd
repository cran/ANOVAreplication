\name{Fbar.exact}
\alias{Fbar.exact}
\title{F-bar for exact constraints}
\description{
The function calculates F-bar for hypotheses constrained with exact values (Type B). See Silvapulle & Sen (2011) for background on the F-bar statistic. The code of Vanbrabant (2017) is the basis to this Fbar function.}
\usage{
Fbar.exact(data,exact)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A dataframe with two variables: (1) a dependent variable, and (2) a grouping variable.}
\item{exact}{A vector of length p, where p is the number of means in the ANOVA model, with the exact values of the constrained hypothesis.}
}
\value{
The value for the F-bar statistic
}
\references{
Silvapulle, M. J., & Sen, P. K. (2011). Constrained statistical inference: Order, inequality, and shape constraints (Vol. 912). John Wiley & Sons. doi: 10.1002/9781118165614.ch1

Vanbrabant, L. (2017). restriktor: Restricted Statistical Estimation and Inference for Linear Models. R package version 0.1-55. https://CRAN.R-project.org/package=restriktor
}
\author{
M. A. J. Zondervan-Zwijnenburg
}
\seealso{
See also \code{\link{runShiny}}, \code{\link{Fbar.ineq}}, and \code{\link{Fbar.dif}}.
}
\examples{
data <- data.frame(y=ChickWeight$weight,g=ChickWeight$Diet)
aggregate(data$y,by=list(data$g),mean)

#specify exact values to be evaluated. Hi: mu1=102,mu2=123,mu3=143,mu4=135.
exact <- c(102,123,143,135)

Fbar.exact(data,exact)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}

