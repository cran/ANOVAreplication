\name{pooled.sd}
\alias{pooled.sd}
\title{Pooled standard deviation calculator}
\description{
Calculates the pooled standard deviation.}
\usage{
pooled.sd(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A dataframe with two variables: the dependent variable in the first column, and the grouping variable in the second column. }
}
\value{
Returns the pooled standard deviation.
}
\author{
M.A.J. Zondervan-Zwijnenburg}
\examples{
data <- data.frame(y=ChickWeight$weight,g=ChickWeight$Diet)
pooled.sd(data)
}
