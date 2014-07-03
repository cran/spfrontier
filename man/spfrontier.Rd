\name{spfrontier}
\alias{spfrontier}
\title{Spatial stochastic frontier model}
\usage{
spfrontier(formula, data, W_y = NULL, W_v = NULL, W_u = NULL,
  inefficiency = "half-normal", initialValues = NULL, logging = c("quiet",
  "info", "debug"), control = NULL, onlyCoef = F)
}
\arguments{
  \item{formula}{an object of class
  "\code{\link{formula}}": a symbolic description of the
  model to be fitted. The details of model specification
  are given under 'Details'.}

  \item{data}{data frame, containing the variables in the
  model}

  \item{W_y}{a spatial weight matrix for spatial lag of the
  dependent variable}

  \item{W_v}{a spatial weight matrix for spatial lag of the
  symmetric error term}

  \item{W_u}{a spatial weight matrix for spatial lag of the
  inefficiency error term}

  \item{initialValues}{an optional vector of initial
  values, used by maximum likelihood estimator. If not
  defined, estimator-specific method of initial values
  estimation is used.}

  \item{logging}{an optional level of logging. Possible
  values are 'quiet','warn','info','debug'. By default set
  to quiet.}

  \item{inefficiency}{sets the distribution for
  inefficiency error component. Possible values are
  'half-normal' (for half-normal distribution) and
  'truncated' (for truncated normal distribution). By
  default set to 'half-normal'. See references for
  explanations}

  \item{onlyCoef}{allows calculating only estimates for
  coefficients (with inefficiencies and other additional
  statistics). Developed generally for testing, to speed up
  the process.}

  \item{control}{an optional list of control parameters,
  passed to \code{\link{optim}} estimator from the
  '\code{\link{stats}} package}
}
\description{
\code{spfrontier} estimates spatial specifications of the
stochastic frontier model.
}
\details{
Models for estimation are specified symbolically, but
without any spatial components. Spatial components are
included implicitly on the base of the \code{model}
argument.
}
\examples{
data( airports )
airports2011 <- subset(airports, Year==2011)
W <- constructW(cbind(airports2011$longitude, airports2011$latitude),airports2011$ICAO)
formula <- log(PAX) ~ log(Population100km) + log(Routes) + log(GDPpc)
ols <- lm(formula , data=airports2011)
summary(ols )
plot(density(stats::residuals(ols)))
skewness(stats::residuals(ols))

# Takes >5 sec, see demo for more examples
# model <- spfrontier(formula , data=airports2011)
# summary(model )

# model <- spfrontier(formula , data=airports2011, W_y=W)
# summary(model )
}
\references{
Kumbhakar, S.C. and Lovell, C.A.K (2000), Stochastic
Frontier Analysis, Cambridge University Press, U.K.
}
\keyword{frontier}
\keyword{spatial}
\keyword{stochastic}

