\name{KIGPO}
\alias{KIGPO}
\alias{dKIGPO}
\alias{pKIGPO}
\alias{qKIGPO}
\alias{rKIGPO}

\title{K-inflated Generalised Poisson distributions for fitting a GAMLSS model
}


\description{

The function \code{KIGPO} defines the K-inflated  Generalised Poisson distribution,a three  parameter distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting  using the function \code{gamlss()}. The functions \code{dKIGPO}, \code{pKIGPO}, \code{qKIGPO} and \code{rKIGPO} define the density, distribution function, quantile function and random generation for the K-inflated Generalised Poisson, \code{KIGPO()}, distribution.
}



\usage{

 KIGPO(mu.link = "log", sigma.link = "log", nu.link = "logit", kinf="K")

dKIGPO(x, mu = 1, sigma = 1, nu = 0.3, kinf=0 ,log = FALSE)

pKIGPO(q, mu = 1, sigma = 1, nu = 0.3, kinf=0, lower.tail = TRUE,
    log.p = FALSE)

qKIGPO(p, mu = 1, sigma = 1, nu = 0.3, kinf=0, lower.tail = TRUE,
    log.p = FALSE)

rKIGPO(n, mu = 1, sigma = 1, nu = 0.3, kinf=0)
}



\arguments{

\item{mu.link}{ Defines the \code{mu.link}, with "log" link as the default for the mu parameter}

\item{sigma.link}{Defines the  \code{sigma.link}, with "log" link as the default for the sigma parameter}

\item{nu.link}{Defines the \code{nu.link}, with "logit" link as the default for the nu parameter}

\item{x}{ vector of (non-negative integer) quantiles}

\item{mu}{vector of positive means}

\item{sigma}{vector of positive despersion parameter}

\item{nu}{vector of inflated point probability }


\item{p}{vector of probabilities}

\item{q}{vector of quantiles}

\item{n}{number of random values to return}

\item{kinf}{ defines inflated point in generating K-inflated distribution}

\item{log,log.p}{logical; if TRUE, probabilities p are given as log(p)}

\item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x],
          otherwise, P[X > x]}
}



\details{
The  definition for the K-inflated Generalised Poisson distribution.
}



\value{
The functions \code{KIGPO} return a \code{gamlss.family} object which can be used to fit K-inflated Generalised Poisson distribution in the \code{gamlss()} function.
}



\references{

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion),\emph{Appl. Statist.},\bold{54}, part 3, pp 507-554.



Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R. Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).



Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.



Rigby, R. A. and Stasinopoulos D. M. (2010) The gamlss.family distributions, (distributed with this package or see\url{http://www.gamlss.org/})



Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)\emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC.



Najafabadi, A. T. P. and MohammadPour, S. (2017). A k-Inflated Negative Binomial Mixture Regression Model: Application to Rate-Making Systems. Asia-Pacific Journal of Risk and Insurance, 12.
}



\author{
Saeed Mohammadpour <\email{s.mohammadpour1111@gamlil.com}>, Mikis Stasinopoulos <\email{d.stasinopoulos@londonmet.ac.uk}>
}


\seealso{
 \code{\link{gamlss.family}}, \code{\link{KIGPO}} }


\examples{
#--------------------------------------------------------------------------------

# gives information about the default links for the  Generalised Poisson distribution
KIGPO()
#--------------------------------------------------------------------------------

# generate zero inflated Generalised Poisson distribution
gen.Kinf(family=GPO, kinf=0)

# generate random sample from zero inflated Generalised Poisson distribution
x<-rinf0GPO(1000,mu=1, sigma=.5, nu=.2)

# fit the zero inflated Generalised Poisson distribution using gamlss
data<-data.frame(x=x)
\dontrun{
gamlss(x~1, family=inf0GPO, data=data)
histDist(x, family=inf0GPO)}
#--------------------------------------------------------------------------------

# generated one inflated Generalised Poisson distribution
gen.Kinf(family=GPO, kinf=1)

# generate random sample from one inflated Generalised Poisson distribution
x<-rinf1GPO(1000,mu=1, sigma=.5, nu=.2)

# fit the one inflated Generalised Poisson distribution using gamlss
data<-data.frame(x=x)
\dontrun{
gamlss(x~1, family=inf1GPO, data=data)
histDist(x, family=inf1GPO)}
#--------------------------------------------------------------------------------

mu=4; sigma=.5; nu=.2;
par(mgp=c(2,1,0),mar=c(4,4,4,1)+0.1)

#plot the pdf using plot
plot(function(x) dinf1GPO(x, mu=mu, sigma=sigma, nu=nu), from=0, to=20,
n=20+1, type="h",xlab="x",ylab="f(x)",cex.lab=1.5)
#--------------------------------------------------------------------------------

#plot the cdf using plot
cdf <- stepfun(0:19, c(0,pinf1GPO(0:19, mu=mu, sigma=sigma, nu=nu)), f = 0)
plot(cdf, xlab="x", ylab="F(x)", verticals=FALSE, cex.points=.8, pch=16, main="",cex.lab=1.5)
#--------------------------------------------------------------------------------

#plot the qdf using plot
invcdf <- stepfun(seq(0.01,.99,length=19), qinf1GPO(seq(0.1,.99,length=20),mu,        sigma), f = 0)
plot(invcdf, ylab=expression(x[p]==F^{-1}(p)), do.points=FALSE,verticals=TRUE,
     cex.points=.8, pch=16, main="",cex.lab=1.5, xlab="p")
#--------------------------------------------------------------------------------

# generate random sample
Ni <- rinf1GPO(1000, mu=mu, sigma=sigma, nu=nu)
 hist(Ni,breaks=seq(min(Ni)-0.5,max(Ni)+0.5,by=1),col="lightgray",main="",cex.lab=2)
barplot(table(Ni))
#--------------------------------------------------------------------------------
}

\keyword{distribution}
\keyword{regression}%


