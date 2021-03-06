
#___________________________________________________________________________________________________________________
dKILG=
  function (x, mu = .1, sigma = 0.1, kinf=0, log = FALSE)
  {
    if (any(mu <= 0) | any(mu >= 1))
      stop(paste("mu must be greater than 0 and less than 1",
                 "\n", ""))
    if (any(sigma <= 0) | any(sigma >= 1))
      stop(paste("sigma must be between 0 and 1", "\n", ""))
    if (any(x <= 0))
      stop(paste("x must be 0 or greater than 0", "\n", ""))
    ly <- max(length(x), length(mu))
    x <- rep(x, length = ly)
    sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)
    inf<- rep(kinf, length = ly)
    fy <- dLG(x, mu = mu,log = T)
    logfy <- rep(kinf, ly)
    logfy <- ifelse((x == inf), log(sigma + (1 - sigma) * exp(fy)), (log(1 -
                                                                           sigma) + fy))
    if (log == FALSE)
      fy <- exp(logfy)
    else fy <- logfy
    fy
  }
#___________________________________________________________________________________________________________________

pKILG=
  function (q, mu = .1, sigma = 0.1, kinf=0, lower.tail = TRUE, log.p = FALSE)
  {
    if (any(mu <= 0) | any(mu >= 1))
      stop(paste("mu must be greater than 0 and less than 1",
                 "\n", ""))
    if (any(sigma <= 0) | any(sigma >= 1))
      stop(paste("sigma must be between 0 and 1", "\n", ""))
    if (any(q <= 0))
      stop(paste("y must be greater than 0", "\n", ""))
    ly <- max(length(q), length(mu), length(sigma))
    q <- rep(q, length = ly)
    mu <- rep(mu, length = ly)
    sigma <- rep(sigma, length = ly)
    inf<- rep(kinf, length = ly)
    cdf <- rep(kinf, ly)
    cdf <- pLG(q, mu = mu, lower.tail = TRUE, log.p = FALSE)
    cdf <- ifelse((q< inf), (1 - sigma) * cdf, sigma + (1 - sigma) * cdf)
    if (lower.tail == TRUE)
      cdf <- cdf
    else cdf <- 1 - cdf
    if (log.p == FALSE)
      cdf <- cdf
    else cdf <- log(cdf)
    cdf
  }

#___________________________________________________________________________________________________________________
qKILG=
  function (p, mu = 1, sigma = 0.1, kinf=0, lower.tail = TRUE, log.p = FALSE)
  {
    if (any(mu <= 0) | any(mu >= 1))
      stop(paste("mu must be greater than 0 and less than 1",
                 "\n", ""))
    if (any(sigma <= 0))
      stop(paste("sigma must be greater than 0", "\n", ""))
    if (any(p <= 0) | any(p >= 1))
      stop(paste("p must be between 0 and 1", "\n", ""))
    if (log.p == TRUE)
      p <- exp(p)
    else p <- p
    if (lower.tail == TRUE)
      p <- p
    else p <- 1 - p
    ly <- max(length(p), length(mu), length(sigma))
    p <- rep(p, length = ly)
    sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)
    inf<- rep(kinf, length = ly)
    pnew <-ifelse(p>pKILG(max(inf-1,1), mu=mu, sigma=sigma, kinf=inf),
                  (p - sigma)/(1 - sigma)-(1e-07), p /(1 - sigma)-(1e-07))
    pnew <-ifelse( inf==1,(p - sigma)/(1 - sigma)-(1e-07), pnew)
    pnew <- ifelse((pnew > 0), pnew, 0)
    q <- qLG(pnew, mu = mu,lower.tail = TRUE, log.p = FALSE, max.value = 10000)
    q
  }

#___________________________________________________________________________________________________________________
rKILG=
  function (n, mu = 1, sigma = 0.1, kinf=0)
  {
    if (any(mu <= 0) | any(mu >= 1))
      stop(paste("mu must be greater than 0 and less than 1",
                 "\n", ""))
    if (any(sigma <= 0))
      stop(paste("sigma must greated than 0", "\n", ""))
    if (any(n <= 0))
      stop(paste("n must be a positive integer", "\n", ""))
    n <- ceiling(n)
    p <- runif(n)
    r=c()
    for(i in 1:n){
      if(p[i]<sigma){r[i]=kinf}
      else{r[i] = rLG(1, mu = mu)}
    }
    r
  }

#___________________________________________________________________________________________________________________
KILG=
  function (mu.link = "logit", sigma.link = "logit", kinf="K")
  {
    mstats <- checklink("mu.link", "KILG", substitute(mu.link),
                        c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    dstats <- checklink("sigma.link", "KILG", substitute(sigma.link),
                        c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    structure(list(family = c(paste("inf",kinf,"LG", sep = ""), paste(kinf,"-inflated Logarithmic", sep = "") ),
                   parameters = list(mu = TRUE, sigma = TRUE), nopar = 2,
                   type = "Discrete", mu.link = as.character(substitute(mu.link)),
                   sigma.link = as.character(substitute(sigma.link)), mu.linkfun = mstats$linkfun,
                   sigma.linkfun = dstats$linkfun, mu.linkinv = mstats$linkinv,
                   sigma.linkinv = dstats$linkinv, mu.dr = mstats$mu.eta,
                   sigma.dr = dstats$mu.eta, dldm = function(y, mu, sigma) {
                     dldm0 <- (1 - sigma) * ((sigma + (1 - sigma) * dLG(kinf, mu))
                                             ^(-1)) * dLG(kinf, mu) * LG()$dldm(kinf, mu)
                     dldm <- ifelse(y == kinf, dldm0, LG()$dldm(y, mu))
                     dldm
                   }, d2ldm2 = function(y, mu, sigma) {
                     dldm0 <- (1 - sigma) * ((sigma + (1 - sigma) * dLG(kinf, mu))
                                             ^(-1)) * dLG(kinf, mu) * LG()$dldm(kinf, mu)
                     dldm <- ifelse(y == kinf, dldm0, LG()$dldm(y, mu))
                     d2ldm2 <- -dldm * dldm
                     d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                     d2ldm2
                   }, dldd = function(y, mu, sigma) {
                     dldd0 <- ((sigma + (1 - sigma) * dLG(kinf, mu))^(-1)) *
                       (1 - dLG(kinf, mu))
                     dldd <- ifelse(y == kinf, dldd0, -1/(1 - sigma))
                     dldd
                   }, d2ldd2 = function(y, mu, sigma) {
                     dldd0 <- ((sigma + (1 - sigma) * dLG(kinf, mu))^(-1)) *
                       (1 - dLG(kinf, mu))
                     dldd <- ifelse(y == kinf, dldd0, -1/(1 - sigma))
                     d2ldd2 <- -dldd * dldd
                     d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                     d2ldd2
                   }, d2ldmdd = function(y, mu, sigma) {
                     dldm0 <- (1 - sigma) * ((sigma + (1 - sigma) * dLG(kinf, mu))
                                             ^(-1)) * dLG(kinf, mu) * LG()$dldm(kinf, mu)
                     dldm <- ifelse(y == kinf, dldm0, LG()$dldm(y, mu))
                     dldd0 <- ((sigma + (1 - sigma) * dLG(kinf, mu))^(-1)) *
                       (1 - dLG(kinf, mu))
                     dldd <- ifelse(y == kinf, dldd0, -1/(1 - sigma))
                     d2ldmdd <- -dldm * dldd
                     d2ldmdd
                   }, G.dev.incr = function(y, mu, sigma, kinf, ...) -2 * dKILG(y,
                      mu, sigma, kinf=kinf, log = TRUE), rqres = expression(rqres(pfun = "pKILG",
                      type = "Discrete", ymin = 1, y = y, mu = mu, sigma = sigma, kinf=kinf)),
                      mu.initial = expression({mu <- 0.9}), sigma.initial = expression(sigma <- rep(0.1,length(y))),
                      mu.valid = function(mu) all(mu > 0 & mu < 1),sigma.valid = function(sigma) all(sigma > 0 & sigma <1),
                      y.valid = function(y) all(y > 0)),
              class = c("gamlss.family","family"))
  }
#___________________________________________________________________________________________________________________
