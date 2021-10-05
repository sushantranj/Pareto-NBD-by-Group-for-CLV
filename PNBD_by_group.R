
# Main model params are :
# beta unobserved shape parameter for dropout process
# s unobserved scale parameter for dropout process
# r unobserved shape parameter for NBD transaction
# alpha unobserved scale parameter for NBD transaction

############################################
## LOAD NEEDED PACKAGES
############################################

library(ggplot2)
library(BTYD)
library(reshape2)
library(plyr)
library(lubridate)
library(hypergeo)
library(gsl)
library(optimx)
library(stats)
library(ucminf)
library(minqa)
library(Rcgmin)
library(Rvmmin)
library(BB)

############################################
## LOAD DATA
#############################################

# Clear the workspace

rm(list = ls(all = TRUE))

# Choose a physical location for the log file export

sink("C:/Users/susha/Downloads/clv/logonlycalibration.txt" ,
     append = FALSE ,
     split = TRUE)

start.time <- Sys.time()

# Choose an end date for calibrating the input period

end.of.cal.period <- as.Date("2021-12-31")

options(scipen = 999)
options (warn = -1)

# Import the data 

my_data <- read.csv("marketing.csv", stringsAsFactors = FALSE)
my_data$date <- as.Date(my_data$date, format = "%d/%m/%Y")
my_data$cust <- as.character(my_data$cust)
my_data$region <- as.character(my_data$region)

# Check and view the data snapshot

head(my_data,3)


#Start of the functions used for PNBD

h2f1 <- function(a, b, c, z) {
  lenz <- length(z)
  j = 0
  uj <- 1:lenz
  uj <- uj / uj
  y <- uj
  lteps <- 0
  while (lteps < lenz) {
    lasty <- y
    j <- j + 1
    uj <- uj * (a + j - 1) * (b + j - 1) / (c + j - 1) * z / j
    y <- y + uj
    lteps <- sum(y == lasty)
  }
  y
}


pnbd.cbs.LL <- function(params, cal.cbs) {
  dc.check.model.params(c("r", "alpha", "s", "beta"), params, "pnbd.cbs.LL")
  
  tryCatch(
    x <-
      cal.cbs[, "x"],
    error = function(e)
      stop(
        "Error in pnbd.cbs.LL: cal.cbs must have a frequency column labelled \"x\""
      )
  )
  tryCatch(
    t.x <-
      cal.cbs[, "t.x"],
    error = function(e)
      stop(
        "Error in pnbd.cbs.LL: cal.cbs must have a recency column labelled \"t.x\""
      )
  )
  tryCatch(
    T.cal <-
      cal.cbs[, "T.cal"],
    error = function(e)
      stop(
        "Error in pnbd.cbs.LL: cal.cbs must have a column for length of time observed labelled \"T.cal\""
      )
  )
  
  if ("custs" %in% colnames(cal.cbs)) {
    custs <- cal.cbs[, "custs"]
  } else {
    custs <- rep(1, length(x))
  }
  return(sum(custs * pnbd.LL(params, x, t.x, T.cal)))
}


pnbd.LL = function (params, x, t.x, T.cal)
{
  max.length <- max(length(x), length(t.x), length(T.cal))
  if (max.length %% length(x))
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length %% length(t.x))
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length %% length(T.cal))
    warning("Maximum vector length not a multiple of the length of T.cal")
  dc.check.model.params(c("r", "alpha", "s", "beta"), params,
                        "pnbd.LL")
  if (any(x < 0) || !is.numeric(x))
    stop("x must be numeric and may not contain negative numbers.")
  if (any(t.x < 0) || !is.numeric(t.x))
    stop("t.x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal))
    stop("T.cal must be numeric and may not contain negative numbers.")
  x <- rep(x, length.out = max.length)
  t.x <- rep(t.x, length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  r <- params[1]
  alpha <- params[2]
  s <- params[3]
  beta <- params[4]
  maxab <- max(alpha, beta)
  absab <- abs(alpha - beta)
  param2 <- s + 1
  if (alpha < beta) {
    param2 <- r + x
  }
  part1 <-
    r * log(alpha) + s * log(beta) - lgamma(r) + lgamma(r + x)
  part2 <- -(r + x) * log(alpha + T.cal) - s * log(beta + T.cal)
  if (absab == 0) {
    part2_times_F1_min_F2 <-
      ((alpha + T.cal) / (maxab + t.x)) ^ (r + x) * (beta + T.cal) ^ s /
      ((maxab + t.x) ^ s) -
      ((alpha + T.cal) / (maxab + T.cal)) ^ (r + x) * (beta + T.cal) ^ s /
      ((maxab + t.x) ^ s)
  }
  else {
    part2_times_F1_min_F2 = h2f1(r + s + x, param2, r + s + x + 1, absab / (maxab +
                                                                              t.x)) *
      ((alpha + T.cal) / (maxab + t.x)) ^ (r + x) * (beta + T.cal) ^ s /
      ((maxab + t.x) ^ s) -
      h2f1(r + s + x, param2, r + s + x + 1, absab / (maxab + T.cal)) *
      ((alpha + T.cal) / (maxab + T.cal)) ^ (r + x) * (beta + T.cal) ^ s /
      ((maxab + t.x) ^ s)
  }
  return(part1 + part2 + log(1 + (s / (r + s + x)) * part2_times_F1_min_F2))
}


pnbd.compress.cbs <- function(cbs, rounding = 3) {
  if (!("x" %in% colnames(cbs)))
    stop("Error in pnbd.compress.cbs: cbs must have a frequency column labelled \"x\"")
  if (!("t.x" %in% colnames(cbs)))
    stop("Error in pnbd.compress.cbs: cbs must have a recency column labelled \"t.x\"")
  if (!("T.cal" %in% colnames(cbs)))
    stop(
      "Error in pnbd.compress.cbs: cbs must have a column for length of time observed labelled \"T.cal\""
    )
  
  orig.rows <- nrow(cbs)
  
  if (!("custs" %in% colnames(cbs))) {
    custs <- rep(1, nrow(cbs))
    cbs <- cbind(cbs, custs)
  }
  
  other.colnames <-
    colnames(cbs)[!(colnames(cbs) %in% c("x", "t.x", "T.cal"))]
  
  ## Round x, t.x and T.cal to the desired level
  cbs[, c("x", "t.x", "T.cal")] <-
    round(cbs[, c("x", "t.x", "T.cal")], rounding)
  
  ## Aggregate every column that is not x, t.x or T.cal by those columns. Do this
  ## by summing entries which have the same x, t.x and T.cal.
  cbs <-
    as.matrix(aggregate(cbs[,!(colnames(cbs) %in% c("x", "t.x", "T.cal"))],
                        by = list(
                          x = cbs[, "x"],
                          t.x = cbs[, "t.x"],
                          T.cal = cbs[, "T.cal"]
                        ), sum))
  
  colnames(cbs) <- c("x", "t.x", "T.cal", other.colnames)
  final.rows <- nrow(cbs)
  message("Data reduced from ", orig.rows, " rows to ", final.rows, " rows.")
  return(cbs)
}


pnbd.EstimateParameters <-
  function(cal.cbs,
           par.start = c(2, 2, 2, 2),
           max.param.value = 10000) {
    dc.check.model.params(c("r", "alpha", "s", "beta"),
                          par.start,
                          "pnbd.EstimateParameters")
    
    
    ## helper function to be optimized
    pnbd.eLL <- function(params, cal.cbs, max.param.value) {
      params <- exp(params)
      params[params > max.param.value] <- max.param.value
      return(-1 * pnbd.cbs.LL(params, cal.cbs))
    }
    logparams <- log(par.start)
    #logparams=c(0,0,0,0)
    results <-
      optimx(
        logparams,
        pnbd.eLL,
        cal.cbs = cal.cbs,
        max.param.value = max.param.value,
        method = c("nlminb","spg"),
        hessian = FALSE ,
        control = list(
          trace = 2,
          kkt = TRUE,
          usenumDeriv = TRUE ,
          save.failures = TRUE,
          dowarn = FALSE
        )
      )
    # results <- optimx(logparams, pnbd.eLL,cal.cbs=cal.cbs, max.param.value=max.param.value, method=c("spg"), hessian=FALSE ,control=list(trace=2,kkt=TRUE,usenumDeriv=TRUE ,save.failures=TRUE,dowarn=FALSE))
    results <- summary(results , order = list(convcode, value))[1,]
    r <- results[1]
    alpha <- results[2]
    s <- results[3]
    beta <- results[4]
    results <- as.numeric(c(r, alpha, s, beta))
    estimated.params <- exp(results)
    estimated.params[estimated.params > max.param.value] <-
      max.param.value
    return(estimated.params)
    
  }


pnbd.Estimateconverge <-
  function(cal.cbs,
           par.start = c(2, 2, 2, 2),
           max.param.value = 10000) {
    dc.check.model.params(c("r", "alpha", "s", "beta"),
                          par.start,
                          "pnbd.EstimateParameters")
    
    
    ## helper function to be optimized
    pnbd.eLL <- function(params, cal.cbs, max.param.value) {
      params <- exp(params)
      params[params > max.param.value] <- max.param.value
      return(-1 * pnbd.cbs.LL(params, cal.cbs))
    }
    logparams <- log(par.start)
    #logparams=c(0,0,0,0)
    results <-
      optimx(
        logparams,
        pnbd.eLL,
        cal.cbs = cal.cbs,
        max.param.value = max.param.value,
        method = c("nlminb"),
        hessian = FALSE ,
        control = list(
          trace = 2,
          kkt = TRUE,
          usenumDeriv = TRUE ,
          save.failures = TRUE,
          dowarn = FALSE
        )
      )
    # results <- optimx(logparams, pnbd.eLL,cal.cbs=cal.cbs, max.param.value=max.param.value, method=c("spg"), hessian=FALSE ,control=list(trace=2,kkt=TRUE,usenumDeriv=TRUE ,save.failures=TRUE,dowarn=FALSE))
    results <- summary(results , order = list(convcode, value))[1,]
    conv1 <- as.character(results[10])
    conv2 <- as.character(results[11])
    processtime  <- as.numeric(results[12])
    return(list(conv1, conv2, processtime))
    
  }



pnbd.pmf <- function(params, t, x) {
  max.length <- max(length(t), length(x))
  if (max.length %% length(t))
    warning("Maximum vector length not a multiple of the length of t")
  if (max.length %% length(x))
    warning("Maximum vector length not a multiple of the length of x")
  
  dc.check.model.params(c("r", "alpha", "s", "beta"), params, "pnbd.pmf")
  
  if (any(t < 0) || !is.numeric(t))
    stop("t must be numeric and may not contain negative numbers.")
  if (any(x < 0) || !is.numeric(x))
    stop("x must be numeric and may not contain negative numbers.")
  
  t. <- rep(t, length.out = max.length)
  x <- rep(x, length.out = max.length)
  
  return(pnbd.pmf.General(params, 0, t, x))
}


pnbd.pmf.General <- function(params, t.start, t.end, x) {
  max.length <- max(length(t.start), length(t.end), length(x))
  
  if (max.length %% length(t.start))
    warning("Maximum vector length not a multiple of the length of t.start")
  if (max.length %% length(t.end))
    warning("Maximum vector length not a multiple of the length of t.end")
  if (max.length %% length(x))
    warning("Maximum vector length not a multiple of the length of x")
  
  dc.check.model.params(c("r", "alpha", "s", "beta"), params, "pnbd.pmf.General")
  
  if (any(t.start < 0) || !is.numeric(t.start))
    stop("t.start must be numeric and may not contain negative numbers.")
  if (any(t.end < 0) || !is.numeric(t.end))
    stop("t.end must be numeric and may not contain negative numbers.")
  if (any(x < 0) || !is.numeric(x))
    stop("x must be numeric and may not contain negative numbers.")
  
  
  t.start <- rep(t.start, length.out = max.length)
  t.end <- rep(t.end, length.out = max.length)
  x <- rep(x, length.out = max.length)
  
  if (any(t.start > t.end)) {
    stop("Error in pnbd.pmf.General: t.start > t.end.")
  }
  
  r <- params[1]
  alpha <- params[2]
  s <- params[3]
  beta <- params[4]
  
  equation.part.0 <- rep(0, max.length)
  equation.part.0[x == 0] <-
    1 - exp(s * log(beta) - s * log(beta + t.start))
  
  ## (t.end - t.start)^x is left outside the exp() for this reason: exp(0 *
  ## log(0))=NaN; 0^0=1.
  equation.part.1 <-
    exp(
      lgamma(r + x) - lgamma(r) - lfactorial(x) + r * log(alpha) -
        r * log(alpha + t.end - t.start) - x * log(alpha + t.end - t.start) + s *
        log(beta) - s * log(beta + t.end)
    ) * (t.end - t.start) ^ x
  
  equation.part.2 <-
    r * log(alpha) + s * log(beta) + lbeta(r + x, s + 1) - lbeta(r,
                                                                 s)
  
  B.1 <- rep(NA, max.length)
  B.1[alpha > beta + t.start] <-
    hyperg_2F1(r + s, s + 1, r + s + x + 1, (alpha -
                                               beta - t.start) /
                 alpha) / (alpha ^ (r + s))
  B.1[alpha <= beta + t.start] <-
    hyperg_2F1(r + s, r + x, r + s + x + 1, (beta +
                                               t.start - alpha) /
                 (beta + t.start)) / ((beta + t.start) ^ (r + s))
  
  B.2 <- function(r, alpha, s, beta, t.start, t.end, x, ii) {
    if (alpha > (beta + t.start)) {
      hyperg_2F1(r + s + ii,
                 s + 1,
                 r + s + x + 1,
                 (alpha - beta - t.start) / (alpha +
                                               t.end - t.start)) /
        ((alpha + t.end - t.start) ^ (r + s + ii))
    } else {
      hyperg_2F1(r + s + ii, r + x, r + s + x + 1, (beta + t.start - alpha) /
                   (beta +
                      t.end)) /
        ((beta + t.end) ^ (r + s + ii))
    }
  }
  
  equation.part.2.summation <- rep(NA, max.length)
  
  for (i in 1:max.length) {
    ii <- c(1:x[i])
    equation.part.2.summation[i] <-
      B.2(r, alpha, s, beta, t.start[i], t.end[i],
          x[i], 0)
    if (x[i] > 0) {
      equation.part.2.summation[i] = equation.part.2.summation[i] + sum((t.end[i] -
                                                                           t.start[i]) ^
                                                                          ii / (ii * beta(r + s, ii)) * B.2(r, alpha, s, beta, t.start[i],
                                                                                                            t.end[i], x[i], ii)
      )
    }
  }
  return(equation.part.0 + equation.part.1 + exp(equation.part.2 + log(B.1 - equation.part.2.summation)))
}




pnbd.ConditionalExpectedTransactions <-
  function(params, T.star, x, t.x,
           T.cal) {
    max.length <-
      max(length(T.star), length(x), length(t.x), length(T.cal))
    
    if (max.length %% length(T.star))
      warning("Maximum vector length not a multiple of the length of T.star")
    if (max.length %% length(x))
      warning("Maximum vector length not a multiple of the length of x")
    if (max.length %% length(t.x))
      warning("Maximum vector length not a multiple of the length of t.x")
    if (max.length %% length(T.cal))
      warning("Maximum vector length not a multiple of the length of T.cal")
    
    dc.check.model.params(c("r", "alpha", "s", "beta"),
                          params,
                          "pnbd.ConditionalExpectedTransactions")
    
    if (any(T.star < 0) || !is.numeric(T.star))
      stop("T.star must be numeric and may not contain negative numbers.")
    if (any(x < 0) || !is.numeric(x))
      stop("x must be numeric and may not contain negative numbers.")
    if (any(t.x < 0) || !is.numeric(t.x))
      stop("t.x must be numeric and may not contain negative numbers.")
    if (any(T.cal < 0) || !is.numeric(T.cal))
      stop("T.cal must be numeric and may not contain negative numbers.")
    
    
    T.star <- rep(T.star, length.out = max.length)
    x <- rep(x, length.out = max.length)
    t.x <- rep(t.x, length.out = max.length)
    T.cal <- rep(T.cal, length.out = max.length)
    
    r <- params[1]
    alpha <- params[2]
    s <- params[3]
    beta <- params[4]
    
    P1 <- (r + x) * (beta + T.cal) / ((alpha + T.cal) * (s - 1))
    P2 <- (1 - ((beta + T.cal) / (beta + T.cal + T.star)) ^ (s - 1))
    P3 <- pnbd.PAlive(params, x, t.x, T.cal)
    return(P1 * P2 * P3)
    
  }


# pnbd.PAlive
pnbd.PAlive <- function(params, x, t.x, T.cal) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  
  if (max.length %% length(x))
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length %% length(t.x))
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length %% length(T.cal))
    warning("Maximum vector length not a multiple of the length of T.cal")
  
  dc.check.model.params(c("r", "alpha", "s", "beta"), params, "pnbd.PAlive")
  
  if (any(x < 0) || !is.numeric(x))
    stop("x must be numeric and may not contain negative numbers.")
  if (any(t.x < 0) || !is.numeric(t.x))
    stop("t.x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal))
    stop("T.cal must be numeric and may not contain negative numbers.")
  
  
  x <- rep(x, length.out = max.length)
  t.x <- rep(t.x, length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  
  r <- params[1]
  alpha <- params[2]
  s <- params[3]
  beta <- params[4]
  
  A0 <- 0
  if (alpha >= beta) {
    F1 <-
      hyperg_2F1(r + s + x, s + 1, r + s + x + 1, (alpha - beta) / (alpha +
                                                                      t.x))
    F2 <-
      hyperg_2F1(r + s + x, s + 1, r + s + x + 1, (alpha - beta) / (alpha +
                                                                      T.cal))
    # A0 <- F1/((alpha + t.x)^(r + s + x)) - F2/((alpha + T.cal)^(r + s + x))
    X1 <-
      F1 * ((alpha + T.cal) / (alpha + t.x)) ^ (r + x) * ((beta + T.cal) / (alpha +
                                                                              t.x)) ^ s
    X2 <- F2 * ((beta + T.cal) / (alpha + T.cal)) ^ s
    
  } else {
    F1 <-
      hyperg_2F1(r + s + x, r + x, r + s + x + 1, (beta - alpha) / (beta +
                                                                      t.x))
    F2 <-
      hyperg_2F1(r + s + x, r + x, r + s + x + 1, (beta - alpha) / (beta +
                                                                      T.cal))
    # A0 <- F1/((beta + t.x)^(r + s + x)) - F2/((beta + T.cal)^(r + s + x))
    X1 <-
      F1 * ((alpha + T.cal) / (beta + t.x)) ^ (r + x) * ((beta + T.cal) / (beta +
                                                                             t.x)) ^ s
    X2 <- F2 * ((alpha + T.cal) / (beta + T.cal)) ^ (r + x)
    return((1 + s / (r + s + x) * (X1 - X2)) ^ (-1))
  }
  # return((1 + s/(r + s + x) * (alpha + T.cal)^(r + x) * (beta + T.cal)^s * A0)^(-1))
  return((1 + s / (r + s + x) * (X1 - X2)) ^ (-1))
}


pnbd.Expectation <- function(params, t) {
  dc.check.model.params(c("r", "alpha", "s", "beta"), params, "pnbd.Expectation")
  
  if (any(t < 0) || !is.numeric(t))
    stop("t must be numeric and may not contain negative numbers.")
  
  r = params[1]
  alpha = params[2]
  s = params[3]
  beta = params[4]
  
  return((r * beta) / (alpha * (s - 1)) * (1 - (beta / (beta + t)) ^ (s - 1)))
}


pnbd.ExpectedCumulativeTransactions <-
  function(params, T.cal, T.tot,
           n.periods.final) {
    dc.check.model.params(c("r", "alpha", "s", "beta"),
                          params,
                          "pnbd.ExpectedCumulativeTransactions")
    
    if (any(T.cal < 0) || !is.numeric(T.cal))
      stop("T.cal must be numeric and may not contain negative numbers.")
    
    if (length(T.tot) > 1 || T.tot < 0 || !is.numeric(T.tot))
      stop("T.cal must be a single numeric value and may not be negative.")
    if (length(n.periods.final) > 1 ||
        n.periods.final < 0 || !is.numeric(n.periods.final))
      stop("n.periods.final must be a single numeric value and may not be negative.")
    
    ## Divide up time into equal intervals
    intervals <-
      seq(T.tot / n.periods.final, T.tot, length.out = n.periods.final)
    
    cust.birth.periods <- max(T.cal) - T.cal
    
    expected.transactions <- sapply(intervals, function(interval) {
      if (interval <= min(cust.birth.periods))
        return(0)
      sum(pnbd.Expectation(params, interval - cust.birth.periods[cust.birth.periods <=
                                                                   interval]))
    })
    
    return(expected.transactions)
  }



pnbd.DERT <- function(params, x, t.x, T.cal, d) {
  max.length <- max(length(x), length(t.x), length(T.cal))
  
  if (max.length %% length(x))
    warning("Maximum vector length not a multiple of the length of x")
  if (max.length %% length(t.x))
    warning("Maximum vector length not a multiple of the length of t.x")
  if (max.length %% length(T.cal))
    warning("Maximum vector length not a multiple of the length of T.cal")
  
  dc.check.model.params(c("r", "alpha", "s", "beta"), params, "pnbd.DERT")
  
  if (any(x < 0) || !is.numeric(x))
    stop("x must be numeric and may not contain negative numbers.")
  if (any(t.x < 0) || !is.numeric(t.x))
    stop("t.x must be numeric and may not contain negative numbers.")
  if (any(T.cal < 0) || !is.numeric(T.cal))
    stop("T.cal must be numeric and may not contain negative numbers.")
  
  
  x <- rep(x, length.out = max.length)
  t.x <- rep(t.x, length.out = max.length)
  T.cal <- rep(T.cal, length.out = max.length)
  
  r <- params[1]
  alpha <- params[2]
  s <- params[3]
  beta <- params[4]
  
  maxab = max(alpha, beta)
  absab = abs(alpha - beta)
  param2 = s + 1
  if (alpha < beta) {
    param2 = r + x
  }
  part1 <- (alpha ^ r * beta ^ s / gamma(r)) * gamma(r + x)
  part2 <- 1 / ((alpha + T.cal) ^ (r + x) * (beta + T.cal) ^ s)
  if (absab == 0) {
    F1 <- 1 / ((maxab + t.x) ^ (r + s + x))
    F2 <- 1 / ((maxab + T.cal) ^ (r + s + x))
  } else {
    F1 <-
      hyperg_2F1(r + s + x, param2, r + s + x + 1, absab / (maxab + t.x)) / ((maxab +
                                                                                t.x) ^
                                                                               (r + s + x))
    F2 <-
      hyperg_2F1(r + s + x, param2, r + s + x + 1, absab / (maxab + T.cal)) /
      ((maxab +
          T.cal) ^
         (r + s + x))
  }
  
  likelihood = part1 * (part2 + (s / (r + s + x)) * (F1 - F2))
  
  z <- d * (beta + T.cal)
  
  tricomi.part.1 = ((z) ^ (1 - s)) / (s - 1) * hyperg_1F1(1, 2 - s, z)
  tricomi.part.2 = gamma(1 - s) * hyperg_1F1(s, s, z)
  tricomi = tricomi.part.1 + tricomi.part.2
  
  result <-
    exp(
      r * log(alpha) + s * log(beta) + (s - 1) * log(d) + lgamma(r +
                                                                   x + 1) + log(tricomi) - lgamma(r) - (r + x + 1) * log(alpha + T.cal) - log(likelihood)
    )
  
  return(result)
}


##############################################
#end of function writing

#create an empty frame to append the calibration output

combinedcal <-
  data.frame (
    cust = character(),
    region = character(),
    T.cal = numeric(),
    t.x = numeric(),
    x = numeric(),
    expectedtranscustcal = integer(),
    probcustcal = double(),
    residcustcal = integer(),
    newpurchase = numeric()
    
  )


# get the unique groups for loop
uniqreg <- unique(my_data$region)

#Testing for one ( Remove the line for running for multiple groupps )
#please make sure that you have enough data points by group to run an optimization 
# Change the optimizer as per your needs 

#uniqreg <- "House Ads"

for (i in 1:length(uniqreg)) {
  print(uniqreg[i])
  dfori <-subset(my_data, my_data$region == uniqreg[i])
  dfori_c <-aggregate(dfori$date, by = list(dfori$cust, dfori$region),frequency )
  dfori_c$x <- NULL
  names(dfori_c) <-   c("cust", "region")
  df <-  dfori[c("cust","date")]

  ####################creation of calibration matrix #####
  
  y <- subset(df, df$date <= end.of.cal.period)
  min <- aggregate(y$date, by = list(y$cust), min)
  rownames(min) <- min$Group.1
  min$Group.1 <- NULL
  
  max <- aggregate(y$date, by = list(y$cust), max)
  rownames(max) <- max$Group.1
  max$Group.1 <- NULL
  
  n <- as.data.frame(table(y$cust))
  colnames(n)[1] <- "cust"
  colnames(n)[2] <- "x"
  
  custdata <- merge(min, max, by = "row.names")
 
  custdata$T.cal <-
    as.numeric((end.of.cal.period - custdata$x.x + 1) / 365.24)
  custdata$t.x <-
    as.numeric((custdata$x.y - custdata$x.x + 1) / 365.24)
  
  custdata$x.x <- NULL
  custdata$x.y <- NULL
  
  colnames(custdata)[1] <- "cust"
  custdata$cust <- as.character(custdata$cust)

  cal.cbs <- merge (custdata, n, by = "cust")
  min(cal.cbs$T.cal)

  cal1.cbs <- cal.cbs
  
  ######################################
  # getting the output for full calibration set
  caloptmix.cbs <- cal.cbs

  (params7 <-
      pnbd.EstimateParameters(caloptmix.cbs))
  
  (LL <- pnbd.cbs.LL(params7, caloptmix.cbs))
  print(params7)
  
  p.matrix <- c(params7, LL)
  print(p.matrix)
  
  
  # use final set of values
  (params2 <- p.matrix[1:4])
  print(params2)
  
  #########################################
  #Output for the full Calibration set
  caltot.cbs <-   cal.cbs
  rownames(caltot.cbs) <- caltot.cbs$cust
  caltot.cbs$cust <- NULL
  class(caltot.cbs)
  caltot.cbs$x <- as.numeric(caltot.cbs$x)
  caltotm.cbs <- as.matrix(caltot.cbs)
    probcustcal <-
    pnbd.PAlive(params2, caltotm.cbs[, "x"], caltotm.cbs[, "t.x"], caltotm.cbs[, "T.cal"])
  
  #conditional expected transaction in a given time period ( below is 1 years)
  
  expectedtranscustcal <-
    pnbd.ConditionalExpectedTransactions(params2, T.star = 1, caltotm.cbs[, "x"], caltotm.cbs[, "t.x"], caltotm.cbs[, "T.cal"])
  
  # change the discounting of 0.l for residual lifetime decay
  residcustcal <-
    pnbd.DERT(params2, (caltotm.cbs[, "x"] / 5), caltotm.cbs[, "t.x"], caltotm.cbs[, "T.cal"], 0.10)
  residcustcal.df <- as.data.frame(residcustcal)
  residcustcal.df$cust <- rownames(residcustcal.df)
  residcustcal.df$row.names <- NULL
  residcustcal.df$residcustcal <- 5 * (residcustcal.df$residcustcal)
 
  
  probcustcal.df <- as.data.frame(probcustcal)
  
  custtempcal <-
    merge(probcustcal.df, residcustcal.df, by = 0)
  rownames(custtempcal) = custtempcal$Row.names
  custtempcal$Row.names <- NULL
 
  predictholdcal.df <-
    as.data.frame(expectedtranscustcal)
  
  custtranscompcal <-
    merge(predictholdcal.df, custtempcal, by = 0)
  rownames(custtranscompcal) <-
    custtranscompcal$Row.names
  custtranscompcal$Row.names <- NULL
  
  
  finalcustcal <-
    merge(caltot.cbs, custtranscompcal, by = 0)
  
  finalcustcal$cust <- NULL

  finalcustcal$residcustcal <-
    ceiling(finalcustcal$residcustcal)
  finalcustcal$expectedtranscustcal <-
    ceiling(finalcustcal$expectedtranscustcal)
  head(finalcustcal)
  
  
  #number of purchases by a new individual
  newpurchase <- pnbd.Expectation(params2, 1)
  colnames(finalcustcal)[1] <- "cust"
  #colnames(dfori_c)[1] <- "cust"
  Outputcal <-  merge(dfori_c, finalcustcal, by = "cust")
  
  Outputcal$newpurchase <- round(newpurchase)
  head(Outputcal)
  
  combinedcal <- rbind.fill(combinedcal, Outputcal)
  
}

combinedcal$T.cal <- NULL
combinedcal$t.x <- NULL
combinedcal$x <- NULL
names(combinedcal) <-   c("cust", "region","expected transactions","probability of activity","residual transactions","expected purchase from a new member bygroup ")

write.csv (combinedcal,
           paste(
             "Output",end.of.cal.period,
             Sys.Date(),
             ".csv"
           ))


end.time <- Sys.time()
time.taken <- end.time - start.time

print(time.taken)

sink()

head(combinedcal) # calibration period output



