library(reshape2)
library(ggplot2)
library(plyr)
library(cit)

#'
#' Calculate fitted values directly
#'
#' @param y Response variable
#' @param  x Predictor
#' @export
#' @return Array
getFittedVals <- function(y, x)
{
  n <- length(x)
  bhat <- cov(y, x) / var(x)
  ahat <- mean(y) - bhat * mean(x)
  fitted <- ahat + x * bhat
  return(fitted)
}

#' Get pval directly
#' @param y Response variable
#' @param x Predictor
#' @export
#' @return p value
getPval <- function(y, x)
{
  n <- length(x)
  fitted <- getFittedVals(y, x)
  residuals <- y - fitted
  SSR <- sum((residuals - mean(residuals))^2)
  SSF <- sum((fitted - mean(fitted))^2)
  Fval <- (SSF) / (SSR/(n-2))
  pval <- pf(Fval, 1, n-2, lowe=F)
  return(pval)
}

makeProxy <- function(x, noise, bias)
{
  y <- x*bias + rnorm(length(x), 0, sd=noise)
  return(y)
}

make_geno <- function(n, p)
{
  stopifnot(all(p >= 0 & p <= 1))
  geno <- drop(sapply(p, function(x) rbinom(n, 2, x)))
  return(geno)
}


#' Create a phenotype using arbitrary number of known causal inputs
#'
#' For a set of \code{x} variables and effect sizes for each variable (\code{b})
#' \code{y} is constructed such that
#' \code{y = Xb + e}
#' Given that the variance explained in \code{y} by \code{X}, is
#' \code{r^2 = sum(b * var(x) / (sqrt(x)*sqrt(y)))}
#' we can model \code{e ~ N(0, 1 - r^2)}
#'
#' @param effs Array of beta values for each input. Leave the vx and vy values to default (=1) to allow effs to be equal to the correlation between y and each x
#' @param indep Matrix of independent variables corresponding to effs
#' @param vy=1 The output variance of y
#' @param vx=rep(1, length(effs)) The desired scaled variance of x
#'
#' @export
#' @return Numeric array, simulated phenotype
#'
#' @examples \dontrun{
#' g1 <- make_geno(1000, 0.5)
#' g2 <- make_geno(1000, 0.3)
#' x1 <- rnorm(1000)
#' x2 <- rnorm(1000)
#' y <- make_phen(effs=c(0.2, 0.1, 0.15, 0.4), cbind(g1, g2, x1, x2))
#'
#'}
#'rahul added this line: https://www.probabilitycourse.com/chapter5/5_3_2_bivariate_normal_dist.php generating bivariate normal, given individual mean, var, and cor
make_phen <- function(effs, indep, vy=1, vx=rep(1, length(effs)))
{
  if(is.null(dim(indep))) indep <- cbind(indep)
  stopifnot(ncol(indep) == length(effs))
  stopifnot(length(vx) == length(effs))
  cors <- effs * vx / sqrt(vx) / sqrt(vy)
  stopifnot(sum(cors) <= 1)
  cors <- c(cors, sqrt(1-sum(cors^2)))
  indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))
  y <- drop(scale(rowSums(indep)) * sqrt(vy))
  return(y)
}


make_system <- function(n, p, r_ab, r_za, noisea, noiseb)
{
  L <- make_geno(n, p)
  A <- make_phen(r_za, L)
  B <- make_phen(r_ab, A)
  Ap <- makeProxy(A, noisea, 1)
  Bp <- makeProxy(B, noiseb, 1)
  L1 <- rep(0, length(L))
  L1[L == 1] <- 1
  L2 <- rep(0, length(L))
  L2[L == 2] <- 1
  return(data.frame(L, A, B, Ap, Bp, L1, L2))
}

# parameters <- expand.grid(
#   n = c(1000),
#   p = 0.5,
#   r_za = c(sqrt(0.6)),
#   r_ab = sqrt(c(0.6)),
#   noisea = sqrt(0.8),
#   noiseb = c(sqrt(0.4)),
#   nsim = 1
# )
# 
# i <- 1
# dat <- with(parameters[i,], make_system(n, p, r_ab, r_za, noisea, noiseb))
# 
# # n <- 10000
# # L <- rbinom(n,2,.5)
# b0 <- 1
# b1 <- 1
# b2 <- 2*b1
# 
# b3 <- 1
# b4 <- 1
# 
# v_eA <- 0.2
# v_eB <- 0.0
# 
# L1 <- rep(0, length(L))
# L1[L == 1] <- 1
# L2 <- rep(0, length(L))
# L2[L == 2] <- 1
# 
# 
# A <- b0 + b1*L1 + b2*L2 + rnorm(n)
# Ap <- A + rnorm(n, sd=sqrt(v_eA))
# B <- b3 + b4*A + b1*L1 + b2*L2 + rnorm(n)
# Bp <- B + rnorm(n, sd = sqrt(v_eB))
# 
# dat <- data.frame(L, L1, L2, A, B, Ap, Bp)




# data_generate_causal_for_yeast ------------------------------------------


make_system_yeast <- function(n, p, r_ab, r_za, noisea, noiseb) {
  make_geno <- function(n, p)
  {
    stopifnot(all(p >= 0 & p <= 1))
    geno <- drop(sapply(p, function(x) rbinom(n, 1, x)))
    return(geno)
  }
  
  L <- make_geno(n, p)
  A <- make_phen(r_za, L)
  B <- make_phen(r_ab, A)
  Ap <- makeProxy(A, noisea, 1)
  Bp <- makeProxy(B, noiseb, 1)
  return(data.frame(L, A, B, Ap, Bp))
  
}
