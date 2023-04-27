# Since I detected when adjcit1.3 4th test run alone, it's working as expected but when run along with
# other tests 1,2,3 it behaves unexpected. So to resolve the issue, this is the updated version of 1.3 to 1.4.
# It's logic is completely same as 1.3, except each test is defined into separate functions.
# library(MASS)
##Least change null model try
##Ideally:input(L,A,B), find (L, At,B) where At = A - coef*L. Now for H0, bootstrap from (L, At, B)
##Reality:input(L,Ap, Bp), find(L, Ap_t, Bp) where Ap_t = At - adjCoef*L. Now H0 , bootstrap from (L, Ap_t, Bp)

## y = b0 + b1x1 + b2L1 + b3L2
# continuos  and error prone :x1, y. But error in y don't affect parameters estimations
# L1, L2: discrete 
get_adj_coeff <- function(y, x1, L1=NA, L2=NA,  v_ex1) {
  
  # #just to make sure not all L1 and L2 are all 0 else singular error happens
  # if(!missing(L1) & all(L1 == 0)) {
  #   #set randomly any L1 to 1 whose L2 is 0.
  #   L1[match(0, L2)] <- 1
  # }
  # 
  # if(!missing(L2) & all(L2 == 0)) {
  #   #set randomly any L2 to 1 whose L1 is 0
  #   L2[match(0, L1)] <- 1
  # }
  
  data <- cbind( rep(1, length(x1)), x1, L1, L2, y)
  mat_x <- matrix(nrow = 4, ncol = 4)
  
  for (i in c(1:nrow(mat_x))) {
    for (j in c(i:ncol(mat_x))) {
      mat_x[i,j] <- mean(data[,i]*data[,j])
      mat_x[j,i] <- mat_x[i,j]
    }
  }
  #only need to adjust mat_x[2,2]
  mat_x[2,2] <- mat_x[2,2] - v_ex1
  
  
  
  mat_y <- matrix(nrow = 4, ncol = 1)
  for (i in c(1:nrow(mat_y))) {
    mat_y[i,1] <- mean(data[,i]*data[,5])
  }
  
  if (missing(L2)) {
    mat_x <- mat_x[1:3, 1:3]
    mat_y <- mat_y[1:3]
  }
  
  if(missing(L1)) {
    mat_x <- mat_x[1:2, 1:2]
    mat_y <- mat_y[1:2]
  }
  
  coef_h <- MASS::ginv(mat_x) %*% mat_y
  
  # print(mat_x)
  # print("------------")
  # print(mat_y)
  return(coef_h)
  
}

get_adj_cittest4_f <- function(L1, L2, Gp, Tp, v_eG, v_eT) {
  #---------------------------------------------------------------------------------------------
  #cit test4 (L indp T|G). Note for independence test, v_h0: complex model. v_h1 represents simpler model. therefore v_h1 > v_h0
  #procedure:
  #s1- estimate f_obs by estimating v_h0 and v_h1.
  #s2- semiparametric (bootstrap + moment matching param est of f*. Then z transformation)
  
  #H0: 
  adj_coef <- get_adj_coeff(y = Tp, x1 = Gp, L1 = L1, L2 = L2,  v_ex1 = v_eG)
  adj_resd <- Tp - (adj_coef[1] + adj_coef[2]*Gp + adj_coef[3]*L1 + adj_coef[4]*L2)
  v_h0 <- var(adj_resd) - adj_coef[2]^2*v_eG - v_eT
  
  #H1:
  adj_coef <- get_adj_coeff(y = Tp, x1 = Gp,  v_ex1 = v_eG)
  adj_res <- Tp - (adj_coef[1] + adj_coef[2]*Gp)
  v_h1 <- var(adj_res) - adj_coef[2]^2*v_eG - v_eT
  
  q <- 2
  p <- 4
  
  f <-  ((v_h1 - v_h0)/v_h0)*(length(Tp) - p)/(p - q)
  
  return(f)
  
}

get_adj_cittest3_p <- function(L1, L2, Gp, Tp, v_eG, v_eT) {
  #cit test 3(L,G|T)
  #H0: Gp = b0 + b1Tp + (-)
  adj_coef <- get_adj_coeff(y = Gp, x1 = Tp,  v_ex1 = v_eT)
  adj_resd <- Gp - (adj_coef[1] + adj_coef[2]*Tp)
  v_h0 <- var(adj_resd) - adj_coef[2]^2*v_eT - v_eG
  
  #H1: Gp = b0 + b1Tp + b2L1 + b3L2 + (-b1*e_t + e_g + e1)
  adj_coef <- get_adj_coeff(y = Gp, x1 = Tp, L1 = L1, L2 = L2,  v_ex1 = v_eT)
  adj_resd <- Gp - (adj_coef[1] + adj_coef[2]*Tp + adj_coef[3]*L1 + adj_coef[4]*L2)
  v_h1 <- var(adj_resd) - adj_coef[2]^2*v_eT - v_eG
  
  q <- 2
  p <- 4
  f_obs <- ((v_h0 - v_h1)/v_h1)*(length(Tp) - p)/(p - q)
  return(pf(f_obs, df1 = p - q, df2 = length(Tp) - p, lower.tail = FALSE))
}

test1_p <- function(L1, L2, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed) {
  
  if(!is.null(rseed))set.seed(rseed)
  
  # cit test 1(L,T)
  #H0 : T = b0 + e0, H1: T = b0 + b1L1 + b2L2 + e1. (Remind: b0, b1 in H1 are diff from H0. bcz both are diff lin reg)
  v_h0 <- var(Tp) - v_eT
  v_h1 <- var(residuals(lm(Tp ~ L1+L2))) - v_eT

  q <- 1
  p <- 3
  f_obs <- ((v_h0 - v_h1)/v_h1)*(length(Tp) - p)/(p - q)
  p_val <- pf(f_obs, df1 = p - q, df2 = length(Tp) - p, lower.tail = FALSE)
  return(p_val)

}

test2_p <- function(L1, L2, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed) {
  
  if(!is.null(rseed))set.seed(rseed)
  
  # #cit test2(G,T|L)
  #
  v_h0 <- var(residuals(lm(Tp ~ L1 + L2))) - v_eT

  adj_coef <- get_adj_coeff(y = Tp, x1 = Gp, L1 = L1, L2 = L2, v_ex1 = v_eG)
  adj_resd <- Tp - (adj_coef[1] + adj_coef[2]*Gp + adj_coef[3]*L1 + adj_coef[4]*L2)
  v_h1 <- var(adj_resd) - adj_coef[2]^2*v_eG - v_eT

  q <- 3
  p <- 4
  f_obs <- ((v_h0 - v_h1)/v_h1)*(length(Tp) - p)/(p - q)
  p_val <- pf(f_obs, df1 = p - q, df2 = length(Tp) - p, lower.tail = FALSE)
  
  return(p_val)

}

test3_p <- function(L1, L2, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed) {
  
  if(!is.null(rseed))set.seed(rseed)
  
  #test (L, G|T)
  p_h1 <- get_adj_cittest3_p(L1 = L1, L2 = L2, Gp = Gp, Tp = Tp, v_eG = v_eG, v_eT = v_eT)
  n <- length(Gp)

  #steps to get Gp_t i'e G tilda prime
  adj_coef <- get_adj_coeff(y = Gp, x1 = Tp, L1 = L1, L2 = L2, v_ex1 = v_eT)
  Gp_t <- Gp - adj_coef[3]*L1 - adj_coef[4]*L2

  p_bootstrap <- rep(NA, bootstrap)
  for (i in 1:bootstrap) {
    indices <- sample(1:n, replace = T)
    L1_b <- L1[indices]
    L2_b <- L2[indices]
    Tp_b <- Tp[indices]
    Gp_t_b <- Gp_t[indices]

    p_bootstrap[i] <- get_adj_cittest3_p(L1 = L1_b, L2 = L2_b, Gp = Gp_t_b, Tp = Tp_b, v_eG = v_eG, v_eT = v_eT)
  }


  p_val <- sum(p_bootstrap < p_h1)/length(p_bootstrap)

  return(p_val)
  
}

test4_p <- function(L1, L2, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed) {
  
  if(!is.null(rseed))set.seed(rseed)
  
  f_obs <- get_adj_cittest4_f(L1, L2, Gp, Tp, v_eG, v_eT)
  
  bootstrap <- resampl
  #steps to generate G*. In code G* is represented as G_.
  lmG_ <- lm(Gp ~ L1+L2)
  # alt = summary(lmG_)$coefficients["(Intercept)",1]
  # blt1 = summary(lmG_)$coefficients["L1",1]
  # blt2 = summary(lmG_)$coefficients["L2",1]
  
  alt <- coefficients(lmG_)['(Intercept)']
  alt <- ifelse(is.na(alt), 0, alt)
  blt1 <- coefficients(lmG_)['L1']
  blt1 <- ifelse(is.na(blt1), 0, blt1)
  blt2 <- coefficients(lmG_)['L2']
  blt2 <- ifelse(is.na(blt2), 0, blt2)
  
  resd_G_ <- residuals(lmG_)
  f_bs <- array(0, length(bootstrap))
  for (i in c(1:bootstrap)) {
    G_ <- alt + blt1*L1 + blt2*L2 + sample(resd_G_)
    f_bs[i] <- get_adj_cittest4_f(L1 = L1, L2 = L2, Gp = G_, Tp = Tp, v_eG =  v_eG, v_eT = v_eT)
  }
  
  
  #####F Method
  f_bs = f_bs[!is.na(f_bs)]
  q <- 2
  p <- 4
  df1 = p - q
  df2 = length(Gp) - p
  fncp = mean(f_bs,na.rm=TRUE)*(df1/df2)*(df2-df1)-df1
  if(fncp < 0) fncp = 0
  
  ######### Transform F to normal
  npvals = pf(f_bs,df1,df2,ncp=fncp,lower.tail=TRUE)
  nfvecr = qnorm(npvals)
  
  npf = pf(f_obs,df1,df2,ncp=fncp,lower.tail=TRUE) #Transform observed F
  zf = qnorm(npf)
  p_val = pnorm(zf,mean=mean(nfvecr),sd=sd(nfvecr))
  
  return(p_val)
  
}



##adj_cit_pvals using f-stats
get_adj_cit_pvals <- function(L1, L2, Gp, Tp, v_eG, v_eT, bootstrap = 300, resampl=50, rseed=NULL) {
  
  out <- tryCatch({
    if(!is.null(rseed))set.seed(rseed)
    
    # #just to make sure not all L1 and L2 are all 0 else singular error happens
    # if(!missing(L1) & all(L1 == 0)) {
    #   #set randomly any L1 to 1 whose L2 is 0.
    #   L1[match(0, L2)] <- 1
    # }
    # 
    # if(!missing(L2) & all(L2 == 0)) {
    #   #set randomly any L2 to 1 whose L1 is 0
    #   L2[match(0, L1)] <- 1
    # }
    
    
    p_vals <- rep(NA, 5)
    
    p_vals[2] <- test1_p(L1, L2, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed)
  
    p_vals[3] <- test2_p(L1, L2, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed)
  
    p_vals[4] <- test3_p(L1, L2, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed)
  
    p_vals[5] = test4_p(L1, L2, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed)
    
    
    p_vals[1] <- max(p_vals[2:5])
    names(p_vals) = c("adj_p_cit", "adj_p_TassocL", "adj_p_TassocGgvnL", "adj_p_GassocLgvnT", 
                      "adj_p_LindTgvnG")
    
    # return(list(p_vals, p_h1, p_s_vec)[1])
    return(p_vals)
    
  }, error = function(cond) {
    p_vals <- rep(NA, 5)
    return(p_vals)
  })
  
  return(out)
  
}

get_cit_direction <- function(p1, p2, thresh)
{
  res <- array(0, length(p1))
  # 1:causal, 2: reactive, 3: indp, -5 No call
  res[p1 < thresh & p2 > thresh] <- 1
  res[p1 > thresh & p2 < thresh] <- 2
  res[p1 > thresh & p2 > thresh] <- 3
  res[p1 <= thresh & p2 <= thresh] <- -5
  return(res)
}

