## y = b0 + b1x1 + b2L
get_adj_coeff <- function(y, x1, L=NA, v_ex1) {
  
  # aa = length(y) == length(x1)
  # if (!aa) 
  #   stop("Error: length of y must equal rows of x1.")
  # 
  # if(!is.na(L))
  # aa = length(y) == length(L)
  # if(!aa)
  #   stop('Error: length of y must equal of L.')
  
  y <- as.numeric(y)
  x1 <- as.numeric(x1)
  L <- as.integer(L)
  
  design_matrix <- cbind(rep(1, length(y)), x1, L)
  
  xtp_x <- matrix(nrow = 3, ncol = 3)
  for (i in seq(1, 3)) {
    for (j in seq(1, 3)) {
      xtp_x[i, j] <- mean(design_matrix[, i] * design_matrix[, j])
    }
  }
  
  xtp_x[2, 2] <- xtp_x[2, 2] - v_ex1
  
  xtp_y <- matrix(nrow = 3, ncol = 1)
  for (i in seq(1, 3)) {
    xtp_y[i, 1] <- mean(design_matrix[, i] * y)
  }
  
  if(is.na(L[1])) {
    xtp_x <- xtp_x[1:2, 1:2]
    xtp_y <- xtp_y[1:2, ]
  }
  
  b_hats <- MASS::ginv(xtp_x)%*%xtp_y
  
  return(c(b_hats))
  
  testing <- function() {
    n <- 10000
    L <- rbinom(n, 1, 0.5)
    A <- rnorm(n, mean = 5, sd = sqrt(4))
    B <- 0.0 + 2*A + L + rnorm(n)
    
    v_ea <- 0.5
    v_eb <- 0.5
    Ap <- A + rnorm(n, mean = 0, sd = sqrt(v_ea))
    Bp <- B + rnorm(n, mean = 0, sd = sqrt(v_eb))
    
    get_adj_coeff(y = Bp, x1 = Ap, L = L, v_ex1 = v_ea)
    lm(B ~ A + L)
    mean(A); mean(Ap)
    mean(A*L); mean(Ap*L)
    mean(B); mean(Bp)
    mean(A*B); mean(Ap*Bp)
    mean(B*L); mean(Bp*L)
    
    #'bivariate normal data generation
    z1 <- rnorm(500)
    z2 <- rnorm(500)
    x <- sqrt(4)*z1 + 5
    y <- sqrt(6)*(0.5*z1 + sqrt(1 - 0.5)*z2) + 3
    hist(x)
    hist(y)
    var(x)
    var(y)
    ggplot() + geom_point(aes(x = x, y = y))
    lm(y~x)
    x <- x + rnorm(500, sd = sqrt(0.4))
    get_adj_coeff(y = y, x1 = x, v_ex1 = 0.4)
    
  }
  
}

test1_p <- function(L, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed = NULL) {
  
  if(!is.null(rseed))set.seed(rseed)
  
  # cit test 1(L,T)
  #H0 : T = b0 + e0, H1: T = b0 + b1L1 + b2L2 + e1. (Remind: b0, b1 in H1 are diff from H0. bcz both are diff lin reg)
  v_h0 <- var(Tp) - v_eT
  v_h1 <- var(residuals(lm(Tp ~ L))) - v_eT
  
  q <- 1
  p <- 2
  f_obs <- ((v_h0 - v_h1)/v_h1)*(length(Tp) - p)/(p - q)
  p_val <- pf(f_obs, df1 = p - q, df2 = length(Tp) - p, lower.tail = FALSE)
  return(p_val)
  
  testing <- function() {
    source('/data/users/cs18s008/projects/extract_trio/code/simulation/data_generate_causal.R')
    # source('/data/users/cs18s008/projects/extract_trio/code/simulation/data_generate_confounding.R')
    parameters <- expand.grid(
      n = c(500),
      p = 0.5,
      r_za = sqrt(c(0.10, 0.50, 0.75)),
      r_ab = sqrt(c(0.10, 0.50, 0.75)),
      noisea = sqrt(seq(0, 1, by=0.2)),
      noiseb = sqrt(c(0.2)),
      nsim = 1:1
    )
    
    parameters$p_val <- sapply(1:nrow(parameters), FUN = function(i) {
      
      data <- make_system_yeast(n = parameters$n[i], p = parameters$p[i], r_ab = parameters$r_za[i], r_za = parameters$r_ab[i],
                                noisea = sqrt(parameters$noisea), noiseb = sqrt(parameters$noiseb[i]))
      test1_p(L = data$L, Gp = data$Ap, Tp = data$Bp, v_eG = parameters$noisea[i]^2, v_eT = parameters$noiseb[i]^2, bootstrap = 500, rseed = 5)
      
    })
    
    ggplot() + geom_histogram(aes(x = parameters$p_val)) 
    
    }
  
}

test2_p <- function(L, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed = NULL) {
  
  if(!is.null(rseed))set.seed(rseed)
  
  # #cit test2(G,T|L)
  #
  v_h0 <- var(residuals(lm(Tp ~ L))) - v_eT
  
  adj_coef <- get_adj_coeff(y = Tp, x1 = Gp, L = L, v_ex1 = v_eG)
  adj_resd <- Tp - (adj_coef[1] + adj_coef[2]*Gp + adj_coef[3]*L)
  v_h1 <- var(adj_resd) - adj_coef[2]^2*v_eG - v_eT
  
  q <- 2
  p <- 3
  f_obs <- ((v_h0 - v_h1)/v_h1)*(length(Tp) - p)/(p - q)
  p_val <- pf(f_obs, df1 = p - q, df2 = length(Tp) - p, lower.tail = FALSE)
  
  return(p_val)
  
  testing <- function() {
    # source('/data/users/cs18s008/projects/extract_trio/code/simulation/data_generate_causal.R')
    source('/data/users/cs18s008/projects/extract_trio/code/simulation/data_generate_confounding.R')
    parameters <- expand.grid(
      n = c(500),
      p = 0.5,
      r_za = sqrt(c(0.10, 0.50, 0.75)),
      r_ab = sqrt(c(0.10, 0.50, 0.75)),
      noisea = sqrt(seq(0, 1, by=0.2)),
      noiseb = sqrt(c(0.2)),
      nsim = 1:1
    )
    
    parameters$p_val <- sapply(1:nrow(parameters), FUN = function(i) {
      
      data <- make_system_yeast(n = parameters$n[i], p = parameters$p[i], r_ab = parameters$r_za[i], r_za = parameters$r_ab[i],
                                noisea = sqrt(parameters$noisea), noiseb = sqrt(parameters$noiseb[i]))
      test2_p(L = data$L, Gp = data$Ap, Tp = data$Bp, v_eG = parameters$noisea[i]^2, v_eT = parameters$noiseb[i]^2, bootstrap = 500, rseed = 5)
      
    })
    
    ggplot() + geom_histogram(aes(x = parameters$p_val)) 
    
  }
  
}




test3_p <- function(L, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed = NULL) {
  
  if(!is.null(rseed))set.seed(rseed)
  
  get_adj_cittest3_p <- function(L, Gp, Tp, v_eG, v_eT) {
    #cit test 3(L,G|T)
    #H0: Gp = b0 + b1Tp + (-)
    adj_coef <- get_adj_coeff(y = Gp, x1 = Tp,  v_ex1 = v_eT)
    adj_resd <- Gp - (adj_coef[1] + adj_coef[2]*Tp)
    v_h0 <- var(adj_resd) - (adj_coef[2]^2)*v_eT - v_eG
    
    #H1: Gp = b0 + b1Tp + b2L + (-b1*e_t + e_g + e1)
    adj_coef <- get_adj_coeff(y = Gp, x1 = Tp, L = L, v_ex1 = v_eT)
    adj_resd <- Gp - (adj_coef[1] + adj_coef[2]*Tp + adj_coef[3]*L)
    v_h1 <- var(adj_resd) - (adj_coef[2]^2)*v_eT - v_eG
    
    q <- 2
    p <- 3
    f_obs <- ((v_h0 - v_h1)/v_h1)*(length(Tp) - p)/(p - q)
    return(pf(f_obs, df1 = p - q, df2 = length(Tp) - p, lower.tail = FALSE))
    # return(f_obs)
  }
  
  #test (L, G|T)
  p_h1 <- get_adj_cittest3_p(L = L, Gp = Gp, Tp = Tp, v_eG = v_eG, v_eT = v_eT)
  n <- length(Gp)
  
  #steps to get Gp_t i'e G tilda prime
  adj_coef <- get_adj_coeff(y = Gp, x1 = Tp, L = L, v_ex1 = v_eT)
  Gp_t <- Gp - adj_coef[3]*L
  
  p_bootstrap <- rep(NA, bootstrap)
  for (i in 1:bootstrap) {
    indices <- sample(1:n, replace = T)
    # L1_b <- L1[indices]
    # L2_b <- L2[indices]
    L_b <- L[indices]
    Tp_b <- Tp[indices]
    Gp_t_b <- Gp_t[indices]
    
    p_bootstrap[i] <- get_adj_cittest3_p(L = L_b, Gp = Gp_t_b, Tp = Tp_b, v_eG = v_eG, v_eT = v_eT)
  }
  
  # p_uniform <- sapply(p_bootstrap, FUN = function(p_val) {return(sum(p_bootstrap < p_val)/length(p_bootstrap))})
  # plot(ggplot() + geom_histogram(aes(p_uniform)))
  
  p_val <- sum(p_bootstrap < p_h1)/length(p_bootstrap)
  
  return(p_val)
  
  testing <- function() {
    # source('/data/users/cs18s008/projects/extract_trio/code/simulation/data_generate_causal.R')
    source('/data/users/cs18s008/projects/extract_trio/code/simulation/data_generate_confounding.R')
    parameters <- expand.grid(
      n = c(500),
      p = 0.5,
      r_za = sqrt(c(0.10, 0.50, 0.75)),
      r_ab = sqrt(c(0.10, 0.50, 0.75)),
      noisea = sqrt(seq(0, 1, by=0.2)),
      noiseb = sqrt(c(0.2)),
      nsim = 1:1
    )
    
    parameters$p_val <- sapply(1:nrow(parameters), FUN = function(i) {
      
      data <- make_system_yeast(n = parameters$n[i], p = parameters$p[i], r_ab = parameters$r_za[i], r_za = parameters$r_ab[i],
                                noisea = sqrt(parameters$noisea), noiseb = sqrt(parameters$noiseb[i]))
      test3_p(L = data$L, Gp = data$Ap, Tp = data$Bp, v_eG = parameters$noisea[i]^2, v_eT = parameters$noiseb[i]^2, bootstrap = 500, rseed = 5)
      
    })
    
    ggplot() + geom_histogram(aes(x = parameters$p_val)) + scale_x_log10()
  }
  
}


test4_p <- function(L, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed = NULL) {
  
  if(!is.null(rseed))set.seed(rseed)
  
  get_adj_cittest4_f <- function(L, Gp, Tp, v_eG, v_eT) {
    #---------------------------------------------------------------------------------------------
    #cit test4 (L indp T|G). Note for independence test, v_h0: complex model. v_h1 represents simpler model. therefore v_h1 > v_h0
    #procedure:
    #s1- estimate f_obs by estimating v_h0 and v_h1.
    #s2- semiparametric (bootstrap + moment matching param est of f*. Then z transformation)
    
    #H0: 
    adj_coef <- get_adj_coeff(y = Tp, x1 = Gp, L = L,  v_ex1 = v_eG)
    adj_resd <- Tp - (adj_coef[1] + adj_coef[2]*Gp + adj_coef[3]*L)
    v_h0 <- var(adj_resd) - adj_coef[2]^2*v_eG - v_eT
    
    #H1:
    adj_coef <- get_adj_coeff(y = Tp, x1 = Gp,  v_ex1 = v_eG)
    adj_res <- Tp - (adj_coef[1] + adj_coef[2]*Gp)
    v_h1 <- var(adj_res) - adj_coef[2]^2*v_eG - v_eT
    
    q <- 2
    p <- 3
    
    f <-  ((v_h1 - v_h0)/v_h0)*(length(Tp) - p)/(p - q)
    
    return(f)
    
  }
  
  f_obs <- get_adj_cittest4_f(L, Gp, Tp, v_eG, v_eT)
  
  bootstrap <- resampl
  #steps to generate G*. In code G* is represented as G_.
  lmG_ <- lm(Gp ~ L)
  # alt = summary(lmG_)$coefficients["(Intercept)",1]
  # blt1 = summary(lmG_)$coefficients["L1",1]
  # blt2 = summary(lmG_)$coefficients["L2",1]
  
  alt <- coefficients(lmG_)['(Intercept)']
  alt <- ifelse(is.na(alt), 0, alt)
  blt1 <- coefficients(lmG_)['L']
  blt1 <- ifelse(is.na(blt1), 0, blt1)

  resd_G_ <- residuals(lmG_)
  f_bs <- array(0, length(bootstrap))
  for (i in c(1:bootstrap)) {
    # since here sample(replace=FALSE), so it acts as permutation
    G_ <- alt + blt1*L + sample(resd_G_)
    f_bs[i] <- get_adj_cittest4_f(L = L, Gp = G_, Tp = Tp, v_eG =  v_eG, v_eT = v_eT)
  }
  
  
  #####F Method
  f_bs = f_bs[!is.na(f_bs)]
  q <- 2
  p <- 3
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
  
  testing <- function() {
    rm(list = ls())
    source('/data/users/cs18s008/projects/extract_trio/code/simulation/data_generate_confounding.R')
    source('/data/users/cs18s008/projects/extract_trio/code/simulation/data_generate_causal.R')
    parameters <- expand.grid(
      n = c(500),
      p = 0.5,
      r_za = sqrt(c(0.10, 0.50, 0.75)),
      r_ab = sqrt(c(0.10, 0.50, 0.75)),
      noisea = sqrt(seq(0, 1, by=0.2)),
      noiseb = sqrt(c(0.2)),
      nsim = 1:1
    )
    
    parameters$p_val <- sapply(1:nrow(parameters), FUN = function(i) {
      
      data <- make_system_yeast(n = parameters$n[i], p = parameters$p[i], r_ab = parameters$r_za[i], r_za = parameters$r_ab[i],
                                noisea = sqrt(parameters$noisea), noiseb = sqrt(parameters$noiseb[i]))
      test4_p(L = data$L, Gp = data$Ap, Tp = data$Bp, v_eG = parameters$noisea[i]^2, v_eT = parameters$noiseb[i]^2, bootstrap = 500, resampl= 100, rseed = 5)
      
    })
    
    ggplot() + geom_histogram(aes(parameters$p_val)) + scale_x_log10()
  }
  
}



##adj_cit_pvals using f-stats
get_adj_cit_pvals <- function(L, Gp, Tp, v_eG, v_eT, bootstrap = 300, resampl=50, rseed=NULL) {
  
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
    
    p_vals[2] <- test1_p(L, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed)
    
    p_vals[3] <- test2_p(L, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed)
    
    p_vals[4] <- test3_p(L, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed)
    
    p_vals[5] = test4_p(L, Gp, Tp, v_eG, v_eT, bootstrap, resampl, rseed)
    
    
    p_vals[1] <- max(p_vals[2:5])
    names(p_vals) = c("adj_p_cit", "adj_p_TassocL", "adj_p_TassocGgvnL", "adj_p_GassocLgvnT", 
                      "adj_p_LindTgvnG")
    
    # return(list(p_vals, p_h1, p_s_vec)[1])
    return(p_vals)
    
    testing <- function() {
      source('/data/users/cs18s008/projects/extract_trio/code/simulation/data_generate_causal.R')
      # source('/data/users/cs18s008/projects/extract_trio/code/simulation/data_generate_confounding.R')
      parameters <- expand.grid(
        n = c(500),
        p = 0.5,
        r_za = sqrt(c(0.10, 0.50, 0.75)),
        r_ab = sqrt(c(0.10, 0.50, 0.75)),
        noisea = sqrt(seq(0, 1, by=0.2)),
        noiseb = sqrt(c(0.2)),
        nsim = 1:1
      )
      
      result <- mapply(1:nrow(parameters), FUN = function(i) {
        
        data <- make_system_yeast(n = parameters$n[i], p = parameters$p[i], r_ab = parameters$r_za[i], r_za = parameters$r_ab[i],
                                  noisea = sqrt(parameters$noisea), noiseb = sqrt(parameters$noiseb[i]))
        get_adj_cit_pvals(L = data$L, Gp = data$Ap, Tp = data$Bp, v_eG = parameters$noisea[i]^2, v_eT = parameters$noiseb[i]^2, bootstrap = 500, resampl= 100, rseed = 5)
        
      })
      
      # ggplot() + geom_histogram(aes(x = parameters$p_val)) 
      
    }
    
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
# test4_p(L = data$L, Gp = data$Ap, Tp = data$Bp, v_eG = parameters$noisea[i]^2, v_eT = parameters$noiseb[i]^2, bootstrap = 500, resampl= 500, rseed = 5)



