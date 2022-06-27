library(mvtnorm)
# read data and initialization
# df = read.table("data.txt", header = TRUE)
idx0 = df[, 1] == 0
X = cbind(1, as.matrix(df[, 2:3]))
Y = c(df[, 1])
sd = 1
beta0 = c(0, 0, 0)
sigmaB = 1000* diag(c(1, 1, 1))
gamma0 = c(0, 0, 0)
sigmaG = 1000*diag(c(1, 1, 1))

# density fucntion of GP
GP_f = function(y, x, beta, alpha = 0) {
  mu = x %*% beta
  return((mu/(1 + alpha*mu))^y * (1 + alpha*y)^(y-1) * exp(-mu*(1 + alpha*y) / (1 + alpha*mu)) / factorial(y))
}

# density function of Possion
P_f = function(y, x, beta) {
  mu = x %*% beta
  return(dpois(y, mu))
}

# log prior density of beta
logp_beta = function(x) {
  return(dmvnorm(x, beta0, sigmaB, log = TRUE))
}

# log prior density of gamma
logp_gamma = function(x) {
  return(dmvnorm(x, gamma0, sigmaG, log = TRUE))
}

# log likelihood

log_lik = function(alpha, beta, gamma) {
  rho = exp(X %*% gamma) / (1 + exp(X %*% gamma))
  ret = sum(log(rho[idx0] + (1 - rho[idx0]) * GP_f(Y[idx0], X[idx0, ], beta, alpha))) + sum(log(1 - rho[!idx0]) + log(GP_f(Y[!idx0], X[!idx0, ], beta, alpha)))
  return(ret)
}

# probability of one sample
prob = function(y, x, alpha, beta, gamma) {
  rho = exp(x %*% gamma) / (1 + exp(x %*% gamma))
  if (y == 0) {
    return(rho + (1 - rho) * GP_f(y, x, beta, alpha))
  }
  else {
    return((1 - rho) * GP_f(y, x, beta, alpha))
  }
}


#full conditional in GP regression
logfull_alpha = function(alpha, beta, gamma) {
  rho = exp(X %*% gamma) / (1 + exp(X %*% gamma))
  ret = sum(log(rho[idx0] + (1 - rho[idx0]) * GP_f(Y[idx0], X[idx0, ], beta, alpha))) + sum(log(GP_f(Y[!idx0], X[!idx0, ], beta, alpha)))
  return(ret)
}

logfull_beta = function(alpha, beta, gamma) {
  ret = logfull_alpha(alpha, beta, gamma) - 0.5*t(beta - beta0) %*% solve(sigmaB) %*% (beta - beta0)
  return(ret)
}

logfull_gamma = function(alpha, beta, gamma) {
  rho = exp(X %*% gamma) / (1 + exp(X %*% gamma))
  ret = sum(log(rho[idx0] + (1 - rho[idx0]) * GP_f(Y[idx0], X[idx0, ], beta, alpha))) + sum(log(1 - rho[!idx0]))
  ret = ret - 0.5*t(gamma - gamma0) %*% solve(sigmaG) %*% (gamma - gamma0)
  return(ret)
}


# fit GP regression

check = function(alpha, beta, gamma) {
  mu = X %*% beta
  if (any((mu/(1 + alpha*mu)) < 0) | any((1 + alpha*Y) < 0)) {
    return(FALSE)
  }
  return(TRUE)
}

update_alpha = function(alpha, beta, gamma) {
  alpha_star = rnorm(1, alpha, sd)
  while (TRUE) {
    if (check(alpha_star, beta, gamma)){
      break
    }
    alpha_star = rnorm(1, alpha, sd)
  }
  log_r = logfull_alpha(alpha_star, beta, gamma) - logfull_alpha(alpha, beta, gamma)
  if (log(runif(1)) < log_r) {
    alpha = alpha_star
  }
  return(alpha)
}

update_beta = function(alpha, beta, gamma) {
  for (i in 1:3) {
    new_beta = beta
    new_beta[i] = rnorm(1, beta[i], sd)
    while (TRUE) {
      if (check(alpha, new_beta, gamma)){
        break
      }
      new_beta[i] = rnorm(1, beta[i], sd)
    }
    
    log_r = logfull_beta(alpha, new_beta, gamma) + dnorm(new_beta[i], 0, sqrt(1000), log = TRUE) - logfull_beta(alpha, beta, gamma) - dnorm(beta[i], 0, sqrt(1000), log = TRUE)
    if (log(runif(1)) < log_r) {
      beta[i] = new_beta[i]
    }
  }
  return(beta)
}

update_gamma = function(alpha, beta, gamma) {
  for (i in 1:3) {
    gammai_star = rnorm(1, gamma[i], sd)
    new_gamma = gamma
    new_gamma[i] = gammai_star
    while (TRUE) {
      if (check(alpha, beta, new_gamma)){
        break
      }
      new_gamma[i] = rnorm(1, gamma[i], sd)
    }
    
    log_r = logfull_gamma(alpha, beta, new_gamma) + dnorm(new_gamma[i], 0, sqrt(1000), log = TRUE) - logfull_gamma(alpha, beta, gamma) - dnorm(gamma[i], 0, sqrt(1000), log = TRUE)
    if (log(runif(1)) < log_r) {
      gamma[i] = new_gamma[i]
    }
  }
  return(gamma)
}


fit_GP = function(N){
  ret = NULL
  alpha = runif(1)
  beta = c(rmvnorm(1, beta0, sigmaB/1000))
  gamma = c(rmvnorm(1, gamma0, sigmaG/1000))
  
  while (TRUE) {
    if (check(alpha, beta, gamma)) {
      break
    }
    alpha = runif(1)
    beta = c(rmvnorm(1, beta0, sigmaB/1000))
    gamma = c(rmvnorm(1, gamma0, sigmaG/1000))
  }
  for (i in 1:N) {
    alpha = update_alpha(alpha, beta, gamma)
    beta = update_beta(alpha, beta, gamma)
    gamma = update_gamma(alpha, beta, gamma)
    if (i %% 20 == 0) {
      ret = rbind(ret, c(alpha, beta, gamma))
    }
  }
  return(ret)
}


# test 
ret = fit_GP(50000)
ret = ret[501:2500, ]
a = colMeans(ret)
log_lik(a[1], a[2:4], a[5:7])
log_lik(0.37, c(1.8,-0.53, 0.23), c(0.38, -0.4, -0.52))


# fit Possion regression
fit_P = function(N) {
  ret = NULL
  beta = c(rmvnorm(1, beta0, sigmaB/1000))
  gamma = c(rmvnorm(1, gamma0, sigmaG/1000))
  
  while (TRUE) {
    if (check(0, beta, gamma)) {
      break
    }
    beta = c(rmvnorm(1, beta0, sigmaB/1000))
    gamma = c(rmvnorm(1, gamma0, sigmaG/1000))
  }
  for (i in 1:N) {
    beta = update_beta(0, beta, gamma)
    gamma = update_gamma(0, beta, gamma)
    if (i %% 20 == 0) {
      ret = rbind(ret, c(beta, gamma))
    }
  }
  return(ret)
  
}
#ret = fit_P(1000)

# DIC of GP regression
DIC_GP = function(samples) {
  mean_sample = colMeans(samples)
  D_hat = -2*log_lik(mean_sample[1], mean_sample[2:4], mean_sample[5:7])
  res = rep(0, nrow(samples)) 
  for (i in 1: nrow(samples)) {
    res[i] = -2*log_lik(samples[i, 1], samples[i, 2:4], samples[i, 5:7])
  }
  D_bar = mean(res)
  return(0.5*var(res) + mean((res)))
}

#DIC_GP(samples)

# DIC of Possion regression
DIC_P = function(samples) {
  mean_sample = colMeans(samples)
  D_hat = -2*log_lik(0, mean_sample[1:3], mean_sample[4:6])
  res = rep(0, nrow(samples)) 
  for (i in 1: nrow(samples)) {
    res[i] = -2*log_lik(0, samples[i, 1:3], samples[i, 4:6])
  }
  D_bar = mean(res)
  return(0.5*var(res) + mean((res)))
}

# B stats of GP regression

B_GP = function(samples) {
  cpo = rep(0, nrow(df))
  n = nrow(samples)
  for (i in 1: nrow(df)) {
    curr = rep(0, n)
    for (j in 1:n) {
      curr[j] = 1/ prob(Y[i], X[i, ], samples[j, 1], samples[j, 2:4], samples[j, 5:7])
    }
    cpo[i] = (mean(curr)) ^ (-1)
  }
  return(mean(log(cpo)))
}
B_GP(samples)


# B stats of Possion regression
B_P = function(samples) {  cpo = rep(0, nrow(df))
  n = nrow(samples)
  for (i in 1: nrow(df)) {
    curr = rep(0, n)
    for (j in 1:n) {
      curr[j] = 1/ prob(Y[i], X[i, ], 0, samples[j, 1:3], samples[j, 4:6])
    }
    cpo[i] = (mean(curr)) ^ (-1)
  }
  return(mean(log(cpo)))
}

# K stats of GP regression
K_GP = function(samples, i) {
  n = nrow(samples)
  curr = rep(0, n)
  p = rep(0, n)
  for (j in 1:n) {
    curr[j] = 1/ prob(Y[i], X[i, ], samples[j, 1], samples[j, 2:4], samples[j, 5:7])
    p[j] = log(prob(Y[i], X[i, ], samples[j, 1], samples[j, 2:4], samples[j, 5:7]))
  }
  ret = log(mean(curr)) + mean(p)
  return(ret)
}

# K stats of Possion regression
K_P = function(samples, i) {
  n = nrow(samples)
  curr = rep(0, n)
  p = rep(0, n)
  for (j in 1:n) {
    curr[j] = 1/ prob(Y[i], X[i, ], 0, samples[j, 1:3], samples[j, 4:6])
    p[j] = log(prob(Y[i], X[i, ], 0, samples[j, 1:3], samples[j, 4:6]))
  }
  ret = log(mean(curr)) + mean(p)
  return(ret)
}





