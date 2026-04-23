
# Discrete Wrapped Stable PMF 

wrapped_stable_pmf <- function( rho, alpha, mu) {
  r_values <- 0:(m-1)
  pmf <- sapply(r_values, function(r) {
    sum_term <- sum(sapply(1:infinity, function(k) {
      rho^(k^alpha) * cos(k * (2 * pi * r / m - mu))
    }))
    pw <- (1 + 2 * sum_term) / m
    return(pw)
  })
  # Normalize due to finite P
  pmf <- pmf / sum(pmf)
  return(list(r = r_values, prob = pmf))
}

# Set up layout for 3x3 grid of PMF plots
par(mfrow = c(3, 3), mar = c(3.5, 3.5, 1))

param_list <- list(
  list(rho = 0.2, alpha = 0.3,  mu = 0.7),
  list(rho = 0.8, alpha = 0.5,  mu = 1.1),
  list(rho = 0.6, alpha = 0.4,  mu = 0.5),
  list(rho = 0.5, alpha = 0.5, mu = 3.0),
  list(rho = 0.5, alpha = 1.5, mu = 1.0),
  list(rho = 0.8, alpha = 0.3, mu = 2.0),
  list(rho = 0.9, alpha = 2.0,  mu = 1.2),
  list(rho = 0.4, alpha = 1.8, mu = 0.5),
  list(rho = 0.7, alpha = 0.7, mu = 1.5)
)

# Loop through each parameter set and plot the PMF
for (params in param_list) {
  res <- wrapped_stable_pmf( 
                            rho = params$rho, 
                            alpha = params$alpha, 
                            mu = params$mu)
  barplot(res$prob, names.arg = res$r, col = "blue",
          xlab = "y", ylab = "pmf",
          main = bquote(rho == .(params$rho) ~ 
                        alpha == .(params$alpha) ~ 
                        mu == .(params$mu)))
}




# Function to compute DWS CDF

dws_cdf <- function(m, q, rho, alpha, mu, P = 50) {
  term1 <- (q + 1) / m
  term2 <- 0
  for (k in 1:P) {
    term2 <- term2 + rho^(k^alpha) * sum(cos(k * (2 * pi * (0:q) / m - mu)))
  }
  cdf_val <- term1 + (2 / m) * term2
  return(cdf_val)
}

# Parameters
m <- 36
P <- 50  # truncation for k-sum
param_list <- list(
  list(rho=0.5, alpha=0.4, mu=0.2),
  list(rho=0.5, alpha=0.6, mu=0.4),
  list(rho=0.5, alpha=0.4, mu=0.3),
  list(rho=0.5, alpha=0.5, mu=1.1)
)

# Values of q (0 to m-1)
q_vals <- 0:(m-1)

# Compute CDFs for each parameter set
cdf_results <- lapply(param_list, function(par) {
  sapply(q_vals, function(q) dws_cdf(m, q, par$rho, par$alpha, par$mu, P))
})

# Plot
plot(q_vals, cdf_results[[1]], type="s", col=4, ylim=c(0,1),
     xlab="q", ylab="CDF F_w(theta)",
     main="CDF of Discrete Wrapped Stable Distribution")
lines(q_vals, cdf_results[[2]], type="s", col=2)
lines(q_vals, cdf_results[[3]], type="s", col=3)
lines(q_vals, cdf_results[[4]], type="s", col=1)
legend("bottomright", legend=c(
  expression(rho==0.5~alpha==0.4~mu==0.2),
  expression(rho==0.5~alpha==0.6~mu==0.4),
  expression(rho==0.5~alpha==0.4~mu==0.3),
  expression(rho==0.5~alpha==0.5~mu==1.1)
), col=1:4, lty=3)










#Simulations 

# PMF function with truncation at K
dws_pmf <- function(r, m, rho, alpha, mu, K = 50) {
  k_seq <- 1:K
  term_sum <- sum(rho^(k_seq^alpha) * cos(k_seq * (2 * pi * r / m - mu)))
  prob <- (1 + 2 * term_sum) / m
  prob
}

# Function to sample from DWS
rdws <- function(n, m, rho, alpha, mu, K = 50) {
  r_vals <- 0:(m - 1)
  pmf_vals <- sapply(r_vals, dws_pmf, m = m, rho = rho, alpha = alpha, mu = mu, K = K)
  pmf_vals <- pmf_vals / sum(pmf_vals)
  sample(r_vals, size = n, replace = TRUE, prob = pmf_vals)
}

# Log-likelihood for MLE
loglik_dws <- function(params, data, m, K = 50) {
  rho <- params[1]
  alpha <- params[2]
  mu <- params[3]
  
  if (rho <= 0 || rho >= 1 || alpha <= 0 || alpha > 5) return(1e10)  # penalty
  # wrap mu into [0, 2pi)
  mu <- mu %% (2 * pi)
  
  pmf_vals <- sapply(data, dws_pmf, m = m, rho = rho, alpha = alpha, mu = mu, K = K)
  pmf_vals[pmf_vals < 1e-12] <- 1e-12
  -sum(log(pmf_vals)) # return negative log-lik for minimization
}

# Method of Moments Estimation (MME)
mme_dws <- function(data, m) {
  theta <- 2 * pi * data / m
  C1 <- mean(cos(theta))
  S1 <- mean(sin(theta))
  rho1 <- sqrt(C1^2 + S1^2)
  mu_hat <- atan2(S1, C1)
  if (mu_hat < 0) mu_hat <- mu_hat + 2 * pi
  
  C2 <- mean(cos(2 * theta))
  S2 <- mean(sin(2 * theta))
  rho2 <- sqrt(C2^2 + S2^2)
  
  alpha_hat <- log(log(rho2) / log(rho1)) / log(2)
  rho_hat <- rho1^(1 / (1^alpha_hat))
  
  c(rho_hat, alpha_hat, mu_hat)
}

set.seed(123)
m <- 36
rho_true <- 0.5
alpha_true <- 0.5
mu_true <- 3
params_true <- c(rho_true, alpha_true, mu_true)
n_vals <- c(50, 100, 200)
R <- 3000
K <- 50

results <- list()

for (n in n_vals) {
  mle_est <- matrix(NA, R, 3)
  mme_est <- matrix(NA, R, 3)
  
  for (i in 1:R) {
    data <- rdws(n, m, rho_true, alpha_true, mu_true, K)
    
    # MLE
    fit_mle <- optim(
      par = c(0.5, 1, pi / 3),
      fn = loglik_dws, data = data, m = m, K = K,
      method = "L-BFGS-B",
      lower = c(0.01, 0.01, -pi),
      upper = c(0.99, 5, 3 * pi)
    )
    mle_est[i, ] <- fit_mle$par
    
    # MME
    mme_est[i, ] <- tryCatch(mme_dws(data, m), error = function(e) c(NA, NA, NA))
  }
  
  bias_mle <- colMeans(mle_est, na.rm = TRUE) - params_true
  mse_mle  <- colMeans((mle_est - matrix(params_true, R, 3, byrow = TRUE))^2, na.rm = TRUE)
  
  bias_mme <- colMeans(mme_est, na.rm = TRUE) - params_true
  mse_mme  <- colMeans((mme_est - matrix(params_true, R, 3, byrow = TRUE))^2, na.rm = TRUE)
  
  results[[paste0("n=", n)]] <- list(
    bias_mle = bias_mle,
    mse_mle = mse_mle,
    bias_mme = bias_mme,
    mse_mme = mse_mme
  )
}

results


# Rose Diagram

library(circular)

data <- c()

angles <- circular(data, units = "degrees", template = "geographics", modulo = "2pi")

rose.diag(angles, bins = 36, col = "lightblue", prop = 1.5,
          main = "Rose Diagram of Data")



# Fitted PMF

dws_pmf <- function(r, m, rho, alpha, mu, K = 50) {
  # r = 0,...,m-1
  pmf <- numeric(length(r))
  for (i in seq_along(r)) {
    sum_term <- 0
    for (k in 1:K) {
      sum_term <- sum_term + rho^(k^alpha) * cos(k * (2*pi*r[i]/m - mu))
    }
    pmf[i] <- (1 + 2*sum_term) / m
  }
  return(pmf)
}

data <- c()
m <- 36   
r_data <- floor(data / 10) %% m

obs_freq <- table(factor(r_data, levels = 0:(m-1)))
obs_prob <- obs_freq / sum(obs_freq)

rho   <- 0.7
alpha <- 1.2
mu    <- 0

r_vals <- 0:(m-1)
pmf_vals <- dws_pmf(r_vals, m, rho, alpha, mu)

barplot(obs_prob, names.arg = r_vals, col = "lightblue",
        main = "Observed Histogram vs DWS PMF",
        xlab = "Bin (r)", ylab = "Probability")

points(r_vals, pmf_vals, col = "red", pch = 19, type = "b", lwd = 2)
legend("topright", legend = c("Observed", "DWS Pmf"),
       fill = c("lightblue", NA), border = c("black", NA),
       pch = c(NA,19), col = c("black","red"), lty = c(NA,1),cex=0.6)

data <- read.csv("C:\\Users\\smrit\\OneDrive\\Desktop\\tezpur.csv")
data_col5 <- data[[5]]
cat("c(", paste(data_col5, collapse = ", "), ")", sep = "")

mat <- rbind(obs_prob, pmf_vals)
barplot(mat, beside = TRUE, names.arg = r_vals,
        col = c("lightgreen", "Purple"),
        main = "Observed vs Fitted DWS (Barplot)",
        xlab = "Bin (r)", ylab = "Probability", cex.names=0.7)

legend("topright", legend = c("Observed", "DWS PMF"),
       fill = c("lightblue", "red"), border = "black", cex=0.7)


dws_pmf <- function(r, m, rho, alpha, mu, K = 50) {
  # r = 0,...,m-1
  pmf <- numeric(length(r))
  for (i in seq_along(r)) {
    sum_term <- 0
    for (k in 1:K) {
      sum_term <- sum_term + rho^(k^alpha) * cos(k * (2*pi*r[i]/m - mu))
    }
    pmf[i] <- (1 + 2*sum_term) / m
  }
  return(pmf)
}

data <- c()
# Convert degrees into discrete bins (0–360 in steps of 10)
m <- 36   # since bins of 10 degrees
r_data <- floor(data / 10) %% m

obs_freq <- table(factor(r_data, levels = 0:(m-1)))
obs_prob <- obs_freq / sum(obs_freq)

rho   <- 0.7
alpha <- 1.2
mu    <- 0

# Theoretical pmf
r_vals <- 0:(m-1)
pmf_vals <- dws_pmf(r_vals, m, rho, alpha, mu)

barplot(obs_prob, names.arg = r_vals, col = "lightblue",
        main = "Observed Histogram vs DWS PMF",
        xlab = "Bin (r)", ylab = "Probability")

points(r_vals, pmf_vals, col = "red", pch = 19, type = "b", lwd = 2)
legend("topright", legend = c("Observed", "DWS Pmf"),
       fill = c("lightblue", NA), border = c("black", NA),
       pch = c(NA,19), col = c("black","red"), lty = c(NA,1),cex=0.6)


# Read CSV
data <- read.csv("C:\\Users\\smrit\\OneDrive\\Desktop\\tezpur.csv")

# Select the 5th column as a vector
data_col5 <- data[[5]]

cat("c(", paste(data_col5, collapse = ", "), ")", sep = "")


mat <- rbind(obs_prob, pmf_vals)

# Side-by-side barplot
barplot(mat, beside = TRUE, names.arg = r_vals,
        col = c("lightgreen", "Purple"),
        main = "Observed vs Fitted DWS (Barplot)",
        xlab = "Bin (r)", ylab = "Probability", cex.names=0.7)

legend("topright", legend = c("Observed", "DWS PMF"),
       fill = c("lightblue", "red"), border = "black", cex=0.7)





#AIC, BIC

rm(list=ls())
library(MASS)

safe_log <- function(x) log(pmax(x, 1e-12))

get_se <- function(res) {
  if (is.null(res$hessian)) return(rep(NA, length(res$par)))

  H <- res$hessian
  H <- H + diag(1e-6, nrow(H))

  invH <- tryCatch(MASS::ginv(H), error=function(e) NULL)
  if (is.null(invH)) return(rep(NA, length(res$par)))

  d <- diag(invH)
  d[d < 0] <- NA

  sqrt(d)
}

dws_pmf <- function(r, m, rho, alpha, mu, K=100) {
  series <- sum(sapply(1:K, function(k)
    rho^(k^alpha) * cos(k*(2*pi*r/m - mu))
  ))
  pmax((1 + 2*series)/m, 1e-12)
}

wg_pmf <- function(r,m,p){
  pmax((p*(1-p)^r)/(1-(1-p)^m),1e-12)
}

wp_pmf <- function(r,m,lambda,K=100){
  pmax(sum(dpois(r + (0:K)*m, lambda)),1e-12)
}

wb_pmf <- function(r,m,n,p){
  ks <- 0:floor(n/m)
  pmax(sum(dbinom(r+ks*m, n, p)),1e-12)
}

wnb_pmf <- function(x,m,r,p,K=120){
  pmax(sum(dnbinom(x + (0:K)*m, size=r, prob=p)),1e-12)
}

wnorm_pmf <- function(r,m,mu,sigma,K=30){
  pmax(sum(dnorm(r + (-K:K)*m, mu, sigma)),1e-12)
}

dws_loglik <- function(params,data,m,K=100){
  rho<-params[1]; alpha<-params[2]; mu<-params[3]
  if(rho<=0||rho>=1||alpha<=0||alpha>2||mu<0||mu>2*pi) return(-1e10)
  sum(safe_log(sapply(data,function(r)dws_pmf(r,m,rho,alpha,mu,K))))
}

wg_loglik <- function(params,data,m){
  p<-params[1]
  if(p<=0||p>=1) return(-1e10)
  sum(safe_log(sapply(data,function(r)wg_pmf(r,m,p))))
}

wp_loglik <- function(params,data,m){
  lambda<-params[1]
  if(lambda<=0) return(-1e10)
  sum(safe_log(sapply(data,function(r)wp_pmf(r,m,lambda))))
}

wb_loglik <- function(params,data,m,n){
  p<-params[1]
  if(p<=0||p>=1) return(-1e10)
  sum(safe_log(sapply(data,function(r)wb_pmf(r,m,n,p))))
}

wnb_loglik <- function(params,data,m){
  r<-params[1]; p<-params[2]
  if(r<=0||p<=0||p>=1) return(-1e10)
  sum(safe_log(sapply(data,function(x)wnb_pmf(x,m,r,p))))
}

wnorm_loglik <- function(params,data,m){
  mu<-params[1]; sigma<-params[2]
  if(sigma<=0) return(-1e10)
  sum(safe_log(sapply(data,function(r)wnorm_pmf(r,m,mu,sigma))))
}

fit_model <- function(loglik_fun,start,lower,upper,data,m,extra_args=list()){

  fn <- function(p){
    val <- tryCatch({
      ll <- do.call(loglik_fun,c(list(p,data,m),extra_args))
      if(!is.finite(ll)) return(1e10)
      -ll
    }, error=function(e) 1e10)
    val
  }

  res <- optim(start, fn, method="L-BFGS-B",
               lower=lower, upper=upper, hessian=TRUE)

  loglik <- -res$value
  k <- length(start)

  list(
    params=res$par,
    SE=get_se(res),
    loglik=loglik,
    AIC=-2*loglik+2*k,
    BIC=-2*loglik+k*log(length(data))
  )
}

data <- c()

m <- 24

fit_dws <- fit_model(dws_loglik,c(0.3,1.2,1),c(1e-4,0.2,0),c(0.95,2,2*pi),data,m)
fit_wg  <- fit_model(wg_loglik,c(0.4),c(1e-4),c(0.95),data,m)
fit_wp  <- fit_model(wp_loglik,c(mean(data)),c(1e-3),c(50),data,m)
fit_wb  <- fit_model(wb_loglik,c(0.5),c(1e-4),c(0.95),data,m,list(n=30))
fit_wnb <- fit_model(wnb_loglik,c(5,0.4),c(1e-3,1e-4),c(50,0.95),data,m)

print(list(DWS=fit_dws, WG=fit_wg, WP=fit_wp,
           WB=fit_wb, WNB=fit_wnb)