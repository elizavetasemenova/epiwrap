library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#==============================================
# wrapper functions
#==============================================

#----------- SIS generate ----------------

sis_generate <- function(T, y0, t0, ts,lambda, beta, gamma, mu, N) {
  standata <-  list (T  = T, y0 = y0, t0 = t0, ts = t,lambda =lambda, beta = beta, gamma =gamma, mu =mu, N=N)
  model <- stan_model("stan_models/sis_generate.stan")
  out <- rstan::sampling(model, data = standata, algorithm="Fixed_param", chains = 1, iter =1)
  return(out)
}

#----------- SIS infer ----------------

sis_infer <- function(T, y0, t0, ts,lambda, gamma, mu, N, y, beta_mu, beta_sigma, y_sigma) {
  standata <-  list (T  = T, y0 = y0, t0 = t0, ts = ts,lambda =lambda, gamma =gamma, mu =mu, N=N, y2_data=y,
                     beta_mu =beta_mu, beta_sigma =beta_sigma, y_sigma=y_sigma)
  model <- stan_model("stan_models/sis_infer.stan")
  out <- rstan::sampling(model, 
                         data = standata, 
                         iter=2000,
                         chains=2,
                         refresh=20,
                         control=list(stepsize=0.001))
  return(out)
}


#==============================================
# tests of wrapper functions
#==============================================

#----------- SIS generate ----------------

#times
n_years <- 6;  n_t <- n_years*12; t	<- seq(0, n_t, by = 1); t0 = t[1]; t <- t[-1]; T <- length(t);

# parameters
lambda <- 0.0001; beta <- 0.5; gamma <-0.01; mu <- 14/(1000*12);

#initial conditions
N <- 1000; i0 <- 1; s0 <-  N-i0; y0 = c( S=s0, I=i0);

# test sis_generate
out <- sis_generate(T, y0, t0, ts,lambda, beta, gamma, mu, N)

s <- rstan::extract(out)
str(s$y)
y <- data.frame(s$y[1,1:T,1:2])
names(y) <- c('S','I')
head(y)

plot(1, xlim=c(min(t), max(t)), ylim=c(min(y),max(y)), type='n', main='STAN Generated noisy SIS' )
lines(t, y[,1], col="red")
lines(t, y[,2], col="blue")

y2_data= y[,2]


#----------- SIS infer ----------------

#ODE parameters
lambda <- 0.0001; gamma <-0.01; mu <- 14/(1000*12);


#times
n_years <- 6;  n_t <- n_years*12; t	<- seq(0, n_t, by = 1); t0 = t[1]; ts <- t[-1]; T <- length(ts);

#initial conditions
N <- 1000
i0 <- 1; s0 <-  N-i0; y0 = c( S=s0, I=i0);

# priors
beta_mu <-0; beta_sigma <-1; y_sigma <-10

# test sis_infer
out <- sis_infer(T, y0, t0, ts,lambda, gamma, mu, N, y2_data, beta_mu, beta_sigma, y_sigma)


