functions {
  //theta[1] = beta
  
  //x_r[1] = lambda, total population
  //x_r[2] = gamma, total population
  //x_r[3] = mu, total population
  
  //x_i[1] = N, total population
  
  real[] sis(real t,
             real[] y,
             real[] theta,
             real[] x_r,
             int[] x_i) {
    real dydt[2];
    dydt[1] =    x_r[1]*x_i[1] -  theta[1]* y[1]*y[2]/x_i[1] +  x_r[2] *y[2] - x_r[3]*y[1];
    dydt[2] =    theta[1]* y[1]*y[2]/x_i[1] - x_r[2] *y[2] - x_r[3]*y[2];
    return dydt;
    
  }
}

data {
  int<lower=1> T;
  real y0[2];
  real t0;
  real ts[T];
  real lambda;
  real gamma;
  real mu;
  int N;
  int cases[T];
  real beta_mu;
  real<lower=0> beta_sigma;
  real<lower=0> a_phi;
  real<lower=0> b_phi;
}

transformed data {
  real x_r[3];
  int x_i[1];
  
  x_i[1]=N;
  
  x_r[1]=lambda;
  x_r[2]=gamma;
  x_r[3]=mu;
}

parameters {
  real<lower=0> beta;
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[T,2];
  real phi = 1. / phi_inv;
  {
    real theta[1];
    theta[1] = beta;
    
    //y = integrate_ode(sis, y0, t0, ts, theta, x_r, x_i);
    y = integrate_ode_bdf(sis, y0, t0, ts, theta, x_r, x_i);
  }
} 



model {
  beta ~ normal(beta_mu, beta_sigma);
  phi_inv ~ cauchy(a_phi, b_phi);
  //col(matrix x, int n) - The n-th column of matrix x

  cases ~ neg_binomial_2(col(to_matrix(y),2), phi);
}

generated quantities {
  real pred_cases[T];
  pred_cases = neg_binomial_2_rng(col(to_matrix(y),2), phi);
}


