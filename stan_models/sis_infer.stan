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
  vector[T] y2_data;
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
}

transformed parameters{
  real y[T,2];
  {
    real theta[1];
    theta[1] = beta;
    
    //y = integrate_ode(sis, y0, t0, ts, theta, x_r, x_i);
    y = integrate_ode_bdf(sis, y0, t0, ts, theta, x_r, x_i);
  }
} 



model {
  beta ~ normal(0, 1);
  //col(matrix x, int n) - The n-th column of matrix x
  y2_data ~ normal(col(to_matrix(y),2), 10);
  
}

generated quantities {
}
