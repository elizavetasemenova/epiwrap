functions {
  //x_r[1] = lambda
  //x_r[2] = gamma
  //x_r[3] = mu
  //x_r[4] = beta
  
  //x_i[1] = N
  
         real[] sis(real t,
         real[] y,
         real[] theta,
         real[] x_r,
         int[] x_i) {
                  real dydt[2];
                  //y[1] - S
                  //y[2] - I
                  dydt[1] =    x_r[1]*x_i[1] -  x_r[4]* y[1]*y[2]/x_i[1] +  x_r[2] *y[2] - x_r[3]*y[1];
                  dydt[2] =    x_r[4]* y[1]*y[2]/x_i[1] - x_r[2] *y[2] - x_r[3]*y[2];
                  return dydt;

    }
}

data {
   int<lower=1> T;
   real y0[2];
   real t0;
   real ts[T];
   real lambda;
   real beta;
   real gamma;
   real mu;
   int N;
   real<lower=0> phi_inv;
}

transformed data {
           real x_r[4];
           int x_i[1];
           real theta[0];
           real y[T,2];
           real<lower=0> phi = 1. / phi_inv;
           
           x_i[1]=N;
           
           x_r[1]=lambda;
           x_r[2]=gamma;
           x_r[3]=mu;
           x_r[4]=beta;
           
           y = integrate_ode_bdf(sis, y0, t0, ts, theta, x_r, x_i);
}

model {
}

generated quantities {
  real pred_cases[T];
  pred_cases = neg_binomial_2_rng(col(to_matrix(y),2), phi);
}



    
 