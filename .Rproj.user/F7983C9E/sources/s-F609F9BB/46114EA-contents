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
}

transformed data {
           real x_r[4];
           int x_i[1];
           real theta[0];
           
           x_i[1]=N;
           
           x_r[1]=lambda;
           x_r[2]=gamma;
           x_r[3]=mu;
           x_r[4]=beta;
}

model {
}

generated quantities {
    real y[T,2];
    real n1;
    real n2;

    //y = integrate_ode(sis, y0, t0, ts, theta, x_r, x_i);
    y = integrate_ode_bdf(sis, y0, t0, ts, theta, x_r, x_i);
    
   // add noise
  for (t in 1:T) {
     n1 = normal_rng(0, 10);
     n2 = normal_rng(0, 10);
     
     if (y[t,1]+  n1>0)
       y[t,1] =  y[t,1]+n1;
     else 
      y[t,1] = 0;

     if (y[t,2]+  n2>0)
       y[t,2] =  y[t,2]+n2;
     else 
      y[t,2] = 0;
  }
}

