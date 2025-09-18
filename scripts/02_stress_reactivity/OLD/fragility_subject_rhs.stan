data {
  int<lower=1> I;                 // soggetti
  int<lower=1> D_base;            // 5 (PID-5 baseline)
  int<lower=1> D_ema;             // 4 (EMA: det, ant, dis, psy)
  vector[I] y;                    // fragilità osservata: NA_pre - NA_post

  matrix[I, D_base] X_base;       // z-score per colonna
  matrix[I, D_ema]  X_ema;        // z-score per colonna

  array[I] int<lower=0, upper=1> female;
  int<lower=0, upper=1> use_gender;
  int<lower=0, upper=1> use_ema;  // 0 = Modello A, 1 = Modello B
}

parameters {
  real a;

  // Regularized Horseshoe - baseline
  vector[D_base] z_base;
  real<lower=0> hs_tau_base;
  vector<lower=0>[D_base] hs_lambda_base;
  real<lower=0> hs_c0_base;

  // Regularized Horseshoe - EMA
  vector[D_ema] z_ema;
  real<lower=0> hs_tau_ema;
  vector<lower=0>[D_ema] hs_lambda_ema;
  real<lower=0> hs_c0_ema;

  real b_female;

  real<lower=0> sigma;
  real<lower=1> nu;
}

transformed parameters {
  vector[D_base] b_base;
  vector[D_ema]  c_ema;
  {
    vector[D_base] lambda_tilde_base;
    vector[D_ema]  lambda_tilde_ema;

    for (d in 1:D_base) {
      lambda_tilde_base[d] = sqrt( square(hs_c0_base) * square(hs_lambda_base[d]) /
                                   ( square(hs_c0_base) + square(hs_tau_base) * square(hs_lambda_base[d]) ) );
      b_base[d] = z_base[d] * lambda_tilde_base[d] * hs_tau_base;
    }
    for (d in 1:D_ema) {
      lambda_tilde_ema[d] = sqrt( square(hs_c0_ema) * square(hs_lambda_ema[d]) /
                                  ( square(hs_c0_ema) + square(hs_tau_ema) * square(hs_lambda_ema[d]) ) );
      c_ema[d] = z_ema[d] * lambda_tilde_ema[d] * hs_tau_ema;
    }
  }
}

model {
  // priors
  a ~ normal(0, 1);

  hs_tau_base   ~ cauchy(0, 0.5);
  hs_lambda_base ~ cauchy(0, 1);
  hs_c0_base    ~ inv_gamma(2, 8);
  z_base        ~ normal(0, 1);

  hs_tau_ema    ~ cauchy(0, 0.5);
  hs_lambda_ema ~ cauchy(0, 1);
  hs_c0_ema     ~ inv_gamma(2, 8);
  z_ema         ~ normal(0, 1);

  b_female ~ normal(0, 0.3);

  sigma ~ exponential(1.5);
  nu    ~ gamma(2, 0.3);   // code robuste

  // likelihood
  for (i in 1:I) {
    real mu;
    mu = a + dot_product( row(X_base, i), b_base );
    if (use_ema == 1)    mu += dot_product( row(X_ema, i),  c_ema );
    if (use_gender == 1) mu += female[i] * b_female;
    target += student_t_lpdf(y[i] | nu, mu, sigma);
  }
}

generated quantities {
  vector[I] log_lik_subj;
  for (i in 1:I) {
    real mu;
    mu = a + dot_product( row(X_base, i), b_base );
    if (use_ema == 1)    mu += dot_product( row(X_ema, i),  c_ema );
    if (use_gender == 1) mu += female[i] * b_female;
    log_lik_subj[i] = student_t_lpdf(y[i] | nu, mu, sigma);
  }
}
