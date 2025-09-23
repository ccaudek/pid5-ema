data {
  int<lower=1> N;                 // observations
  int<lower=1> S;                 // subjects
  int<lower=1> K;                 // fixed effects
  int<lower=1> R;                 // random effects per subject
  matrix[N, K] X;                 // fixed-effect design
  matrix[N, R] Z;                 // random-effect design
  array[N] int<lower=1, upper=S> sid; // subject id per row
  vector[N] y;                    // outcome (preferably standardized)
}

parameters {
  // Fixed effects
  vector[K] beta;

  // Random effects (non-centered, correlated)
  vector<lower=0>[R] tau;                 // RE standard deviations (half-N prior)
  cholesky_factor_corr[R] L_Omega;        // RE correlation (Cholesky)
  matrix[R, S] z_b;                       // standard normal REs

  // Heteroscedastic residual scale (subject-specific)
  real alpha_sigma;                       // global log sigma
  real<lower=0> sigma_delta;              // SD of log-sigma deviations
  vector[S] delta;                        // subject log-sigma deviations

  // Student-t residual df
  real<lower=2> nu;                       // degrees of freedom (>2 ensures finite variance)
}

transformed parameters {
  // Subject-specific random effects and residual sigmas
  matrix[R, S] b;                         // random effects by subject
  vector[S] sigma_s;                      // residual SD per subject

  b = diag_pre_multiply(tau, L_Omega) * z_b;    // R x S
  for (s in 1:S) {
    sigma_s[s] = exp(alpha_sigma + delta[s]);
  }
}

model {
  // Priors
  beta ~ normal(0, 1);                    // standardized predictors → N(0,1)

  tau ~ normal(0, 0.5);                   // slightly stronger shrinkage on RE scales
  L_Omega ~ lkj_corr_cholesky(4);         // stronger shrinkage toward independence
  to_vector(z_b) ~ normal(0, 1);          // non-centered REs

  alpha_sigma ~ normal(0, 1);
  sigma_delta ~ normal(0, 1);             // half-N via <lower=0>
  delta ~ normal(0, sigma_delta);

  nu ~ gamma(2, 0.1);                     // mildly regularizing heavy tails

  // Likelihood (conditional on b_s and sigma_s)
  {
    vector[N] mu;
    for (n in 1:N) {
      int s = sid[n];
      mu[n] = dot_product(row(X, n), beta) + dot_product(row(Z, n), col(b, s));
    }
    for (n in 1:N) {
      target += student_t_lpdf(y[n] | nu, mu[n], sigma_s[sid[n]]);
    }
  }
}

generated quantities {
  // Pointwise log-likelihood per observation (conditional), plus subject-summed version.
  vector[N] log_lik_obs;
  vector[S] subj_loglik_sum;
  {
    vector[N] mu;
    for (n in 1:N) {
      int s = sid[n];
      mu[n] = dot_product(row(X, n), beta) + dot_product(row(Z, n), col(b, s));
    }
    for (n in 1:N) {
      log_lik_obs[n] = student_t_lpdf(y[n] | nu, mu[n], sigma_s[sid[n]]);
    }
    // group-sum by subject
    for (s in 1:S) subj_loglik_sum[s] = 0;
    for (n in 1:N) {
      subj_loglik_sum[sid[n]] += log_lik_obs[n];
    }
  }
}

