data {
  int<lower=1> N;                 // observations
  int<lower=1> S;                 // subjects
  int<lower=1> K;                 // fixed effects
  int<lower=1> R;                 // random effects per subject
  matrix[N, K] X;                 // fixed-effect design
  matrix[N, R] Z;                 // random-effect design
  array[N] int<lower=1, upper=S> sid; // subject id per row
  vector[N] y;                    // outcome
}

parameters {
  // Fixed effects
  vector[K] beta;

  // Random effects (non-centered, correlated)
  vector<lower=0>[R] tau;                 // RE standard deviations
  cholesky_factor_corr[R] L_Omega;        // RE correlation (Cholesky)
  matrix[R, S] z_b;                       // standard normal REs

  // Heteroscedastic residual scale (subject-specific)
  real alpha_sigma;                       // global log sigma
  real<lower=0> sigma_delta;              // SD of log-sigma deviations
  vector[S] delta;                        // subject log-sigma deviations

  // (Optional) jitter or small ridge for numerical stability if needed
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
  tau ~ normal(0, 1);                     // half-N implied by <0 constraint not set, so:
                                          // use normal(0,1) with <lower=0> gives half-N
  L_Omega ~ lkj_corr_cholesky(2);         // mild shrinkage to identity

  to_vector(z_b) ~ normal(0, 1);          // non-centered REs

  alpha_sigma ~ normal(0, 1);
  sigma_delta ~ normal(0, 1);             // half-N(0,1) via <lower=0>
  delta ~ normal(0, sigma_delta);

  // Likelihood (conditional on b_s and sigma_s)
  // μ_n = x_n'β + z_n' b_{sid[n]}
  {
    vector[N] mu;
    for (n in 1:N) {
      int s = sid[n];
      mu[n] = dot_product(row(X, n), beta) + dot_product(row(Z, n), col(b, s));
    }
    for (n in 1:N) {
      target += normal_lpdf(y[n] | mu[n], sigma_s[sid[n]]);
    }
  }
}

generated quantities {
  // ---- Subject-wise *marginal* log-likelihood (integrating out b_s) ----
  // For each subject s:
  // y_s ~ Normal( X_s beta,   Z_s diag(tau^2) Z_s' + sigma_s^2 I )
  // We compute the exact Gaussian log-lik using Cholesky.
  vector[S] subj_loglik;
  {
    // Precompute per-subject blocks
    for (s in 1:S) {
      // gather indices for subject s
      int ns = 0;
      for (n in 1:N) if (sid[n] == s) ns += 1;

      // allocate blocks
      matrix[ns, K] Xs;
      matrix[ns, R] Zs;
      vector[ns] ys;
      int j = 1;
      for (n in 1:N) {
        if (sid[n] == s) {
          Xs[j] = row(X, n);
          Zs[j] = row(Z, n);
          ys[j] = y[n];
          j += 1;
        }
      }

      // Σ_s = Zs * diag_matrix(tau .* tau) * Zs' + sigma_s^2 * I
      matrix[ns, ns] Sigma = Zs * (diag_matrix(square(tau))) * Zs';
      for (i in 1:ns) Sigma[i, i] += square(sigma_s[s]);

      // μ_s = Xs * beta
      vector[ns] mu_s = Xs * beta;

      // log |Σ_s| and quadratic form via Cholesky
      matrix[ns, ns] L = cholesky_decompose(Sigma);
      vector[ns] r = mdivide_left_tri_low(L, ys - mu_s);
      real quad = dot_self(r);
      real logdet = 2 * sum(log(diagonal(L)));

      subj_loglik[s] = -0.5 * (ns * log(2 * pi()) + logdet + quad);
    }
  }
}
