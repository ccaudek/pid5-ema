functions {
  real loglik_subject_marginal_normal_sigmaS(
      vector r_s,    // residui y_s - Xs*beta
      matrix Zs,     // design random (n_s x R)
      matrix L,      // Cholesky dei RE (R x R)
      real sigma_s   // sd residua del soggetto s
  ) {
    int ns = rows(r_s);
    matrix[ns, ns] Sigma = diag_matrix(rep_vector(square(sigma_s), ns))
                           + (Zs * L) * (Zs * L)';
    matrix[ns, ns] Ls = cholesky_decompose(Sigma);
    return multi_normal_cholesky_lpdf(r_s | rep_vector(0, ns), Ls);
  }
}

data {
  int<lower=1> N; int<lower=1> S; int<lower=1> K; int<lower=1> R;
  vector[N] y; matrix[N, K] X; matrix[N, R] Z;
  array[N] int<lower=1, upper=S> subj;
  array[N] int<lower=1, upper=N> obs_by_subj;
  array[S] int<lower=1, upper=N> start_subj;
  array[S] int<lower=0, upper=N> n_per_subj;
}

parameters {
  vector[K] beta;
  matrix[S, R] z_b;
  vector<lower=0>[R] tau;
  cholesky_factor_corr[R] L_Omega;
  real<lower=0> sigma;       // scala media
  vector[S] delta_raw;       // etero per soggetto (log-scale)
  real<lower=0> sd_delta;    // iper-sd su delta
}

transformed parameters {
  matrix[R, R] L = diag_pre_multiply(tau, L_Omega);
  matrix[S, R] b_subj;
  vector<lower=0>[S] sigma_s;

  for (s in 1:S)
    b_subj[s] = (L * to_vector(z_b[s]'))';

  for (s in 1:S)
    sigma_s[s] = sigma * exp(sd_delta * delta_raw[s]);
}

model {
  // priori
  beta ~ normal(0, 1.5);
  to_vector(z_b) ~ normal(0, 1);
  tau ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(1.5);
  sigma ~ normal(0, 1);
  delta_raw ~ normal(0, 1);
  sd_delta ~ normal(0, 0.5);

  // likelihood condizionale
  for (n in 1:N) {
    real mu_n = row(X, n) * beta
              + dot_product( to_vector(row(Z, n)'), to_vector(b_subj[subj[n]]') );
    target += normal_lpdf(y[n] | mu_n, sigma_s[subj[n]]);
  }
}

generated quantities {
  vector[S] log_lik_subj;
  vector[N] resid = y - X * beta; // residui rispetto alla parte fissa

  for (s in 1:S) {
    int ns = n_per_subj[s];
    if (ns == 0) { log_lik_subj[s] = 0; }
    else {
      int start = start_subj[s];
      vector[ns] r_s;
      matrix[ns, R] Zs;
      for (j in 1:ns) {
        int i = obs_by_subj[start + j - 1];
        r_s[j]  = resid[i];
        Zs[j, ] = row(Z, i);
      }
      log_lik_subj[s] = loglik_subject_marginal_normal_sigmaS(r_s, Zs, L, sigma_s[s]);
    }
  }
}
