functions {
  // Log-lik marginale (già marginale su b_s) per dati del soggetto s
  // con sigma_s fissata:
  real ll_s_normal_given_sigma(
      vector r_s,             // residui y_s - X_s*beta  (n_s)
      matrix Zs,              // design random per s (n_s x R)
      matrix L,               // Cholesky dei RE (R x R)
      real sigma_s            // sd residua del soggetto s
  ) {
    int ns = rows(r_s);
    matrix[ns, ns] Sigma = diag_matrix(rep_vector(square(sigma_s), ns))
                           + (Zs * L) * (Zs * L)'; // sigma_s^2 I + (ZL)(ZL)'
    matrix[ns, ns] Ls = cholesky_decompose(Sigma);
    return multi_normal_cholesky_lpdf(r_s | rep_vector(0, ns), Ls);
  }
}

data {
  int<lower=1> N;
  int<lower=1> S;
  int<lower=1> K;
  int<lower=1> R;
  vector[N] y;
  matrix[N, K] X;
  matrix[N, R] Z;
  array[N] int<lower=1, upper=S> subj;
  array[N] int<lower=1, upper=N> obs_by_subj;
  array[S] int<lower=1, upper=N> start_subj;
  array[S] int<lower=0, upper=N> n_per_subj;

  // --- NEW: Gauss–Hermite per integrare delta_s (sigma_s soggetto) ---
  int<lower=5> Q;                 // es. 15 o 20
  array[Q] real gh_nodes;         // nodi Hermite standard
  array[Q] real gh_weights;       // pesi Hermite standard
}

parameters {
  vector[K] beta;
  matrix[S, R] z_b;
  vector<lower=0>[R] tau;
  cholesky_factor_corr[R] L_Omega;
  real<lower=0> sigma;
  vector[S] delta_raw;
  real<lower=0> sd_delta;
}

transformed parameters {
  matrix[R, R] L = diag_pre_multiply(tau, L_Omega);
  matrix[S, R] b_subj;
  for (s in 1:S)
    b_subj[s] = (L * to_vector(z_b[s]'))';
}

model {
  // priors
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
              + dot_product(to_vector(row(Z, n)'), to_vector(b_subj[subj[n]]'));
    // usa sigma_s soggetto-specifica in sampling (come prima)
    real sigma_s = sigma * exp(sd_delta * delta_raw[subj[n]]);
    target += normal_lpdf(y[n] | mu_n, sigma_s);
  }
}

generated quantities {
  vector[S] log_lik_subj;      // come prima: NON usata per LOO ora
  vector[S] log_lik_subj_gh;   // NEW: LOO-READY, con sigma_s marginalizzata
  vector[N] resid = y - X * beta;

  // log_lik_subj (con sigma_s puntuale) - utile solo per diagnostica
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
      // sigma_s puntuale (data delta_raw[s]) - solo per confronto interno
      real sigma_s_point = sigma * exp(sd_delta * delta_raw[s]);
      log_lik_subj[s] = ll_s_normal_given_sigma(r_s, Zs, L, sigma_s_point);
    }
  }

  // NEW: log_lik_subj_gh = ∫ ll_s(r_s | sigma_s(delta)) * φ(delta) d delta
  //     ≈ (1/sqrt(pi)) * Σ w_j * exp( ll_s( delta = sqrt(2)*x_j ) ), in log-sum-exp
  {
    real log_const = -0.5 * log(pi());
    for (s in 1:S) {
      int ns = n_per_subj[s];
      if (ns == 0) { log_lik_subj_gh[s] = 0; }
      else {
        int start = start_subj[s];
        vector[ns] r_s;
        matrix[ns, R] Zs;
        vector[Q] logs;  // log terms per nodo GH
        for (j in 1:ns) {
          int i = obs_by_subj[start + j - 1];
          r_s[j]  = resid[i];
          Zs[j, ] = row(Z, i);
        }
        for (q in 1:Q) {
          real delta = sqrt(2) * gh_nodes[q];
          real sigma_s_int = sigma * exp(sd_delta * delta);
          logs[q] = log(gh_weights[q]) + ll_s_normal_given_sigma(r_s, Zs, L, sigma_s_int);
        }
        // log-sum-exp per stabilità numerica
        real m = max(logs);
        log_lik_subj_gh[s] = log_const + m + log(sum(exp(logs - m)));
      }
    }
  }
}
