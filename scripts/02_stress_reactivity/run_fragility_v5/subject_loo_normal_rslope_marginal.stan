functions {
  // Log-lik marginale per soggetto con Normale e random effects multivariati:
  // y_s | beta ~ N( Xs*beta , Sigma_s ), con Sigma_s = sigma^2 I + Zs G Zs'
  // G = L L' dove L = diag_pre_multiply(tau, L_Omega) è la Cholesky dei RE
  real loglik_subject_marginal_normal(
      vector r_s,              // residui y_s - Xs*beta  (n_s)
      matrix Zs,               // design random effects   (n_s x R)
      matrix L,                // Cholesky dei RE         (R x R)  (diag_pre * L_Omega)
      real sigma               // sd residua
  ) {
    int ns = rows(r_s);
    matrix[ns, ns] Sigma = diag_matrix(rep_vector(square(sigma), ns))
                           + (Zs * L) * (Zs * L)';   // sigma^2 I + (ZL)(ZL)'
    matrix[ns, ns] Ls = cholesky_decompose(Sigma);
    return multi_normal_cholesky_lpdf(r_s | rep_vector(0, ns), Ls);
  }
}

data {
  int<lower=1> N;                          // osservazioni
  int<lower=1> S;                          // soggetti
  int<lower=1> K;                          // predittori fissi (in X)
  int<lower=1> R;                          // n. random effects (in Z)
  vector[N] y;                             
  matrix[N, K] X;                          // design fisso
  matrix[N, R] Z;                          // design random (es. [1, is_pre, is_post, EMA...])
  array[N] int<lower=1, upper=S> subj;     // soggetto per riga

  // Indici per costruire blocchi per soggetto
  array[N] int<lower=1, upper=N> obs_by_subj;
  array[S] int<lower=1, upper=N> start_subj;
  array[S] int<lower=0, upper=N> n_per_subj;
}

parameters {
  vector[K] beta;                               // effetti fissi
  matrix[S, R] z_b;                             // non-centered RE (standard normal)
  vector<lower=0>[R] tau;                       // scale RE
  cholesky_factor_corr[R] L_Omega;              // correlazioni RE (Cholesky)
  real<lower=0> sigma;                          // sd residua
}

transformed parameters {
  matrix[R, R] L = diag_pre_multiply(tau, L_Omega); // Cholesky dei RE
  matrix[S, R] b_subj;                              // RE in scala
  for (s in 1:S) {
    // z_b[s] è row_vector -> converti a vector, moltiplica e torna row_vector
    b_subj[s] = (L * to_vector(z_b[s]'))';
  }
}

model {
  // Priori debolmente informative
  beta     ~ normal(0, 1.5);
  to_vector(z_b) ~ normal(0, 1);
  tau      ~ normal(0, 1);
  L_Omega  ~ lkj_corr_cholesky(1.5);
  sigma    ~ normal(0, 1);

  // Likelihood condizionale
  for (n in 1:N) {
    real mu_n = row(X, n) * beta
              + dot_product( to_vector(row(Z, n)'), to_vector(b_subj[subj[n]]') );
    target += normal_lpdf(y[n] | mu_n, sigma);
  }
}

generated quantities {
  vector[S] log_lik_subj;     // LOO per soggetto (MARGINALIZZATO sui RE)
  vector[N] resid;            // residui rispetto alla parte fissa (per diagnostica)

  // residui rispetto a X*beta (senza RE), usati per la marginale
  {
    vector[N] mu_fixed = X * beta;
    resid = y - mu_fixed;
  }

  for (s in 1:S) {
    int ns = n_per_subj[s];
    if (ns == 0) {
      log_lik_subj[s] = 0;
    } else {
      int start = start_subj[s];
      vector[ns] r_s;
      matrix[ns, R] Zs;
      for (j in 1:ns) {
        int i = obs_by_subj[start + j - 1];
        r_s[j]   = resid[i];
        Zs[j, ]  = row(Z, i);
      }
      log_lik_subj[s] = loglik_subject_marginal_normal(r_s, Zs, L, sigma);
    }
  }
}
