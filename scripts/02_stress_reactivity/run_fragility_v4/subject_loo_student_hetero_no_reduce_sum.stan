functions {
  // LOO per soggetto: marginalizza b_s ~ Normal(0, sigma_b) con residui Student-t
  // e scala eteroschedastica per osservazione. Integrazione 1D via Gauss–Hermite.
  // log ∫ N(b|0,sigma_b^2) ∏_i t(r_i | nu, b, sigma_i) db
  real loglik_subject_marginal_student_t_hetero(
      vector r,                    // residui y_s - X_s*beta (lunghezza ns)
      vector sigma_obs,            // scala per ciascuna osservazione del soggetto
      real sigma_b,                // sd random intercept
      real nu,                     // df t
      array[] real gh_nodes,       // nodi GH
      array[] real gh_weights      // pesi GH
  ) {
    int Q = size(gh_nodes);
    int ns = rows(r);
    vector[Q] lterms;
    for (j in 1:Q) {
      real b = sqrt(2) * sigma_b * gh_nodes[j];
      real acc = 0;
      for (i in 1:ns) {
        acc += student_t_lpdf(r[i] | nu, b, sigma_obs[i]);
      }
      lterms[j] = acc + log(gh_weights[j]);
    }
    return log_sum_exp(lterms) - 0.5 * log(pi()); // log(1/sqrt(pi))
  }
}

data {
  int<lower=1> N;                         // osservazioni
  int<lower=1> S;                         // soggetti
  int<lower=1> K;                         // predittori fissi (no intercetta)
  vector[N] y;                            // outcome
  matrix[N, K] X;                         // design matrix
  array[N] int<lower=1, upper=S> subj;    // soggetto per riga

  // Dummies periodo per eteroschedasticità
  array[N] int<lower=0, upper=1> is_pre;
  array[N] int<lower=0, upper=1> is_post;

  // Indici per GQ marginale a livello soggetto
  array[N] int<lower=1, upper=N> obs_by_subj;
  array[S] int<lower=1, upper=N> start_subj;
  array[S] int<lower=0, upper=N> n_per_subj;

  // Gauss-Hermite
  int<lower=1> Q;
  array[Q] real gh_nodes;
  array[Q] real gh_weights;
}

parameters {
  vector[K] beta;            // effetti fissi (intercetta è nel random)
  vector[S] z_b;             // non-centered
  real<lower=0> sigma_b;     // sd random intercept
  real<lower=0> sigma;       // scala base Student-t
  real<lower=2> nu;          // df Student-t (>=2)
  real gamma_pre;            // effetti log-scala per periodo
  real gamma_post;
}

transformed parameters {
  vector[S] b_subj = sigma_b * z_b;
}

model {
  // Priori debolmente informative
  beta    ~ normal(0, 1.5);
  z_b     ~ normal(0, 1);
  sigma_b ~ normal(0, 1);
  sigma   ~ normal(0, 1);
  target += exponential_lpdf(nu - 2 | 1.0 / 15.0); // df ~ 2 + Exp(15)
  gamma_pre  ~ normal(0, 0.5);
  gamma_post ~ normal(0, 0.5);

  // Likelihood Student-t eteroschedastica
  for (n in 1:N) {
    real mu_n    = row(X, n) * beta + b_subj[subj[n]];
    real sig_n   = sigma * exp( (is_pre[n] ? gamma_pre : 0)
                              + (is_post[n] ? gamma_post : 0) );
    target += student_t_lpdf(y[n] | nu, mu_n, sig_n);
  }
}

generated quantities {
  vector[S] log_lik_subj;   // LOO per soggetto (marginalizzato su b_s)
  vector[N] resid;          // residui rispetto alla parte fissa
  vector[N] sigma_obs;      // scala per osservazione (utile per PPC/diagnostica)

  {
    vector[N] mu_fixed = X * beta;
    for (n in 1:N) {
      resid[n] = y[n] - mu_fixed[n];
      sigma_obs[n] = sigma * exp( (is_pre[n] ? gamma_pre : 0)
                                + (is_post[n] ? gamma_post : 0) );
    }
  }

  for (s in 1:S) {
    int ns = n_per_subj[s];
    if (ns == 0) {
      log_lik_subj[s] = 0;
    } else {
      int start = start_subj[s];
      vector[ns] r_s;
      vector[ns] sig_s;
      for (j in 1:ns) {
        int i = obs_by_subj[start + j - 1];
        r_s[j]   = resid[i];
        sig_s[j] = sigma_obs[i];
      }
      log_lik_subj[s] = loglik_subject_marginal_student_t_hetero(
        r_s, sig_s, sigma_b, nu, gh_nodes, gh_weights
      );
    }
  }
}
