functions {
  // LOO per soggetto: marginalizza l’intercetta casuale b_s ~ Normal(0, sigma_b)
  // con residui Student-t. Integrazione 1D via Gauss–Hermite.
  // ∫ N(b|0,sigma_b^2) ∏_i t(r_i | nu, b, sigma) db
  // Con b = sqrt(2)*sigma_b*z, GH: (1/sqrt(pi)) * sum_j w_j * ∏_i t(r_i | nu, b_j, sigma)
  // dove b_j = sqrt(2)*sigma_b*nodes_j
  real loglik_subject_marginal_student_t(
      vector r,                    // residui y_s - X_s*beta
      real sigma,                  // scale t
      real sigma_b,                // sd random intercept
      real nu,                     // df t
      array[] real gh_nodes,       // nodi GH
      array[] real gh_weights      // pesi GH
  ) {
    int Q = size(gh_nodes);
    vector[Q] lterms;
    for (j in 1:Q) {                                 // <-- FIX: sintassi Stan
      real b = sqrt(2) * sigma_b * gh_nodes[j];
      lterms[j] = student_t_lpdf(r | nu, b, sigma)    // già somma sui r
                  + log(gh_weights[j]);
    }
    return log_sum_exp(lterms) - 0.5 * log(pi());     // log(1/sqrt(pi))
  }
}

data {
  int<lower=1> N;                         // osservazioni
  int<lower=1> S;                         // soggetti
  int<lower=1> K;                         // predittori fissi (no intercetta)
  vector[N] y;                            // outcome
  matrix[N, K] X;                         // design matrix
  array[N] int<lower=1, upper=S> subj;    // soggetto per riga

  // Indici per GQ marginale a livello soggetto
  array[N] int<lower=1, upper=N> obs_by_subj;
  array[S] int<lower=1, upper=N> start_subj;
  array[S] int<lower=0, upper=N> n_per_subj;

  // Gauss-Hermite (es. Q=20)
  int<lower=1> Q;
  array[Q] real gh_nodes;
  array[Q] real gh_weights;
}

parameters {
  vector[K] beta;            // effetti fissi (intercetta nel random)
  vector[S] z_b;             // non-centered
  real<lower=0> sigma_b;     // sd random intercept
  real<lower=0> sigma;       // scale Student-t
  real<lower=2> nu;          // df Student-t (>=2 → varianza finita)
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
  target += exponential_lpdf(nu - 2 | 1.0 / 15.0); // df ~ 2 + Exp(1/15) ⇒ E[nu]≈17

  // Likelihood Student-t con intercetta casuale
  for (n in 1:N) {
    real mu_n = row(X, n) * beta + b_subj[subj[n]];
    target += student_t_lpdf(y[n] | nu, mu_n, sigma);
  }
}

generated quantities {
  vector[S] log_lik_subj;   // LOO per soggetto (marginalizzato su b_s)

  // residui rispetto alla parte fissa (b_s integrato fuori)
  vector[N] resid;
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
      for (j in 1:ns) {
        int i = obs_by_subj[start + j - 1];
        r_s[j] = resid[i];
      }
      log_lik_subj[s] = loglik_subject_marginal_student_t(
        r_s, sigma, sigma_b, nu, gh_nodes, gh_weights
      );
    }
  }
}
