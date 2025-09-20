functions {
  // Log-lik marginale per un soggetto (intercetta casuale integrata esattamente)
  real loglik_subject_marginal_normal(
      vector r,    // residui y_s - X_s*beta
      real sigma,  // sd residua (>0)
      real sigma_b // sd random intercept (>0)
  ) {
    int ns = rows(r);
    real s2  = sum(r);
    real s1  = dot_self(r);
    real ss  = square(sigma);
    real sb2 = square(sigma_b);
    real log_det = ns * log(ss) + log1p( (ns * sb2) / ss );
    real quad    = (s1 / ss) - (sb2 / (ss * (ss + ns * sb2))) * square(s2);
    return -0.5 * ( ns * log(2 * pi()) + log_det + quad );
  }
}

data {
  int<lower=1> N;                         // osservazioni
  int<lower=1> S;                         // soggetti
  int<lower=1> K;                         // predittori fissi (no intercetta)
  vector[N] y;                             // outcome
  matrix[N, K] X;                          // design matrix
  array[N] int<lower=1, upper=S> subj;     // soggetto per riga

  // Indici per GQ marginale a livello soggetto
  array[N] int<lower=1, upper=N> obs_by_subj;
  array[S] int<lower=1, upper=N> start_subj;
  array[S] int<lower=0, upper=N> n_per_subj;
}

parameters {
  vector[K] beta;            // effetti fissi (intercetta nel random)
  vector[S] z_b;             // non-centered
  real<lower=0> sigma_b;     // sd random intercept
  real<lower=0> sigma;       // sd residua
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

  // Likelihood (loop semplice, senza reduce_sum)
  for (n in 1:N) {
    real mu_n = row(X, n) * beta + b_subj[subj[n]];
    target += normal_lpdf(y[n] | mu_n, sigma);
  }
}

generated quantities {
  vector[S] log_lik_subj;   // LOO per soggetto (marginalizzato su b_s)

  // residui rispetto alla parte fissa
  vector[N] resid;
  {
    vector[N] mu_fixed = X * beta; // solo parte fissa; b_s integrato fuori
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
      log_lik_subj[s] = loglik_subject_marginal_normal(r_s, sigma, sigma_b);
    }
  }
}
