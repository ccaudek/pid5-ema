functions {
  // Somma parziale della log-likelihood osservazionale per reduce_sum
  // (modello gaussiano con intercetta casuale b[subj])
  real partial_sum_normal_ll(
      array[] real y,            // <-- array "sliceable"
      int start, int end,
      matrix X,                  // N x K
      vector beta,               // K
      vector b_subj,             // S
      array[] int subj,          // N (1..S)
      real sigma                 // >0
  ) {
    real lps = 0;
    for (n in start:end) {
      real mu_n = row(X, n) * beta + b_subj[subj[n]];
      lps += normal_lpdf(y[n] | mu_n, sigma);
    }
    return lps;
  }

  // Log-lik marginale per un soggetto (integrazione esatta dell'intercetta casuale)
  // y_s ~ Normal( X_s * beta + 1 * b_s , sigma ), b_s ~ Normal(0, sigma_b)
  // => y_s | beta ~ MVN( X_s * beta ,  sigma^2 I + sigma_b^2 J )
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
    real quad = (s1 / ss) - (sb2 / (ss * (ss + ns * sb2))) * square(s2);

    return -0.5 * ( ns * log(2 * pi()) + log_det + quad );
  }
}

data {
  int<lower=1> N;            // numero osservazioni
  int<lower=1> S;            // numero soggetti
  int<lower=1> K;            // numero predittori fissi
  array[N] real y;           // <-- esito come array
  matrix[N, K] X;            // design matrix fissa (senza intercetta)
  array[N] int<lower=1, upper=S> subj; // soggetto 1..S

  // Per reduce_sum
  int<lower=1> grainsize;

  // Indici per GQ marginale a livello soggetto
  array[N] int<lower=1, upper=N> obs_by_subj;
  array[S] int<lower=1, upper=N> start_subj;
  array[S] int<lower=0, upper=N> n_per_subj;
}

parameters {
  vector[K] beta;            // effetti fissi (no intercetta: è nel random)
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

  // Likelihood con parallelizzazione within-chain (sliceable = y)
  target += reduce_sum(
    partial_sum_normal_ll, y, grainsize,
    X, beta, b_subj, subj, sigma
  );
}

generated quantities {
  // LOO per soggetto: log-likelihood marginalizzata su b_s
  vector[S] log_lik_subj;

  // residui rispetto alla parte fissa
  vector[N] resid;
  {
    vector[N] mu_fixed = X * beta;
    for (n in 1:N) resid[n] = y[n] - mu_fixed[n];  // y è array → loop
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
