
// reactivity_latent_eta.stan (modern Stan syntax)
// Hierarchical latent reactivity model with period-specific factor per subject.
// - Baseline latent state fixed to 0; estimate pre and post with MVN pooling and predictors.
// - Measurement model links ping-level outcomes to the latent state with loadings.
//   Scale set by lambda[1] = 1.
//
// Period coding: 1=baseline, 2=pre_exam, 3=post_exam

data {
  int<lower=1> I;               // subjects
  int<lower=2> P;               // periods (here P=3)
  int<lower=1> K;               // number of outcomes

  // Outcome-specific long-format data (K=2 implemented below)
  array[K] int<lower=0> N;      // number of observations per outcome
  vector[N[1]] y1;
  array[N[1]] int<lower=1, upper=I> subj1;
  array[N[1]] int<lower=1, upper=P> period1;

  // If K >= 2:
  vector[N[2]] y2;
  array[N[2]] int<lower=1, upper=I> subj2;
  array[N[2]] int<lower=1, upper=P> period2;

  int<lower=1> Q;               // number of subject-level predictors (force >=1)
  matrix[I, Q] W;               // predictors (z-scored; if no real predictors, pass a zero column)
}

parameters {
  // Measurement model
  vector[K] alpha;                          // outcome intercepts
  vector<lower=0>[K] sigma;                 // outcome residual SDs
  vector<lower=0>[K-1] lambda_free;         // loadings for outcomes 2..K

  // Latent reactivity for [pre, post] (size P-1)
  vector[P-1] mu0;                          // population intercepts
  matrix[Q, P-1] B;                         // predictor effects
  vector<lower=0>[P-1] tau_eta;             // SDs
  cholesky_factor_corr[P-1] L_Omega;        // Cholesky corr
  matrix[I, P-1] z_eta;                     // standard normals (non-centered)
}

transformed parameters {
  vector[K] lambda;
  lambda[1] = 1.0;
  for (k in 2:K) lambda[k] = lambda_free[k-1];

  // Subject means for [pre, post]
  matrix[I, P-1] mu_eta;
  for (i in 1:I) {
    row_vector[P-1] mu_row = mu0';
    mu_row += W[i] * B;           // W[i] is row_vector[Q]; result is row_vector[P-1]
    mu_eta[i] = mu_row;
  }
}

model {
  // Priors
  alpha ~ student_t(3, 0, 1);
  sigma ~ exponential(1);
  lambda_free ~ normal(1, 0.5);

  mu0 ~ normal(0, 1);
  to_vector(B) ~ normal(0, 0.5);
  tau_eta ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(2);

  // standard normals for non-centered latents
  to_vector(z_eta) ~ normal(0, 1);

  // Likelihood
  {
    matrix[P-1, P-1] L_Sigma = diag_pre_multiply(tau_eta, L_Omega);

    // outcome 1
    for (n in 1:N[1]) {
      int i = subj1[n];
      int p = period1[n];
      real eta_ip;
      if (p == 1) {
        eta_ip = 0.0;
      } else {
        row_vector[P-1] zi = z_eta[i];
        vector[P-1] dev = L_Sigma * to_vector(zi');  // deviation for [pre, post]
        eta_ip = mu_eta[i, p-1] + dev[p-1];
      }
      y1[n] ~ normal(alpha[1] + lambda[1] * eta_ip, sigma[1]);
    }

    // outcome 2
    if (K >= 2) {
      for (n in 1:N[2]) {
        int i = subj2[n];
        int p = period2[n];
        real eta_ip;
        if (p == 1) {
          eta_ip = 0.0;
        } else {
          row_vector[P-1] zi = z_eta[i];
          vector[P-1] dev = L_Sigma * to_vector(zi');
          eta_ip = mu_eta[i, p-1] + dev[p-1];
        }
        y2[n] ~ normal(alpha[2] + lambda[2] * eta_ip, sigma[2]);
      }
    }
  }
}

generated quantities {
  corr_matrix[P-1] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
