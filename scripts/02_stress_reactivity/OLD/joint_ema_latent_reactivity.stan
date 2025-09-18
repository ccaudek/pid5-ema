// Joint measurement-structural model for EMA negative affect reactivity
// - Ordinal measurement model for 4 affect items (1=happy, 2=sad, 3=satisfied, 4=angry)
// - Higher eta = more negative affect
// - Correlated uniquenesses via per-observation method factors u_pos (happy,satisfied)
//   and u_neg (sad,angry)
// - Random intercept + AR(1) dynamics over irregular time gaps
// - Reactivity to exam periods + between-subject predictors
// - sigma_eta is the stationary SD of the AR(1) process

data {
  // Sample sizes and dimensions
  int<lower=1> I;                           // subjects
  int<lower=1> N_obs;                       // EMA observations
  int<lower=1> K;                           // items (=4)
  int<lower=1> P;                           // periods (=3: 1=baseline,2=pre,3=post)
  int<lower=1> Q;                           // between-subject predictors

  // Item responses (ordinal 1..7), flattened long
  int<lower=1> N_items;
  array[N_items] int<lower=1, upper=7> y_item;
  array[N_items] int<lower=1, upper=K> item_id;       // 1=happy,2=sad,3=satisfied,4=angry
  array[N_items] int<lower=1, upper=N_obs> obs_id;    // observation index for each item

  // Observation-level data
  array[N_obs] int<lower=1, upper=I> subject;         // subject per obs
  array[N_obs] int<lower=1, upper=P> period;          // period per obs
  vector[N_obs] time_days;                             // time in days from start

  // Between-subject predictors (standardized)
  matrix[I, Q] W;

  // DASS validation (optimized with precomputed obs index)
  int<lower=0> N_dass;
  vector[N_dass] dass_stress;
  array[N_dass] int<lower=1, upper=N_obs> dass_obs;   // obs index aligned to each DASS row

  // (Optional) fast lookup of first EMA obs by (subject, period); 0=missing
  array[I, P] int<lower=0, upper=N_obs> first_obs_by_ip;
}

parameters {
  // Measurement model
  vector[K] nu;                                        // item intercepts
  vector[3] lambda_raw;                                // raw params for free loadings
  array[K] ordered[6] tau;                             // thresholds per item

  // Correlated uniquenesses via method factors (per observation)
  real<lower=0> sigma_pos;                             // SD of u_pos
  real<lower=0> sigma_neg;                             // SD of u_neg
  vector[N_obs] u_pos;                                 // method factor for [happy,satisfied]
  vector[N_obs] u_neg;                                 // method factor for [sad,angry]

  // Random intercepts (between)
  vector[I] mu_i;                                      // standardized RE for subjects
  real mu_0;                                           // population mean
  real<lower=0> sigma_mu;                              // SD of RE

  // AR(1) dynamics (stationary SD parametrization)
  real<lower=0, upper=1> phi;                          // AR coefficient per day
  real<lower=0> sigma_eta;                             // stationary SD of eta

  // Period + between-subject effects
  vector[P-1] gamma;                                   // period effects (baseline ref)
  vector[Q] beta;                                      // between-subject effects (on mu)

  // State innovations (standardized)
  vector[N_obs] eta_raw;

  // DASS validation parameters
  real alpha_dass;
  real lambda_dass;
  real<lower=0> sigma_dass;
}

transformed parameters {
  vector[K] lambda;                                    // factor loadings
  vector[N_obs] eta;                                   // latent states
  vector[I] mu_subj;                                   // subject-specific means

  // Loadings mapping:
  // sad fixed to +1; angry positive; happy/satisfied negative (higher eta = worse affect)
  lambda[2] = 1.0;                     // sad
  lambda[4] =  exp(lambda_raw[3]);     // angry  ( >0 )
  lambda[1] = -exp(lambda_raw[1]);     // happy  ( <0 )
  lambda[3] = -exp(lambda_raw[2]);     // satisfied ( <0 )

  // Subject intercepts
  mu_subj = mu_0 + W * beta + sigma_mu * mu_i;

  // AR(1) with irregular spacing; sigma_eta is stationary SD
  for (n in 1:N_obs) {
    int i = subject[n];
    int p = period[n];
    real period_effect = (p == 1) ? 0.0 : gamma[p-1];

    if (n == 1 || subject[n] != subject[n-1]) {
      // first obs for this subject
      eta[n] = mu_subj[i] + period_effect + sigma_eta * eta_raw[n];
    } else {
      real dt = time_days[n] - time_days[n-1];
      real phi_adj = pow(phi, dt);
      eta[n] = mu_subj[i]
             + period_effect
             + phi_adj * ( eta[n-1]
                           - mu_subj[ subject[n-1] ]
                           - ((period[n-1] == 1) ? 0.0 : gamma[period[n-1]-1]) )
             + sigma_eta * sqrt(1 - phi_adj * phi_adj) * eta_raw[n];
    }
  }
}

model {
  // Priors
  nu ~ normal(0, 1);
  lambda_raw ~ normal(0, 0.5);
  for (k in 1:K) tau[k] ~ normal(0, 2);                // ordered[] enforces monotonicity

  // Method factors
  sigma_pos ~ exponential(1);
  sigma_neg ~ exponential(1);
  u_pos ~ normal(0, sigma_pos);
  u_neg ~ normal(0, sigma_neg);

  // Between / dynamics
  mu_0 ~ normal(0, 1);
  sigma_mu ~ exponential(1);
  mu_i ~ std_normal();

  phi ~ beta(2, 1);                                    // favors persistence
  sigma_eta ~ exponential(1);                          // stationary SD

  gamma ~ normal(0, 0.5);
  beta  ~ normal(0, 0.5);

  eta_raw ~ std_normal();

  // DASS
  alpha_dass  ~ normal(0, 1);
  lambda_dass ~ normal(1, 0.5);
  sigma_dass  ~ exponential(1);

  // Measurement likelihood: add method factor by pair
  {
    vector[N_items] linpred;
    for (n in 1:N_items) {
      real m = (item_id[n] == 1 || item_id[n] == 3)      // happy or satisfied
               ? u_pos[ obs_id[n] ]
               : u_neg[ obs_id[n] ];                     // sad or angry
      linpred[n] = nu[item_id[n]] + lambda[item_id[n]] * eta[obs_id[n]] + m;
    }
    for (n in 1:N_items)
      y_item[n] ~ ordered_logistic(linpred[n], tau[item_id[n]]);
  }

  // DASS validation (optimized)
  if (N_dass > 0) {
    vector[N_dass] mu_dass = alpha_dass + lambda_dass * to_vector(eta[dass_obs]);
    dass_stress ~ normal(mu_dass, sigma_dass);
  }
}

generated quantities {
  // --- Fit on DASS
  real R2_dass;
  vector[N_dass] dass_pred;
  if (N_dass > 0) {
    vector[N_dass] mu_dass = alpha_dass + lambda_dass * to_vector(eta[dass_obs]);
    vector[N_dass] dass_resid = dass_stress - mu_dass;
    dass_pred = mu_dass;
    R2_dass   = 1 - variance(dass_resid) / variance(dass_stress);
  } else {
    R2_dass = negative_infinity();
    // dass_pred stays default-initialized length 0
  }

  // --- Period contrasts per subject (using fast lookup if available)
  matrix[I, P] eta_by_period;
  vector[I] delta_pre_baseline;
  vector[I] delta_post_baseline;
  vector[I] delta_post_pre;           // post - pre (positive = worse at post)
  vector[I] reactivity_pre_post;      // max(pre, post) - baseline

  for (i in 1:I) {
    for (p in 1:P) {
      if (first_obs_by_ip[i, p] > 0) {
        eta_by_period[i, p] = eta[first_obs_by_ip[i, p]];
      } else {
        eta_by_period[i, p] = negative_infinity();
      }
    }

    // Δ pre - baseline
    if (eta_by_period[i, 1] != negative_infinity() &&
        eta_by_period[i, 2] != negative_infinity()) {
      delta_pre_baseline[i]  = eta_by_period[i, 2] - eta_by_period[i, 1];
    } else delta_pre_baseline[i] = negative_infinity();

    // Δ post - baseline
    if (eta_by_period[i, 1] != negative_infinity() &&
        eta_by_period[i, 3] != negative_infinity()) {
      delta_post_baseline[i] = eta_by_period[i, 3] - eta_by_period[i, 1];
    } else delta_post_baseline[i] = negative_infinity();

    // Δ post - pre
    if (eta_by_period[i, 2] != negative_infinity() &&
        eta_by_period[i, 3] != negative_infinity()) {
      delta_post_pre[i] = eta_by_period[i, 3] - eta_by_period[i, 2];
    } else delta_post_pre[i] = negative_infinity();

    // Reactivity = max(pre, post) - baseline
    if (eta_by_period[i, 1] != negative_infinity()) {
      real mx = negative_infinity();
      if (eta_by_period[i, 2] != negative_infinity()) mx = fmax(mx, eta_by_period[i, 2]);
      if (eta_by_period[i, 3] != negative_infinity()) mx = fmax(mx, eta_by_period[i, 3]);
      reactivity_pre_post[i] =
        (mx != negative_infinity()) ? (mx - eta_by_period[i, 1]) : negative_infinity();
    } else {
      reactivity_pre_post[i] = negative_infinity();
    }
  }
}
