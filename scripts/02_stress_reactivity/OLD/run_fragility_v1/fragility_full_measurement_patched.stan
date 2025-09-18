data {
  int<lower=1> I;
  int<lower=1> N_obs;
  int<lower=1> K;
  int<lower=1> P;
  int<lower=1> D;

  int<lower=1> N_items;
  array[N_items] int<lower=1, upper=7> y_item;
  array[N_items] int<lower=1, upper=K> item_id;
  array[N_items] int<lower=1, upper=N_obs> obs_id;

  array[N_obs] int<lower=1, upper=I> subject;
  array[N_obs] int<lower=1, upper=P> period;

  int<lower=0> M_ema;
  vector[M_ema] ema_val;
  array[M_ema] int<lower=2, upper=D> ema_dim;
  array[M_ema] int<lower=1, upper=N_obs> ema_obs;

  matrix[I, D] X_base;
  array[I] int<lower=0, upper=1> female;

  int<lower=0> N_interact;

  real<lower=0, upper=1> use_ema;
  int<lower=0, upper=1> use_interact;
}

parameters {
  vector[K] nu;
  vector[3] lambda_raw;
  array[K] ordered[6] tau;

  real<lower=0> sigma_pos;
  real<lower=0> sigma_neg;
  vector[N_obs] u_pos;
  vector[N_obs] u_neg;

  real gamma_pre;
  real gamma_post;
  real<lower=0> tau_pre;
  real<lower=0> tau_post;
  vector[I] delta_pre;
  vector[I] delta_post;

  real<lower=0> sigma_state;
  real<lower=1> nu_state;
  vector[N_obs] eps_obs;

  matrix[I, D] theta;
  vector<lower=0>[D] sigma_ema;

  real a_frag;
  vector[D] b_base;
  vector[D] c_ema;
  real b_female;
  vector[D] b_base_female;
  vector[D] c_ema_female;

  vector[max(1, N_interact)] b_interact;

  real b_var_base;
  real b_var_ema;
  real b_var_female;

  real<lower=0> tau_base;
  real<lower=0> tau_ema;
  real<lower=0> tau_interact;

  real<lower=0> sigma_diff;
  real<lower=1> nu_diff;
}

transformed parameters {
  vector[K] lambda;
  vector[N_obs] eta;
  vector[I] diff_frag;
  vector[I] var_emotional_base;

  lambda[2] = 1.0;
  lambda[4] =  exp(lambda_raw[3]);
  lambda[1] = -exp(lambda_raw[1]);
  lambda[3] = -exp(lambda_raw[2]);

  for (n in 1:N_obs) {
    int i = subject[n];
    int p = period[n];
    real m_pre  = (p == 2) ? (gamma_pre  + delta_pre[i])  : 0.0;
    real m_post = (p == 3) ? (gamma_post + delta_post[i]) : 0.0;
    eta[n] = theta[i, 1] + m_pre + m_post + sigma_state * eps_obs[n];
  }

  for (i in 1:I)
    diff_frag[i] = (gamma_pre + delta_pre[i]) - (gamma_post + delta_post[i]);

  for (i in 1:I) {
    int cnt = 0;
    real mu = 0;
    for (n in 1:N_obs)
      if (subject[n] == i && period[n] == 1) { mu += eta[n]; cnt += 1; }
    mu = (cnt > 0) ? mu / cnt : 0;
    real ss = 0;
    int cnt2 = 0;
    for (n in 1:N_obs)
      if (subject[n] == i && period[n] == 1) { ss += square(eta[n] - mu); cnt2 += 1; }
    var_emotional_base[i] = (cnt2 > 1) ? ss / (cnt2 - 1) : 0;
  }
}

model {
  nu ~ normal(0, 1);
  lambda_raw ~ normal(0, 0.3);
  for (k in 1:K) tau[k] ~ normal(0, 1);

  sigma_pos ~ exponential(3);
  sigma_neg ~ exponential(3);
  u_pos ~ normal(0, sigma_pos);
  u_neg ~ normal(0, sigma_neg);

  gamma_pre  ~ normal(0, 1);
  gamma_post ~ normal(0, 1);
  tau_pre  ~ exponential(3);
  tau_post ~ exponential(3);
  delta_pre  ~ normal(0, tau_pre);
  delta_post ~ normal(0, tau_post);

  sigma_state ~ exponential(3);
  nu_state ~ gamma(2, 0.3);
  eps_obs ~ student_t(nu_state, 0, 1);

  to_vector(theta) ~ normal(0, 1);
  sigma_ema ~ exponential(3);

  for (m in 1:M_ema) {
    int i = subject[ema_obs[m]];
    int d = ema_dim[m];
    ema_val[m] ~ normal(theta[i, d], sigma_ema[d]);
  }

  {
    vector[N_items] linpred;
    for (n in 1:N_items) {
      real m = (item_id[n] == 1 || item_id[n] == 3) ? u_pos[obs_id[n]] : u_neg[obs_id[n]];
      linpred[n] = nu[item_id[n]] + lambda[item_id[n]] * eta[obs_id[n]] + m;
    }
    for (n in 1:N_items)
      y_item[n] ~ ordered_logistic(linpred[n], tau[item_id[n]]);
  }

  a_frag ~ normal(0, 1.5);
  b_female ~ normal(0, 0.5);

  tau_base ~ exponential(2);
  tau_ema  ~ exponential(2);
  b_base ~ normal(0, tau_base);
  c_ema  ~ normal(0, tau_ema);
  b_base_female ~ normal(0, tau_base);
  c_ema_female  ~ normal(0, tau_ema);

  if (use_interact == 1 && N_interact > 0) {
    tau_interact ~ exponential(2);
    b_interact[1:N_interact] ~ normal(0, tau_interact);
  }

  b_var_base   ~ normal(0, 0.2);
  b_var_ema    ~ normal(0, 0.2);
  b_var_female ~ normal(0, 0.2);

  sigma_diff ~ exponential(3);
  nu_diff ~ gamma(2, 0.3);

  for (i in 1:I) {
    real interact_term = 0.0;
    if (use_interact == 1 && N_interact > 0) {
      int idx = 1;
      for (d1 in 1:(D-1))
        for (d2 in (d1+1):D) {
          interact_term += b_interact[idx] * X_base[i, d1] * theta[i, d2];
          idx += 1;
        }
    }
    real mu_frag =
      a_frag +
      (X_base[i] * b_base) +
      use_ema * (theta[i] * c_ema) +
      female[i] * b_female +
      female[i] * (X_base[i] * b_base_female) +
      female[i] * use_ema * (theta[i] * c_ema_female) +
      use_ema * interact_term +
      var_emotional_base[i] * (b_var_base + female[i] * b_var_female) +
      use_ema * var_emotional_base[i] * b_var_ema;

    target += student_t_lpdf(diff_frag[i] | nu_diff, mu_frag, sigma_diff);
  }
}

generated quantities {
  vector[I] mu_frag_hat;
  real R2_frag;

  vector[N_items] log_lik_item;
  vector[N_obs] log_lik_obs;
  for (n in 1:N_obs) log_lik_obs[n] = 0;

  {
    vector[N_items] linpred_gq;
    for (n in 1:N_items) {
      real m = (item_id[n] == 1 || item_id[n] == 3)
               ? u_pos[obs_id[n]] : u_neg[obs_id[n]];
      linpred_gq[n] = nu[item_id[n]] + lambda[item_id[n]] * eta[obs_id[n]] + m;
      log_lik_item[n] = ordered_logistic_lpmf(y_item[n] | linpred_gq[n], tau[item_id[n]]);
      log_lik_obs[obs_id[n]] += log_lik_item[n];
    }
  }

  for (i in 1:I) {
    real interact_term = 0.0;
    if (use_interact == 1 && N_interact > 0) {
      int idx = 1;
      for (d1 in 1:(D-1))
        for (d2 in (d1+1):D) {
          interact_term += b_interact[idx] * X_base[i, d1] * theta[i, d2];
          idx += 1;
        }
    }
    mu_frag_hat[i] =
      a_frag +
      (X_base[i] * b_base) +
      use_ema * (theta[i] * c_ema) +
      female[i] * b_female +
      female[i] * (X_base[i] * b_base_female) +
      female[i] * use_ema * (theta[i] * c_ema_female) +
      use_ema * interact_term +
      var_emotional_base[i] * (b_var_base + female[i] * b_var_female) +
      use_ema * var_emotional_base[i] * b_var_ema;
  }

  {
    real v_fit = variance(mu_frag_hat);
    real v_resid = (nu_diff > 2)
                   ? (nu_diff / (nu_diff - 2)) * square(sigma_diff)
                   : positive_infinity();
    R2_frag = v_fit / (v_fit + v_resid);
  }

  vector[I] log_lik_subj;
  for (i in 1:I)
    log_lik_subj[i] = student_t_lpdf(diff_frag[i] | nu_diff, mu_frag_hat[i], sigma_diff);
}
