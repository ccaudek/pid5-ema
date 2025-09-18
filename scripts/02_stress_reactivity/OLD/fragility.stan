data {
  int<lower=1> I;                 // soggetti
  int<lower=1> N_obs;             // osservazioni EMA
  int<lower=1> K;                 // items = 4 (happy, sad, satisfied, angry)
  int<lower=1> P;                 // periodi = 3 (1=base,2=pre,3=post)
  int<lower=1> D;                 // dimensioni EMA a livello soggetto = 5

  // Item responses ordinali 1..7 (solo non-mancanti, in long)
  int<lower=1> N_items;
  array[N_items] int<lower=1, upper=7> y_item;
  array[N_items] int<lower=1, upper=K> item_id;
  array[N_items] int<lower=1, upper=N_obs> obs_id;

  // Mappa osservazione -> soggetto/periodo
  array[N_obs] int<lower=1, upper=I> subject;
  array[N_obs] int<lower=1, upper=P> period;

  // EMA osservate (SOLO per dimensioni 2..5) già z-score in R
  int<lower=0> M_ema;
  vector[M_ema] ema_val;
  array[M_ema] int<lower=2, upper=D> ema_dim;
  array[M_ema] int<lower=1, upper=N_obs> ema_obs;

  // Predittori baseline (I x 5), già z-score
  matrix[I, D] X_base;
  
  // Genere: 1=femmina, 0=maschio
  array[I] int<lower=0, upper=1> female;
  
  // Numero di interazioni baseline x EMA (triangolare superiore)
  int<lower=0> N_interact;

  // Interruttore per includere le EMA nel modello di fragilità
  real<lower=0, upper=1> use_ema;
}

parameters {
  // Misurazione ordinale per i 4 item
  vector[K] nu;
  vector[3] lambda_raw;
  array[K] ordered[6] tau;

  // Fattori di metodo (correlated uniqueness per valenza)
  real<lower=0> sigma_pos;
  real<lower=0> sigma_neg;
  vector[N_obs] u_pos;
  vector[N_obs] u_neg;

  // Effetti di periodo e variazioni soggetto-specifiche
  real gamma_pre;
  real gamma_post;
  real<lower=0> tau_pre;
  real<lower=0> tau_post;
  vector[I] delta_pre;
  vector[I] delta_post;

  // Variabilità intra-periodo (stato)
  real<lower=0> sigma_state;
  real<lower=1> nu_state;
  vector[N_obs] eps_obs;

  // Tratti EMA per soggetto (5 dimensioni)
  matrix[I, D] theta;
  vector<lower=0>[D] sigma_ema;

  // === Modello di fragilità ===
  real a_frag;              // intercetta
  vector[D] b_base;         // effetti principali baseline
  vector[D] c_ema;          // effetti principali EMA
  real b_female;            // effetto principale genere
  vector[D] b_base_female;  // interazioni genere x baseline
  vector[D] c_ema_female;   // interazioni genere x EMA

  // !!! QUI: niente ridichiarazione di N_interact !!!
  vector[N_interact] b_interact;   // coefficienti interazioni baseline x EMA

  // === Variabilità emotiva ===
  real b_var_base;
  real b_var_ema;
  real b_var_female;

  // Parametri di shrinkage
  real<lower=0> tau_base;
  real<lower=0> tau_ema;
  real<lower=0> tau_interact;
  
  // Residui robusti
  real<lower=0> sigma_diff;
  real<lower=1> nu_diff;
}

transformed parameters {
  vector[K] lambda;
  vector[N_obs] eta;
  vector[I] diff_frag;        // fragilità = (pre) - (post)
  vector[I] var_emotional;    // variabilità emotiva per soggetto

  // Carichi degli item
  lambda[2] = 1.0;                     // sad (riferimento)
  lambda[4] = exp(lambda_raw[3]);      // angry > 0
  lambda[1] = -exp(lambda_raw[1]);     // happy < 0
  lambda[3] = -exp(lambda_raw[2]);     // satisfied < 0

  // Stati latenti di negative affect
  for (n in 1:N_obs) {
    int i = subject[n];
    int p = period[n];
    real m_pre  = (p == 2) ? (gamma_pre + delta_pre[i])  : 0.0;
    real m_post = (p == 3) ? (gamma_post + delta_post[i]) : 0.0;
    eta[n] = theta[i,1] + m_pre + m_post + sigma_state * eps_obs[n];
  }

  // Fragilità come differenza pre-post
  for (i in 1:I)
    diff_frag[i] = (gamma_pre + delta_pre[i]) - (gamma_post + delta_post[i]);

  // Variabilità emotiva per soggetto (varianza empirica degli eta)
  for (i in 1:I) {
    vector[N_obs] is_subj_i;
    vector[N_obs] eta_subj;
    int n_subj_obs = 0;
    for (n in 1:N_obs) {
      if (subject[n] == i) {
        n_subj_obs += 1;
        is_subj_i[n] = 1.0;
        eta_subj[n] = eta[n];
      } else {
        is_subj_i[n] = 0.0;
        eta_subj[n] = 0.0;
      }
    }
    if (n_subj_obs > 1) {
      real mean_eta = sum(eta_subj .* is_subj_i) / n_subj_obs;
      real sum_sq_dev = 0.0;
      for (n in 1:N_obs)
        if (subject[n] == i)
          sum_sq_dev += square(eta[n] - mean_eta);
      var_emotional[i] = sum_sq_dev / (n_subj_obs - 1);
    } else {
      var_emotional[i] = 0.0;
    }
  }
}

model {
  // --- Sanity check sulla coerenza di N_interact ---
  if (N_interact != D * (D - 1) / 2)
    reject("N_interact (", N_interact, ") non uguale a D*(D-1)/2 (", D * (D - 1) / 2, ").");

  // Priors misurazione
  nu ~ normal(0, 1);
  lambda_raw ~ normal(0, 0.3);
  for (k in 1:K) tau[k] ~ normal(0, 1);

  sigma_pos ~ exponential(1);
  sigma_neg ~ exponential(1);
  u_pos ~ normal(0, sigma_pos);
  u_neg ~ normal(0, sigma_neg);

  gamma_pre ~ normal(0, 1);
  gamma_post ~ normal(0, 1);
  tau_pre ~ exponential(1);
  tau_post ~ exponential(1);
  delta_pre ~ normal(0, tau_pre);
  delta_post ~ normal(0, tau_post);

  sigma_state ~ exponential(1);
  nu_state ~ gamma(2, 0.1);
  eps_obs ~ student_t(nu_state, 0, 1);

  to_vector(theta) ~ normal(0, 1);
  sigma_ema ~ exponential(1);

  // Misurazione EMA (dimensioni 2..5)
  for (m in 1:M_ema) {
    int i = subject[ema_obs[m]];
    int d = ema_dim[m];
    ema_val[m] ~ normal(theta[i, d], sigma_ema[d]);
  }

  // Likelihood item ordinali
  {
    vector[N_items] linpred;
    for (n in 1:N_items) {
      real m = (item_id[n] == 1 || item_id[n] == 3) 
               ? u_pos[obs_id[n]] : u_neg[obs_id[n]];
      linpred[n] = nu[item_id[n]] + lambda[item_id[n]] * eta[obs_id[n]] + m;
    }
    for (n in 1:N_items)
      y_item[n] ~ ordered_logistic(linpred[n], tau[item_id[n]]);
  }

  // Priors per il modello di fragilità
  a_frag ~ normal(0, 2);
  b_female ~ normal(0, 1);
  b_var_base ~ normal(0, 0.5);
  b_var_ema ~ normal(0, 0.5);
  b_var_female ~ normal(0, 0.5);
  
  tau_base ~ exponential(1);
  tau_ema ~ exponential(1);
  tau_interact ~ exponential(2);
  
  b_base ~ normal(0, tau_base);
  c_ema ~ normal(0, tau_ema);
  b_base_female ~ normal(0, tau_base);
  c_ema_female ~ normal(0, tau_ema);
  b_interact ~ normal(0, tau_interact);
  
  sigma_diff ~ exponential(1);
  nu_diff ~ gamma(2, 0.1);

  // Regressione della fragilità
  for (i in 1:I) {
    real mu_frag;
    real interact_term = 0.0;
    int idx = 1;
    for (d1 in 1:(D-1)) {
      for (d2 in (d1+1):D) {
        interact_term += b_interact[idx] * X_base[i, d1] * theta[i, d2];
        idx += 1;
      }
    }
    mu_frag = a_frag +
              (X_base[i] * b_base) +
              use_ema * (theta[i] * c_ema) +
              female[i] * b_female +
              female[i] * (X_base[i] * b_base_female) +
              female[i] * use_ema * (theta[i] * c_ema_female) +
              use_ema * interact_term +
              var_emotional[i] * (b_var_base + female[i] * b_var_female) +
              use_ema * var_emotional[i] * b_var_ema;
    target += student_t_lpdf(diff_frag[i] | nu_diff, mu_frag, sigma_diff);
  }
}

generated quantities {
  vector[I] mu_frag_hat;
  real R2_frag;
  vector[I] y_rep;   // nuove fragilità replicate
  vector[N_items] y_item_rep;

  
  // Log-likelihood per osservazione (per LOO)
  vector[N_items] log_lik_item;
  vector[N_obs] log_lik_obs;
  for (n in 1:N_obs) log_lik_obs[n] = 0;

  // Likelihood items
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

  // Predizioni medie e R²
  for (i in 1:I) {
    real interact_term = 0.0;
    int idx = 1;
    for (d1 in 1:(D-1)) {
      for (d2 in (d1+1):D) {
        interact_term += b_interact[idx] * X_base[i, d1] * theta[i, d2];
        idx += 1;
      }
    }
    mu_frag_hat[i] = a_frag +
                     (X_base[i] * b_base) +
                     use_ema * (theta[i] * c_ema) +
                     female[i] * b_female +
                     female[i] * (X_base[i] * b_base_female) +
                     female[i] * use_ema * (theta[i] * c_ema_female) +
                     use_ema * interact_term +
                     var_emotional[i] * (b_var_base + female[i] * b_var_female) +
                     use_ema * var_emotional[i] * b_var_ema;
  }

  {
    real v_fit = variance(mu_frag_hat);
    real v_resid = (nu_diff > 2) 
                   ? (nu_diff / (nu_diff - 2)) * square(sigma_diff)
                   : positive_infinity();
    R2_frag = v_fit / (v_fit + v_resid);
  }

  // Genera dati replicati
  for (i in 1:I) {
    y_rep[i] = student_t_rng(nu_diff, mu_frag_hat[i], sigma_diff);
  }

  for (n in 1:N_items) {
    real m = (item_id[n] == 1 || item_id[n] == 3)
           ? u_pos[obs_id[n]] : u_neg[obs_id[n]];
    real lp = nu[item_id[n]] + lambda[item_id[n]] * eta[obs_id[n]] + m;
    y_item_rep[n] = ordered_logistic_rng(lp, tau[item_id[n]]);
  }

}
