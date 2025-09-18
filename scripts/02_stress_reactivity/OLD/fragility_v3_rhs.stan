data {
  int<lower=1> I;                 // soggetti
  int<lower=1> N_obs;             // osservazioni EMA
  int<lower=1> K;                 // items ordinali (4: happy, sad, satisfied, angry)
  int<lower=1> P;                 // periodi (1=base,2=pre,3=post)
  int<lower=1> D;                 // dimensioni EMA/trait (es. 5)

  // --- Item ordinali (long, non-mancanti)
  int<lower=1> N_items;
  array[N_items] int<lower=1, upper=7> y_item;
  array[N_items] int<lower=1, upper=K> item_id;
  array[N_items] int<lower=1, upper=N_obs> obs_id;

  // --- Mappa osservazione -> soggetto/periodo
  array[N_obs] int<lower=1, upper=I> subject;
  array[N_obs] int<lower=1, upper=P> period;

  // --- Misure EMA continue (subset, per ancorare theta)
  int<lower=0> M_ema;
  array[M_ema] real ema_val;
  array[M_ema] int<lower=2, upper=D> ema_dim;      // 2..D
  array[M_ema] int<lower=1, upper=N_obs> ema_obs;

  // --- Covariate baseline (z-score) per soggetto
  matrix[I, D] X_base;

  // Genere (0/1) e interruttore di uso
  array[I] int<lower=0, upper=1> female;
  int<lower=0, upper=1> use_gender;

  // Interazioni baseline x EMA
  int<lower=0, upper=1> use_interact;
  // N_interact deve essere D*(D-1)/2 se use_interact=1, altrimenti può essere 0
  int<lower=0> N_interact;

  // Interruttore per includere EMA come predittori nel modello di fragilità
  int<lower=0, upper=1> use_ema;
}

parameters {
  // ---- misurazione item ordinale ----
  vector[K] nu;
  vector[3] lambda_raw;
  array[K] ordered[6] tau;

  // fattori di metodo (non-centered)
  real<lower=0> sigma_pos;
  real<lower=0> sigma_neg;
  vector[N_obs] z_pos;
  vector[N_obs] z_neg;

  // effetti di periodo (non-centered)
  real gamma_pre;
  real gamma_post;
  real<lower=0> tau_pre;
  real<lower=0> tau_post;
  vector[I] z_delta_pre;
  vector[I] z_delta_post;

  // stato intra-osservazione (robusto)
  real<lower=0> sigma_state;
  real<lower=1> nu_state;
  vector[N_obs] eps_obs;

  // tratti EMA per soggetto (non-centered)
  matrix[I, D] z_theta;
  vector<lower=0>[D] sigma_ema;

  // ---- modello di fragilità (diff pre-post) ----
  real a_frag;

  // blocco baseline: Regularized Horseshoe
  // ----------------
  vector[D] b_base_raw;
  real<lower=0> hs_tau_base;
  vector<lower=0>[D] hs_lambda_base;
  real<lower=0> hs_c0_base;

  // blocco EMA: Regularized Horseshoe
  // ----------------
  vector[D] c_ema_raw;
  real<lower=0> hs_tau_ema;
  vector<lower=0>[D] hs_lambda_ema;
  real<lower=0> hs_c0_ema;

  // blocco interazioni baseline×EMA (opzionale)
  // ----------------
  vector[max(1, N_interact)] b_interact_raw;
  real<lower=0> hs_tau_int;
  vector<lower=0>[max(1, N_interact)] hs_lambda_int;
  real<lower=0> hs_c0_int;

  // moderazioni genere (opzionali) con shrinkage semplice
  real b_female;
  vector[D] b_base_female;
  vector[D] c_ema_female;

  // effetti della variabilità emotiva (shrinkage semplice)
  real b_var_base;
  real b_var_ema;
  real b_var_female;

  // residui robusti
  real<lower=0> sigma_diff;
  real<lower=1> nu_diff;
}

transformed parameters {
  // carichi
  vector[K] lambda;
  // fattori di metodo ricostruiti
  vector[N_obs] u_pos = sigma_pos * z_pos;
  vector[N_obs] u_neg = sigma_neg * z_neg;

  // period effects per soggetto
  vector[I] delta_pre  = tau_pre  * z_delta_pre;
  vector[I] delta_post = tau_post * z_delta_post;

  // tratti EMA per soggetto
  matrix[I, D] theta;
  for (d in 1:D)
    for (i in 1:I)
      theta[i, d] = sigma_ema[d] * z_theta[i, d];

  // stato latente per osservazione
  vector[N_obs] eta;
  for (n in 1:N_obs) {
    int i = subject[n];
    int p = period[n];
    real m_pre  = (p == 2) ? (gamma_pre  + delta_pre[i])  : 0.0;
    real m_post = (p == 3) ? (gamma_post + delta_post[i]) : 0.0;
    eta[n] = theta[i, 1] + m_pre + m_post + sigma_state * eps_obs[n];
  }

  // diff di fragilità (pre - post) osservabile
  vector[I] diff_frag;
  for (i in 1:I)
    diff_frag[i] = (gamma_pre + delta_pre[i]) - (gamma_post + delta_post[i]);

  // carichi con vincoli di segno (sad positivo, happy/satisfied negativi)
  lambda[2] = 1.0;
  lambda[4] =  exp(lambda_raw[3]);
  lambda[1] = -exp(lambda_raw[1]);
  lambda[3] = -exp(lambda_raw[2]);
}

model {
  // --- sanity interazioni ---
  if (use_interact == 1 && N_interact != D*(D-1)/2)
    reject("N_interact must be D*(D-1)/2 when use_interact=1.");

  // ---- priors: misurazione ----
  nu ~ normal(0, 1);
  lambda_raw ~ normal(0, 0.3);
  for (k in 1:K) tau[k] ~ normal(0, 1);

  sigma_pos ~ exponential(4);   // più stretti di v2
  sigma_neg ~ exponential(4);
  z_pos ~ normal(0, 1);
  z_neg ~ normal(0, 1);

  gamma_pre  ~ normal(0, 1);
  gamma_post ~ normal(0, 1);
  tau_pre  ~ exponential(3);
  tau_post ~ exponential(3);
  z_delta_pre  ~ normal(0, 1);
  z_delta_post ~ normal(0, 1);

  sigma_state ~ exponential(3);
  nu_state ~ gamma(2, 0.15);      // media ≈ 13
  eps_obs ~ student_t(nu_state, 0, 1);

  to_vector(z_theta) ~ normal(0, 1);
  sigma_ema ~ exponential(3);

  // --- misure EMA continue (anchor)
  for (m in 1:M_ema) {
    int i = subject[ema_obs[m]];
    int d = ema_dim[m];
    ema_val[m] ~ normal(theta[i, d], sigma_ema[d]);
  }

  // --- item ordinali
  {
    vector[N_items] linpred;
    for (n in 1:N_items) {
      real m = (item_id[n] == 1 || item_id[n] == 3) ? u_pos[obs_id[n]] : u_neg[obs_id[n]];
      linpred[n] = nu[item_id[n]] + lambda[item_id[n]] * eta[obs_id[n]] + m;
    }
    y_item ~ ordered_logistic(linpred, tau[item_id]);
  }

  // ---- priors: fragilità con Regularized Horseshoe ----
  // baseline
  {
    // Piironen & Vehtari (2017): global tau più stretta della classica
    hs_tau_base ~ cauchy(0, 0.5);
    hs_lambda_base ~ cauchy(0, 1);
    hs_c0_base ~ inv_gamma(2, 8);
    b_base_raw ~ normal(0, 1);
  }
  // ema
  {
    hs_tau_ema ~ cauchy(0, 0.5);
    hs_lambda_ema ~ cauchy(0, 1);
    hs_c0_ema ~ inv_gamma(2, 8);
    c_ema_raw ~ normal(0, 1);
  }
  // interazioni (se usate)
  {
    if (use_interact == 1) {
      hs_tau_int ~ cauchy(0, 0.5);
      hs_lambda_int[1:N_interact] ~ cauchy(0, 1);
      hs_c0_int ~ inv_gamma(2, 8);
      b_interact_raw[1:N_interact] ~ normal(0, 1);
    } else {
      // tieni il parametro definito ma non usato
      hs_tau_int ~ cauchy(0, 0.5);
      hs_lambda_int ~ cauchy(0, 1);
      hs_c0_int ~ inv_gamma(2, 8);
      b_interact_raw ~ normal(0, 1);
    }
  }

  // moderazioni e varianza emotiva (shrunk)
  a_frag ~ normal(0, 1.2);
  b_female ~ normal(0, 0.4);
  b_base_female ~ normal(0, 0.3);
  c_ema_female ~ normal(0, 0.3);
  b_var_base ~ normal(0, 0.3);
  b_var_ema ~ normal(0, 0.3);
  b_var_female ~ normal(0, 0.3);

  sigma_diff ~ exponential(3);
  nu_diff ~ gamma(2, 0.3);   // media ~6.7 -> più robusto di v2
}

generated quantities {
  // --- ricostruzione horseshoe (lambda tilde) e coefficienti effettivi ---
  vector[D] b_base;
  vector[D] c_ema;
  vector[max(1, N_interact)] b_interact;
  {
    vector[D] lambda_tilde_base;
    vector[D] lambda_tilde_ema;
    for (d in 1:D) {
      lambda_tilde_base[d] = sqrt( square(hs_c0_base) * square(hs_lambda_base[d]) /
                                   ( square(hs_c0_base) + square(hs_tau_base) * square(hs_lambda_base[d]) ) );
      lambda_tilde_ema[d]  = sqrt( square(hs_c0_ema) * square(hs_lambda_ema[d]) /
                                   ( square(hs_c0_ema) + square(hs_tau_ema) * square(hs_lambda_ema[d]) ) );
      b_base[d] = b_base_raw[d] * lambda_tilde_base[d] * hs_tau_base;
      c_ema[d]  = c_ema_raw[d]  * lambda_tilde_ema[d]  * hs_tau_ema;
    }
    if (use_interact == 1) {
      for (l in 1:N_interact) {
        real lt = sqrt( square(hs_c0_int) * square(hs_lambda_int[l]) /
                        ( square(hs_c0_int) + square(hs_tau_int) * square(hs_lambda_int[l]) ) );
        b_interact[l] = b_interact_raw[l] * lt * hs_tau_int;
      }
    } else {
      b_interact[1] = 0; // placeholder non usato
    }
  }

  // --- varianza emotiva intra-soggetto ---
  vector[I] var_emotional;
  {
    // media per soggetto
    vector[I] mu_i;
    for (i in 1:I) mu_i[i] = 0;
    for (n in 1:N_obs) mu_i[ subject[n] ] += eta[n];
    for (i in 1:I) {
      int count_i = 0;
      for (n in 1:N_obs) if (subject[n]==i) count_i += 1;
      if (count_i > 0) mu_i[i] /= count_i; else mu_i[i] = 0;
    }
    // varianza
    for (i in 1:I) {
      int count_i = 0;
      real ss = 0;
      for (n in 1:N_obs) if (subject[n]==i) { ss += square(eta[n] - mu_i[i]); count_i += 1; }
      var_emotional[i] = (count_i > 1) ? ss / (count_i - 1) : 0.0;
    }
  }

  // --- predizione media fragilità e log-lik per soggetto ---
  vector[I] mu_frag_hat;
  vector[I] log_lik_subj;

  // --- (facoltativo) log-lik per osservazione su item, utile solo a fini diagnostici ---
  vector[N_obs] log_lik_obs;
  for (n in 1:N_obs) log_lik_obs[n] = 0;

  {
    // helper per prodotto incrociato baseline×EMA (solo se attivo)
    int idx;
    for (i in 1:I) {
      real interact_term = 0.0;
      if (use_interact == 1) {
        idx = 1;
        for (d1 in 1:(D-1))
          for (d2 in (d1+1):D) {
            interact_term += b_interact[idx] * X_base[i, d1] * theta[i, d2];
            idx += 1;
          }
      }
      real mu =
        a_frag +
        dot_product( row(X_base, i), b_base ) +
        (use_ema == 1 ? dot_product( row(theta, i), c_ema ) : 0) +
        (use_gender == 1 ? female[i] * b_female : 0) +
        (use_gender == 1 ? female[i] * dot_product( row(X_base, i), b_base_female ) : 0) +
        (use_gender == 1 && use_ema == 1 ? female[i] * dot_product( row(theta, i), c_ema_female ) : 0) +
        (use_interact == 1 && use_ema == 1 ? interact_term : 0) +
        var_emotional[i] * ( b_var_base + (use_gender == 1 ? female[i] * b_var_female : 0) ) +
        (use_ema == 1 ? var_emotional[i] * b_var_ema : 0);

      mu_frag_hat[i] = mu;
      log_lik_subj[i] = student_t_lpdf( (gamma_pre + delta_pre[i]) - (gamma_post + delta_post[i])
                                        | nu_diff, mu, sigma_diff );
    }

    // log-lik per osservazione (sommando item): diagnostico, NON usarlo per il confronto
    for (n in 1:N_items) {
      real m = (item_id[n]==1 || item_id[n]==3) ? u_pos[obs_id[n]] : u_neg[obs_id[n]];
      real lp = nu[item_id[n]] + lambda[item_id[n]] * eta[obs_id[n]] + m;
      log_lik_obs[ obs_id[n] ] += ordered_logistic_lpmf(y_item[n] | lp, tau[item_id[n]]);
    }
  }
  
    // ---- LOO marginalizzato per soggetto (integra i random-effects) ----
  int S_prior = 8;                 // n. campioni prior-condizionati
  vector[I] log_lik_subj_marg;
  for (i in 1:I) log_lik_subj_marg[i] = negative_infinity();

  for (i in 1:I) {
    // accumulator sullo spazio log
    real lse;
    lse = negative_infinity();

    for (s in 1:S_prior) {
      // 1) campiona effetti i dal loro prior condizionale (non centrato)
      //    NB: NON usare i delta_pre/post e theta[i,*] "allenati"!
      real delta_pre_i;
      real delta_post_i;
      row_vector[D] theta_i;

      delta_pre_i  = tau_pre  * normal_rng(0, 1);
      delta_post_i = tau_post * normal_rng(0, 1);
      for (d in 1:D) theta_i[d] = sigma_ema[d] * normal_rng(0, 1);

      // 2) costruisci mu_frag per questo draw prior-condizionato
      real interact_term;
      interact_term = 0.0;
      if (use_interact == 1) {
        int idx;
        idx = 1;
        for (d1 in 1:(D-1)) {
          for (d2 in (d1+1):D) {
            interact_term += b_interact[idx] * X_base[i, d1] * theta_i[d2];
            idx += 1;
          }
        }
      }

      real mu;
      mu =
        a_frag +
        dot_product( row(X_base, i), b_base ) +
        (use_ema == 1 ? dot_product( theta_i, c_ema ) : 0) +
        (use_gender == 1 ? female[i] * b_female : 0) +
        (use_gender == 1 ? female[i] * dot_product( row(X_base, i), b_base_female ) : 0) +
        (use_gender == 1 && use_ema == 1 ? female[i] * dot_product( theta_i, c_ema_female ) : 0) +
        (use_interact == 1 && use_ema == 1 ? interact_term : 0);

      // 3) var_emotional marginalizzata: per ora usiamo 0 (prudente).
      //    (Se vuoi, puoi sostituire con una piccola MC per eta_i e calcolare la varianza.)
      real mu_var_i;
      mu_var_i = 0;

      real mu_i;
      mu_i = mu +
             mu_var_i * ( b_var_base + (use_gender == 1 ? female[i] * b_var_female : 0) ) +
             (use_ema == 1 ? mu_var_i * b_var_ema : 0);

      // 4) fragilità osservabile come (pre - post): usa i delta_i campionati qui
      real diff_i;
      diff_i = (gamma_pre + delta_pre_i) - (gamma_post + delta_post_i);

      // 5) contribuzione di likelihood (student-t robusta) per questo draw
      {
        real ll_i;
        ll_i = student_t_lpdf(diff_i | nu_diff, mu_i, sigma_diff);
        lse = log_sum_exp(lse, ll_i);
      }
    }

    // media sui S campioni
    log_lik_subj_marg[i] = lse - log(S_prior);
  }

}
