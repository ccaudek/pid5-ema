data {
  int<lower=1> I;                 // soggetti
  int<lower=1> N_obs;             // osservazioni EMA
  int<lower=1> K;                 // items = 4 (happy, sad, satisfied, angry)
  int<lower=1> P;                 // periodi = 3 (1=base,2=pre,3=post)
  int<lower=1> D;                 // dimensioni EMA a livello soggetto = 5
                                  // 1 = Negative Affect (latente dagli item)
                                  // 2..5 = det/ant/dis/psy (osservate con errore)

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
  array[M_ema] int<lower=2, upper=D> ema_dim;  // NB: 2..5
  array[M_ema] int<lower=1, upper=N_obs> ema_obs;

  // Predittori baseline (I x 5), già z-score
  matrix[I, D] X_base;

  // Interruttore per includere le EMA nel modello di fragilità
  real<lower=0, upper=1> use_ema;
}

parameters {
  // Misurazione ordinale per i 4 item
  vector[K] nu;
  vector[3] lambda_raw;                // happy, satisfied, angry (sad fissato a 1)
  array[K] ordered[6] tau;             // 7 categorie -> 6 soglie/item

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
  vector[I] delta_pre;                 // ~ normal(0, tau_pre)
  vector[I] delta_post;                // ~ normal(0, tau_post)

  // Variabilità intra-periodo (stato)
  real<lower=0> sigma_state;
  real<lower=1> nu_state;              // robusto
  vector[N_obs] eps_obs;

  // Tratti EMA per soggetto (5 dimensioni)
  matrix[I, D] theta;                  // theta[,1] = NA (latent dagli item)
  vector<lower=0>[D] sigma_ema;        // sd di misura per le EMA osservate (2..5)

  // Regressione della fragilità (pre - post)
  real a_frag;
  vector[D] b_base;                    // coeff. baseline (5D)
  vector[D] c_ema;                     // coeff. EMA (5D)
  real<lower=0> sigma_diff;            // sd residuo
  real<lower=1> nu_diff;               // df t-robusto

  // (opzionale) shrinkage globale sui c_ema
  real<lower=0> tau_ema;
}

transformed parameters {
  vector[K] lambda;
  vector[N_obs] eta;                   // stato latente NA per osservazioni
  vector[I] diff_frag;                 // fragilità = (pre) - (post)

  // carichi: sad positivo, happy/satisfied negativi, angry positivo
  lambda[2] = 1.0;                     // sad
  lambda[4] =  exp(lambda_raw[3]);     // angry  > 0
  lambda[1] = -exp(lambda_raw[1]);     // happy  < 0
  lambda[3] = -exp(lambda_raw[2]);     // satisfied < 0

  // NA osservazionale: soggetto-NA + effetti di periodo + stato
  for (n in 1:N_obs) {
    int i = subject[n];
    int p = period[n];
    real m_pre  = (p == 2) ? (gamma_pre  + delta_pre[i])  : 0.0;
    real m_post = (p == 3) ? (gamma_post + delta_post[i]) : 0.0;
    // theta[i,1] è il tratto NA del soggetto
    eta[n] = theta[i,1] + m_pre + m_post + sigma_state * eps_obs[n];
  }

  // Fragilità per soggetto sulla scala dei mezzi latenti
  for (i in 1:I)
    diff_frag[i] = (gamma_pre + delta_pre[i]) - (gamma_post + delta_post[i]);
}

model {
  // Priors misurazione
  nu ~ normal(0, 1);
  lambda_raw ~ normal(0, 0.3);
  for (k in 1:K) tau[k] ~ normal(0, 1);

  // Method effects
  sigma_pos ~ exponential(1);
  sigma_neg ~ exponential(1);
  u_pos ~ normal(0, sigma_pos);
  u_neg ~ normal(0, sigma_neg);

  // Periodo e dinamiche
  gamma_pre  ~ normal(0, 1);
  gamma_post ~ normal(0, 1);
  tau_pre  ~ exponential(1);
  tau_post ~ exponential(1);
  delta_pre  ~ normal(0, tau_pre);
  delta_post ~ normal(0, tau_post);

  sigma_state ~ exponential(1);
  nu_state ~ gamma(2, 0.1);
  eps_obs ~ student_t(nu_state, 0, 1);

  // Tratti EMA soggetto: priori moderati
  to_vector(theta) ~ normal(0, 1);
  sigma_ema ~ exponential(1);

  // Misurazione per EMA osservate (SOLO d in 2..5)
  for (m in 1:M_ema) {
    int i = subject[ ema_obs[m] ];
    int d = ema_dim[m];                 // 2..5
    ema_val[m] ~ normal(theta[i, d], sigma_ema[d]);
  }

  // Likelihood per item ordinali (NA)
  {
    vector[N_items] linpred;
    for (n in 1:N_items) {
      real m = (item_id[n] == 1 || item_id[n] == 3)
               ? u_pos[ obs_id[n] ] : u_neg[ obs_id[n] ];
      linpred[n] = nu[item_id[n]] + lambda[item_id[n]] * eta[ obs_id[n] ] + m;
    }
    for (n in 1:N_items)
      y_item[n] ~ ordered_logistic(linpred[n], tau[item_id[n]]);
  }

  // Regressione fragilità (una riga per soggetto)
  a_frag  ~ normal(0, 2);
  b_base  ~ normal(0, 0.5);
  tau_ema ~ exponential(1);
  c_ema   ~ normal(0, tau_ema);
  sigma_diff ~ exponential(1);
  nu_diff    ~ gamma(2, 0.1);

  for (i in 1:I) {
    real mu_frag = a_frag
                   + (X_base[i] * b_base)          // baseline (5D)
                   + use_ema * (theta[i] * c_ema); // EMA (5D: NA + altre 4)
    target += student_t_lpdf(diff_frag[i] | nu_diff, mu_frag, sigma_diff);
  }
}

generated quantities {
  vector[I] mu_frag_hat;
  real R2_frag;

  // log-lik al livello item e osservazione (per LOO a livello beep)
  vector[N_items] log_lik_item;
  vector[N_obs]   log_lik_obs;
  for (n in 1:N_obs) log_lik_obs[n] = 0;

  {
    vector[N_items] linpred_gq;
    for (n in 1:N_items) {
      real m = (item_id[n] == 1 || item_id[n] == 3)
               ? u_pos[ obs_id[n] ] : u_neg[ obs_id[n] ];
      linpred_gq[n] = nu[item_id[n]] + lambda[item_id[n]] * eta[ obs_id[n] ] + m;
      log_lik_item[n] = ordered_logistic_lpmf(y_item[n] | linpred_gq[n], tau[item_id[n]]);
      log_lik_obs[ obs_id[n] ] += log_lik_item[n];
    }
  }

  // Pred media per soggetto & R2 t-robusto
  for (i in 1:I)
    mu_frag_hat[i] = a_frag + (X_base[i] * b_base) + use_ema * (theta[i] * c_ema);

  {
    real v_fit   = variance(mu_frag_hat);
    real v_resid = (nu_diff > 2)
                   ? (nu_diff / (nu_diff - 2)) * square(sigma_diff)
                   : positive_infinity();
    R2_frag = v_fit / (v_fit + v_resid);
  }
}

