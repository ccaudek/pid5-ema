// reactivity_latent_eta.stan
// Modello latente di reattività (soggetto x periodo) con predittori tra-soggetto W
// - Gestione dei missing in W dentro Stan (imputazione con prior std_normal su scala z)
// - Eta[i,p] = mu0[p] + W[i,*] * B[,p] + L_Sigma * z_i  (z_i ~ N(0, I_P))
// - Due outcome osservati (y1, y2) ~ Normal(alpha_k + lambda_k * eta, sigma_k)

data {
  int<lower=1> I;                                // numero soggetti
  int<lower=1> P;                                // numero periodi (es. 3: baseline, pre, post)
  int<lower=1> K;                                // numero outcome (qui atteso 2)
  int<lower=1> Q;                                // numero predittori tra-soggetto (colonne di W)

  // Osservazioni outcome 1
  array[2] int<lower=1> N;                       // N[1]=n_y1, N[2]=n_y2
  vector[N[1]] y1;
  array[N[1]] int<lower=1, upper=I> subj1;
  array[N[1]] int<lower=1, upper=P> period1;

  // Osservazioni outcome 2
  vector[N[2]] y2;
  array[N[2]] int<lower=1, upper=I> subj2;
  array[N[2]] int<lower=1, upper=P> period2;

  // ---- Predittori W con missing gestiti in Stan ----
  int<lower=0> NW_obs;                           // numero celle osservate in W
  array[NW_obs] int<lower=1, upper=I> w_obs_i;   // riga osservata
  array[NW_obs] int<lower=1, upper=Q> w_obs_j;   // colonna osservata
  vector[NW_obs] W_obs;                          // valore osservato

  int<lower=0> NW_mis;                           // numero celle mancanti
  array[NW_mis] int<lower=1, upper=I> w_mis_i;   // riga mancante
  array[NW_mis] int<lower=1, upper=Q> w_mis_j;   // colonna mancante
}

parameters {
  // Celle mancanti in W (assunte su scala z; prior std_normal)
  vector[NW_mis] W_mis;

  // Struttura latente di eta per periodo
  vector[P] mu0;                                 // intercette per periodo (media di eta)
  matrix[Q, P] B;                                // effetti dei predittori W su eta (Q x P)
  vector<lower=0>[P] tau_eta;                    // deviazioni standard per periodo
  cholesky_factor_corr[P] L_Omega;               // correlazione tra periodi (Cholesky)
  matrix[I, P] z_eta;                            // effetti soggetto (standardizzati) I x P

  // Modello di misura per gli outcome (K=2)
  vector[K] alpha;                               // intercette outcome
  vector[K] lambda;                              // loading outcome -> eta
  vector<lower=0>[K] sigma;                      // sd residuali outcome
}

transformed parameters {
  // Ricostruzione completa di W e calcolo di eta
  matrix[I, Q] W;
  matrix[I, P] eta;

  // Ricostruisci W (osservati + imputati)
  W = rep_matrix(0, I, Q);
  for (n in 1:NW_obs)
    W[w_obs_i[n], w_obs_j[n]] = W_obs[n];
  for (m in 1:NW_mis)
    W[w_mis_i[m], w_mis_j[m]] = W_mis[m];

  // Eta = mu0 + W * B + L_Sigma * z_i
  {
    matrix[P, P] L_Sigma = diag_pre_multiply(tau_eta, L_Omega);
    for (i in 1:I) {
      row_vector[Q] Wi = W[i, ];
      // mean per periodo: mu0' + Wi * B (QxP) -> row_vector[P]
      eta[i] = to_row_vector(mu0) + (Wi * B);
      // aggiungi deviazione soggetto (correlata tra periodi)
      eta[i] += to_row_vector( L_Sigma * to_vector(z_eta[i]) );
    }
  }
}

model {
  // ---- Priors ----
  // Celle mancanti di W (colonne di W sono z-standardizzate)
  W_mis ~ std_normal();

  // Regressione W -> eta
  to_vector(B) ~ normal(0, 0.5);
  mu0 ~ normal(0, 1);

  // Struttura random effects su periodi
  tau_eta ~ student_t(3, 0, 1);                  // half-t via vincolo <lower=0>
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(z_eta) ~ std_normal();               // i.i.d. N(0,1)

  // Misura
  alpha ~ normal(0, 1);
  lambda ~ normal(1, 0.5);                       // loading attorno a 1 (puoi usare normal(0,1) se preferisci)
  sigma ~ student_t(3, 0, 2);                    // half-t via <lower=0>

  // ---- Likelihood ----
  // Outcome 1
  for (n in 1:N[1]) {
    int i = subj1[n];
    int p = period1[n];
    y1[n] ~ normal( alpha[1] + lambda[1] * eta[i, p], sigma[1] );
  }
  // Outcome 2
  for (n in 1:N[2]) {
    int i = subj2[n];
    int p = period2[n];
    y2[n] ~ normal( alpha[2] + lambda[2] * eta[i, p], sigma[2] );
  }
}

generated quantities {
  // Correlazione tra periodi (per diagnosi/interpretazione)
  corr_matrix[P] Omega = multiply_lower_tri_self_transpose(L_Omega);

  // R^2 out-of-model (empirico) per ciascun outcome (versione semplice)
  real R2_y1;
  real R2_y2;
  {
    vector[N[1]] y1_hat;
    vector[N[2]] y2_hat;

    for (n in 1:N[1]) {
      int i = subj1[n];
      int p = period1[n];
      y1_hat[n] = alpha[1] + lambda[1] * eta[i, p];
    }
    for (n in 1:N[2]) {
      int i = subj2[n];
      int p = period2[n];
      y2_hat[n] = alpha[2] + lambda[2] * eta[i, p];
    }

    // Nota: varianza dei dati "as-is" (dipende dai dati passati; qui è una metrica utile per PPC veloci)
    real var_y1 = variance(y1);
    real var_y2 = variance(y2);
    R2_y1 = (var_y1 - variance(y1 - y1_hat)) / var_y1;
    R2_y2 = (var_y2 - variance(y2 - y2_hat)) / var_y2;
  }
}
