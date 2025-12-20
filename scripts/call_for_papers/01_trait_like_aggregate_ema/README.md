# pid5-ema

In questa cartella sono contenuti gli script per l'analisi dei dati relativi al progetto PID-5 EMA che riproducono i risultati dell'abstract che è stato accettato.

version: Sat Dec 20 10:01:40 2025

------

Lo scopo di questo progetto è di valutare se le misure "dinamiche" PID-5 somministrate ripetutamente nel tempo con la metodologia EMA sono più efficaci nell'individuare partecipanti che hanno un atteggiamento/comportamento disadattato/disfunzionale rispetto alla somministrazione "statica" del PID-5 una volta soltanto. La somministrazione "statica" del PID-5 include il test completo di 220 item (in realtà dovrei escludere gli item che uso per l'EMA); la somministrazione "dinamica" include un subset di 15 item del PID-5 (3 item per ciascuno dei 5 domini) -- tali item sono stati individuati con un'analisi fattoriale su un campione diverso (ma sempre studenti universitari) di partecipanti. 

I dati di cui dispongo sono i seguenti: 

[1] "user_id" "esi_bf" "domain_negative_affect" 
[4] "domain_detachment" "domain_antagonism" "domain_disinhibition" 
[7] "domain_psychoticism" "happy" "sad" 
[10] "satisfied" "angry" "pid5_13" 
[13] "pid5_15" "pid5_11" "pid5_3" 
[16] "pid5_2" "pid5_7" "pid5_14" 
[19] "pid5_6" "pid5_4" "pid5_12" 
[22] "pid5_1" "pid5_9" "pid5_5" 
[25] "pid5_8" "pid5_10" "tripm_1" 
[28] "tripm_3" "tripm_2" "tripm_4" 
[31] "dass21_2" "dass21_5" "dass21_3" 
[34] "dass21_4" "dass21_6" "dass21_1" 
[37] "scs3_pos" "scs5_neg" "scs7_pos" 
[40] "scs4_neg" "scs2_neg" "scs1_pos" 
[43] "scs8_neg" "scs6_pos" "vq_2" 
[46] "vq_1" "vq_4" "vq_3" 
[49] "cope_nvi_2" "cope_nvi_7" "cope_nvi_5" 
[52] "cope_nvi_8" "cope_nvi_9" "cope_nvi_10" 
[55] "cope_nvi_1" "cope_nvi_3" "cope_nvi_4" 
[58] "cope_nvi_6" "day" "hour" 
[61] "calendar_day" "bysubj_day" "context_quality" 
[64] "context_control" "context_support" "context_threat" 
[67] "pid5_sum" "pid5_negative_affectivity" "pid5_detachment" 
[70] "pid5_antagonism" "pid5_disinhibition" "pid5_psychoticism" 
[73] "ipv_sum" "tripm_4_rev" "tripm_sum" 
[76] "tripm_boldness" "tripm_meanness" "dass_sum" 
[79] "dass_stress" "dass_depression" "dass_anxiety" 
[82] "cope_10_rev" "cope_avoid" "cope_prob_or" 
[85] "cope_social_support" "cope_positive_att" "cope_trascendent_or" 
[88] "cs_pos" "ucs_neg" "neg_aff_ema" 
[91] "neg_aff_ema_c" 

Ci sono "domain_negative_affect"   "domain_detachment"         "domain_antagonism"         "domain_disinhibition"   "domain_psychoticism" che misurano i punteggi di ciascun individuo sui 5 domini PID-5 (e rimangono invarianti in ciascuna notifica EMA).

Le variabili "pid5_sum" "pid5_negative_affectivity" "pid5_detachment" "pid5_antagonism" "pid5_disinhibition" "pid5_psychoticism" rappresentano la somma di tutti i 15 PID-5 item usanti nell'EMA (pid5_sum), oppure la somma dei tre item di ciascun dominio ( "pid5_negative_affectivity" "pid5_detachment" "pid5_antagonism" "pid5_disinhibition" "pid5_psychoticism" ). 

Abbiamo 10 item del COPE, anche questi combinati nelle variabili cope_10_rev" "cope_avoid" "cope_prob_or" "cope_social_support" "cope_positive_att" "cope_trascendent_or" . 

Ci sono 4 item relativi al contesto: "context_quality"   "context_control"           "context_support"           "context_threat" 

Ci sono le variabili che misurano, negli EMA, il TRI-PM, somma "tripm_sum"   oppure  distinguendo tra  "tripm_boldness"            "tripm_meanness" 

Ci sono 4 item per l’umore istantaneo EMA: "happy"     "sad"   "satisfied"                 "angry"

Ci sono le misure EMA sul DASS-21 (solo alcuni item): "dass_sum"  oppure   "dass_stress"               "dass_depression"           "dass_anxiety" 

Ci sono le misure EMA del COPE (alcuni item): "cope_avoid"                "cope_prob_or"     "cope_social_support"       "cope_positive_att"         "cope_trascendent_or"

C’è la self compassion di stato, nelle sue due componenti: "cs_pos"                    "ucs_neg" 

Le misure si possono separare in 3 "finestre contestuali" che suscitano una emotività/reattività diversa. I partecipanti sono tutti studenti universitari. Il primo contesto temporale è lontano da un esame universitario. Il secondo contesto è "il giorno prima dell'esame". Il terzo contesto temporale è subito dopo l'esame (avendo terminato l'esame, ovvero il giorno dopo).

Essendo uno studio EMA, il data set contiene dati mancanti, quando un partecipante non ha risposto alla notifica. 



ESI-BF ---------------------------------------------------------------------

The Externalizing Spectrum Inventory-Brief Form (ESI-BF) is a 160-item questionnaire designed to assess the externalizing spectrum of psychopathology. It includes three subfactors: Disinhibition (ESI-BFDIS), Callous-Aggression (ESI-BFAGG), and Substance Abuse (ESI-BFSUB). The ESI-BF is used to evaluate tendencies towards impulsive, aggressive, and substance-related behaviors. Here's a more detailed breakdown: 

- ESI-BFDIS (Disinhibition): This subscale consists of 20 items assessing general externalizing tendencies, including problematic impulsivity, irresponsibility, and proneness to boredom.
- ESI-BFAGG (Callous-Aggression): This subscale contains 19 items that focus on traits like deficient empathy, relational aggression, and destructiveness.
- ESI-BFSUB (Substance Abuse): This subscale includes 18 items specifically targeting problematic and recreational substance use, including alcohol and other drugs.

The questionnaire uses a 4-point Likert scale for responses, ranging from "Not true at all" to "Completely true". The ESI-BF is valuable for research on the correlates and mechanisms of externalizing behaviors like risk-taking, delinquency, aggression, and substance abuse, according to research papers from the American Psychological Association. 

Nello studio, esi_bf è ESI-BFSUB.

---------------------------------------------------------------------

















Un modo per esaminare un atteggiamento disfunzionale è quello di considerare la componente negativa della self-compassion. 

La domanda che motiva il progetto di ricerca è se le misure “dinamiche” EMA (somministrate in molteplici occasioni) aggiungono potere predittivo rispetto a dimensioni disfunzionali degli atteggiamenti/comportamenti degli individui, al di là di ciò che viene offerto dalle misure “statiche” (somministrate una volta soltanto) corrispondenti ai 5 dominii del PID-5. Le misure basali del PID-5 sono state calcolate *escludendo* i 15 item (3 per dominio) usati nell’EMA.

Per testare l’ipotesi possiamo solo usare le misure a disposizione. Abbiamo un questionario Externalizing Spectrum Inventory BF (esi_bf), che è stato somministrato una volta soltanto. Questa sarebbe la variabile di esito più ovvia, ma sembra che il valore aggiunto delle componenti “dinamiche” del PID-5, rispetto al questionario somministrato una volta soltanto, emerge più chiaramente quando si vuole prevedere qualcosa che varia nel tempo. 

Le misure che variano nel tempo, misurate con l’EMA, sono elencate sopra. 

Una prima analisi è stata fatta confrontando questi due modelli. 

model_base <- brm(
  ucs_neg ~ neg_aff_ema + 
    domain_negative_affect + domain_detachment +
    domain_antagonism + domain_disinhibition + domain_psychoticism + 
    (1 + neg_aff_ema | user_id),
  data = df_self_comp_ema_scaled,
  family = skew_normal(),
  prior = c(
    prior(normal(0, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4,
  cores = 4,
  iter = 2000,
  seed = 123,
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE)
)

model_alt <- brm(
  ucs_neg ~
    (neg_aff_ema + domain_negative_affect + domain_detachment + 
       domain_antagonism + domain_disinhibition + domain_psychoticism) *
      (pid5_negative_affectivity + pid5_detachment + pid5_antagonism +
         pid5_disinhibition + pid5_psychoticism) +
    (1 + neg_aff_ema | user_id),
  data = df_self_comp_ema_scaled,
  family = skew_normal(),
  prior = c(
    prior(normal(0, 1), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4,
  cores = 4,
  iter = 2000,
  seed = 123,
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE)
)

La componente negativa della self-compassion di stato viene usata qui come proxy di uno stato disfunzionale del soggetto.

Una seconda analisi riguarda la correlazione tra la componente positiva e negativa della self-compassion di stato. Si assume che le persone ben adattate mostrano un equilibrio tra le due componenti (quando una componente è forte, l’altra è debole, quindi forte correlazione negativa). Invece, un’assenza di correlazione o correlazione di segno opposto potrebbe essere un indice di disadattamento. Anche in questo caso, userei due modelli, come sopra, con questa nuova variabile dipendente.

Una terza analisi riguarda la reattività allo stress, ovvero l’interazione tra il giudizio sul contesto (soprattutto il carattere di minaccia) e il negative affect istantaneo, calcolato mediante i 4 item relativi all’umore usati negli EMA. Anche in questo caso, confronterei due modelli, come sopra.

