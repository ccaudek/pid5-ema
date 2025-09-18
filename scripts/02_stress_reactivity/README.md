# PID-5 EMA

In questo progetto, abbiamo raccolto intensive longitudinal data are by using ecological momentary assessment (EMA) method. 

Il progetto è durato circa 2 mesi e mezzo, e le notifiche EMA venivano fornite due volte la settimana.

Abbiamo anche somministrato dei questionari di baseline all'inizio del progetto. Il PID-5 è stato somministrato inegralmente a inizio progetto e i punteggi sulle sue sottoscale sono: "pid5_negative_affect_baseline" "pid5_detachment_baseline"  "pid5_antagonism_baseline"      "pid5_disinhibition_baseline"  "pid5_psychoticism_baseline" 

Il DASS-21 è stato somministrato integralmente a inizio progetto:
"dass_stress_baseline"          "dass_anxiety_baseline"        "dass_depression_baseline"

Il ESI-BF (uso si sostanze) è stato somministrato integralmente a inizio progetto:
esi_bf_baseline

Il giorno della somministrazione EMA è nella colonna:
day

L'ora della somministrazione è nella colonna
hour

Se utile, la colonna bysubj_day riporta il numero sequenziale delle notiriche all'interno di ciascun soggetto (i soggetti non hanno risposto sempre alle notifiche)

Per alcuni soggetti, in questo arco temporale, ci sono stati degli eventi stressogeni: i soggetti erano studenti universitari e, per alcuni soggetti, l'arco temporale del progetto, includeva delle giornate in cui si svolgeva un esame. La variabile exam_period descrive questo aspetto e ha tre modalità:

- baseline -- non ci sono eventi stressanti (esami) in prossimità della rilevazione EMA;
- pre_exam: la rilevazione EMA è avvenuta il giorno prima di un esame; questo dovrebbe essere il momento di maggiore stress;
- post_exam: la rilevazione EMA è avvenuta la sera del giorno di un esame (a punteggi esame conosciuti, almento per la maggioranza dei soggetti); questo dovrebbe essere un momento in cui lo stress è diminuito, rispetto a pre_exam, ma non è detto che sia ritornato a livello baseline.

È importante notare che gli effetti dello stress (esame) risultano più marcato per le femmine che per i maschi. Il genere è codificato nella variabile: 
sex

Sono state misurate varie dimensioni psicologiche mediante le notifiche EMA:

cs_pos, ucs_neg: rilevazioe EMA (ovvero, ripetuta nel tempo) delle due componenti della self-compassion scale di stato.

Componenti del COPE in versione EMA (misurazione ripetuta nel tempo). Sono stati usati alcuni item del COPE per presentarli sul cellulare durante le notifiche EMA:
cope_avoid      
cope_prob_or      
cope_social_support         
cope_positive_att           
cope_trascendent_or      

Componenti della DASS-21 in versione EMA: sono stati somministrati alcuni item per ciascuna sottoscala DASS-21 per misurare:
dass_stress, dass_depression, dass_anxiety 

Componenti del Tri-PM in versione EMA: Sono stati selezionati alcuni item del Tri-PM (psicopatia) per misurare le dimensioni:
tripm_boldness
tripm_meanness

Intimate Partner Violence in versione EMA:
ipv_sum 

Componenti del PID-5 in versione EMA: Dal PID-5 sono stati selezionati alcuni item rappresentativi per ciascuna dimensione e, in base a questi sono state calcolati i punteggi EMA seguenti:
pid5_negative_affectivity, pid5_detachment, pid5_antagonism, pid5_disinhibition,
pid5_psychoticism

Sono stati misurati, nelle notifiche EMA, 4 item relativi alla valutazione del contesto:
context_quality, context_control, context_support, context_threat

Sono stati utilizzati 4 item (2 positivi e 2 negativi) per misurare lo stato emotivo (negative affect) in ciascuna rilevazione EMA: happy, sad, satisfied, e angry.

Nel data frame dei dati grezzi, le variabili hanno questi nomi:
> glimpse(d)
Rows: 4,825
Columns: 111
$ happy                         <int> 100, 85, 89, 71, 93, 88, 83, 100, 85, 100, 70, 82, 26, 92, …
$ sad                           <int> 0, 11, 15, 16, 17, 0, 0, 0, 8, 0, 1, 0, 79, 2, 0, 0, 59, 5,…
$ satisfied                     <int> 64, 82, 32, 54, 40, 72, 73, 87, 49, 76, 57, 66, 39, 63, 73,…
$ angry                         <int> 16, 0, 100, 0, 67, 61, 50, 29, 66, 25, 57, 0, 100, 0, 0, 0,…
$ pid5_13                       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0,…
$ pid5_15                       <int> 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 2, 0, 0,…
$ pid5_11                       <int> 0, 1, 2, 1, 2, 2, 3, 2, 3, 2, 3, 2, 2, 1, 1, 3, 2, 2, 3, 3,…
$ pid5_3                        <int> 0, 2, 3, 2, 3, 2, 3, 3, 3, 3, 3, 2, 2, 3, 0, 2, 2, 2, 2, 2,…
$ pid5_2                        <int> 0, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 2, 2, 3, 0, 3, 2, 2, 2,…
$ pid5_7                        <int> 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,…
$ pid5_14                       <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1,…
$ pid5_6                        <int> 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0,…
$ pid5_4                        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ pid5_12                       <int> 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,…
$ pid5_1                        <int> 0, 0, 0, 0, 1, 0, 0, 1, 2, 2, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1,…
$ pid5_9                        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,…
$ pid5_5                        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0,…
$ pid5_8                        <int> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,…
$ pid5_10                       <int> 2, 2, 2, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 2, 1, 0, 2, 0, 3, 0,…
$ tripm_1                       <int> 3, 3, 3, 2, 1, 1, 3, 3, 1, 3, 3, 3, 1, 3, 3, 3, 2, 3, 3, 4,…
$ tripm_3                       <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1,…
$ tripm_2                       <int> 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1,…
$ tripm_4                       <int> 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 1, 4, 3, 4, 4, 3,…
$ dass21_2                      <int> 0, NA, 1, NA, NA, NA, 1, 1, NA, NA, 0, 1, 2, NA, 1, NA, 2, …
$ dass21_5                      <int> 0, NA, 0, NA, NA, NA, 0, 0, NA, NA, 0, 3, 0, NA, 0, NA, 1, …
$ dass21_3                      <int> 0, NA, 0, NA, NA, NA, 1, 1, NA, NA, 0, 0, 1, NA, 2, NA, 1, …
$ dass21_4                      <int> 0, NA, 0, NA, NA, NA, 1, 0, NA, NA, 2, 0, 0, NA, 1, NA, 0, …
$ dass21_6                      <int> 0, NA, 0, NA, NA, NA, 1, 0, NA, NA, 0, 0, 0, NA, 1, NA, 0, …
$ dass21_1                      <int> 0, NA, 2, NA, NA, NA, 1, 0, NA, NA, 0, 0, 0, NA, 2, NA, 0, …
$ scs3_pos                      <int> 1, NA, 3, NA, NA, NA, 2, 2, NA, NA, 2, 2, 2, NA, -2, NA, 2,…
$ scs5_neg                      <int> -3, NA, -3, NA, NA, NA, 2, -3, NA, NA, -2, -3, -3, NA, -1, …
$ scs7_pos                      <int> 2, NA, 2, NA, NA, NA, 1, 2, NA, NA, 1, 2, 2, NA, -3, NA, 1,…
$ scs4_neg                      <int> 1, NA, 1, NA, NA, NA, 3, 3, NA, NA, 2, 2, 2, NA, 2, NA, 2, …
$ scs2_neg                      <int> -3, NA, -2, NA, NA, NA, -1, -1, NA, NA, -2, -2, -1, NA, 2, …
$ scs1_pos                      <int> 2, NA, 2, NA, NA, NA, 1, 2, NA, NA, 2, 2, 1, NA, -3, NA, 1,…
$ scs8_neg                      <int> 2, NA, 2, NA, NA, NA, 3, 3, NA, NA, 2, 2, 2, NA, 2, NA, 2, …
$ scs6_pos                      <int> -3, NA, 1, NA, NA, NA, 2, -1, NA, NA, 2, 2, -1, NA, -2, NA,…
$ vq_2                          <int> NA, 0, NA, 0, 0, 0, NA, NA, 0, 0, 0, 0, NA, 0, NA, 0, NA, N…
$ vq_1                          <int> NA, 0, NA, 0, 0, 0, NA, NA, 0, 0, 0, 0, NA, 0, NA, 0, NA, N…
$ vq_4                          <int> NA, 0, NA, 0, 0, 0, NA, NA, 0, 0, 0, 0, NA, 0, NA, 0, NA, N…
$ vq_3                          <int> NA, 0, NA, 0, 0, 0, NA, NA, 0, 0, 0, 0, NA, 0, NA, 0, NA, N…
$ cope_nvi_2                    <int> NA, 1, NA, 2, 2, 2, NA, NA, 2, 1, 1, 2, NA, 1, NA, 1, NA, N…
$ cope_nvi_7                    <int> NA, 3, NA, 2, 1, 2, NA, NA, 3, 2, 2, 2, NA, 2, NA, 2, NA, N…
$ cope_nvi_5                    <int> NA, 3, NA, 4, 3, 2, NA, NA, 4, 4, 4, 4, NA, 4, NA, 3, NA, N…
$ cope_nvi_8                    <int> NA, 3, NA, 3, 3, 2, NA, NA, 3, 3, 3, 3, NA, 3, NA, 3, NA, N…
$ cope_nvi_9                    <int> NA, 3, NA, 2, 1, 2, NA, NA, 2, 2, 2, 2, NA, 2, NA, 2, NA, N…
$ cope_nvi_10                   <int> NA, 2, NA, 2, 2, 3, NA, NA, 2, 3, 3, 3, NA, 3, NA, 3, NA, N…
$ cope_nvi_1                    <int> NA, 2, NA, 1, 4, 2, NA, NA, 1, 2, 2, 2, NA, 2, NA, 2, NA, N…
$ cope_nvi_3                    <int> NA, 3, NA, 2, 2, 2, NA, NA, 2, 2, 3, 3, NA, 3, NA, 3, NA, N…
$ cope_nvi_4                    <int> NA, 3, NA, 3, 3, 2, NA, NA, 2, 2, 3, 2, NA, 2, NA, 3, NA, N…
$ cope_nvi_6                    <int> NA, 3, NA, 3, 3, 4, NA, NA, 3, 4, 4, 3, NA, 4, NA, 4, NA, N…
$ user_id                       <chr> "ad_pa_2006_01_27_630", "ad_pa_2006_01_27_630", "ad_pa_2006…
$ day                           <IDate> 2025-03-12, 2025-03-15, 2025-03-19, 2025-03-22, 2025-03-2…
$ hour                          <int> 18, 9, 18, 21, 18, 9, 20, 12, 18, 10, 19, 21, 10, 9, 19, 19…
$ calendar_day                  <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, …
$ bysubj_day                    <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, …
$ context_quality               <int> 1, -1, -1, 0, 0, 2, 1, 2, 0, 2, 1, 1, 1, 2, 2, 1, -1, 1, 2,…
$ context_control               <int> 4, 4, 0, 3, 3, 4, 5, 5, 4, 4, 4, 4, 4, 5, 5, 4, 1, 5, 5, 5,…
$ context_support               <int> 4, 4, 5, 3, 4, 4, 5, 5, 4, NA, 5, NA, NA, NA, NA, 3, 1, NA,…
$ context_threat                <int> 1, 2, 4, 2, 3, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,…
$ pid5_sum                      <int> 3, 8, 10, 6, 11, 6, 8, 8, 12, 9, 9, 9, 6, 12, 16, 8, 9, 8, …
$ pid5_negative_affectivity     <int> 0, 4, 6, 5, 7, 4, 5, 6, 7, 7, 5, 4, 4, 6, 4, 2, 5, 4, 4, 5,…
$ pid5_detachment               <int> 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 1, 0, 1, 4, 0, 0, 0, 0, 0,…
$ pid5_antagonism               <int> 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0,…
$ pid5_disinhibition            <int> 2, 3, 4, 1, 4, 2, 3, 2, 3, 2, 3, 4, 2, 5, 2, 3, 4, 2, 6, 3,…
$ pid5_psychoticism             <int> 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 3, 0, 2, 1, 1,…
$ ipv_sum                       <int> NA, 0, NA, 0, 0, 0, NA, NA, 0, 0, 0, 0, NA, 0, NA, 0, NA, N…
$ tripm_4_rev                   <int> 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 4, 1, 2, 1, 1, 2,…
$ tripm_sum                     <int> 9, 9, 9, 7, 7, 7, 9, 9, 7, 9, 8, 8, 6, 9, 9, 9, 7, 9, 9, 9,…
$ tripm_boldness                <int> 5, 4, 4, 3, 2, 2, 4, 4, 2, 4, 4, 4, 2, 4, 6, 4, 3, 4, 4, 5,…
$ tripm_meanness                <int> 3, 2, 2, 3, 2, 2, 2, 2, 2, 2, 3, 3, 3, 2, 6, 2, 3, 2, 2, 3,…
$ dass_sum                      <int> 0, NA, 3, NA, NA, NA, 5, 2, NA, NA, 2, 4, 3, NA, 7, NA, 4, …
$ dass_stress                   <int> 0, NA, 3, NA, NA, NA, 2, 1, NA, NA, 0, 1, 2, NA, 3, NA, 2, …
$ dass_depression               <int> 0, NA, 0, NA, NA, NA, 2, 1, NA, NA, 2, 0, 1, NA, 3, NA, 1, …
$ dass_anxiety                  <int> 0, NA, 0, NA, NA, NA, 1, 0, NA, NA, 0, 3, 0, NA, 1, NA, 1, …
$ cope_10_rev                   <int> NA, 3, NA, 3, 3, 2, NA, NA, 3, 2, 2, 2, NA, 2, NA, 2, NA, N…
$ cope_avoid                    <int> NA, 3, NA, 3, 6, 4, NA, NA, 3, 3, 3, 4, NA, 3, NA, 3, NA, N…
$ cope_prob_or                  <int> NA, 6, NA, 5, 5, 4, NA, NA, 4, 4, 6, 5, NA, 5, NA, 6, NA, N…
$ cope_social_support           <int> NA, 6, NA, 7, 6, 6, NA, NA, 7, 8, 8, 7, NA, 8, NA, 7, NA, N…
$ cope_positive_att             <int> NA, 6, NA, 5, 4, 4, NA, NA, 6, 5, 5, 5, NA, 5, NA, 5, NA, N…
$ cope_trascendent_or           <int> NA, 6, NA, 5, 4, 4, NA, NA, 5, 4, 4, 4, NA, 4, NA, 4, NA, N…
$ cs_pos                        <int> 2, NA, 8, NA, NA, NA, 6, 5, NA, NA, 7, 8, 4, NA, -10, NA, 3…
$ ucs_neg                       <int> -3, NA, -2, NA, NA, NA, 7, 2, NA, NA, 0, -1, 0, NA, 5, NA, …
$ rosenberg_score_baseline      <int> 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,…
$ self_kindness_baseline        <int> 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,…
$ common_humanity_baseline      <int> 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,…
$ mindfulness_baseline          <int> 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,…
$ self_judgment_baseline        <int> 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,…
$ isolation_baseline            <int> 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,…
$ over_identification_baseline  <int> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,…
$ scs_total_score_baseline      <int> 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62,…
$ boldness_baseline             <int> 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27,…
$ meanness_baseline             <int> 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,…
$ dass_stress_baseline          <int> 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,…
$ dass_anxiety_baseline         <int> 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,…
$ dass_depression_baseline      <int> 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,…
$ pid5_negative_affect_baseline <int> 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19,…
$ pid5_detachment_baseline      <int> 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,…
$ pid5_antagonism_baseline      <int> 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,…
$ pid5_disinhibition_baseline   <int> 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27,…
$ pid5_psychoticism_baseline    <int> 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38,…
$ esi_bf_baseline               <int> 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,…
$ inhibitory_ius_baseline       <int> 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,…
$ prospective_ius_baseline      <int> 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,…
$ ius_tot_baseline              <int> 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,…
$ .is_ema                       <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,…
$ n_ema                         <int> 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31,…
$ course                        <chr> "Psicometria", "Psicometria", "Psicometria", "Psicometria",…
$ sex                           <chr> "Femmina", "Femmina", "Femmina", "Femmina", "Femmina", "Fem…
$ date_meta                     <IDate> 2025-04-01, 2025-04-01, 2025-04-01, 2025-04-01, 2025-04-0…
$ exam_period                   <chr> "baseline", "baseline", "baseline", "baseline", "baseline",…

La domanda della ricerca è di stabilire se la "fragilità psicologica" degli studenti, operazionalizzata dalla loro reazione emotiva all'esame, sia meglio misurata dalle caratteristiche di tratto, fornite dal PID-5 somministrato integralmente una volta soltanto, all'inizio del progetto, o da un piccolo sottoinsieme di item del PID-5 somministrate molte volte nel tempo. I compositi basati sugli item PID-5 somministrati durane le notifiche EMA sono:

pid5_negative_affectivity"     "pid5_detachment"   "pid5_antagonism"   "pid5_disinhibition"  "pid5_psychoticism"

Per operazionalizzare la fragilità psicologica (maggiore differenza nel negative affect tra pre-test e post-test) utilizzerò i 4 item utilizzati nelle notifiche EMA che misurano lo stato emotivo: happy, sad, satisfied, e angry.

In un’EMA (circa 170 soggetti, ≈10 timepoint ciascuno, molti NA), forse è opportuno passare dal composito “somma di 4 item” a un punteggio latente stimato congiuntamente al modello strutturale, farlo *jointly* nello stesso modello Stan, non in due fasi. 

I dati di input al modello Stan sono stati creati con questo script:

# 02_clean_after_merge.R
# Data cleaning post-merge EMA + baseline
# - Esclusione careless responding (lista a priori)
# - Filtro numerosità EMA: 5 <= n_ema <= 40
# - Import metadati (genere) + tag esami (exam_period)
# - QA essenziale + export dataset finale e log esclusioni

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(lubridate)
  library(readr)
  library(rio)
})

# ============================
# 0) PERCORSI / INPUT
# ============================

# File "raw totale" già salvato dallo script precedente (RDS/CSV)
# -> Se esistono entrambi, preferisco RDS
INPUT_RDS <- here::here("data", "processed", "ema_plus_scales_merged.RDS")
INPUT_CSV <- here::here("data", "processed", "ema_plus_scales_merged.csv")

# Metadati (almeno: user_id, genere) — adatta il percorso/nome file
META_PATH <- here::here("data", "raw", "meta", "all_combined_sex_NEW_1.xlsx")
# deve contenere colonne: user_id, genere

# Tag esami (finestre) — adatta il percorso/nome file
# formato atteso: start_date, end_date, exam_period (YYYY-MM-DD)
EXAM_TAGS_PATH <- here::here("data", "raw", "meta", "exam_periods.csv")

# Output
OUT_CLEAN_RDS <- here::here("data", "processed", "ema_plus_scales_cleaned.rds")
OUT_CLEAN_CSV <- here::here("data", "processed", "ema_plus_scales_cleaned.csv")
OUT_LOG_EXCL <- here::here("data", "processed", "ema_exclusion_log.csv")
OUT_QA_TXT <- here::here("data", "processed", "ema_cleaning_QA.txt")

# ============================
# 1) LETTURA RAW TOTALE
# ============================

dat_raw <- if (file.exists(INPUT_RDS)) {
  readRDS(INPUT_RDS)
} else if (file.exists(INPUT_CSV)) {
  suppressMessages(readr::read_csv(INPUT_CSV, show_col_types = FALSE))
} else {
  stop(
    "Non trovo il file raw totale (né RDS né CSV). Aggiorna INPUT_RDS / INPUT_CSV."
  )
}

# Normalizzo alcuni nomi tipici (se presenti)
dat_raw <- dat_raw %>%
  rename(
    user_id = any_of(c("user_id", "UserID", "id_anon")),
    date = any_of(c("date", "Date", "baseline_date")),
    ema_time = any_of(c("ema_time", "timestamp", "time", "datetime"))
  )

# ============================
# 2) LISTA CARELESS A PRIORI
#    (presa dallo script allegato)
# ============================

careless_ids <- c(
  "ma_se_2005_11_14_490",
  "reve20041021036",
  "di_ma_2005_10_20_756",
  "pa_sc_2005_09_10_468",
  "il_re_2006_01_18_645",
  "so_ma_2003_10_13_804",
  "lo_ca_2005_05_07_05_437",
  "va_ma_2005_05_31_567",
  "no_un_2005_06_29_880",
  "an_bo_1988_08_24_166",
  "st_ma_2004_04_21_426",
  "an_st_2005_10_16_052",
  "vi_de_2002_12_30_067",
  "gi_ru_2005_03_08_033",
  "al_mi_2005_03_05_844",
  "la_ma_2006_01_31_787",
  "gi_lo_2004_06_27_237",
  "ch_bi_2001_01_28_407",
  "al_pe_2001_04_20_079",
  "le_de_2003_09_05_067",
  "fe_gr_2002_02_19_434",
  "ma_ba_2002_09_09_052",
  "ca_gi_2003_09_16_737",
  "an_to_2003_08_06_114",
  "al_se_2003_07_28_277",
  "ja_tr_2002_10_06_487",
  "el_ci_2002_02_15_057",
  "se_ti_2000_03_04_975",
  "co_ga_2003_10_29_614",
  "al_ba_2003_18_07_905",
  "bi_ro_2003_09_07_934",
  "an_va_2004_04_08_527",
  "ev_cr_2003_01_27_573"
)

# ============================
# 3) CALCOLO n_ema (robusto)
# ============================

# ============================
# 3) CALCOLO n_ema (robusto)
# ============================

if ("n_ema" %in% names(dat_raw)) {
  n_ema_tbl <- dat_raw %>% dplyr::distinct(user_id, n_ema)
} else {
  # colonne candidate a indicare una riga EMA
  ind_candidates <- c("ema_time", "day", "beep", "ema_id", "ema_wave")
  present_inds <- intersect(ind_candidates, names(dat_raw))

  # matrice logica: per ogni colonna presente TRUE se non NA
  flag_mat <- NULL
  if (length(present_inds) > 0) {
    flag_mat <- sapply(present_inds, function(nm) !is.na(dat_raw[[nm]]))
    # sapply con 1 sola colonna restituisce un vettore: forziamo matrice
    if (!is.matrix(flag_mat)) flag_mat <- matrix(flag_mat, ncol = 1)
  }

  # flag da "source" se esiste ed è testuale con "ema" dentro
  source_flag <- if ("source" %in% names(dat_raw)) {
    grepl("ema", dat_raw[["source"]], ignore.case = TRUE)
  } else {
    rep(FALSE, nrow(dat_raw))
  }

  # combina gli indicatori senza mai referenziare colonne mancanti
  if (is.null(flag_mat)) {
    is_ema <- source_flag
  } else {
    is_ema <- (rowSums(flag_mat) > 0) | source_flag
  }

  dat_raw$.is_ema <- is_ema

  n_ema_tbl <- dat_raw %>%
    dplyr::group_by(user_id) %>%
    dplyr::summarise(n_ema = sum(.is_ema, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(n_ema = as.integer(n_ema))
}

# Attacca n_ema al dataset principale (se non c'è già)
if (!"n_ema" %in% names(dat_raw)) {
  dat_raw <- dat_raw %>% dplyr::left_join(n_ema_tbl, by = "user_id")
}

# ============================
# 4) METADATI (genere) + TAG ESAMI (exam_period)
# ============================

# ============================
# 4) METADATI (genere, corso) + TAG ESAMI (exam_period)
# ============================

# -- 4a) Metadati: user_id, course, sex (data in meta è "dd_mm_yyyy", qui NON usata per il tag)
meta_df <- suppressMessages(rio::import(META_PATH)) %>%
  dplyr::rename(
    user_id = any_of(c("subj_code", "user_id", "id_anon")),
    date_meta = any_of(c("date", "Date", "data"))
  ) %>%
  dplyr::mutate(
    # parsiamo la data se presente (non è necessaria per il tag esami, ma utile per QA)
    date_meta = suppressWarnings(as.Date(date_meta, format = "%d_%m_%Y")),
    # normalizziamo il nome corso
    course = dplyr::case_when(
      stringr::str_detect(course, regex("^psico", ignore_case = TRUE)) ~
        "Psicometria",
      stringr::str_detect(course, regex("^test", ignore_case = TRUE)) ~
        "Testing",
      stringr::str_detect(course, regex("^inter", ignore_case = TRUE)) ~
        "Interventi",
      stringr::str_detect(course, regex("^clin", ignore_case = TRUE)) ~
        "Clinica",
      TRUE ~ course
    ),
    sex = dplyr::case_when(
      stringr::str_detect(sex, regex("^f", ignore_case = TRUE)) ~ "Femmina",
      stringr::str_detect(sex, regex("^m", ignore_case = TRUE)) ~ "Maschio",
      TRUE ~ sex
    )
  ) %>%
  dplyr::select(user_id, course, sex, date_meta) %>%
  dplyr::distinct()

# -- 4b) Join metadati
# NB: dat_final deve esistere già dal passo precedente (dopo esclusioni careless + n_ema)
dat_final <- dat_raw %>%
  dplyr::left_join(meta_df, by = "user_id")

# -- 4c) Definizione date d’esame per corso (baseline altrimenti)
psico_pre <- as.Date(c("2025-04-14", "2025-05-21"))
psico_post <- as.Date(c("2025-04-15", "2025-05-22"))

test_pre <- as.Date(c("2025-04-14", "2025-05-25"))
test_post <- as.Date(c("2025-04-15", "2025-05-26"))

interv_pre <- as.Date("2025-05-12")
interv_post <- as.Date("2025-05-13")

# -- 4d) Tag exam_period in base al corso e alla data EMA (colonna 'day')
dat_final <- dat_final %>%
  dplyr::mutate(
    day = as.Date(day),
    exam_period = dplyr::case_when(
      course == "Clinica" ~ "baseline",

      course == "Psicometria" & day %in% psico_pre ~ "pre_exam",
      course == "Psicometria" & day %in% psico_post ~ "post_exam",

      course == "Testing" & day %in% test_pre ~ "pre_exam",
      course == "Testing" & day %in% test_post ~ "post_exam",

      course == "Interventi" & day %in% interv_pre ~ "pre_exam",
      course == "Interventi" & day %in% interv_post ~ "post_exam",

      TRUE ~ "baseline"
    ),
    exam_period = factor(
      exam_period,
      levels = c("baseline", "pre_exam", "post_exam")
    )
  )

# -- 4e) QA rapido (opzionale)
dplyr::count(dat_final, course, exam_period)
table(dat_final$exam_period)

# ============================
# 5) ESCLUSIONI
#   A) careless responding (lista a priori)
#   B) numerosità EMA: 5 <= n_ema <= 40
# ============================

# (A) careless
excl_careless <- dat_final %>%
  dplyr::distinct(user_id) %>%
  dplyr::filter(user_id %in% careless_ids) %>%
  mutate(reason = "careless_responding (a priori)")

# (B) n_ema fuori range
excl_nema <- n_ema_tbl %>%
  dplyr::filter(is.na(n_ema) | n_ema < 5 | n_ema > 40) %>%
  mutate(reason = "n_ema_out_of_range (must be 5..40)")

# Log esclusioni unificato
excl_log <- bind_rows(excl_careless, excl_nema) %>%
  dplyr::distinct(user_id, reason)

# Soggetti da mantenere
keep_ids <- dat_final %>%
  dplyr::distinct(user_id) %>%
  dplyr::filter(!(user_id %in% excl_log$user_id)) %>%
  pull(user_id)

dat_clean <- dat_final %>%
  dplyr::filter(user_id %in% keep_ids)

# ============================
# 6) QA ESSENZIALE (stampata e salvata su file)
# ============================

qa_lines <- c()

n_subj_raw <- n_distinct(dat_raw$user_id)
n_subj_clean <- n_distinct(dat_clean$user_id)
n_excl <- n_distinct(excl_log$user_id)

qa_lines <- c(
  qa_lines,
  sprintf("Totale soggetti (raw): %d", n_subj_raw),
  sprintf("Esclusi (tot): %d", n_excl),
  sprintf("   - Careless: %d", n_distinct(excl_careless$user_id)),
  sprintf("   - n_ema fuori [5,40]: %d", n_distinct(excl_nema$user_id)),
  sprintf("Soggetti inclusi (finale): %d", n_subj_clean),
  ""
)

# Distribuzione n_ema
nema_keep <- n_ema_tbl %>% filter(user_id %in% keep_ids)
qa_lines <- c(
  qa_lines,
  sprintf(
    "n_ema (included) — min: %d, q1: %d, median: %d, mean: %.2f, q3: %d, max: %d",
    min(nema_keep$n_ema, na.rm = TRUE),
    quantile(nema_keep$n_ema, 0.25, na.rm = TRUE) %>% as.integer(),
    median(nema_keep$n_ema, na.rm = TRUE) %>% as.integer(),
    mean(nema_keep$n_ema, na.rm = TRUE),
    quantile(nema_keep$n_ema, 0.75, na.rm = TRUE) %>% as.integer(),
    max(nema_keep$n_ema, na.rm = TRUE)
  )
)

# Genere (se presente)
if ("sex" %in% names(dat_clean)) {
  gen_tab <- dat_clean %>%
    distinct(user_id, sex) %>%
    count(sex, sort = TRUE)
  qa_lines <- c(qa_lines, "", "Distribuzione genere (soggetti unici):")
  qa_lines <- c(qa_lines, paste0(" - ", gen_tab$sex, ": ", gen_tab$n))
}

# Exam period (se disponibile)
if ("exam_period" %in% names(dat_clean)) {
  ex_tab <- dat_clean %>%
    distinct(user_id, exam_period) %>%
    count(exam_period, sort = TRUE)
  qa_lines <- c(qa_lines, "", "Distribuzione exam_period (soggetti unici):")
  qa_lines <- c(
    qa_lines,
    paste0(
      " - ",
      ifelse(is.na(ex_tab$exam_period), "NA", ex_tab$exam_period),
      ": ",
      ex_tab$n
    )
  )
}

# (A) careless
excl_careless <- dat_final %>%
  dplyr::distinct(user_id) %>%
  dplyr::filter(user_id %in% careless_ids) %>%
  dplyr::left_join(n_ema_tbl, by = "user_id") %>%
  dplyr::mutate(reason = "careless_responding (a priori)")

# (B) n_ema fuori range
excl_nema <- n_ema_tbl %>%
  dplyr::filter(is.na(n_ema) | n_ema < 5 | n_ema > 40) %>%
  dplyr::mutate(reason = "n_ema_out_of_range (must be 5..40)")

# Log
excl_log <- dplyr::bind_rows(excl_careless, excl_nema) %>%
  dplyr::distinct(user_id, n_ema, reason)

# Stampa a console
cat(paste(qa_lines, collapse = "\n"), "\n")

# Salva QA su file
writeLines(qa_lines, con = OUT_QA_TXT)

# ============================
# 7) EXPORT: dataset finale + log esclusioni
# ============================

# Dataset finale
rio::export(dat_clean, OUT_CLEAN_CSV)
saveRDS(dat_clean, OUT_CLEAN_RDS)

# Log esclusioni (uno per riga con motivazione)
excl_log %>%
  arrange(reason, user_id) %>%
  rio::export(OUT_LOG_EXCL)

message(
  "Pulizia completata.\n- Dataset: ",
  OUT_CLEAN_CSV,
  " | ",
  OUT_CLEAN_RDS,
  "\n- Log esclusioni: ",
  OUT_LOG_EXCL,
  "\n- QA: ",
  OUT_QA_TXT
)

# eof ---

Il modello Stan usato è

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

Lo script usato per l'analisi è

> set.seed(20250917)
> # z-score "sicuro"
> z_ <- function(x) {
+   m <- mean(x, na.rm = TRUE)
+   s <- sd(x, na.rm = TRUE)
+   if (is.na(s) || s == 0) return(x - m)
+   (x - m) / s
+ }
> # ---------- Carica il dataset "clean" ----------
> d <- rio::import(
+   here::here("data", "processed", "ema_plus_scales_cleaned.csv")
+ ) |>
+   dplyr::rename(
+     # rinomina per coerenza con lo Stan (baseline a 5 dimensioni)
+     pid5_negative_affect_baseline = domain_negative_affect_baseline,
+     pid5_detachment_baseline = domain_detachment_baseline,
+     pid5_antagonism_baseline = domain_antagonism_baseline,
+     pid5_disinhibition_baseline = domain_disinhibition_baseline,
+     pid5_psychoticism_baseline = domain_psychoticism_baseline
+   )
> # ---------- Funzione di preparazione ----------
> prepare_fragility_data <- function(d, min_obs_per_period = 1) {
+   d <- d %>% mutate(day = as.Date(day))
+ 
+   # 1) Item ordinali (happy/sad/satisfied/angry) -> 1..7
+   item_data <- d %>%
+     filter(!is.na(exam_period)) %>%
+     select(
+       user_id,
+       day,
+       hour,
+       exam_period,
+       happy,
+       sad,
+       satisfied,
+       angry,
+       # EMA osservate (4 dimensioni: det/ant/dis/psy)
+       pid5_detachment,
+       pid5_antagonism,
+       pid5_disinhibition,
+       pid5_psychoticism,
+       # baseline 5D (una per soggetto)
+       pid5_negative_affect_baseline,
+       pid5_detachment_baseline,
+       pid5_antagonism_baseline,
+       pid5_disinhibition_baseline,
+       pid5_psychoticism_baseline
+     ) %>%
+     mutate(
+       # mappa 0..100 -> 1..7
+       happy_ord = pmax(1L, pmin(7L, as.integer(round(happy / 100 * 6) + 1))),
+       sad_ord = pmax(1L, pmin(7L, as.integer(round(sad / 100 * 6) + 1))),
+       satisfied_ord = pmax(
+         1L,
+         pmin(7L, as.integer(round(satisfied / 100 * 6) + 1))
+       ),
+       angry_ord = pmax(1L, pmin(7L, as.integer(round(angry / 100 * 6) + 1))),
+       period_numeric = case_when(
+         exam_period == "baseline" ~ 1L,
+         exam_period == "pre_exam" ~ 2L,
+         exam_period == "post_exam" ~ 3L,
+         TRUE ~ NA_integer_
+       )
+     ) %>%
+     filter(!is.na(period_numeric)) %>%
+     arrange(user_id, day, hour)
+ 
+   # 2) indicizza soggetti, tieni solo chi ha almeno un pre e un post
+   subj0 <- item_data %>%
+     distinct(user_id) %>%
+     arrange(user_id) %>%
+     mutate(subject_numeric0 = row_number())
+ 
+   obs0 <- item_data %>%
+     left_join(subj0, by = "user_id") %>%
+     arrange(subject_numeric0, day, hour)
+ 
+   keep_ids <- obs0 %>%
+     group_by(user_id) %>%
+     summarise(
+       has_pre = any(period_numeric == 2L),
+       has_post = any(period_numeric == 3L),
+       .groups = "drop"
+     ) %>%
+     filter(has_pre & has_post) %>%
+     pull(user_id)
+ 
+   obs1 <- obs0 %>% filter(user_id %in% keep_ids)
+ 
+   # 3) re-index finale soggetti
+   subj_map <- obs1 %>%
+     distinct(user_id) %>%
+     arrange(user_id) %>%
+     mutate(subject_numeric = row_number())
+ 
+   obs2 <- obs1 %>%
+     select(-subject_numeric0) %>%
+     left_join(subj_map, by = "user_id") %>%
+     arrange(subject_numeric, day, hour)
+ 
+   # 4) filtra osservazioni con almeno 1 item presente
+   obs_items <- obs2 %>%
+     filter(if_any(
+       c(happy_ord, sad_ord, satisfied_ord, angry_ord),
+       ~ !is.na(.)
+     )) %>%
+     arrange(subject_numeric, day, hour)
+ 
+   # 5) vincolo minimo di osservazioni per PRE e POST
+   keep_subj <- obs_items %>%
+     count(subject_numeric, period_numeric, name = "n_obs") %>%
+     pivot_wider(
+       names_from = period_numeric,
+       values_from = n_obs,
+       names_prefix = "p",
+       values_fill = 0
+     ) %>%
+     mutate(keep = (p2 >= min_obs_per_period) & (p3 >= min_obs_per_period)) %>%
+     filter(keep) %>%
+     pull(subject_numeric)
+ 
+   obs_data <- obs_items %>%
+     filter(subject_numeric %in% keep_subj) %>%
+     arrange(subject_numeric, day, hour) %>%
+     mutate(obs_id = row_number())
+ 
+   # 6) long degli item (solo non-NA)
+   items_long <- obs_data %>%
+     select(obs_id, happy_ord, sad_ord, satisfied_ord, angry_ord) %>%
+     pivot_longer(
+       ends_with("_ord"),
+       names_to = "item_name",
+       values_to = "response"
+     ) %>%
+     filter(!is.na(response)) %>%
+     mutate(
+       item_id = case_when(
+         item_name == "happy_ord" ~ 1L,
+         item_name == "sad_ord" ~ 2L,
+         item_name == "satisfied_ord" ~ 3L,
+         item_name == "angry_ord" ~ 4L
+       ),
+       response = as.integer(response)
+     ) %>%
+     arrange(obs_id, item_id)
+ 
+   # 7) Baseline 5D per soggetto, z-score + imputazione a 0 se NA
+   base_dims <- obs_data %>%
+     group_by(subject_numeric, user_id) %>%
+     summarise(
+       naff_b = first(pid5_negative_affect_baseline),
+       det_b = first(pid5_detachment_baseline),
+       ant_b = first(pid5_antagonism_baseline),
+       dis_b = first(pid5_disinhibition_baseline),
+       psy_b = first(pid5_psychoticism_baseline),
+       .groups = "drop"
+     ) %>%
+     arrange(subject_numeric)
+ 
+   X_base_df <- base_dims %>%
+     mutate(across(
+       c(naff_b, det_b, ant_b, dis_b, psy_b),
+       z_,
+       .names = "z_{.col}"
+     )) %>%
+     mutate(across(starts_with("z_"), ~ coalesce(.x, 0))) %>%
+     select(starts_with("z_"))
+ 
+   X_base <- as.matrix(X_base_df)
+   colnames(X_base) <- c("z_naff_b", "z_det_b", "z_ant_b", "z_dis_b", "z_psy_b")
+ 
+   # 8) EMA osservate (SOLO 4 dimensioni: det/ant/dis/psy) per osservazione
+   ema_mat <- obs_data %>%
+     transmute(
+       obs_id,
+       subj = subject_numeric,
+       det_e = z_(pid5_detachment),
+       ant_e = z_(pid5_antagonism),
+       dis_e = z_(pid5_disinhibition),
+       psy_e = z_(pid5_psychoticism)
+     )
+ 
+   ema_long <- ema_mat %>%
+     pivot_longer(
+       cols = c(det_e, ant_e, dis_e, psy_e),
+       names_to = "dim",
+       values_to = "value"
+     ) %>%
+     filter(is.finite(value)) %>%
+     mutate(
+       # Mappa alle dimensioni Stan: 1=NA (latente dagli item), 2..5 = det/ant/dis/psy
+       dim_id = dplyr::recode(
+         dim,
+         det_e = 2L,
+         ant_e = 3L,
+         dis_e = 4L,
+         psy_e = 5L
+       ),
+       dim_id = as.integer(dim_id)
+     ) %>%
+     arrange(obs_id, dim_id)
+ 
+   # 9) controlli
+   I <- n_distinct(obs_data$subject_numeric)
+   N_obs <- nrow(obs_data)
+   stopifnot(identical(sort(unique(obs_data$subject_numeric)), seq_len(I)))
+   stopifnot(max(obs_data$subject_numeric) == I)
+   stopifnot(nrow(X_base) == I, ncol(X_base) == 5)
+ 
+   # 10) dati per Stan
+   stan_data <- list(
+     I = I,
+     N_obs = N_obs,
+     K = 4L, # happy, sad, satisfied, angry
+     P = 3L, # baseline, pre, post
+     D = 5L, # 5 dimensioni EMA a livello soggetto (1=NA latente, 2..5 osservate)
+     N_items = nrow(items_long),
+     y_item = as.integer(items_long$response),
+     item_id = as.integer(items_long$item_id),
+     obs_id = as.integer(items_long$obs_id),
+     subject = as.integer(obs_data$subject_numeric),
+     period = as.integer(obs_data$period_numeric),
+     M_ema = nrow(ema_long),
+     ema_val = as.numeric(ema_long$value),
+     ema_dim = as.integer(ema_long$dim_id), # solo {2,3,4,5}
+     ema_obs = as.integer(ema_long$obs_id),
+     X_base = X_base,
+     use_ema = 1.0
+   )
+ 
+   list(
+     stan_data = stan_data,
+     obs_data = obs_data,
+     subject_map = subj_map,
+     items_long = items_long
+   )
+ }
> # ---------- Prepara e fit ----------
> frag_data <- prepare_fragility_data(d, min_obs_per_period = 1)
> with(frag_data$stan_data, {
+   cat("I =", I, "| N_obs =", N_obs, "| N_items =", N_items, "\n")
+ })
I = 172 | N_obs = 4825 | N_items = 19300 
> stan_file <- here::here(
+   "scripts",
+   "02_stress_reactivity",
+   "fragility_ema_vs_baseline.stan"
+ )
> model_frag <- cmdstan_model(stan_file)
Model executable is up to date!
> # Fit A: baseline-only
> data_A <- frag_data$stan_data
> data_A$use_ema <- 0.0
> fit_A <- model_frag$variational(
+   data = data_A,
+   seed = 20250912,
+   algorithm = "meanfield",
+   elbo_samples = 100,
+   adapt_engaged = TRUE,
+   tol_rel_obj = 0.001,
+   eval_elbo = 100,
+   output_samples = 2000,
+   refresh = 200
+ )
------------------------------------------------------------ 
EXPERIMENTAL ALGORITHM: 
  This procedure has not been thoroughly tested and may be unstable 
  or buggy. The interface is subject to change. 
------------------------------------------------------------ 
Gradient evaluation took 0.008218 seconds 
1000 transitions using 10 leapfrog steps per transition would take 82.18 seconds. 
MEDIAN ELBO CONVERGED 
Drawing a sample of size 2000 from the approximate posterior...  
COMPLETED. 
Finished in  40.9 seconds.
> # Fit B: baseline + EMA (inclusa NA latente dagli item)
> data_B <- frag_data$stan_data
> data_B$use_ema <- 1.0
> fit_B <- model_frag$variational(
+   data = data_B,
+   seed = 20250912,
+   algorithm = "meanfield",
+   elbo_samples = 100,
+   adapt_engaged = TRUE,
+   tol_rel_obj = 0.001,
+   eval_elbo = 100,
+   output_samples = 2000,
+   refresh = 200
+ )
------------------------------------------------------------ 
EXPERIMENTAL ALGORITHM: 
  This procedure has not been thoroughly tested and may be unstable 
  or buggy. The interface is subject to change. 
------------------------------------------------------------ 
Gradient evaluation took 0.008205 seconds 
1000 transitions using 10 leapfrog steps per transition would take 82.05 seconds. 
MEDIAN ELBO CONVERGED 
Drawing a sample of size 2000 from the approximate posterior...  
COMPLETED. 
Finished in  41.3 seconds.
> # 5) Estrai Bayesian R2 della regressione della fragilità
> R2_A <- as_draws_df(fit_A$draws("R2_frag"))$R2_frag
> R2_B <- as_draws_df(fit_B$draws("R2_frag"))$R2_frag
> cat("Bayesian R2 – baseline-only:     ", round(mean(R2_A), 3), "\n")
Bayesian R2 – baseline-only:      0.173 
> cat("Bayesian R2 – baseline + EMA:    ", round(mean(R2_B), 3), "\n")
Bayesian R2 – baseline + EMA:     0.234 
> cat(
+   "ΔR2 (B − A):                     ",
+   round(mean(R2_B - R2_A), 3),
+   "  (95% CrI ≈ [",
+   round(quantile(R2_B - R2_A, c(.025, .975)), 3),
+   "] )\n",
+   sep = ""
+ )
ΔR2 (B − A):                     0.061  (95% CrI ≈ [-0.3030.424] )
> # =======================
> # LOO: test cruciale EMA
> # =======================
> 
> # 1) Estrarre log-lik per osservazione (beep)
> ll_A <- fit_A$draws("log_lik_obs", format = "matrix") # draws x N_obs
> ll_B <- fit_B$draws("log_lik_obs", format = "matrix")
> 
> # controllo dimensioni coerenti
> stopifnot(ncol(ll_A) == ncol(ll_B))
> N_obs <- ncol(ll_A)
> 
> # 2) LOO con moment matching (riduce impatto k alti)
> loo_A <- loo::loo(ll_A, moment_match = TRUE)
Warning message:
Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
 
> loo_B <- loo::loo(ll_B, moment_match = TRUE)
Warning message:
Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
 
> 
> cat("\n=== Pareto-k (beep level) ===\n")

=== Pareto-k (beep level) ===
> print(loo::pareto_k_table(loo_A))
Pareto k diagnostic values:
                         Count Pct.    Min. ESS
(-Inf, 0.7]   (good)      897  18.6%   48      
   (0.7, 1]   (bad)      2570  53.3%   <NA>    
   (1, Inf)   (very bad) 1358  28.1%   <NA>    
> print(loo::pareto_k_table(loo_B))
Pareto k diagnostic values:
                         Count Pct.    Min. ESS
(-Inf, 0.7]   (good)      888  18.4%   48      
   (0.7, 1]   (bad)      2583  53.5%   <NA>    
   (1, Inf)   (very bad) 1354  28.1%   <NA>    
> 
> # 3) Confronto modelli (beep level)
> cmp_obs <- loo::loo_compare(list(B = loo_B, A = loo_A))
> cat("\n=== LOO comparison (beep level) ===\n")

=== LOO comparison (beep level) ===
> print(cmp_obs)
  elpd_diff se_diff
A  0.0       0.0   
B -7.2       0.7   
> 
> # Δelpd = elpd(B) - elpd(A)
> delta_elpd_obs <- as.numeric(cmp_obs["A", "elpd_diff"]) * -1
> per_obs <- delta_elpd_obs / N_obs
> pct_gain <- (exp(per_obs) - 1) * 100
> 
> cat(sprintf(
+   "\nΔelpd (B−A): %.1f  | per-beep: %.4f  | mean %% gain per beep: %.2f%%\n",
+   delta_elpd_obs,
+   per_obs,
+   pct_gain
+ ))

Δelpd (B−A): -0.0  | per-beep: -0.0000  | mean % gain per beep: 0.00%
> 
> # 4) (Opzionale) LOO aggregato per soggetto: somma log-lik di tutti i beep dello stesso soggetto
> #    È utile come analisi di robustezza (meno punti = minore varianza del PSIS)
> subj_idx <- frag_data$stan_data$subject # length = N_obs, valori 1..I
> I <- max(subj_idx)
> 
> # Matrice di raggruppamento JxN (J=I soggetti). Sparse per efficienza.
> G <- sparseMatrix(i = subj_idx, j = seq_len(N_obs), x = 1, dims = c(I, N_obs)) # I x N_obs
> 
> # draws x N_obs  %*%  t(G)  ->  draws x I (somma per soggetto)
> ll_A_subj <- ll_A %*% t(as.matrix(G))
> ll_B_subj <- ll_B %*% t(as.matrix(G))
> 
> loo_A_subj <- loo::loo(ll_A_subj, moment_match = TRUE)
Warning message:
Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
 
> loo_B_subj <- loo::loo(ll_B_subj, moment_match = TRUE)
Warning message:
Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
 
> 
> cat("\n=== Pareto-k (subject level) ===\n")

=== Pareto-k (subject level) ===
> print(loo::pareto_k_table(loo_A_subj))
Pareto k diagnostic values:
                         Count Pct.    Min. ESS
(-Inf, 0.7]   (good)       0     0.0%  <NA>    
   (0.7, 1]   (bad)        0     0.0%  <NA>    
   (1, Inf)   (very bad) 172   100.0%  <NA>    
> print(loo::pareto_k_table(loo_B_subj))
Pareto k diagnostic values:
                         Count Pct.    Min. ESS
(-Inf, 0.7]   (good)       0     0.0%  <NA>    
   (0.7, 1]   (bad)        0     0.0%  <NA>    
   (1, Inf)   (very bad) 172   100.0%  <NA>    
> 
> cmp_subj <- loo::loo_compare(list(B = loo_B_subj, A = loo_A_subj))
> cat("\n=== LOO comparison (subject level) ===\n")

=== LOO comparison (subject level) ===
> print(cmp_subj)
  elpd_diff se_diff
A  0.0       0.0   
B -5.1       1.1   
> 
> delta_elpd_subj <- as.numeric(cmp_subj["A", "elpd_diff"]) * -1
> per_subj <- delta_elpd_subj / I
> pct_gain_subj <- (exp(per_subj) - 1) * 100
> 
> cat(sprintf(
+   "\nΔelpd (B−A): %.1f  | per-subject: %.4f  | mean %% gain per subject: %.2f%%\n",
+   delta_elpd_subj,
+   per_subj,
+   pct_gain_subj
+ ))

Δelpd (B−A): -0.0  | per-subject: -0.0000  | mean % gain per subject: 0.00%
> 
> # 5) Riassunto finale — “decisione” per la domanda di ricerca
> winner_obs <- rownames(cmp_obs)[1]
> winner_subj <- rownames(cmp_subj)[1]
> cat("\n=== Verdict ===\n")

=== Verdict ===
> cat(sprintf(
+   "Beep level:   %s ha migliore ELPD (Δelpd = %.1f).\n",
+   winner_obs,
+   delta_elpd_obs
+ ))
Beep level:   A ha migliore ELPD (Δelpd = -0.0).
> cat(sprintf(
+   "Subject level:%s ha migliore ELPD (Δelpd = %.1f).\n",
+   winner_subj,
+   delta_elpd_subj
+ ))
Subject level:A ha migliore ELPD (Δelpd = -0.0).

Come vedi non c'è alcuna evidenza che le informazioni EMA del PID-5 aggiungano qulcosa alle informazioni di baseline. Ho provato diverse alternative relative alla selezione iniziale dei soggetti (da 172 come qui fino a 232), ma il risultato non cambia.

C'è un modo migliore per fare le analisi sperando di trovare evidenze a favore dell'ipotesi della ricerca (se la "fragilità psicologica" degli studenti, operazionalizzata dalla loro reazione emotiva all'esame, sia meglio misurata dalle caratteristiche di tratto, fornite dal PID-5 somministrato integralmente una volta soltanto, all'inizio del progetto, o da un piccolo sottoinsieme di item del PID-5 somministrate molte volte nel tempo)?