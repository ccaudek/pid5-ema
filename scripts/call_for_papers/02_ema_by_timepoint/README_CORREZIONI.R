# ==============================================================================
# RIEPILOGO: Problemi e Correzioni per Replicare l'Abstract
# ==============================================================================

# PROBLEMA PRINCIPALE IDENTIFICATO
# ================================

# Nello script 02 originale, linee 187-208, i dati PID-5 EMA vengono
# elaborati in modo errato:

# 1. AGGREGAZIONE (linee 187-195):
#    pid5_agg <- pid5_df %>%
#      group_by(ID, timepoint) %>%   # Aggrega per ID E timepoint
#      summarise(...)

# 2. ESTRAZIONE (linee 205-208):
#    pid5_baseline <- df_clean %>%
#      filter(timepoint == "baseline") %>%   # ← PROBLEMA: solo baseline!
#      select(ID, starts_with("pid5"))

# Questo significa che si usano SOLO i valori EMA raccolti durante il
# periodo "baseline" (lontano dall'esame), ignorando tutti i dati EMA
# raccolti nei periodi pre e post esame.

# L'ABSTRACT DICE:
# ================
# "PID-5 at baseline, then conducted intensive EMA over 2.5 months
# (twice-weekly assessments using 15 PID-5 items; 3 per domain)"
#
# "each EMA PID-5 domain moderated stress reactivity"
#
# "EMA measures added predictive value beyond baseline assessment"

# Questo suggerisce che i PID-5 EMA dovrebbero essere aggregati su TUTTO
# il periodo di 2.5 mesi come misure trait-like stabili, non solo sul
# periodo baseline.

# SOLUZIONE
# =========

# Negli script corretti (02_CORRECTED e 03_CORRECTED), la variabile
# PID5_VERSION controlla quale aggregazione usare:

# PID5_VERSION = "ema_all"       ← RACCOMANDATO
#   - Media di TUTTI i valori EMA per soggetto
#   - Coerente con "intensive EMA over 2.5 months"
#   - Crea predittori trait-like stabili

# PID5_VERSION = "ema_baseline"  ← ERRATO (versione attuale)
#   - Solo valori EMA del periodo baseline
#   - Ignora ~2/3 dei dati EMA

# PID5_VERSION = "pid5_full"     ← Alternativa
#   - PID-5 completo (220 item) somministrato al baseline
#   - Utile per confronto con EMA

# PID5_VERSION = "ema_timevarying" ← Sperimentale
#   - Valori EMA diversi per ogni timepoint vocale
#   - Può aggiungere rumore (solo 3 misure per soggetto)

# RISULTATI ATTESI DALL'ABSTRACT
# ==============================

# L'abstract riporta questi effetti chiave (su vocale /a/ primariamente):

# 1. Negative Affectivity × Stress → F0 Mean
#    β = 5.03 Hz [1.26, 8.76]
#    Interpretazione: +1 SD in NA amplifica l'aumento del pitch sotto stress

# 2. Detachment × Recovery → F0 Mean
#    β = -4.37 [-8.30, -0.46]
#    Interpretazione: +1 SD in Det impair il recupero del pitch post-stress

# 3. Antagonism × Recovery
#    β = 3.44-4.17 across vowels
#    Interpretazione: +1 SD in Ant facilita il recupero

# 4. Psychoticism
#    "pervasive pitch elevation and voice quality degradation"
#    Effetti su f0_mean e jitter/nne

# 5. Disinhibition × F2
#    β = 63 Hz [26, 99]
#    Interpretazione: Alterazioni articolatorie (secondo formante)

# 6. Formant effects
#    β = 63-117 Hz per moderazione formanti

# COME PROCEDERE
# ==============

# 1. Esegui 00_diagnostic_pid5_versions.R per ispezionare i dati
# 2. Esegui 02_voice_personality_analysis_CORRECTED.R con PID5_VERSION = "ema_all"
# 3. Esegui 03_moderation_analysis_CORRECTED.R
# 4. Confronta i risultati con l'abstract

# Se i risultati non sono coerenti con "ema_all", prova:
# - PID5_VERSION = "pid5_full" (PID-5 completo baseline)
# - Verifica la qualità dei dati EMA (missing, outlier)
# - Considera che il campione è leggermente diverso (111 vs N attuale)

# NOTA IMPORTANTE
# ===============
# I valori esatti nell'abstract potrebbero differire perché:
# - Ora hai incluso alcuni soggetti in più
# - Potrebbero esserci state piccole modifiche nel preprocessing
# - I seed casuali dell'imputazione possono variare leggermente

# L'obiettivo è replicare il PATTERN QUALITATIVO:
# - NA amplifica stress response
# - Detachment impairs recovery
# - Antagonism facilitates recovery
# - Psychoticism shows pervasive effects
# - Disinhibition affects F2
