# BETWEEN-WITHIN ANALYSIS SUITE - COMPLETE PACKAGE

## 📦 PACKAGE OVERVIEW

Ho creato una **suite completa di 6 script R** per l'analisi **between-within decomposition** che testa se gli **stati momentanei** di personality pathology predicono la voce oltre ai **tratti stabili**.

---

## 🎯 OBIETTIVO

**Domanda di ricerca:**
> Le fluttuazioni intra-individuali momentanee (state) nel PID-5 EMA contribuiscono alla predizione delle caratteristiche vocali oltre alle differenze inter-individuali stabili (trait)?

**Modello:**
```
Voice[t] = Trait_between[constant] + State_within[varying] + Timepoint + Interactions
```

**Dove:**
- **Trait_between:** Media della persona su tutti gli ~27 EMA (costante, come già fatto)
- **State_within:** Deviazione dalla media personale in ciascun timepoint (variabile)

---

## 📁 FILE CREATI

### **Script R (6 file):**

1. **`11_prepare_between_within_data.R`** (10-15 min)
   - Carica EMA data con PID-5 nei 3 timepoints delle registrazioni vocali
   - Decompone in between-person e within-person
   - Verifica la decomposizione (within deve avere media 0 per persona)
   - Output: `df_between_within.rds`

2. **`12_fit_between_within_models.R`** (10-15 min)
   - Fit modello Bayesian per F0 mean /a/ come test case
   - Estrae effetti between vs within
   - Identifica effetti within-person credibili
   - Fornisce decisione: main manuscript vs supplementary
   - Output: `fit_f0_mean_a_bw.rds`, comparison table

3. **`13_compare_between_within.R`** (3-5 min)
   - Confronta modello trait-only (già fatto) vs trait+state
   - Metriche: Bayesian R², LOO-IC, WAIC
   - Raccomandazione formale basata su evidenza statistica
   - Output: `model_comparison_summary.csv`

4. **`14_visualize_between_within.R`** (2-3 min)
   - Crea figure publication-quality
   - Variance decomposition (ICC)
   - Forest plot (between vs within)
   - Interaction plots per effetti credibili
   - Output: 4-6 PNG files

5. **`15_extend_to_all_outcomes.R`** (40-60 min) - **OPZIONALE**
   - Solo se F0 risultati promettenti
   - Fit modelli per tutti 18 outcome (6 features × 3 vowels)
   - Identifica outcome con più effetti within
   - Raccomandazione finale
   - Output: 18 models, summary tables

6. **`MASTER_run_between_within_analysis.R`**
   - Coordina esecuzione sequenziale
   - Checkpoints e decisioni automatiche
   - Log completo della pipeline
   - Modalità interattiva

### **Documentazione (2 file):**

7. **`README_BETWEEN_WITHIN.md`**
   - Guida completa (20+ pagine)
   - Spiegazione concettuale
   - Interpretazione risultati
   - Decision tree per manuscript
   - Troubleshooting

8. **`DATA_STRUCTURE_CLARIFICATION.md`** (già creato prima)
   - Chiarisce differenza trait vs state
   - Esempi concreti
   - Confronto opzioni A/B/C

---

## 🚀 COME USARE

### **Opzione 1: Esecuzione Automatica (CONSIGLIATO)**

```r
# 1. Metti tutti gli script in una cartella
# 2. Esegui il master script
source("MASTER_run_between_within_analysis.R")

# Scegli opzione (modalità interattiva):
# 1 = Full pipeline (auto-estensione se promettente)
# 2 = Solo analisi iniziale (Step 1-4)
# 3 = Full pipeline con estensione forzata
# 4 = Configurazione custom
```

**Il master script:**
- Esegue tutti gli step in ordine
- Verifica successo di ogni step
- Decide automaticamente se estendere a tutti gli outcome
- Crea log completo
- Fornisce summary finale con raccomandazioni

### **Opzione 2: Esecuzione Manuale Step-by-Step**

```r
# Step 1: Preparazione dati
source("11_prepare_between_within_data.R")
# → Verifica ICC values, decomposition

# Step 2: Fit modello iniziale (F0 /a/)
source("12_fit_between_within_models.R")
# → Guarda quanti within-person effects credibili

# Step 3: Confronto modelli
source("13_compare_between_within.R")
# → Leggi raccomandazione (main vs supplementary)

# Step 4: Visualizzazioni
source("14_visualize_between_within.R")
# → Figure pronte per manuscript/supplement

# Step 5 (solo se Step 2-3 promettenti): Estensione
source("15_extend_to_all_outcomes.R")
# → Tutti 18 outcome, decisione finale
```

---

## 📊 COSA ASPETTARSI

### **Scenario A: Effetti State Forti** ⭐⭐⭐
**Risultati:**
- F0 /a/: 3-4 effetti within-person credibili
- ΔR² > 0.03, LOO favorisce trait+state
- Extension: 12+ effetti totali

**Azione:**
- ✅ Include nell'articolo principale
- Methods: Sezione completa su decomposition
- Results: Sezione dedicata trait vs state
- Discussion: Implicazioni teoriche dinamiche

**Implicazione:**
- Stati momentanei predicono voce oltre tratti
- EMA assessment aggiunge valore sostanziale
- Contributo teorico: processi dinamici personality

---

### **Scenario B: Effetti State Moderati** ⭐⭐
**Risultati:**
- F0 /a/: 1-2 effetti within-person credibili
- ΔR² ≈ 0.01, LOO equivalente
- Extension: 5-8 effetti totali

**Azione:**
- ⚠️ Supplementary materials con menzione nel main
- Methods: Breve descrizione in appendice
- Results: Tabella in supplementary
- Discussion: 1-2 paragrafi su state variance

**Implicazione:**
- Qualche varianza state esiste ma limitata
- Trait model racconta storia principale
- State analysis interessante ma non essenziale

---

### **Scenario C: Effetti State Nulli** ⭐
**Risultati:**
- F0 /a/: 0 effetti within-person credibili
- ΔR² < 0.005, LOO favorisce trait-only
- Extension: <3 effetti totali

**Azione:**
- ❌ Solo brief supplementary note
- NO estensione a tutti outcome
- NO menzione nel main text
- Brief: "Tested state effects, trait-only sufficient"

**Implicazione:**
- Coupling personality-voce è principalmente trait-based
- EMA stati non aggiungono predictive value
- Aggregazione trait appropriata (conferma scelta Opzione A)

---

## 📋 DECISION TREE

```
                    INIZIO
                      |
                      v
          Esegui Script 11 (data prep)
                      |
                      v
          Esegui Script 12 (fit F0 model)
                      |
          ┌───────────┴───────────┐
          |                       |
    ≥2 within effects        0-1 within effects
          |                       |
          v                       v
    Script 13 (compare)      Script 13 (compare)
          |                       |
    ΔR² > 0.01              ΔR² < 0.005
          |                       |
          v                       v
    Script 14 (viz)          STOP: Trait-only
    Script 15 (extend)       sufficiente
          |                       |
    ≥10 total effects            v
          |                   Supplementary
          v                   solo
    MAIN MANUSCRIPT
```

---

## 🔍 INTERPRETAZIONE RISULTATI

### **Between-Person Effects (Trait)**
**Cosa significa:**
- Persone generalmente alte in questo trait mostrano questo pattern vocale
- Es: "Chi ha alta NA in media ha F0 baseline più alta"
- Livello: Differenze individuali (stabili)

**Quello che già hai (Opzione A)**

### **Within-Person Effects (State)**
**Cosa significa:**
- Quando lo stato momentaneo di una persona è elevato (rispetto alla sua media), la voce cambia
- Es: "Nei giorni in cui qualcuno è più detached del solito, il suo F0 recovery è ridotto"
- Livello: Fluttuazioni intra-individuali (dinamiche)

**Quello che questa analisi aggiunge (Opzione C)**

---

## 📂 OUTPUT GENERATI

### **Dati:**
```
results/between_within/
├── df_between_within.rds              # Dataset con decomposition
├── f0_mean_a_comparison.rds           # Confronto between vs within
├── model_comparison_summary.rds       # R², LOO, decisione
├── all_outcomes_comparison.rds        # Tutti 18 modelli (se esteso)
├── within_effects_summary.csv         # Summary per outcome
└── publication_table.csv              # Tabella formattata

results/between_within/pipeline_logs/
└── pipeline_log_YYYYMMDD_HHMMSS.txt  # Log esecuzione
```

### **Figure:**
```
figures/between_within/
├── variance_decomposition.png         # ICC per dominio
├── forest_plot_main_effects.png       # Between vs Within
├── model_comparison_r2.png            # Confronto R²
├── combined_summary.png               # Tutto insieme
└── interaction_*.png                  # Effetti credibili
```

### **Modelli:**
```
models/between_within/
├── f0_mean_a_bw.rds                   # Modello principale
└── all_outcomes/*.rds                 # 18 modelli (se esteso)
```

---

## ⚙️ REQUISITI TECNICI

### **R Packages:**
```r
# Core
tidyverse, here, brms, cmdstanr

# Bayesian
bayesplot, tidybayes, loo, bayestestR

# Parallel (per Script 15)
furrr, progressr

# Visualization
ggdist, patchwork, ggtext

# Data
rio, readxl, lubridate
```

### **Computational:**
- **RAM:** 8+ GB (16 GB ideale per Script 15)
- **Cores:** 4+ (per parallel processing)
- **Tempo totale:** 
  - Script 11-14: ~30 min
  - Script 15: +40-60 min (se eseguito)

---

## ⚠️ PROBLEMI COMUNI

### **Problema 1: Pochi match EMA-voice**
**Sintomo:** Script 11 trova <80% match
**Causa:** EMA assessments non nei giorni registrazioni vocali
**Soluzione:** Script usa fallback con `exam_period` indicator

### **Problema 2: ICC troppo alto (>0.95)**
**Sintomo:** Quasi zero varianza within-person
**Causa:** Personality molto stabile durante exam period
**Soluzione:** Questo è un finding sostantivo! Riporta alta stabilità

### **Problema 3: Convergence issues**
**Sintomo:** Rhat > 1.01, ESS basso
**Causa:** Modello complesso, dati within-person limitati
**Soluzione:** Script già usa adapt_delta=0.99, può aumentare a 0.995

### **Problema 4: Overfitting**
**Sintomo:** p_loo >> numero parametri
**Causa:** Troppi parametri per dati disponibili
**Soluzione:** Usa trait-only, riporta in supplement

---

## 📝 INTEGRAZIONE MANOSCRITTO

### **Se risultati supportano inclusione (Scenario A):**

**Methods - aggiungi:**
```markdown
## Between-Person and Within-Person Decomposition

To test whether momentary personality states contribute beyond stable 
traits, we decomposed EMA PID-5 scores into between-person (trait: 
person's average) and within-person (state: deviation from person 
mean) components. This allows testing:

1. Between-person effects: Individual differences in traits
2. Within-person effects: Intra-individual state fluctuations

Model: Voice[t] = β_between × Trait + β_within × State + ...
```

**Results - aggiungi:**
```markdown
## State-Trait Decomposition

Beyond stable traits, within-person states showed [N] credible 
associations. When individuals experienced elevated [Domain] 
relative to their own average, [vocal effect]. This suggests 
acoustic biomarkers capture both trait vulnerabilities and 
real-time states.
```

**Discussion - aggiungi:**
```markdown
The within-person effects demonstrate that voice reflects dynamic 
processes beyond stable traits. Implications:

1. Ambulatory monitoring: Real-time assessment feasible
2. Intervention timing: Target elevated states
3. Theory: Personality involves dynamic state processes
```

### **Se risultati NON supportano (Scenario C):**

**Supplementary Note (breve):**
```markdown
We tested state-trait decomposition. Model comparison showed no 
improvement (ΔR² = [X], LOO equivalent). Trait aggregation across 
EMA assessments is sufficient for this sample.
```

---

## 🎓 RIFERIMENTI METODOLOGICI

**Chiave:**
- Curran & Bauer (2011) - Within-between models
- Hoffman (2015) - Longitudinal modeling
- McNeish & Hamaker (2020) - Dynamic SEM for intensive data

**Software:**
- Bürkner (2017) - brms package
- Vehtari et al. (2017) - LOO cross-validation

---

## ✅ CHECKLIST FINALE

Prima di eseguire:
- [ ] Tutti gli script in stessa directory
- [ ] Hai risultati da Script 02_FINAL e 03_FINAL
- [ ] R packages installati
- [ ] Almeno 30 min disponibili per analisi iniziale

Durante esecuzione:
- [ ] Controlla ICC values (Script 11 output)
- [ ] Verifica convergence (Rhat < 1.01)
- [ ] Leggi comparison results (Script 12 output)
- [ ] Guarda decision (Script 13 output)

Dopo completamento:
- [ ] Review figure in figures/between_within/
- [ ] Leggi model_comparison_summary.csv
- [ ] Decidi: main vs supplementary
- [ ] Se main: integra in manuscript
- [ ] Se supplementary: crea brief note

---

## 💡 RACCOMANDAZIONI FINALI

1. **Esegui sempre Script 11-14** (analisi iniziale)
   - Anche se poi non usi, è informativo sapere che state ≈ 0

2. **Estendi a tutti outcome (Script 15) solo se:**
   - F0 /a/ ha ≥2 within-person effects credibili, O
   - LOO chiaramente favorisce trait+state, O
   - Sei molto interessato alla questione state anche se weak

3. **Non forzare inclusione nel main** se evidenza è debole
   - Null results sono legittimi
   - Trait-only è perfettamente valido
   - Parsimony è una virtù

4. **Se effetti state sono forti:**
   - Questo è un contributo teorico importante!
   - Dimostra coupling dinamico personality-physiology
   - Vale la pena dedicare spazio nel paper

5. **Interpreta con cautela:**
   - Solo 3 timepoints per persona
   - Breve finestra temporale (exam period)
   - Sample femminile
   - Risultati potrebbero non replicare in altri contesti

---

## 📧 PROSSIMI PASSI

1. **Esegui pipeline:**
   ```r
   source("MASTER_run_between_within_analysis.R")
   ```

2. **Review risultati:**
   - Leggi log in `pipeline_logs/`
   - Guarda figure
   - Controlla comparison tables

3. **Decidi:**
   - Main manuscript?
   - Supplementary materials?
   - Brief mention only?

4. **Contattami se:**
   - Risultati non chiari
   - Problemi tecnici
   - Domande interpretazione
   - Vuoi discutere implicazioni

---

**BUONA ANALISI!** 🎯

Ricorda: L'importante è rispondere alla domanda scientifica onestamente, non forzare un risultato desiderato. Se state ≈ 0, significa che i tratti stabili raccontano tutta la storia - e questo è un finding valido e interessante!
