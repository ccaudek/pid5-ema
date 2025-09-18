# PID-5 EMA — Starter Bundle

This archive contains:
- `01_gen_df_exam_tagged_refactored.qmd`: a refactored, self-contained Quarto notebook
  that builds the analysis dataset by merging baseline (PID-5 domains, ESI-BF, optional
  DASS-21) with EMA data, applies reproducible quality filters, tags exam periods,
  computes state SCS composites and a standardized negative affect composite, runs
  basic QA summaries, and exports the final CSV.
- `00_data_dictionary.R`: script to generate `data/processed/data_dictionary.csv`
  from the processed dataset (or, if not present, from known variable names).
- `USAGE.md`: quick usage notes (this file).

## How to use

1. Place these files in your project root (or adjust paths accordingly).
2. Open `01_gen_df_exam_tagged_refactored.qmd` and verify:
   - `params` YAML (exam dates, `min_n_ema`, export path).
   - Paths to your processed baseline files (`data/processed/esi_bf.csv` and `data/processed/pid5.csv`)
     and EMA RDS (`data/raw/ema/ema_data_scoring.RDS`).
3. Render the QMD (e.g., `quarto render 01_gen_df_exam_tagged_refactored.qmd`).
   This will export `data/processed/ema_plus_baseline_exam_tagged.csv`.
4. Run `00_data_dictionary.R` to create `data/processed/data_dictionary.csv`.

Both files are written with clarity and reproducibility in mind, so you can quote them
in Materials & Methods and rerun easily.