# Mammography AI Equity (Synthetic)

This repository contains a **fully synthetic** dataset and a reference analysis workflow for evaluating
equity of an AI risk score in mammography screening.

**No real patient data is included.**

## Contents

- `data/mammography_ai_equity_data_simulated.csv` — synthetic dataset (5,000 exams)
- `equity-evaluation-mammo-ai.Rmd` — R Markdown analysis (AUC, calibration, ECE, confounding/adjustment)
- `simulate_data.R` / `simulate_data.py` — scripts to regenerate the dataset

## Data columns

- `exam_id`: synthetic exam identifier
- `age`: age in years (approx. 40–74)
- `ethnicity`: synthetic coarse category
- `ses`: socioeconomic status (Low/Medium/High)
- `region`: Swedish region label (for stratification examples)
- `mammography_vendor`: anonymized vendor label **A/B/C/D**
- `breast_density`: A–D (A lowest, D highest)
- `cancer`: 1 = cancer diagnosed, 0 = no cancer
- `ai_score`: AI risk score in [0,1] (synthetic; intentionally imperfect)

## Run the analysis

In RStudio, open `equity-evaluation-mammo-ai.Rmd` and knit to HTML.

Required packages:

```r
install.packages(c("dplyr","readr","knitr","pROC","ggplot2","scales","broom","tidyr","forcats"))
```

## Regenerate the dataset (optional)

**R**

```bash
Rscript simulate_data.R 5000 20260218 data/mammography_ai_equity_data_simulated.csv
```

**Python**

```bash
python simulate_data.py --n 5000 --seed 20260218 --out data/mammography_ai_equity_data_simulated.csv
```

## Important note

This dataset is simulated for educational / demonstration purposes only.
Do not interpret results as medical evidence.

## Notes on simulation

The synthetic generator is tuned to mimic a typical screening-AI homework dataset: cancer prevalence ~0.15, very high discrimination (AUC close to 1), and imperfect calibration (the score behaves like an uncalibrated risk score rather than a well-calibrated probability).
