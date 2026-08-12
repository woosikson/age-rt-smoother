# age-rt-smoother

Reproduction code for

> **An Age-Structured Mathematical Modeling with Particle Smoother for Instantaneous
> Reproduction Number Estimation in Influenza Epidemics**
> Minsoo Kim, Jong-Hoon Kim, Woo-Sik Son (*Scientific Reports*, under revision)

The repository reproduces **every figure and every numeric table of the paper** — the
synthetic validation study, the eleven-season analysis of Korean national influenza
claims data (NHIS), and all sensitivity analyses — from the raw inputs.

---

## What the code does

An age-structured **SEIR** model (six age bands: 0–5, 6–11, 12–17, 18–44, 45–64, 65+)
is combined with a **fixed-window particle smoother** to infer the age-specific
infection probability per contact, `χ_i(t)`, from daily case counts. The age-specific
instantaneous reproduction number is the column sum of the time-dependent
next-generation matrix,

```
K_ij(t) = χ_i(t) · c_ij / γ · S_i(t)/N_i ,      𝓡_j(t) = Σ_i K_ij(t)
```

**Baseline configuration** (all analyses unless stated otherwise)

| Component | Value |
|---|---|
| Model | age-structured SEIR, σ = γ = 1/1.8 day⁻¹ (generation interval: Erlang, mean 3.6 d) |
| Particles | `Np = 1e5` |
| Smoother | fixed window `L = 16` days, per-stratum multinomial resampling, every step |
| Observation model | negative binomial, `r = 30`, evaluated in log space |
| Transmission prior | log-random-walk on χ, `σ_rw = 0.15` |
| Initial χ | non-negative least squares (NNLS) guess targeting 𝓡 ≈ 1.05 |
| Initial state | per particle `E_i(0), I_i(0) ~ iid Poisson(λ_i)`, `R_i(0) = 0`, `S_i(0) = N_i − E_i(0) − I_i(0)`;<br>λ = true value (synthetic) or mean of the 7 days preceding the season / γ (NHIS) |
| ODE solver | `Tsit5`, one fixed step per day |

---

## Repository layout

```
age-rt-smoother/
├── data/        inputs (see data/README.md for sources, licences and the NHIS rename step)
├── code/        00–05 Julia (data preparation and estimation) · make_prem_contact.py
│                figures.R · tables.R
├── results/     CSV produced by the Julia stage — shipped precomputed
├── output/      Figure*.eps and Table*.csv produced by the R stage
└── run_all.sh   the full chain, in order
```

`results/` is shipped **precomputed** so that the figures and tables can be regenerated in
minutes without repeating the ~19 h of estimation. The one exception is `results/NHIS_Flu.csv`,
the tidy copy of the claims data, which is not redistributed: step 0 recreates it from the
downloaded workbook. Delete `results/` and run `run_all.sh` to reproduce everything from the
raw inputs.

---

## Quick start

Regenerate all figures and tables from the shipped results (about 3 minutes):

```bash
Rscript code/figures.R      # → output/Figure2.eps … Figure7.eps, FigureS1.eps … FigureS10.eps
Rscript code/tables.R       # → output/TableS1.csv … TableS7.csv, output/tables.md
```

Reproduce everything from the raw data:

```bash
julia --project=. code/00_prepare_nhis.jl   # data/NHIS.xlsx -> results/NHIS_Flu.csv (20 s)
bash run_all.sh                             # the whole chain (~19 h), including the step above
```

**The NHIS claims data are not distributed with this repository**, in neither the original nor
the tidy form. Step 0 turns the downloaded workbook (`data/NHIS.xlsx`) into
`results/NHIS_Flu.csv`, which step 4 then reads; see [`data/README.md`](data/README.md) for the
download link, the rename step and the fingerprint of the release used in the paper.
The figures and tables do not read either file, so `Rscript code/figures.R` and
`Rscript code/tables.R` work on the shipped results without downloading anything.

### Dependencies

* **Julia** ≥ 1.9 — `DifferentialEquations`, `Distributions`, `StatsBase`, `Statistics`,
  `LinearAlgebra`, `NonNegLeastSquares`, `CSV`, `DataFrames`, `Dates`, `XLSX`, `Random`, `Printf`
  (`julia --project=. -e 'using Pkg; Pkg.instantiate()'`)
* **R** ≥ 4.2 — `readr`, `dplyr`, `tidyr`, `ggplot2`, `cowplot`, `scales`
* **Python** ≥ 3.9 (only to rebuild the Prem contact matrices) — `numpy`, `pandas`, `openpyxl`

Figures are written as **vector EPS** through `grDevices::cairo_ps`. The classic
`postscript()` device is not usable here: it cannot resolve Arial and mishandles the
script 𝓡 (U+211B) used for the reproduction number.

### Runtime and memory

| Stage | Wall clock | Peak memory |
|---|---|---|
| `00_prepare_nhis.jl` | 20 s | < 2 GB |
| `01_synthetic.jl` | seconds | < 1 GB |
| `02_particle_smoother.jl --mode sweep` (×3: clean-Poisson, clean-NB, noisy) | ≈ 3.5 h each | ≈ 30 GB at `Np = 1e6` |
| `02` weight sweeps, `baseline_fit`, `pfps` | ≈ 40 min total | < 4 GB |
| `03_renewal.jl` | seconds | < 1 GB |
| `04_nhis.jl` (11 seasons × 10 variants) | ≈ 3 h | ≈ 7 GB |
| `05_pfps_uniform.jl` | minutes | < 1 GB |
| `figures.R` + `tables.R` | ≈ 3 min | < 2 GB |

The `Np = 1e6` cells exist only to demonstrate convergence (Table S3) and path
degeneracy (Figure S3); the adopted configuration is `Np = 1e5`.

---

## Figure map

Every paper figure is produced by `code/figures.R` from the CSV listed below.
**Figure 1** is a TikZ schematic drawn inside the manuscript and has no code here.

### Main text

| Output | Content | Inputs (`results/`) |
|---|---|---|
| `Figure2.eps` | synthetic epidemic: age-specific cases and prescribed 𝓡_j(t) | `pre_confirm.csv`, `pre_rt.csv` |
| `Figure3.eps` | synthetic: particle smoother (median + 95% CrI) vs renewal equation vs truth | `ps_baseline_curves.csv`, `renewal_rt.csv` |
| `Figure4.eps` | 2018–19: all-age cases and 𝓡_j(t) medians | `nhis_rt_2018.csv` |
| `Figure5.eps` | 2018–19 by age group: model fit and 𝓡_j(t) with 95% CrI, renewal overlay | `nhis_rt_2018.csv` |
| `Figure6.eps` | the other ten seasons, cases and 𝓡_j(t) | `nhis_rt_{2009…2022}.csv` |
| `Figure7.eps` | eleven-season 𝓡_j(t) heat map (season × week × age group) | `nhis_rt_*.csv` |

### Supplement

| Output | Content | Inputs |
|---|---|---|
| `FigureS1.eps` | generation-interval discretisation (Erlang, K = 16) and influenza serial-interval reference | computed in `figures.R` |
| `FigureS2.eps` | particle filter vs particle smoother on step / sinusoidal 𝓡_t | `pfps_uniform_curves.csv` |
| `FigureS3.eps` | path degeneracy: unique ancestors by window length | `ps_degeneracy.csv` |
| `FigureS4.eps` | noisy synthetic data: Poisson vs NB(30) observation model | `ps_weight_curves_r25.csv`, `pre_confirm_noisy_r25.csv` |
| `FigureS5.eps` | NHIS 2018–19: NB(30) vs Poisson, medians only | `nhis_rt_2018.csv`, `nhis_rj_ma_pois_2018.csv` |
| `FigureS6.eps` | 7-day moving average vs raw daily counts, all ages | `nhis_rt_2018.csv`, `nhis_rj_raw_2018.csv` |
| `FigureS7.eps` | same comparison by age group, with 95% CrI | as above |
| `FigureS8.eps` | the three contact matrices | `data/contact_prem{2017,2021}.csv` |
| `FigureS9.eps` | 2018–19 𝓡_j(t) under the three contact matrices | `nhis_rj_prem{2017,2021}_2018.csv` |
| `FigureS10.eps` | 2018–19 with healthcare-seeking and pre-existing-immunity corrections | `nhis_overlay_imhs_2018.csv` |

The manuscript sources exactly these file names, so copying `output/*.eps` next to the
`.tex` file is sufficient.

## Table map

`code/tables.R` regenerates the numeric supplement tables directly from `results/`:

| Output | Content | Input |
|---|---|---|
| `TableS1.csv` | fixed window length L: accuracy, coverage, day-1 unique ancestors | `ps_metrics.csv`, `ps_degeneracy.csv` |
| `TableS2.csv` | observation-likelihood misspecification (r_gen ∉ estimation grid) | `ps_weight_misspec.csv` |
| `TableS3.csv` | `Np` × `σ_rw` joint sweep | `ps_metrics.csv` |
| `TableS4.csv` | resampling schedule: every step vs adaptive | `ps_metrics.csv` |
| `TableS6.csv` | particle smoother vs renewal equation (MAPE, RMSE, lag) | `compare_metrics.csv` |
| `TableS7.csv` | initial estimation bias over days 1–5 | `compare_initbias.csv` |

Table 1 (main text) and Table S5 list configuration choices rather than computed
numbers; the synthetic parameters of Table S5 are defined at the top of
`code/01_synthetic.jl`.

---

## Conventions worth knowing before re-running

* **ISO/KDCA week 53.** Three seasons (2011–12, 2016–17, 2022–23) contain a week 53.
  Weeks 52 and 53 are averaged day-of-week by day-of-week into week 52 and week 53 is
  dropped, so that all seasons align to 364 days (`merge_w53()` in `figures.R`).
* **Time-index convention for the lag.** Estimates refer to the state at the start of a
  day, whereas the truth and the observed counts `Y_t` accumulate to the end of the day.
  Both methods therefore carry the same −1 day offset, and Table S6 reports lags after a
  common +1 day realignment of the raw values in `compare_metrics.csv`.
* **Generation interval.** Continuous Erlang(2, 1.8) → rounded to days → day 0 removed →
  renormalised, truncated at K = 16 days (discrete mean 3.70 days).
* **Heat map colours.** Diverging scale centred on 𝓡 = 1 with limits [0, 3] and colour
  stops at 0 / 0.70 / 0.90 / 0.97 / 1.03 / 1.10 / 1.35 / 3.0; the grey band spans
  0.97–1.03 so that 𝓡 ≈ 1 is visible as a band rather than a single contour.
* **Figure S5** deliberately omits credible intervals: the point of the panel is that the
  Poisson median repeatedly crosses the NB(30) median, which the bands would obscure.

---

## Data

See [`data/README.md`](data/README.md) for the provenance and licence of every input
file, and for how to rebuild the Prem contact matrices from their original releases.

## Citing

If you use this code, please cite the paper and the archived release (see
`CITATION.cff`). Released under the MIT licence — see `LICENSE`.
