#!/usr/bin/env bash
# age-rt-smoother — full reproduction chain (raw data → figures & tables)
#
#   bash run_all.sh              everything (about 19 hours)
#   bash run_all.sh figures      figures and tables only, reusing results/ (about 3 minutes)
#
# Each step reads results/*.csv written by the previous one, so the order matters.
set -euo pipefail
cd "$(dirname "$0")"
JL="julia --threads=auto --project=."

if [ "${1:-all}" != "figures" ]; then

  echo "== 00  NHIS workbook -> tidy CSV =========================================="
  $JL code/00_prepare_nhis.jl

  echo "== 01  synthetic epidemic (pre-defined R_j) =============================="
  $JL code/01_synthetic.jl

  echo "== 02a particle smoother sweep — clean data, Poisson weight (reference) =="
  $JL code/02_particle_smoother.jl --mode sweep
  for f in metrics baseline_curves ess degeneracy adaptive; do
    cp "results/ps_$f.csv" "results/ps_${f}_clean_poisson.csv"      # kept before the NB sweep overwrites it
  done

  echo "== 02b particle smoother sweep — clean data, NB(30) weight (main) ========"
  $JL code/02_particle_smoother.jl --mode sweep --weight-r 30

  echo "== 02c observation-model sweeps =========================================="
  $JL code/02_particle_smoother.jl --mode weight_sweep
  $JL code/02_particle_smoother.jl --mode weight_misspec
  $JL code/02_particle_smoother.jl --mode save_noisy
  $JL code/02_particle_smoother.jl --mode weight_curves_r25
  for rg in 15 50; do
    $JL code/02_particle_smoother.jl --mode baseline_fit --noise-r $rg
  done

  echo "== 02d noisy sweep (supplementary robustness) ============================"
  $JL code/02_particle_smoother.jl --mode sweep --noise-r 25

  echo "== 02e particle filter vs particle smoother =============================="
  $JL code/02_particle_smoother.jl --mode pfps
  $JL code/05_pfps_uniform.jl

  echo "== 03  renewal-equation comparison ======================================="
  $JL code/03_renewal.jl

  echo "== 04  NHIS, 11 seasons + sensitivity analyses ==========================="
  $JL code/04_nhis.jl                                            # baseline + renewal
  $JL code/04_nhis.jl --contact prem2017
  $JL code/04_nhis.jl --contact prem2021
  $JL code/04_nhis.jl --healthcare-seeking
  $JL code/04_nhis.jl --immunity
  $JL code/04_nhis.jl --raw
  $JL code/04_nhis.jl --immunity-hseek --seasons 2018
  $JL code/04_nhis.jl --raw --weight-r 0 --tag raw_pois --seasons 2018
  $JL code/04_nhis.jl --raw --window 0   --tag raw_pf   --seasons 2018
  $JL code/04_nhis.jl        --weight-r 0 --tag ma_pois  --seasons 2018

fi

echo "== figures and tables ===================================================="
Rscript code/figures.R
Rscript code/tables.R
echo "done → output/"
