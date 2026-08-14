## TableS1 — Fixed window length L (Np = 1e5, sigma_rw = 0.15, every-step)

| L | Rj MAPE (%) | 95% coverage (%) | uniq_day1 |
|---|---|---|---|
| 16 | 1.24 | 100 | 2616 |
| 24 | 1.25 | 100 | 2019 |
| 32 | 1.26 | 100 | 1643 |
| whole history | 1.26 | 100 |  454 |

## TableS2 — Observation-likelihood misspecification (MAPE % and 95% coverage %)

| weight | mape_rgen15 | mape_rgen25 | mape_rgen50 | cov_rgen15 | cov_rgen25 | cov_rgen50 |
|---|---|---|---|---|---|---|
| Poisson | 54.74 | 52.58 | 40.67 |   7.6 |   8.2 |  13.8 |
| NB40 |  5.81 |  4.86 |  4.40 |  95.4 |  99.1 |  99.3 |
| NB30 |  5.33 |  4.49 |  4.17 |  98.0 |  99.9 |  99.7 |
| NB20 |  4.93 |  4.00 |  3.85 |  99.9 | 100.0 | 100.0 |
| NB10 |  4.33 |  3.32 |  3.69 | 100.0 | 100.0 | 100.0 |

## TableS3 — Np x sigma_rw joint sweep — cell = Rj MAPE (%) / 95% coverage (%)

| sigma_rw | Np=1e4 | Np=1e5 | Np=1e6 |
|---|---|---|---|
| 0.10 | 2.03 / 96.6 | 1.42 / 98.1 | 1.10 / 99.9 |
| 0.15 | 1.38 / 100.0 | 1.24 / 100.0 | 1.21 / 100.0 |
| 0.20 | 1.48 / 100.0 | 1.45 / 100.0 | 1.43 / 100.0 |
| 0.25 | 1.69 / 100.0 | 1.67 / 100.0 | 1.67 / 100.0 |
| 0.30 | 1.95 / 100.0 | 1.91 / 100.0 | 1.91 / 100.0 |

## TableS4 — Resampling schedule (Np = 1e5, L = 16, sigma_rw = 0.15)

| resample | Rj MAPE (%) | 95% coverage (%) |
|---|---|---|
| every step | 1.24 | 100 |
| adaptive (ESS < Np/2) | 1.83 | 100 |

## TableS6 — Particle smoother vs renewal equation, full period (day 1-149)

| method | MAPE (%) | RMSE | Lag (days) |
|---|---|---|---|
| Renewal equation | 6.45 | 0.402 | 6.8 |
| Particle smoother | 1.47 | 0.024 | 0.0 |

## TableS7 — Initial estimation bias (relative error over days 1-5)

| method | MAPE day 1-5 (%) |
|---|---|
| Renewal equation | 177.04 |
| Particle smoother |   8.08 |

## In-text values (checked against results)

- Baseline (Np = 1e5, L = 16, sigma_rw = 0.15, NB r = 30): Rj MAPE 1.24%, coverage 100%, ESS median 91,561, ESS min 3,073, runtime 49 s

