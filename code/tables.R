#!/usr/bin/env Rscript
# age-rt-smoother / code/tables.R
# ─────────────────────────────────────────────────────────────────────────────
# results/*.csv  →  output/TableS{1,2,3,4,6,7}.csv  +  output/tables.md
#   Regenerates the numeric supplement tables from the estimation output, so that no
#   number is transcribed by hand. Table 1 (main text) and Table S5 list configuration
#   choices rather than computed values; the synthetic parameters of Table S5 are
#   defined at the top of code/01_synthetic.jl.
suppressMessages({library(readr); library(dplyr); library(tidyr)})
.code_dir <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                      error = function(e) ".")
if (length(.code_dir) == 0 || .code_dir == "") .code_dir <- "."
ROOT <- normalizePath(file.path(.code_dir, ".."))
RES  <- file.path(ROOT, "results"); OUT <- file.path(ROOT, "output")
dir.create(OUT, showWarnings = FALSE)
rd <- function(f) read_csv(file.path(RES, f), show_col_types = FALSE)
md <- c()                                                  # accumulates tables.md
emit <- function(tag, title, df, digits = NULL) {
  if (!is.null(digits)) df <- df %>% mutate(across(where(is.numeric), ~ round(.x, digits)))
  write_csv(df, file.path(OUT, paste0(tag, ".csv")))
  md <<- c(md, paste0("## ", tag, " — ", title), "",
           paste0("| ", paste(names(df), collapse = " | "), " |"),
           paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
           apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
  cat(sprintf("wrote: %s.csv (%d rows)\n", tag, nrow(df)))
}

met <- rd("ps_metrics.csv")
BASE <- list(Np = 1e5, srw = 0.15, L = "16")               # adopted baseline

# ── Table S1 — fixed window length L (Np=1e5·σ_rw=0.15·every) ────────────────
#   day-1 unique ancestors are recomputed from ps_degeneracy.csv at Np = 1e6, averaged over age groups
deg1 <- rd("ps_degeneracy.csv") %>% filter(Np == 1e6, day == 1) %>%
  group_by(L) %>% summarise(uniq_day1 = round(mean(uniq)), .groups = "drop")
s1 <- met %>%
  filter(Np == BASE$Np, abs(sigma_rw - BASE$srw) < 1e-9, resample == "every") %>%
  transmute(L, `Rj MAPE (%)` = round(mape_Rj, 2), `95% coverage (%)` = round(cov_Rj, 1)) %>%
  left_join(deg1, by = "L") %>%
  mutate(L = factor(L, levels = c("16", "24", "32", "wh"),
                    labels = c("16", "24", "32", "whole history"))) %>% arrange(L)
emit("TableS1", "Fixed window length L (Np = 1e5, sigma_rw = 0.15, every-step)", s1)

# ── Table S2 — observation-likelihood misspecification ───────────────────────
#   the true overdispersion r_gen (15/25/50) lies outside the estimation grid
#   (Poisson, NB10-40), so the comparison is not circular
ws <- rd("ps_weight_misspec.csv") %>%
  mutate(weight = factor(weight, levels = c("Poisson", "NB40", "NB30", "NB20", "NB10")))
s2 <- ws %>%
  transmute(weight, r_gen, mape = round(mape_Rj, 2), cov = round(cov_Rj, 1)) %>%
  pivot_wider(names_from = r_gen, values_from = c(mape, cov),
              names_glue = "{.value}_rgen{r_gen}") %>% arrange(weight)
emit("TableS2", "Observation-likelihood misspecification (MAPE % and 95% coverage %)", s2)

# ── Table S3 — Np x sigma_rw joint sweep (L=16·every) ────────────────────────
s3 <- met %>% filter(L == BASE$L, resample == "every") %>%
  transmute(sigma_rw, Np = sprintf("Np=1e%d", round(log10(Np))),
            cell = sprintf("%.2f / %.1f", mape_Rj, cov_Rj)) %>%
  pivot_wider(names_from = Np, values_from = cell) %>% arrange(sigma_rw)
emit("TableS3", "Np x sigma_rw joint sweep — cell = Rj MAPE (%) / 95% coverage (%)", s3)

# ── Table S4 — resampling schedule (Np=1e5·L=16·σ_rw=0.15) ───────────────────
s4 <- met %>% filter(Np == BASE$Np, L == BASE$L, abs(sigma_rw - BASE$srw) < 1e-9) %>%
  transmute(resample = ifelse(resample == "every", "every step", "adaptive (ESS < Np/2)"),
            `Rj MAPE (%)` = round(mape_Rj, 2), `95% coverage (%)` = round(cov_Rj, 1))
emit("TableS4", "Resampling schedule (Np = 1e5, L = 16, sigma_rw = 0.15)", s4)

# ── Table S6 — PS vs renewal (full period, 6 age groups averaged) ────────────
#   NOTE lag is reported after the day-index convention correction (+1): in the raw CSV the
#   estimate refers to the state at the start of a day whereas the truth and the observed Y_t
#   accumulate to the end of it, so both methods carry the same -1 day offset.

cm <- rd("compare_metrics.csv") %>% filter(age_group == "all", convention == "full")
s6 <- cm %>% transmute(method = ifelse(method == "renewal", "Renewal equation", "Particle smoother"),
                       `MAPE (%)` = round(mape, 2), RMSE = round(rmse, 3),
                       `Lag (days)` = round(lag + 1, 1))
emit("TableS6", "Particle smoother vs renewal equation, full period (day 1-149)", s6)

# ── Table S7 — initial estimation bias (first 5 days) ────────────────────────
ib <- rd("compare_initbias.csv") %>% filter(age_group == "all")
s7 <- ib %>% transmute(method = ifelse(method == "renewal", "Renewal equation", "Particle smoother"),
                       `MAPE day 1-5 (%)` = round(mape_d1_5, 2))
emit("TableS7", "Initial estimation bias (relative error over days 1-5)", s7)

# ── In-text values quoted outside the tables — printed for checking ──────────
b <- met %>% filter(Np == BASE$Np, L == BASE$L, abs(sigma_rw - BASE$srw) < 1e-9, resample == "every")
md <- c(md, "## In-text values (checked against results)", "",
        sprintf("- Baseline (Np = 1e5, L = 16, sigma_rw = 0.15, NB r = 30): Rj MAPE %.2f%%, coverage %.0f%%, ESS median %s, ESS min %s, runtime %.0f s",
                b$mape_Rj, b$cov_Rj, format(round(b$ess_med), big.mark = ","),
                format(round(b$ess_min), big.mark = ","), b$runtime_s), "")
writeLines(md, file.path(OUT, "tables.md"))
cat("wrote: tables.md\n")
