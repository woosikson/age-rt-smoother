#!/usr/bin/env Rscript
# age-rt-smoother / code/figures.R
# ─────────────────────────────────────────────────────────────────────────────
# results/*.csv  →  output/Figure{2..7}.eps · output/Figure{S1..S10}.eps
#   The estimation is done by code/01_synthetic.jl, 02_particle_smoother.jl, 03_renewal.jl,
#   04_nhis.jl and 05_pfps_uniform.jl; this script only draws.
#
# Figure conventions
#   - Arial 8-12 pt, horizontal grid only, no title inside the figure (titles live in captions)
#   - colour-blind safe Okabe-Ito age colours, grey dashed reference line at R = 1, panels labelled A/B
#   - the reproduction number is set as script R (U+211B) with italic subscripts, rendering exactly
#     like the LaTeX $\mathcal{R}_j(t)$
#   - uncertainty: smoother = 95% credible interval; renewal equation = point estimate
# ─────────────────────────────────────────────────────────────────────────────
suppressMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(cowplot); library(scales); library(ggtext)
})
if (!interactive()) grDevices::pdf(NULL)
.code_dir <- tryCatch(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE))),
          error = function(e) ".")
if (length(.code_dir) == 0 || .code_dir == "") .code_dir <- "."
ROOT   <- normalizePath(file.path(.code_dir, ".."))
here   <- ROOT                                   # the blocks below read from file.path(here, "results")
outdir <- file.path(ROOT, "output"); dir.create(outdir, showWarnings = FALSE)
FONT <- "Arial"

# ── per-figure type scaling ────────────────────────────────────────────────
#   Each figure is inserted into the manuscript at a different width, so a single point size on
#   the canvas prints at a different size in the article. FIG_K rescales the type of one figure
#   (and nothing else: line widths, point sizes and the canvas are untouched) so that every
#   figure prints at the same nominal sizes: axis title 9 pt, ticks and legend 8 pt, panel tag 12 pt.
FIG_K <- 1
#   Per-figure point offsets on top of the shared rule, for the few panels where the common
#   sizes do not fit (set right after use_k(); use_k() resets them).
#   rmult overrides the enlargement of the axis titles (R_LAB_MULT); set it to 1 where the panel
#   subtitles must stay the largest element of the panel.
FIG_D <- list(title = 0, text = 0, strip = 0, legend = 0, rmult = NA)
use_k <- function(k, title = 0, text = 0, strip = 0, legend = 0, rmult = NA) {
  FIG_K <<- k
  FIG_D <<- list(title = title, text = text, strip = strip, legend = legend, rmult = rmult)
  invisible(k)
}
r_mult <- function() if (is.na(FIG_D$rmult)) R_LAB_MULT else FIG_D$rmult

# ── R_j(t) label ──────────────────────────────────────────────────────────
#   The script capital R (U+211B) is drawn small by design, so labels containing it are enlarged.
#   Set with plotmath, not ggtext. ggtext can enlarge the script R on its own (font-size in an
#   inline span), but it cannot control the spacing that follows it: the advance width of the
#   U+211B fallback glyph is narrower than its ink, and every remedy is ignored inside a span —
#   hair/thin/en/non-breaking spaces are all collapsed, and font-family is overridden by the
#   fallback. The residual gap then differs between horizontal and rotated labels, which is worse
#   than the original problem. plotmath places the subscript correctly and identically in both
#   orientations, so the label is set with plotmath and the WHOLE label is enlarged instead
#   (R_LAB_MULT), which is the fallback allowed when per-glyph control is not available.
R_LAB_MULT <- 1.25
LAB_RJT   <- function() expression("ℛ"[italic(j)](italic(t)))
LAB_RJT_H <- function() expression("ℛ"[italic(j)](italic(t)))
LAB_RT    <- function() expression("ℛ"[italic(t)])
r_title   <- function() element_text(size = 9 * FIG_K * R_LAB_MULT, margin = margin(r = 6))

# ── figure numbering: source stem → figure number in the paper ──────────────
#   A stem absent from this table is not used in the paper, so no file is written for it.
PAPER_FIG <- c(
  fig2_true                          = "Figure2",
  fig_ps_combined                    = "Figure3",
  fig_nhis_overlay_2018              = "Figure4",
  fig_nhis_2018                      = "Figure5",
  fig_nhis_overlay_10seasons         = "Figure6",
  fig_nhis_heatmap                   = "Figure7",
  gi_discretization                  = "FigureS1",
  fig_pfps_uniform                   = "FigureS2",
  fig_ps_degeneracy                  = "FigureS3",
  fig_weight_curves_r25              = "FigureS4",
  fig_nhis_ma_ps_poisson_2018        = "FigureS5",
  fig_nhis_raw_overlay_2018          = "FigureS6",
  fig_nhis_raw_cri_2018              = "FigureS7",
  fig_contact_matrices               = "FigureS8",
  fig_nhis_cm_2018                   = "FigureS9",
  fig_nhis_immun_hseek_overlay_2018  = "FigureS10")

# ── axis tick numbers: thousands separated by commas everywhere ─────────────
#   scales::label_number() defaults to a SPACE as big.mark, so an axis without an explicit labels=
#   would print "20 000". The scale constructors are wrapped to change that default; an explicit
#   labels= still wins, and drop0trailing keeps decimal axes unchanged.
lab_comma <- function(x) {
  o <- format(x, big.mark = ",", trim = TRUE, scientific = FALSE, drop0trailing = TRUE)
  o[is.na(x)] <- NA; o
}
scale_y_continuous <- function(..., labels = lab_comma) ggplot2::scale_y_continuous(..., labels = labels)
scale_x_continuous <- function(..., labels = lab_comma) ggplot2::scale_x_continuous(..., labels = labels)
scale_y_log10 <- function(..., labels = lab_comma) ggplot2::scale_y_log10(..., labels = labels)  # same for log axes
scale_x_log10 <- function(..., labels = lab_comma) ggplot2::scale_x_log10(..., labels = labels)

# ── Okabe-Ito palette (colour-blind safe) ───────────────────────────────────
OKABE <- c(orange="#E69F00", skyblue="#56B4E9", green="#009E73", yellow="#F0E442",
           blue="#0072B2", vermillion="#D55E00", purple="#CC79A7", black="#000000")
REF_COL <- "grey50"                          # R = 1 reference line (neutral grey dashed)
PS_COL  <- unname(OKABE["blue"])             # particle smoother
RN_COL  <- unname(OKABE["green"])            # renewal equation
OBS_COL <- "grey35"                          # observed data

AGE  <- c("0-5", "6-11", "12-17", "18-44", "45-64", "65+")
# young groups get warm accent colours, adults grey (de-emphasised), the elderly blue
age_cols <- c("0-5"="#D55E00", "6-11"="#E69F00", "12-17"="#CC79A7",
              "18-44"="#BDBDBD", "45-64"="#737373", "65+"="#0072B2")
age_lw <- c("0-5"=0.7, "6-11"=0.7, "12-17"=0.7, "18-44"=0.45, "45-64"=0.45, "65+"=0.7)

# ── helpers for the manual legend: the smoother median line and the 95% CrI patch are drawn as
#    six age-coloured segments, so a grey key is not mistaken for the grey adult groups
leg_line6 <- function(x0, x1, y, lw = 0.9) {                 # a line made of the six age colours, left to right
  w <- (x1 - x0) / 6
  lapply(seq_len(6), function(i) annotate("segment", x = x0 + (i - 1) * w, xend = x0 + i * w,
    y = y, yend = y, colour = unname(age_cols[AGE[i]]), linewidth = lw))
}
leg_rect6 <- function(x0, x1, ylo, yhi, alpha = 0.55) {      # a rectangle made of the six age colours, left to right
  w <- (x1 - x0) / 6
  lapply(seq_len(6), function(i) annotate("rect", xmin = x0 + (i - 1) * w, xmax = x0 + i * w,
    ymin = ylo, ymax = yhi, fill = unname(age_cols[AGE[i]]), alpha = alpha))
}

# ── manual horizontal legend ───────────────────────────────────────────────
#   Some legends need keys that a normal ggplot guide cannot draw (a six-colour stripe for the
#   per-age median, a six-colour band for the credible interval), so they are drawn by hand.
#   Label widths are measured on the device and the row is laid out in inches, shrinking the key
#   and the separators (and only then the type) until the row fits the figure width.
leg_row <- function(items, size_pt, width_in) {
  cexf <- size_pt / 12                                   # strwidth cex is relative to pointsize 12
  lw   <- vapply(items, function(it) strwidth(it$label, units = "inches", cex = cexf, family = FONT), 0)
  key <- 0.26; gap <- 0.06; sep <- 0.26                  # key length, key-to-label gap, item separation
  while (sum(key + gap + lw + sep) - sep > width_in * 0.98 && key > 0.13) { key <- key - 0.01; sep <- sep - 0.01 }
  while (sum(key + gap + lw + sep) - sep > width_in * 0.98 && size_pt > 7) {
    size_pt <- size_pt - 0.5; cexf <- size_pt / 12
    lw <- vapply(items, function(it) strwidth(it$label, units = "inches", cex = cexf, family = FONT), 0)
  }
  p <- ggplot() + coord_cartesian(xlim = c(0, width_in), ylim = c(0, 1), expand = FALSE) +
    theme_void(base_family = FONT) +
    theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(2, 2, 2, 2))
  x <- (width_in - (sum(key + gap + lw + sep) - sep)) / 2            # centre the row in the figure
  for (i in seq_along(items)) {
    it <- items[[i]]
    p <- p + switch(it$kind,
      line  = annotate("segment", x = x, xend = x + key, y = 0.5, yend = 0.5,
                       colour = it$colour, linewidth = it$lwd,
                       linetype = if (is.null(it$lty)) "solid" else it$lty),
      point = annotate("point", x = x + key / 2, y = 0.5, colour = it$colour, size = 1.2),
      rect  = annotate("rect", xmin = x, xmax = x + key, ymin = 0.34, ymax = 0.66,
                       fill = it$colour, alpha = 0.5),
      line6 = leg_line6(x, x + key, 0.5),
      rect6 = leg_rect6(x, x + key, 0.34, 0.66))
    p <- p + annotate("text", x = x + key + gap, y = 0.5, label = it$label,
                      hjust = 0, size = size_pt / .pt, family = FONT)
    x <- x + key + gap + lw[i] + sep
  }
  p
}

# ── shared theme (Arial, horizontal grid, no title inside the panel) ────────
pub_theme <- function(base = 9)
  theme_minimal_hgrid(font_size = 9 * FIG_K, font_family = FONT) +
  theme(plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.title       = element_text(size = 9 * FIG_K + FIG_D$title),
        # The y title is almost always R_j(t); it is enlarged as a whole (see LAB_RJT above).
        # Plain-text y titles ("Cases (per 100,000)") are enlarged by the same factor, which keeps
        # the two columns of the combined panels visually consistent.
        axis.title.y     = element_text(size = 9 * FIG_K * r_mult() + FIG_D$title, margin = margin(r = 6)),
        axis.text        = element_text(size = 8 * FIG_K + FIG_D$text),
        strip.text       = element_text(size = 9 * FIG_K + FIG_D$strip, face = "bold"),
        legend.title     = element_blank(),
        legend.key.height = unit(9, "pt"),
        legend.text      = element_text(size = 8 * FIG_K + FIG_D$legend),
        legend.position  = "top", legend.justification = "center",
        legend.margin    = margin(0, 0, 0, 0),
        plot.margin      = margin(3, 8, 3, 3))

# ── export: only the figures used in the paper, as vector EPS ──────────────
#   grDevices::cairo_ps is required: the classic postscript() device cannot resolve Arial and
#   mishandles unicode symbols such as the script R (U+211B); cairo_ps handles both.
save_wilke <- function(plot, stem, width, height, tiff = TRUE, png = FALSE) {
  nm <- unname(PAPER_FIG[basename(stem)])
  if (is.na(nm)) return(invisible(NULL))          # not used in the paper: nothing is written
  f <- file.path(outdir, paste0(nm, ".eps"))
  ggplot2::ggsave(f, plot, device = grDevices::cairo_ps, width = width, height = height)
  cat(sprintf("saved: %-10s <- %s (%.1f x %.1f in)\n", paste0(nm, ".eps"), basename(stem), width, height))
}

# ═════════════════════════════════════════════════════════════════════
#  Figure 2 — synthetic data: (A) incidence by age group, (B) true R_j(t)
# ═════════════════════════════════════════════════════════════════════
conf <- read_csv(file.path(here, "results", "pre_confirm.csv"), show_col_types = FALSE) %>%
  mutate(age_group = factor(age_group, levels = AGE))
rt <- read_csv(file.path(here, "results", "pre_rt.csv"), show_col_types = FALSE) %>%
  mutate(age_group = factor(age_group, levels = AGE))

age_scale <- list(
  scale_colour_manual(values = age_cols, breaks = AGE),
  scale_linewidth_manual(values = age_lw, guide = "none"),
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.1))))

use_k(1.60, legend = 1.75)                    # Figure 2
pA <- ggplot(conf, aes(day, incidence_per100k, colour = age_group, linewidth = age_group)) +
  geom_line() + age_scale +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     labels = scales::label_comma()) +
  labs(x = NULL, y = NULL) +
  pub_theme() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

pB <- ggplot(rt, aes(day, Rt, colour = age_group, linewidth = age_group)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.4) +
  geom_line() + age_scale +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, max(3, max(rt$Rt) * 1.05))) +
  labs(x = "Day", y = NULL) +
  pub_theme() + theme(legend.position = "none")

# The two panels have tick labels of very different widths ("4,000" vs "3"), and stacked panels
# are aligned on the panel edge, so a y title left inside each panel would sit at a different x.
# Both titles are therefore drawn in one shared column, which puts them on the same vertical line.
# The spacer under the lower title offsets the x axis of panel B, so each title stays centred on
# its own panel.
y_title <- function(lab) ggdraw() +
  draw_label(lab, angle = 90, fontfamily = FONT, size = 9 * FIG_K * r_mult() + FIG_D$title)
fig2_body <- plot_grid(pA, pB, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1),
                       labels = c("A", "B"), label_fontfamily = FONT, label_size = 12 * FIG_K,
                       label_fontface = "bold")
fig2_titles <- plot_grid(y_title("Cases (per 100,000)"),
                         plot_grid(y_title(LAB_RJT()), NULL, ncol = 1, rel_heights = c(1, 0.16)),
                         ncol = 1, rel_heights = c(1, 1))
# the empty middle column is the gap between the titles and panel A's tick labels
fig2 <- plot_grid(fig2_titles, NULL, fig2_body, nrow = 1, rel_widths = c(0.045, 0.022, 1))
save_wilke(fig2, file.path(outdir, "fig2_true"), width = 6.5, height = 6.9)  # each sub-panel about 2.0:1

# ═════════════════════════════════════════════════════════════════════
#  gi_discretization — (A) Erlang(σ=γ=1/1.8) GI + (B) Weibull(flu)
# ═════════════════════════════════════════════════════════════════════
use_k(1.08)                                   # Figure S1
{
  Kg <- 16; xg <- seq(0, 16, 0.05)
  round_disc <- function(Fc) { w <- numeric(Kg); w[1] <- Fc(1.5) - Fc(0.5)
    for (d in 2:Kg) w[d] <- Fc(d + 0.5) - Fc(d - 0.5); data.frame(d = 1:Kg, w = w / sum(w)) }
  scE <- 1.8
  c1 <- data.frame(x = xg, pdf = dgamma(xg, 2, scale = scE))
  d1 <- round_disc(function(x) ifelse(x <= 0, 0, pgamma(x, 2, scale = scE)))
  mW <- 3.6; sW <- 1.6
  kW <- uniroot(function(k) gamma(1 + 2/k)/gamma(1 + 1/k)^2 - 1 - (sW/mW)^2, c(0.5, 20))$root
  lamW <- mW / gamma(1 + 1/kW)
  c2 <- data.frame(x = xg, pdf = dweibull(xg, kW, lamW))
  ymx <- 1.05 * max(c1$pdf, c2$pdf, d1$w)
  gi_axis <- list(scale_x_continuous(breaks = seq(0, 16, 2), expand = expansion(mult = c(0.01, 0.02))),
                  scale_y_continuous(limits = c(0, ymx), expand = expansion(mult = c(0.05, 0))),  # lower margin so that points near zero stay visible
                  coord_cartesian(xlim = c(0, 16)))
  gA <- ggplot() +
    geom_col(data = d1, aes(d, w), width = 1, fill = OKABE["skyblue"], alpha = 0.30) +
    geom_line(data = c1, aes(x, pdf), colour = "black", linewidth = 0.7) +
    geom_vline(xintercept = 2 * scE, linetype = "dotted", colour = "grey45", linewidth = 0.5) +
    geom_point(data = d1, aes(d, w), colour = PS_COL, size = 1.3) +
    gi_axis + labs(x = NULL, y = "Density / prob.") + pub_theme()
  gB <- ggplot() +
    geom_line(data = c2, aes(x, pdf), colour = "black", linewidth = 0.7) +
    geom_vline(xintercept = mW, linetype = "dotted", colour = "grey45", linewidth = 0.5) +
    gi_axis + labs(x = "Generation interval (days)", y = "Density") + pub_theme()
  gi <- plot_grid(gA, gB, ncol = 1, align = "v", axis = "lr", labels = c("A", "B"),
                  label_fontfamily = FONT, label_size = 12 * FIG_K, label_fontface = "bold")
  save_wilke(gi, file.path(outdir, "gi_discretization"), width = 6.0, height = 4.8)
}

# ═════════════════════════════════════════════════════════════════════
#  Particle-smoother figures (drawn when the ps_*.csv files are present)
# ═════════════════════════════════════════════════════════════════════
psf <- file.path(here, "results", "ps_metrics.csv")
if (file.exists(psf)) {
  met <- read_csv(psf, show_col_types = FALSE)
  cur <- read_csv(file.path(here, "results", "ps_baseline_curves.csv"), show_col_types = FALSE) %>%
    mutate(age_group = factor(age_group, levels = AGE))
  deg <- read_csv(file.path(here, "results", "ps_degeneracy.csv"), show_col_types = FALSE)
  adp <- read_csv(file.path(here, "results", "ps_adaptive.csv"), show_col_types = FALSE) %>%
    mutate(age_group = factor(age_group, levels = AGE))
  bNp <- 1e5                                       # adopted particle count
  cb <- filter(cur, Np == bNp)
  # the R_j and dR fits appear in the combined panel below

  # (3) Np-convergence
  mc <- met %>% filter(L == "16", sigma_rw == 0.15, resample == "every") %>%
    select(Np, mape_Rj, mape_dR) %>% pivot_longer(-Np, names_to = "metric", values_to = "mape")
  p3 <- ggplot(mc, aes(Np, mape, colour = metric)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
    scale_x_log10(breaks = c(1e4, 1e5, 1e6), labels = c("1e4", "1e5", "1e6")) +
    scale_colour_manual(values = c(mape_Rj = PS_COL, mape_dR = unname(OKABE["vermillion"])),
                        labels = c(mape_Rj = expression("ℛ"[italic(j)]~MAPE),
                                   mape_dR = expression(Delta*italic(R)~MAPE))) +
    labs(x = expression("Number of particles"~italic(N)[p]), y = "MAPE (%)") + pub_theme()
  save_wilke(p3, file.path(outdir, "fig_ps_convergence"), width = 4.4, height = 3.2)

  # sigma_rw sensitivity (MAPE and 95% coverage on a dual axis)
  ms <- met %>% filter(Np == bNp, L == "16", resample == "every") %>% select(sigma_rw, mape_Rj, cov_Rj)
  p4 <- ggplot(ms, aes(sigma_rw)) +
    geom_line(aes(y = mape_Rj, colour = "MAPE"), linewidth = 0.6) +
    geom_point(aes(y = mape_Rj, colour = "MAPE"), size = 1.8) +
    geom_line(aes(y = cov_Rj / 20, colour = "95% coverage"), linewidth = 0.6) +
    geom_point(aes(y = cov_Rj / 20, colour = "95% coverage"), size = 1.8) +
    scale_y_continuous(name = expression("ℛ"[italic(j)]~"MAPE (%)"),
                       sec.axis = sec_axis(~ . * 20, name = "Coverage (%)")) +
    scale_colour_manual(values = c("MAPE" = PS_COL, "95% coverage" = unname(OKABE["vermillion"]))) +
    scale_x_continuous(breaks = c(0.1, 0.15, 0.2, 0.25, 0.3)) +
    labs(x = expression(sigma[rw])) + pub_theme()
  save_wilke(p4, file.path(outdir, "fig_ps_sigma"), width = 4.4, height = 3.2)

  # path degeneracy: unique ancestors over time, one line per window length
  degNp <- max(deg$Np)
  dd <- deg %>% filter(Np == degNp) %>% group_by(L, day) %>%
    summarise(uniq = mean(uniq), .groups = "drop") %>%
    mutate(L = factor(L, levels = c("16", "24", "32", "wh")))
  use_k(1.08)                                 # Figure S3
  p5 <- ggplot(dd, aes(day, uniq, colour = L)) + geom_line(linewidth = 0.6) +
    scale_colour_viridis_d(option = "C", end = 0.9) +
    scale_y_log10(labels = lab_comma) +
    labs(x = "Day", y = "Unique ancestors (mean over age)") +
    pub_theme() + theme(legend.position = "right")
  save_wilke(p5, file.path(outdir, "fig_ps_degeneracy"), width = 6.0, height = 3.4)

  # (6) adaptive resampling monitor (SI/diagnostic)
  am <- adp %>% filter(Np == max(adp$Np), age_group == "6-11"); Np_a <- max(adp$Np)
  p6 <- ggplot(am, aes(day, ess)) +
    geom_hline(yintercept = Np_a / 2, linetype = "dashed", colour = unname(OKABE["vermillion"])) +
    annotate("text", x = 5, y = Np_a/2, label = "ESS = N/2", vjust = -0.4, hjust = 0,
             colour = unname(OKABE["vermillion"]), size = 3, family = FONT) +
    geom_line(colour = "grey40", linewidth = 0.5) +
    geom_point(data = filter(am, resampled), aes(day, ess), colour = PS_COL, size = 1.2) +
    labs(x = "Day", y = "ESS") + pub_theme()
  save_wilke(p6, file.path(outdir, "fig_ps_adaptive_monitor"), width = 6.5, height = 3.0)

  # rmult = 1: the grid is dense, so the y title is not enlarged here. Both axis titles are then at
  # the common 9 pt nominal, and the age subtitles are set to the same size (see mk_panel below).
  use_k(1.35, rmult = 1)                      # Figure 3
  # ── combined 4x3 panel: age group x {cases (top), R_j (bottom)}, in two blocks of three ──
  #   age colour = smoother median with 95% CrI band, truth in black; cases share one range
  # age-group population, recovered from the incidence and per-100k columns to avoid drift
  NN <- conf %>% filter(incidence > 0) %>% group_by(age_group) %>%
    summarise(N = median(incidence * 1e5 / incidence_per100k), .groups = "drop")
  AGE_N <- setNames(NN$N, as.character(NN$age_group))       # cases are shown per 100,000 of each age group
  cases_lim <- c(0, max(cb$dR_hi * 1e5 / AGE_N[as.character(cb$age_group)]))  # includes the dR credible interval
  rj_top <- c(min(filter(cb, age_group %in% AGE[1:3])$Rj_lo), max(filter(cb, age_group %in% AGE[1:3])$Rj_hi))
  rj_bot <- c(min(filter(cb, age_group %in% AGE[4:6])$Rj_lo), max(filter(cb, age_group %in% AGE[4:6])$Rj_hi))
  rr_all <- if (file.exists(file.path(here, "results", "renewal_rt.csv")))   # renewal estimate, overlaid
    read_csv(file.path(here, "results", "renewal_rt.csv"), show_col_types = FALSE) %>%
      mutate(age_group = factor(age_group, levels = AGE)) else NULL
  mk_panel <- function(ag, qty, rj_lim, show_ytitle, show_yticks, show_x, title = NULL) {
    col <- unname(age_cols[ag]); d <- filter(cb, age_group == ag)   # smoother: age colour solid; truth: black dotted; renewal: black dashed
    if (qty == "cases") {
      sc <- 1e5 / AGE_N[[as.character(ag)]]                          # per 100,000 of the age group
      p <- ggplot(d, aes(day)) +
        geom_ribbon(aes(ymin = dR_lo * sc, ymax = dR_hi * sc), fill = col, alpha = 0.18) +  # cases 95% CrI(coherent Ih)
        geom_line(aes(y = dR_med * sc), colour = col, linewidth = 0.5) +
        geom_line(aes(y = true_dR * sc), colour = "black", linewidth = 0.4, linetype = "dotted") +
        scale_y_continuous(labels = scales::label_comma()) +
        coord_cartesian(ylim = cases_lim) +
        labs(x = NULL, y = if (show_ytitle) "Cases (per 100,000)" else NULL)
    } else {
      p <- ggplot(d, aes(day)) +
        geom_ribbon(aes(ymin = Rj_lo, ymax = Rj_hi), fill = col, alpha = 0.18) +
        geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.35)
      if (!is.null(rr_all))                                  # renewal 𝓡_j overlay (black dashed)
        p <- p + geom_line(data = filter(rr_all, age_group == ag), aes(y = Rj_est),
                           colour = "black", linewidth = 0.45, linetype = "42")
      p <- p + geom_line(aes(y = Rj_med), colour = col, linewidth = 0.5) +
        geom_line(aes(y = true_Rj), colour = "black", linewidth = 0.4, linetype = "dotted") +
        coord_cartesian(ylim = rj_lim) +
        labs(x = NULL, y = if (show_ytitle) LAB_RJT() else NULL)
    }
    p <- p + scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) + pub_theme(8)
    # the age subtitle is set at the axis-title size (a fixed 9 pt would print far smaller, since
    # the canvas is reduced in the article and only sizes carrying FIG_K follow that reduction)
    if (!is.null(title)) p <- p + ggtitle(title) +
      theme(plot.title = element_text(size = 9 * FIG_K + FIG_D$title, face = "bold", hjust = 0.5,
                                      margin = margin(b = 1)))
    if (!show_yticks) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if (show_x) p <- p + labs(x = "Day")
    else p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                        axis.line.x = element_blank())     # rows 1-3: drop the x axis line, ticks and labels
    p
  }
  b1 <- AGE[1:3]; b2 <- AGE[4:6]
  # The four y titles are drawn once in a column of their own so that they share one baseline: a
  # title attached to its own panel sits next to its own tick labels, which are wider in the cases
  # rows ("4,000") than in the R_j rows ("3"), and the two would not line up.
  # The grid is then built column by column rather than as one 4 x 3 block. Aligning all twelve
  # panels at once makes every cell reserve the axis width of column 1, which leaves columns 2 and 3
  # with a wide empty gutter; per-column assembly aligns the panels within a column and gives the
  # spare width back to the panels (rel_widths compensates for column 1 carrying the tick labels).
  mkcol <- function(k) plot_grid(
    mk_panel(b1[k], "cases", rj_top, FALSE, k == 1, FALSE, b1[k]),
    mk_panel(b1[k], "rj",    rj_top, FALSE, k == 1, FALSE),
    mk_panel(b2[k], "cases", rj_bot, FALSE, k == 1, FALSE, b2[k]),
    mk_panel(b2[k], "rj",    rj_bot, FALSE, k == 1, TRUE),
    ncol = 1, align = "v", axis = "lr")
  gcols <- plot_grid(mkcol(1), mkcol(2), mkcol(3), nrow = 1, rel_widths = c(1.16, 1, 1))
  # spacers put each title beside its own panel: the cases rows carry a bold age title above them,
  # the last row the x axis below it
  ttl <- function(lab, above = 0, below = 0)
    plot_grid(NULL, y_title(lab), NULL, ncol = 1, rel_heights = c(max(above, 1e-6), 1, max(below, 1e-6)))
  gttl <- plot_grid(ttl("Cases (per 100,000)", above = 0.09), ttl(LAB_RJT()),
                    ttl("Cases (per 100,000)", above = 0.09), ttl(LAB_RJT(), below = 0.22), ncol = 1)
  gmain <- plot_grid(gttl, NULL, gcols, nrow = 1, rel_widths = c(0.030, 0.010, 1))
  # manual legend: truth (black dotted), smoother median (six colours), renewal (black dashed), CrI band
  leg <- leg_row(list(
    list(kind = "line",  label = "Truth",     colour = "black", lwd = 0.6, lty = "dotted"),
    list(kind = "line6", label = "PS median"),
    list(kind = "line",  label = "Renewal",   colour = "black", lwd = 0.7, lty = "42"),
    list(kind = "rect6", label = "95% CrI")), size_pt = 11.5, width_in = 7.5)
  fig <- plot_grid(leg, gmain, ncol = 1, rel_heights = c(0.04, 1))
  save_wilke(fig, file.path(outdir, "fig_ps_combined"), width = 7.5, height = 8.5)

  # ── same layout with the renewal overlay replaced by the particle filter (median + 95% CrI) ──
  #   the filter band is visibly wider than the smoother band, which is the point of the comparison
  pff <- file.path(here, "results", "ps_pfps_curves.csv")
  if (file.exists(pff)) {
    pf_r <- read_csv(pff, show_col_types = FALSE) %>%
      filter(regime == "clean", method == "PF") %>% mutate(age_group = factor(age_group, levels = AGE))
    PFCOL <- "grey35"; PFFILL <- "grey65"
    rjlim <- function(ags) range(c(filter(cb, age_group %in% ags)$Rj_lo, filter(cb, age_group %in% ags)$Rj_hi,
                                    filter(pf_r, age_group %in% ags)$Rj_lo, filter(pf_r, age_group %in% ags)$Rj_hi))
    rj_top2 <- rjlim(AGE[1:3]); rj_bot2 <- rjlim(AGE[4:6])
    mk_panel_pf <- function(ag, qty, rj_lim, show_ytitle, show_yticks, show_x, title = NULL) {
      col <- unname(age_cols[ag]); d <- filter(cb, age_group == ag)
      if (qty == "cases") {
        sc <- 1e5 / AGE_N[[as.character(ag)]]
        p <- ggplot(d, aes(day)) +
          geom_ribbon(aes(ymin = dR_lo * sc, ymax = dR_hi * sc), fill = col, alpha = 0.18) +
          geom_line(aes(y = dR_med * sc), colour = col, linewidth = 0.5) +
          geom_line(aes(y = true_dR * sc), colour = "black", linewidth = 0.4, linetype = "dotted") +
          scale_y_continuous(labels = scales::label_comma()) +
          coord_cartesian(ylim = cases_lim) + labs(x = NULL, y = if (show_ytitle) "Cases (per 100,000)" else NULL)
      } else {
        pf <- filter(pf_r, age_group == ag)
        p <- ggplot(d, aes(day)) +
          geom_ribbon(data = pf, aes(day, ymin = Rj_lo, ymax = Rj_hi), fill = PFFILL, alpha = 0.40, inherit.aes = FALSE) +
          geom_ribbon(aes(ymin = Rj_lo, ymax = Rj_hi), fill = col, alpha = 0.22) +
          geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.35) +
          geom_line(data = pf, aes(day, Rj_med), colour = PFCOL, linewidth = 0.45, inherit.aes = FALSE) +
          geom_line(aes(y = Rj_med), colour = col, linewidth = 0.5) +
          geom_line(aes(y = true_Rj), colour = "black", linewidth = 0.4, linetype = "dotted") +
          coord_cartesian(ylim = rj_lim) + labs(x = NULL, y = if (show_ytitle) LAB_RJT() else NULL)
      }
      p <- p + scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) + pub_theme(8)
      if (!is.null(title)) p <- p + ggtitle(title) +
        theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5, margin = margin(b = 1)))
      if (!show_yticks) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      if (show_x) p <- p + labs(x = "Day")
      else p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
      p
    }
    q1 <- lapply(seq_along(b1), function(i) mk_panel_pf(b1[i], "cases", rj_top2, i == 1, i == 1, FALSE, b1[i]))
    q2 <- lapply(seq_along(b1), function(i) mk_panel_pf(b1[i], "rj", rj_top2, i == 1, i == 1, FALSE))
    q3 <- lapply(seq_along(b2), function(i) mk_panel_pf(b2[i], "cases", rj_bot2, i == 1, i == 1, FALSE, b2[i]))
    q4 <- lapply(seq_along(b2), function(i) mk_panel_pf(b2[i], "rj", rj_bot2, i == 1, i == 1, TRUE))
    gmain2 <- plot_grid(plotlist = c(q1, q2, q3, q4), ncol = 3, align = "hv", axis = "tblr")
    leg2 <- ggplot() + coord_cartesian(xlim = c(0, 15), ylim = c(0, 1), expand = FALSE) +
      theme_void(base_family = FONT) +
      theme(plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(2, 2, 2, 2)) +
      annotate("segment", x = 0.2, xend = 0.9, y = 0.5, yend = 0.5, colour = "black", linewidth = 0.6, linetype = "dotted") +
      annotate("text", x = 1.05, y = 0.5, label = "Truth", hjust = 0, size = 2.9, family = FONT) +
      leg_line6(2.4, 3.8, 0.5) +
      annotate("text", x = 3.95, y = 0.5, label = "PS median", hjust = 0, size = 2.9, family = FONT) +
      leg_rect6(6.1, 7.4, 0.34, 0.66) +
      annotate("text", x = 7.55, y = 0.5, label = "PS 95% CrI", hjust = 0, size = 2.9, family = FONT) +
      annotate("segment", x = 10.0, xend = 10.7, y = 0.5, yend = 0.5, colour = PFCOL, linewidth = 0.7) +
      annotate("text", x = 10.85, y = 0.5, label = "PF median", hjust = 0, size = 2.9, family = FONT) +
      annotate("rect", xmin = 13.0, xmax = 13.6, ymin = 0.34, ymax = 0.66, fill = PFFILL, alpha = 0.6) +
      annotate("text", x = 13.75, y = 0.5, label = "PF CrI", hjust = 0, size = 2.9, family = FONT)
    save_wilke(plot_grid(leg2, gmain2, ncol = 1, rel_heights = c(0.04, 1)),
               file.path(outdir, "fig_pfps_combined"), width = 7.5, height = 8.5)
  }
}

# ═════════════════════════════════════════════════════════════════════
#  Observation-model sensitivity on noisy synthetic data (r_gen = 20)
#    Poisson vs NB(r = 10/20/30/40), scored against the true R_j.
#    ps_weight_sweep.csv holds the metrics, ps_weight_curves_r25.csv the representative curves.
# ═════════════════════════════════════════════════════════════════════
wsf <- file.path(here, "results", "ps_weight_sweep.csv")
if (file.exists(wsf)) {
  # (the metrics of this sweep are reported as a supplement table, not as a figure)

    use_k(1.35, title = -1.2, text = -1.2, strip = 1.2, rmult = 1)   # Figure S4
  # representative curves at sigma_rw = 0.15: cases and R_j by age group, Poisson vs NB30 vs truth,
  #   on data generated with NB(r_gen = 25) observation noise
  wcf <- file.path(here, "results", "ps_weight_curves_r25.csv")
  if (file.exists(wcf)) {
    wc <- read_csv(wcf, show_col_types = FALSE) %>%
      mutate(age_group = factor(age_group, levels = AGE),
             weight = factor(weight, levels = c("Poisson", "NB30")))
    wcol <- c("Poisson" = unname(OKABE["vermillion"]), "NB30" = PS_COL)
    NNw <- conf %>% filter(incidence > 0) %>% group_by(age_group) %>%     # age-group population, recovered as in the panel above
      summarise(N = median(incidence * 1e5 / incidence_per100k), .groups = "drop")
    AGE_Nw <- setNames(NNw$N, as.character(NNw$age_group))
    scw <- function(ag) 1e5 / AGE_Nw[[as.character(ag)]]
    wc <- wc %>% mutate(sc = 1e5 / AGE_Nw[as.character(age_group)])
    obs1 <- filter(wc, weight == "Poisson")                              # obs and true_dR do not depend on the weight; de-duplicated
    casesw_lim <- c(0, max(wc$dR_med * wc$sc, wc$true_dR * wc$sc, na.rm = TRUE) * 1.05) # shared range, based on the dR median and the truth
    rjw_top <- c(0, 3); rjw_bot <- c(0, 2.6)                             # shared R_j range for the young and the adult/elderly rows
    mk_wpanel <- function(ag, qty, rj_lim, show_ytitle, show_yticks, show_x, title = NULL) {
      d <- filter(wc, age_group == ag); sc <- scw(ag); do <- filter(obs1, age_group == ag)
      if (qty == "cases") {
        # cases: median fit, observed points and the clean truth only
        p <- ggplot(d, aes(day)) +
          geom_point(data = do, aes(y = obs * sc), colour = "grey65", size = 0.35, alpha = 0.6) +
          geom_line(aes(y = dR_med * sc, colour = weight), linewidth = 0.45) +
          geom_line(data = do, aes(y = true_dR * sc), colour = "black", linewidth = 0.4, linetype = "dotted") +
          scale_y_continuous(labels = scales::label_comma()) +
          coord_cartesian(ylim = casesw_lim) +
          labs(x = NULL, y = if (show_ytitle) "Cases (per 100,000)" else NULL)
      } else {
        p <- ggplot(d, aes(day)) +
          geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.35) +
          geom_ribbon(aes(ymin = Rj_lo, ymax = Rj_hi, fill = weight), alpha = 0.14) +
          geom_line(aes(y = Rj_med, colour = weight), linewidth = 0.45) +
          geom_line(aes(y = true_Rj), colour = "black", linewidth = 0.4, linetype = "dotted") +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +   # no lower expansion, so the x axis sits exactly at y = 0
          coord_cartesian(ylim = rj_lim) +
          labs(x = NULL, y = if (show_ytitle) LAB_RJT() else NULL)
      }
      p <- p + scale_colour_manual(values = wcol) + scale_fill_manual(values = wcol) +
        scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) + pub_theme(8) +
        theme(legend.position = "none")
      # subtitle scaled like the axis titles (as in Figure 3), plus the per-figure offset that keeps
      # it the largest element of the panel here: subtitle > axis title > tick
      if (!is.null(title)) p <- p + ggtitle(title) +
        theme(plot.title = element_text(size = 9 * FIG_K + FIG_D$strip, face = "bold", hjust = 0.5,
                                        margin = margin(b = 1)))
      if (!show_yticks) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      if (show_x) p <- p + labs(x = "Day")
      else p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                          axis.line.x = element_blank())     # rows 1-3: drop the x axis
      p
    }
    b1 <- AGE[1:3]; b2 <- AGE[4:6]
    # assembled as Figure 3: y titles in one shared column, and the grid built column by column so
    # that columns 2 and 3 do not reserve the axis width of column 1
    wcol_k <- function(k) plot_grid(
      mk_wpanel(b1[k], "cases", rjw_top, FALSE, k == 1, FALSE, b1[k]),
      mk_wpanel(b1[k], "rj",    rjw_top, FALSE, k == 1, FALSE),
      mk_wpanel(b2[k], "cases", rjw_bot, FALSE, k == 1, FALSE, b2[k]),
      mk_wpanel(b2[k], "rj",    rjw_bot, FALSE, k == 1, TRUE),
      ncol = 1, align = "v", axis = "lr")
    wcols <- plot_grid(wcol_k(1), wcol_k(2), wcol_k(3), nrow = 1, rel_widths = c(1.16, 1, 1))
    wttl <- function(lab, above = 0, below = 0)
      plot_grid(NULL, y_title(lab), NULL, ncol = 1, rel_heights = c(max(above, 1e-6), 1, max(below, 1e-6)))
    wgttl <- plot_grid(wttl("Cases (per 100,000)", above = 0.09), wttl(LAB_RJT()),
                       wttl("Cases (per 100,000)", above = 0.09), wttl(LAB_RJT(), below = 0.22), ncol = 1)
    wmain <- plot_grid(wgttl, NULL, wcols, nrow = 1, rel_widths = c(0.030, 0.010, 1))
    wleg <- leg_row(list(
      list(kind = "point", label = "Observed (noisy)", colour = "grey65"),
      list(kind = "line",  label = "Truth",   colour = "black", lwd = 0.6, lty = "dotted"),
      list(kind = "line",  label = "Poisson", colour = unname(OKABE["vermillion"]), lwd = 0.7),
      list(kind = "line",  label = "NB30",    colour = PS_COL,  lwd = 0.7),
      list(kind = "rect",  label = "95% CrI", colour = "grey60")), size_pt = 10.2, width_in = 7.5)
    pc <- plot_grid(wleg, wmain, ncol = 1, rel_heights = c(0.04, 1))
    save_wilke(pc, file.path(outdir, "fig_weight_curves_r25"), width = 7.5, height = 8.5)
  }
}

# ── misspecification: true overdispersion outside the estimation grid, by estimated r ──
wmf <- file.path(here, "results", "ps_weight_misspec.csv")
if (file.exists(wmf)) {
  WLEV <- c("Poisson", "NB40", "NB30", "NB20", "NB10")     # increasing overdispersion from left to right (r decreasing)
  wm <- read_csv(wmf, show_col_types = FALSE) %>%
    mutate(weight = factor(weight, levels = WLEV),
           r_gen = factor(r_gen, levels = c(15, 25, 50)))
  wml <- wm %>% select(weight, r_gen, mape_Rj, cov_Rj, width_Rj, roughness) %>%
    pivot_longer(c(mape_Rj, cov_Rj, width_Rj, roughness),
                 names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c("mape_Rj", "cov_Rj", "width_Rj", "roughness"),
      labels = c("'ℛ'[italic(j)]~MAPE~('%')", "'95% coverage (%)'", "'95% CrI width'", "Roughness~(aux.)")))
  pm <- ggplot(wml, aes(weight, value, colour = r_gen, group = r_gen)) +
    geom_line(linewidth = 0.55) + geom_point(size = 1.6) +
    facet_wrap(~ metric, scales = "free_y", labeller = label_parsed, ncol = 2) +
    scale_colour_viridis_d(option = "C", end = 0.85,
                           name = expression(r[gen]~"(true overdispersion)")) +
    labs(x = NULL, y = NULL) + pub_theme(8) +
    theme(legend.position = "top", legend.title = element_text(size = 8))
  # separate data frame so the 95% reference line appears only in the coverage panel
  href <- data.frame(metric = factor("'95% coverage (%)'",
    levels = c("'ℛ'[italic(j)]~MAPE~('%')", "'95% coverage (%)'", "'95% CrI width'", "Roughness~(aux.)")), y = 95)
  pm <- pm + geom_hline(data = href, aes(yintercept = y), linetype = "dashed",
                        colour = REF_COL, linewidth = 0.35)
  save_wilke(pm, file.path(outdir, "fig_weight_misspec"), width = 6.2, height = 4.8)
}

# ═════════════════════════════════════════════════════════════════════
#  NHIS seasons, no ground truth (nhis_rt_<season>.csv)
#    cases = observed (black) and model fit (age colour with 95% CrI); R_j = smoother median with
#    CrI and the renewal estimate overlaid. Same 4x3 layout as the synthetic panel.
# ═════════════════════════════════════════════════════════════════════
POP <- c("0-5" = 1791788, "6-11" = 2707574, "12-17" = 2772098,
         "18-44" = 18023319, "45-64" = 16876969, "65+" = 9267290)   # resident registration population
# smoother = solid age colour, observed = black dotted, renewal = black dashed, 95% CrI = age band
make_nhis <- function(nf, rn_col = "black", rn_lty = "42", tag = "", do_tiff = TRUE) {
  # same type setting as Figure 3: no enlargement of the y title on this dense grid, so that both
  # axis titles and the age subtitles are at the common 9 pt nominal
  use_k(1.35, rmult = 1)                      # Figure 5
  season <- sub(".*nhis_rt_(\\d+)\\.csv", "\\1", basename(nf))
  d0 <- read_csv(nf, show_col_types = FALSE) %>%
    mutate(age_group = factor(age_group, levels = AGE),
           Rj_renewal = suppressWarnings(as.numeric(Rj_renewal)),
           sc = 1e5 / POP[as.character(age_group)],
           obs_pc = obs * sc, dRm_pc = dR_med * sc, dRlo_pc = dR_lo * sc, dRhi_pc = dR_hi * sc)
  cases_lim <- c(0, max(c(d0$obs_pc, d0$dRhi_pc)))            # includes the observed series and the dR credible interval
  rjlim <- function(ags) { s <- filter(d0, age_group %in% ags); c(min(s$Rj_lo), max(s$Rj_hi)) }
  rj_top <- rjlim(AGE[1:3]); rj_bot <- rjlim(AGE[4:6])
  # x axis: day is the coordinate and the KDCA week only the label, because week numbers roll over
  dwk <- distinct(d0, day, wk) %>% arrange(day)
  brk_i <- seq(1, nrow(dwk), by = 56)                        # a tick every 8 weeks (56 days)
  brk_days <- dwk$day[brk_i]; brk_labs <- dwk$wk[brk_i]

  mk_np <- function(ag, qty, rjl, show_yt, show_ytk, show_x, title = NULL) {
    col <- unname(age_cols[ag]); d <- filter(d0, age_group == ag)
    if (qty == "cases") {                                      # model fit: age colour solid with 95% CrI; observed: black dotted on top
      p <- ggplot(d, aes(day)) +
        geom_ribbon(aes(ymin = dRlo_pc, ymax = dRhi_pc), fill = col, alpha = 0.18) +
        geom_line(aes(y = dRm_pc), colour = col, linewidth = 0.5) +
        geom_line(aes(y = obs_pc), colour = "black", linewidth = 0.4, linetype = "dotted") +
        scale_y_continuous(labels = scales::label_comma()) +
        coord_cartesian(ylim = cases_lim) +
        labs(x = NULL, y = if (show_yt) "Cases (per 100,000)" else NULL)
    } else {                                                   # smoother median: age colour solid; renewal: separate colour and dash
      p <- ggplot(d, aes(day)) +
        geom_ribbon(aes(ymin = Rj_lo, ymax = Rj_hi), fill = col, alpha = 0.18) +
        geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.35) +
        geom_line(aes(y = Rj_renewal), colour = rn_col, linewidth = 0.45, linetype = rn_lty) +
        geom_line(aes(y = Rj_med), colour = col, linewidth = 0.5) +
        coord_cartesian(ylim = rjl) +
        labs(x = NULL, y = if (show_yt) LAB_RJT() else NULL)
    }
    p <- p + scale_x_continuous(breaks = brk_days, labels = brk_labs,
                                expand = expansion(mult = c(0.01, 0.02))) + pub_theme(8)
    # age subtitle at the axis-title size (see the note in the Figure 3 panel builder)
    if (!is.null(title)) p <- p + ggtitle(title) +
      theme(plot.title = element_text(size = 9 * FIG_K + FIG_D$title, face = "bold", hjust = 0.5,
                                      margin = margin(b = 1)))
    if (!show_ytk) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if (show_x) p <- p + labs(x = "Week")                     # only the bottom row carries the week axis
    else p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                        axis.line.x = element_blank())        # rows 1-3: drop the x axis line, ticks and labels
    p
  }
  b1 <- AGE[1:3]; b2 <- AGE[4:6]
  # assembled exactly as Figure 3: the four y titles in one shared column so that they line up,
  # and the grid built column by column so that columns 2 and 3 do not reserve the axis width of
  # column 1 (rel_widths compensates for column 1 carrying the tick labels)
  ncol_k <- function(k) plot_grid(
    mk_np(b1[k], "cases", rj_top, FALSE, k == 1, FALSE, b1[k]),
    mk_np(b1[k], "rj",    rj_top, FALSE, k == 1, FALSE),
    mk_np(b2[k], "cases", rj_bot, FALSE, k == 1, FALSE, b2[k]),
    mk_np(b2[k], "rj",    rj_bot, FALSE, k == 1, TRUE),
    ncol = 1, align = "v", axis = "lr")
  # 1.10, not the 1.16 of Figure 3: the tick labels here are shorter ("500" against "4,000")
  ncols <- plot_grid(ncol_k(1), ncol_k(2), ncol_k(3), nrow = 1, rel_widths = c(1.10, 1, 1))
  nttl <- function(lab, above = 0, below = 0)
    plot_grid(NULL, y_title(lab), NULL, ncol = 1, rel_heights = c(max(above, 1e-6), 1, max(below, 1e-6)))
  gttl <- plot_grid(nttl("Cases (per 100,000)", above = 0.09), nttl(LAB_RJT()),
                    nttl("Cases (per 100,000)", above = 0.09), nttl(LAB_RJT(), below = 0.22), ncol = 1)
  gm <- plot_grid(gttl, NULL, ncols, nrow = 1, rel_widths = c(0.030, 0.010, 1))
  # manual legend: observed (black dotted), smoother median (six colours), renewal, 95% CrI band
  leg <- leg_row(list(
    list(kind = "line",  label = "Observed",  colour = "black", lwd = 0.6, lty = "dotted"),
    list(kind = "line6", label = "PS median"),
    list(kind = "line",  label = "Renewal",   colour = rn_col,  lwd = 0.7, lty = rn_lty),
    list(kind = "rect6", label = "95% CrI")), size_pt = 11.5, width_in = 7.5)   # same as Figure 3
  fign <- plot_grid(leg, gm, ncol = 1, rel_heights = c(0.04, 1))
  stem <- if (nchar(tag)) sprintf("fig_nhis_%s_%s", season, tag) else sprintf("fig_nhis_%s", season)
  save_wilke(fign, file.path(outdir, stem), width = 7.5, height = 8.6, tiff = do_tiff)
}
for (nf in Sys.glob(file.path(here, "results", "nhis_rt_*.csv"))) make_nhis(nf)   # black dashed·TIFF

# ── 2018-19 overlay: (A) observed cases for all age groups, (B) smoother R_j medians ──
nf18 <- file.path(here, "results", "nhis_rt_2018.csv")
if (file.exists(nf18)) {
  d18 <- read_csv(nf18, show_col_types = FALSE) %>%
    mutate(age_group = factor(age_group, levels = AGE),
           sc = 1e5 / POP[as.character(age_group)], obs_pc = obs * sc)
  dwk18 <- distinct(d18, day, wk) %>% arrange(day)                 # x axis: day as coordinate, KDCA week as label
  bi <- seq(1, nrow(dwk18), by = 56); bd <- dwk18$day[bi]; bl <- dwk18$wk[bi]
  use_k(1.42, legend = 1.75)                  # Figure 4
  qA <- ggplot(d18, aes(day, obs_pc, colour = age_group, linewidth = age_group)) +
    geom_line() + age_scale +
    scale_x_continuous(breaks = bd, labels = bl, expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = scales::label_comma()) +
    labs(x = NULL, y = NULL) +
    pub_theme() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  qB <- ggplot(d18, aes(day, Rj_med, colour = age_group, linewidth = age_group)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.4) +
    geom_line() + age_scale +
    scale_x_continuous(breaks = bd, labels = bl, expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(0, min(3, max(d18$Rj_med, na.rm = TRUE) * 1.05))) +
    labs(x = "Week", y = NULL) +
    pub_theme() + theme(legend.position = "none")
  # same layout as Figure 2 (identical canvas and panel structure): the two y titles are drawn in
  # one shared column so that they line up, with the spacer column between titles and tick labels
  fig18_body <- plot_grid(qA, qB, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1),
                     labels = c("A", "B"), label_fontfamily = FONT, label_size = 12 * FIG_K, label_fontface = "bold")
  fig18_titles <- plot_grid(y_title("Cases (per 100,000)"),
                            plot_grid(y_title(LAB_RJT()), NULL, ncol = 1, rel_heights = c(1, 0.16)),
                            ncol = 1, rel_heights = c(1, 1))
  fig18 <- plot_grid(fig18_titles, NULL, fig18_body, nrow = 1, rel_widths = c(0.045, 0.022, 1))
  save_wilke(fig18, file.path(outdir, "fig_nhis_overlay_2018"), width = 6.5, height = 6.9)   # each sub-panel about 2.0:1
}

# ═════════════════════════════════════════════════════════════════════
#  Contact-matrix robustness: survey matrix vs Prem 2017 / Prem 2021 (Korea)
#    The Prem matrices are aggregated to the six bands and transposed to the convention used here
#    (code/make_prem_contact.py). (A) the matrices themselves, (B) the resulting R_j medians.
# ═════════════════════════════════════════════════════════════════════
cm_our <- matrix(c(1.76,0.21,0.03,0.20,0.05,0.05, 0.33,3.75,0.30,0.31,0.14,0.13,
                   0.05,0.31,3.65,0.22,0.31,0.12, 2.07,2.11,1.38,1.78,1.26,0.84,
                   0.47,0.88,1.87,1.18,1.97,1.47, 0.30,0.47,0.42,0.45,0.83,2.50),
                 nrow = 6, byrow = TRUE)   # C[i,j]: row i = infectee, column j = infector
cm_f17 <- file.path(ROOT, "data", "contact_prem2017.csv")
cm_f21 <- file.path(ROOT, "data", "contact_prem2021.csv")
if (file.exists(cm_f17) && file.exists(cm_f21)) {
  cm_p17 <- as.matrix(read.csv(cm_f17, header = FALSE))
  cm_p21 <- as.matrix(read.csv(cm_f21, header = FALSE))
  dimnames(cm_our) <- dimnames(cm_p17) <- dimnames(cm_p21) <- list(AGE, AGE)
  # (A) heat map of raw contact numbers on a shared scale; x = participant (ego), y = contact
  cm_tidy <- function(M, lab) {
    data.frame(matrix_lab = lab,
               contact     = factor(rep(AGE, times = 6), levels = AGE),   # row i = contact
               participant = factor(rep(AGE, each = 6),  levels = AGE),   # column j = participant (ego)
               raw = as.vector(M))
  }
  cmdf <- rbind(cm_tidy(cm_our, "Son 2025"), cm_tidy(cm_p17, "Prem 2017"),
                cm_tidy(cm_p21, "Prem 2021")) %>%
    mutate(matrix_lab = factor(matrix_lab, levels = c("Son 2025", "Prem 2017", "Prem 2021")))
  use_k(1.32)                                 # Figure S8
  ph <- ggplot(cmdf, aes(participant, contact, fill = raw)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f", raw)), size = 2.3, family = FONT) +
    facet_wrap(~ matrix_lab, ncol = 3) +
    scale_fill_gradient(low = "white", high = unname(OKABE["vermillion"]),
                        name = "Contacts/day", limits = c(0, max(cmdf$raw))) +
    scale_y_discrete(limits = AGE) +                    # y ascending: 0-5 at the bottom, 65+ at the top
    labs(x = "Participant", y = "Contact") +
    pub_theme(8) + theme(legend.position = "right", panel.grid = element_blank(),
                         axis.text.x = element_text(angle = 45, hjust = 1))
  save_wilke(ph, file.path(outdir, "fig_contact_matrices"), width = 8.2, height = 3.2)

  # (B) R_j(t) medians under the three matrices, one facet per age group
  cm_cols <- c("Son 2025" = "grey25", "Prem 2017" = unname(OKABE["vermillion"]), "Prem 2021" = unname(OKABE["blue"]))
  make_cm_compare <- function(sy) {
    use_k(1.22)                               # Figure S9
    f_our <- file.path(here, "results", sprintf("nhis_rt_%d.csv", sy))
    f17 <- file.path(here, "results", sprintf("nhis_rj_prem2017_%d.csv", sy))
    f21 <- file.path(here, "results", sprintf("nhis_rj_prem2021_%d.csv", sy))
    if (!all(file.exists(c(f_our, f17, f21)))) return(invisible(NULL))
    dO <- read_csv(f_our, show_col_types = FALSE) %>% transmute(day, wk, age_group, Rj_med, cm = "Son 2025")
    d17 <- read_csv(f17, show_col_types = FALSE) %>% transmute(day, wk, age_group, Rj_med, cm = "Prem 2017")
    d21 <- read_csv(f21, show_col_types = FALSE) %>% transmute(day, wk, age_group, Rj_med, cm = "Prem 2021")
    d <- bind_rows(dO, d17, d21) %>%
      mutate(age_group = factor(age_group, levels = AGE),
             cm = factor(cm, levels = c("Son 2025", "Prem 2017", "Prem 2021")))
    dwk <- distinct(dO, day, wk) %>% arrange(day)
    bi <- seq(1, nrow(dwk), by = 56); bd <- dwk$day[bi]; bl <- dwk$wk[bi]
    p <- ggplot(d, aes(day, Rj_med, colour = cm)) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.35) +
      geom_line(linewidth = 0.5) +
      facet_wrap(~ age_group, ncol = 3) +
      scale_colour_manual(values = cm_cols, name = NULL) +
      scale_x_continuous(breaks = bd, labels = bl) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +      # y = 0 is the lower bound of the axis
      coord_cartesian(ylim = c(0, min(3, max(d$Rj_med, na.rm = TRUE) * 1.05))) +
      labs(x = "Week", y = LAB_RJT()) +
      pub_theme(8) + theme(legend.position = "top")
    save_wilke(p, file.path(outdir, sprintf("fig_nhis_cm_%d", sy)), width = 7.2, height = 4.8)
  }
  for (sy in c(2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2022)) make_cm_compare(sy)
}

# ═════════════════════════════════════════════════════════════════════
#  2018-19 with the healthcare-seeking and pre-existing-immunity corrections combined
#    nhis.jl RT_IMHS=1 → nhis_overlay_imhs_2018.csv(cond=baseline/corrected, S_med·cases·Rj_med).
#    (A) susceptible fraction S(t)/N, (B) cases, (C) R_j(t) median. Baseline solid, corrected dashed.
# ═════════════════════════════════════════════════════════════════════
im_f <- file.path(here, "results", "nhis_overlay_imhs_2018.csv")
if (file.exists(im_f)) {
  d <- read_csv(im_f, show_col_types = FALSE) %>%
    mutate(age_group = factor(age_group, levels = AGE),
           cond = factor(ifelse(cond == "baseline", "Baseline", "Immunity + HC-seeking"),
                         levels = c("Baseline", "Immunity + HC-seeking")))
  scv <- function(ag) 1e5 / POP[as.character(ag)]
  dwk <- distinct(d, day, wk) %>% arrange(day); bi <- seq(1, nrow(dwk), by = 56); bd <- dwk$day[bi]; bl <- dwk$wk[bi]
  # "22" rather than "dashed": at this line width the default dash shows only one segment per curve
  im_lty <- c("Baseline" = "solid", "Immunity + HC-seeking" = "22")
  age_c <- scale_colour_manual(values = age_cols, breaks = AGE, guide = guide_legend(nrow = 1, order = 1))
  lty_g <- guide_legend(nrow = 1, order = 2, override.aes = list(colour = "grey30", linewidth = 0.7))
  xsc <- scale_x_continuous(breaks = bd, labels = bl)
  use_k(1.32)                                 # Figure S10
  gA <- ggplot(d, aes(day, S_med / POP[as.character(age_group)], colour = age_group, linetype = cond)) +
    geom_line(linewidth = 0.55) + age_c + scale_linetype_manual(values = im_lty, guide = lty_g) + xsc +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::label_percent(accuracy = 1)) +
    # short title: the rotated long form is taller than the panel once the type is scaled for print
    labs(x = NULL, y = NULL) +
    pub_theme() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                        # stack the age and condition legends instead of placing them side by side
                        legend.box = "vertical", legend.spacing.y = unit(1, "pt"),
                        legend.key.width = unit(26, "pt"))   # room for more than one dash in the key
  gB <- ggplot(d, aes(day, cases * scv(age_group), colour = age_group, linetype = cond)) +
    geom_line(linewidth = 0.55) + age_c + scale_linetype_manual(values = im_lty, guide = "none") + xsc +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::label_comma()) +
    labs(x = NULL, y = NULL) +
    pub_theme() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
  gC <- ggplot(d, aes(day, Rj_med, colour = age_group, linetype = cond)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.4) +
    geom_line(linewidth = 0.55) + age_c + scale_linetype_manual(values = im_lty, guide = "none") + xsc +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(ylim = c(0, min(3, max(d$Rj_med, na.rm = TRUE) * 1.05))) +
    labs(x = "Week", y = NULL) + pub_theme() + theme(legend.position = "none")
  # the three y titles are drawn in one shared column so that they line up; the spacers offset the
  # legend above panel A and the x axis below panel C, so each title stays centred on its own panel
  im_body <- plot_grid(gA, gB, gC, ncol = 1, align = "v", axis = "lr", labels = c("A", "B", "C"),
                     label_fontfamily = FONT, label_size = 12 * FIG_K, label_fontface = "bold")
  im_ttl <- plot_grid(plot_grid(NULL, y_title("Susceptible (%)"), ncol = 1, rel_heights = c(0.13, 1)),
                      y_title("Cases (per 100,000)"),
                      plot_grid(y_title(LAB_RJT()), NULL, ncol = 1, rel_heights = c(1, 0.13)),
                      ncol = 1)
  figim <- plot_grid(im_ttl, NULL, im_body, nrow = 1, rel_widths = c(0.045, 0.022, 1))
  save_wilke(figim, file.path(outdir, "fig_nhis_immun_hseek_overlay_2018"), width = 6.5, height = 9.6)
}

# ═════════════════════════════════════════════════════════════════════
#  2018-19: 7-day moving average (baseline) vs raw daily counts
#    (1) overlay of all age groups, (2) 3x2 per-age medians with 95% CrI; raw de-emphasised in grey
#    baseline=nhis_rt_2018.csv(obs·Rj_med±CrI), raw=nhis_rj_raw_2018.csv(cases·Rj_med±CrI; causal seed).
# ═════════════════════════════════════════════════════════════════════
rf_our <- file.path(here, "results", "nhis_rt_2018.csv")
rf_raw <- file.path(here, "results", "nhis_rj_raw_2018.csv")
if (all(file.exists(c(rf_our, rf_raw)))) {
  dO <- read_csv(rf_our, show_col_types = FALSE) %>% mutate(age_group = factor(age_group, levels = AGE))
  dR <- read_csv(rf_raw, show_col_types = FALSE) %>% mutate(age_group = factor(age_group, levels = AGE))
  scv <- function(ag) 1e5 / POP[as.character(ag)]
  MA <- "7-day MA"; RAW <- "Raw daily"
  RAWLINE <- "grey55"; RAWFILL <- "grey72"
  casesO <- dO %>% transmute(day, age_group, val = obs   * scv(age_group))
  casesR <- dR %>% transmute(day, age_group, val = cases * scv(age_group))
  rjO <- dO %>% transmute(day, age_group, val = Rj_med)
  rjR <- dR %>% transmute(day, age_group, val = Rj_med)
  dwk <- distinct(dO, day, wk) %>% arrange(day); bi <- seq(1, nrow(dwk), by = 56); bd <- dwk$day[bi]; bl <- dwk$wk[bi]
  age_c <- scale_colour_manual(values = age_cols, breaks = AGE, guide = guide_legend(nrow = 1, order = 1))
  age_f <- scale_fill_manual(values = age_cols, breaks = AGE, guide = "none")
  xsc <- scale_x_continuous(breaks = bd, labels = bl)
  ysc_c <- scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::label_comma())
  hl <- geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.4)
  no_x <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  rjmax <- max(rjO$val, rjR$val, na.rm = TRUE)

  # legend for the two conditions: grey raw vs thicker baseline in age colour
  cg_lw <- function(o) guide_legend(nrow = 1, order = o, override.aes = list(colour = c(RAWLINE, "grey25")))
  lw_g <- c("Raw daily" = 0.35, "7-day MA" = 0.85)

  # ── (1) overlay: cases and R_j medians for all age groups ──
  use_k(1.24)                                 # Figure S6
  ovA <- ggplot(mapping = aes(day, val)) +
    geom_line(data = casesR, aes(group = age_group, linewidth = RAW), colour = RAWLINE, alpha = 0.85) +
    geom_line(data = casesO, aes(colour = age_group, linewidth = MA)) + age_c +
    scale_linewidth_manual(values = lw_g, name = NULL, guide = cg_lw(2)) +
    xsc + ysc_c + labs(x = NULL, y = NULL) + pub_theme() + no_x
  ovB <- ggplot(mapping = aes(day, val)) + hl +
    geom_line(data = rjR, aes(group = age_group), colour = RAWLINE, linewidth = 0.35, alpha = 0.85) +
    geom_line(data = rjO, aes(colour = age_group), linewidth = 0.85) + age_c +
    xsc + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + coord_cartesian(ylim = c(0, min(3, rjmax * 1.05))) +
    labs(x = "Week", y = NULL) + pub_theme() + theme(legend.position = "none")
  # y titles in one shared column, as in Figure 2 (same canvas and panel structure)
  ov_body <- plot_grid(ovA, ovB, ncol = 1, align = "v", axis = "lr", labels = c("A", "B"),
                       label_fontfamily = FONT, label_size = 12 * FIG_K, label_fontface = "bold")
  ov_ttl <- plot_grid(plot_grid(NULL, y_title("Cases (per 100,000)"), ncol = 1, rel_heights = c(0.08, 1)),
                      plot_grid(y_title(LAB_RJT()), NULL, ncol = 1, rel_heights = c(1, 0.16)),
                      ncol = 1)
  save_wilke(plot_grid(ov_ttl, NULL, ov_body, nrow = 1, rel_widths = c(0.045, 0.022, 1)),
             file.path(outdir, "fig_nhis_raw_overlay_2018"), width = 6.5, height = 6.9)

  # ── (2) 3x2 per age group: R_j median with 95% CrI. The two intervals are nearly the same width,
  #    so the baseline is filled and the raw interval is drawn as dotted boundaries. ──
  cri_ylim <- coord_cartesian(ylim = c(0, min(3, max(dO$Rj_hi, dR$Rj_hi, na.rm = TRUE) * 1.05)))
  use_k(1.47)                                 # Figure S7
  cri <- ggplot() +
    geom_ribbon(data = dO, aes(day, ymin = Rj_lo, ymax = Rj_hi, fill = age_group), alpha = 0.30) + age_f + hl +
    geom_line(data = dR, aes(day, Rj_lo), colour = RAWLINE, linewidth = 0.28, linetype = "22") +
    geom_line(data = dR, aes(day, Rj_hi), colour = RAWLINE, linewidth = 0.28, linetype = "22") +
    geom_line(data = dR, aes(day, Rj_med, linewidth = RAW), colour = RAWLINE) +
    geom_line(data = dO, aes(day, Rj_med, colour = age_group, linewidth = MA)) + age_c +
    scale_linewidth_manual(values = lw_g, name = NULL, guide = cg_lw(2)) +
    facet_wrap(~ age_group, ncol = 3) + xsc +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + cri_ylim +   # y = 0 is the lower bound of the axis
    labs(x = "Week", y = LAB_RJT()) + pub_theme(9) + theme(legend.position = "top")
  save_wilke(cri, file.path(outdir, "fig_nhis_raw_cri_2018"), width = 8.2, height = 5.4)

  # ── method and weight comparisons for 2018-19: baseline (smoother, NB30) vs a variant ──
  #    baseline = age-coloured median with a filled CrI; variant = grey median with dotted CrI bounds
  #    cri = FALSE draws medians only, which makes repeated crossings of the two medians visible
  make_cri_cmp <- function(varfile, varlab, fname, cri = TRUE) {
    use_k(1.47)                               # Figure S5
    f <- file.path(here, "results", varfile)
    if (!file.exists(f)) return(invisible(NULL))
    dv <- read_csv(f, show_col_types = FALSE) %>% mutate(age_group = factor(age_group, levels = AGE))
    lwv <- setNames(if (cri) c(0.35, 0.85) else c(0.55, 0.85), c(varlab, "Baseline"))
    cg <- guide_legend(nrow = 1, order = 2, override.aes = list(colour = c(RAWLINE, "grey25")))
    ytop <- if (cri) max(dO$Rj_hi, dv$Rj_hi, na.rm = TRUE) else max(dO$Rj_med, dv$Rj_med, na.rm = TRUE)
    p <- ggplot()
    if (cri) p <- p +
      geom_ribbon(data = dO, aes(day, ymin = Rj_lo, ymax = Rj_hi, fill = age_group), alpha = 0.30) + age_f +
      geom_line(data = dv, aes(day, Rj_lo), colour = RAWLINE, linewidth = 0.28, linetype = "22") +
      geom_line(data = dv, aes(day, Rj_hi), colour = RAWLINE, linewidth = 0.28, linetype = "22")
    p <- p + hl +
      geom_line(data = dv, aes(day, Rj_med, linewidth = varlab), colour = RAWLINE) +
      geom_line(data = dO, aes(day, Rj_med, colour = age_group, linewidth = "Baseline")) + age_c +
      scale_linewidth_manual(values = lwv, name = NULL, breaks = c(varlab, "Baseline"), guide = cg) +
      facet_wrap(~ age_group, ncol = 3) + xsc +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +      # y = 0 is the lower bound of the axis
      coord_cartesian(ylim = c(0, min(3, ytop * 1.05))) +
      labs(x = "Week", y = LAB_RJT()) + pub_theme(9) +
      # the two legends are stacked, not placed side by side: at print size the age keys and the
      # long condition labels do not fit on one line
      theme(legend.position = "top", legend.box = "vertical", legend.spacing.y = unit(1, "pt"))
    save_wilke(p, file.path(outdir, fname), width = 8.2, height = 5.4)
  }
  make_cri_cmp("nhis_rj_raw_pois_2018.csv", "PS + Poisson (raw)",      "fig_nhis_raw_ps_poisson_2018")
  make_cri_cmp("nhis_rj_raw_pf_2018.csv",   "PF + NB(30) (raw)",       "fig_nhis_raw_pf_nb_2018")
  make_cri_cmp("nhis_rj_ma_pois_2018.csv",  "PS + Poisson (7-day MA)", "fig_nhis_ma_ps_poisson_2018",
               cri = FALSE)                                             # medians only
}

# ═════════════════════════════════════════════════════════════════════
#  fig_pfps_uniform: uniform SEIR step/sin pre-defined 𝓡_t, PF vs PS (median+95% CrI vs truth). σ=γ=1/1.8.
#    2x2 = {cases, R_t} x {step, sinusoidal}. Smoother blue, filter vermillion, truth dotted.
# ═════════════════════════════════════════════════════════════════════
uf <- file.path(here, "results", "pfps_uniform_curves.csv")
if (file.exists(uf)) {
  u <- read_csv(uf, show_col_types = FALSE) %>% mutate(scenario = factor(scenario, levels = c("step", "sin")))
  PScol <- unname(OKABE["blue"]); PFcol <- unname(OKABE["vermillion"])
  use_k(1.22)                                 # Figure S2
  mk_u <- function(scen, qty, show_x, show_ytitle, title = NULL) {
    d <- filter(u, scenario == scen, quantity == qty) %>% mutate(series = factor(series, levels = c("true", "PS", "PF")))
    band <- filter(d, series %in% c("PS", "PF"))
    p <- ggplot()
    if (qty == "Rt") p <- p + geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COL, linewidth = 0.35)
    p <- p +
      geom_ribbon(data = band, aes(day, ymin = lo, ymax = hi, fill = series), alpha = 0.18) +
      geom_line(data = d, aes(day, mid, colour = series, linetype = series), linewidth = 0.5) +
      scale_fill_manual(values = c(PS = PScol, PF = PFcol), guide = "none") +
      scale_colour_manual(values = c(true = "black", PS = PScol, PF = PFcol), name = NULL,
                          labels = c(true = "Truth", PS = "PS", PF = "PF"),
                          guide = guide_legend(override.aes = list(linetype = c("dotted", "solid", "solid")))) +
      scale_linetype_manual(values = c(true = "dotted", PS = "solid", PF = "solid"), guide = "none") +
      labs(x = if (show_x) "Day" else NULL,
           y = if (show_ytitle) (if (qty == "obs") "Cases" else LAB_RT()) else NULL) +
      pub_theme(9)
    if (qty == "obs") p <- p + scale_y_continuous(labels = scales::label_comma())
    # panel subtitle at the size of the y axis titles (a fixed 10 pt did not follow FIG_K and printed
    # smaller than every other label); this keeps the subtitle >= axis title > tick hierarchy
    if (!is.null(title)) p <- p + ggtitle(title) +
      theme(plot.title = element_text(size = 9 * FIG_K * r_mult(), face = "bold", hjust = 0.5))
    # rows without an x axis also drop the axis LINE, otherwise a stray axis appears under y = 0
    if (!show_x) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                               axis.line.x = element_blank())
    p
  }
  leg <- get_legend(mk_u("step", "Rt", TRUE, TRUE) +
                    theme(legend.position = "top", legend.direction = "horizontal"))
  p11 <- mk_u("step", "obs", FALSE, FALSE, "Step")       + theme(legend.position = "none")
  p12 <- mk_u("sin",  "obs", FALSE, FALSE, "Sinusoidal") + theme(legend.position = "none")
  p21 <- mk_u("step", "Rt",  TRUE,  FALSE) + theme(legend.position = "none")
  p22 <- mk_u("sin",  "Rt",  TRUE,  FALSE) + theme(legend.position = "none")
  body <- plot_grid(p11, p12, p21, p22, ncol = 2, align = "hv", axis = "tblr", labels = c("A", "B", "C", "D"),
                    label_fontfamily = FONT, label_size = 12 * FIG_K, label_fontface = "bold")
  # The two y titles are drawn in one shared column so that they line up (row 1 has the wider tick
  # labels). align = "hv" gives both rows the same padding — the panel title above and the x axis
  # below are reserved in each — so the panels sit the same distance above their cell centre and
  # both titles take the same spacer, rather than one per row.
  ttl_u <- function(lab) plot_grid(y_title(lab), NULL, ncol = 1, rel_heights = c(1, 0.09))
  grid <- plot_grid(plot_grid(ttl_u("Cases"), ttl_u(LAB_RT()), ncol = 1),
                    NULL, body, nrow = 1, rel_widths = c(0.032, 0.020, 1))
  save_wilke(plot_grid(leg, grid, ncol = 1, rel_heights = c(0.06, 1)),
             file.path(outdir, "fig_pfps_uniform"), width = 7.2, height = 5.6)
}

# ═════════════════════════════════════════════════════════════════════
#  fig_nhis_heatmap — season x week heat map of R_j(t) over the eleven NHIS seasons
#    rows = seasons (09-10 at the bottom), columns = time within the season (week 36 to week 35),
#    facets = age group, fill = smoother median.
#    x is the day index within the season and only the tick LABEL is the KDCA week: week numbers
#      cannot serve as a coordinate when a season contains a week 53.
#    Week 53 (only 2011-12, 2016-17 and 2022-23) is merged into week 52 by averaging matching days
#      of the week, so every season is 52 weeks and 364 days long and one set of week labels fits
#      all of them. Week 53 always follows week 52 and both are 7 days long.
#    Colours: grey at R_j = 1, blue below and red above, changing sharply near 1 because the colour
#      stops are compressed around 1 rather than spread uniformly.
hm_files <- Sys.glob(file.path(here, "results", "nhis_rt_*.csv"))
if (length(hm_files) > 0) {
  hm <- bind_rows(lapply(hm_files, function(f) {
    sy <- as.integer(sub(".*nhis_rt_(\\d+)\\.csv", "\\1", basename(f)))
    read_csv(f, show_col_types = FALSE) %>%
      transmute(start_year = sy, day, wk, age_group, Rj = Rj_med)
  }))
  # merge week 53: within one season and age group, average weeks 52 and 53 by day of week, drop 53
  merge_w53 <- function(d) {
    d <- d[order(d$day), ]
    i53 <- which(d$wk == 53)
    if (length(i53) > 0) {
      i52 <- which(d$wk == 52)
      stopifnot(length(i52) == length(i53), all(diff(range(c(i52, i53))) == length(i52) + length(i53) - 1))
      d$Rj[i52] <- (d$Rj[i52] + d$Rj[i53]) / 2                 # average matching days of the week
      d <- d[-i53, ]
    }
    d$day <- seq_len(nrow(d))
    d
  }
  hm <- hm %>% group_by(start_year, age_group) %>% group_modify(~ merge_w53(.x)) %>% ungroup()
  stopifnot(all(table(hm$start_year, hm$age_group) == 364))    # every season must now be 364 days long
  sy_lev <- sort(unique(hm$start_year))                       # 2009 to 2022; the 2019-20 to 2021-22 seasons are absent
  sy_lab <- sprintf("%02d–%02d", sy_lev %% 100, (sy_lev + 1) %% 100)
  hm <- hm %>% mutate(
    season    = factor(sprintf("%02d–%02d", start_year %% 100, (start_year + 1) %% 100),
                       levels = sy_lab),                      # the first level is the bottom row (09-10)
    age_group = factor(age_group, levels = AGE,
                       labels = c("0–5", "6–11", "12–17",
                                  "18–44", "45–64", "65+")))
  # x ticks every 8 weeks (56 days), labelled with the KDCA week of that day. Four-week ticks
  # collide once the type is scaled up for print.
  ref <- filter(hm, start_year == 2018) %>% distinct(day, wk)
  brk <- seq(1, max(hm$day), by = 56)
  lab <- sapply(brk, function(d) { w <- ref$wk[match(d, ref$day)]; if (is.na(w)) "" else as.character(w) })
  # legend range 0-3, which covers the whole data range so nothing is clipped;
  #   grey at 1, changing sharply within about 10%, flattening towards 0 and 3.
  #   The grey must be a BAND (0.97-1.03): a single stop has zero width and would be invisible.
  HM_LIMS <- c(0, 3)
  HM_COLS <- c("#2166AC", "#4393C3", "#9ECAE1", "#BDBDBD", "#BDBDBD", "#FCAE91", "#D6604D", "#B2182B")
  HM_VALS <- scales::rescale(c(HM_LIMS[1], 0.70, 0.90, 0.97, 1.03, 1.10, 1.35, HM_LIMS[2]), from = HM_LIMS)
  cat(sprintf("fig_nhis_heatmap: legend [%g, %g], R_j data range [%.3f, %.3f]\n",
              HM_LIMS[1], HM_LIMS[2], min(hm$Rj), max(hm$Rj)))
  # 11 rows in a short panel: label every other season, and the last three in full (the 2019-20 to
  # 2021-22 gap means the top of the axis is where the labels matter most). Every row keeps its tick.
  sy_keep <- unique(c(sy_lab[seq(1, length(sy_lab), by = 2)], tail(sy_lab, 3)))
  use_k(1.14)                                 # Figure 7
  p_hm <- ggplot(hm, aes(day, season, fill = Rj)) +
    # geom_tile, not geom_raster: cairo embeds a raster layer as a low-resolution bitmap
    # (one pixel per cell), which blurs the cell edges when the EPS is scaled up.
    geom_tile(colour = NA) +
    facet_wrap(~ age_group, ncol = 3) +
    scale_x_continuous(breaks = brk, labels = lab, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0), labels = function(s) ifelse(s %in% sy_keep, s, "")) +
    scale_fill_gradientn(colours = HM_COLS, values = HM_VALS, limits = HM_LIMS,
                         oob = scales::squish, breaks = seq(0, 3, by = 0.5),
                         name = LAB_RJT_H(),        # horizontal legend title
                         guide = guide_colourbar(barwidth = 0.6, barheight = 11,
                                                 frame.colour = "grey40", frame.linewidth = 0.3,
                                                 ticks.colour = "grey40")) +
    labs(x = "Week", y = "Influenza season") + pub_theme(9) +
    theme(panel.grid = element_blank(),
          # 11 season rows in a short panel: the tick labels are set below the shared rule
          axis.text = element_text(size = 7.8),
          panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.3),
          axis.ticks = element_line(colour = "grey30", linewidth = 0.3),
          legend.position = "right", legend.justification = "center",
          legend.title = element_text(size = 9 * FIG_K * R_LAB_MULT), legend.key.height = unit(11, "pt"),
          panel.spacing = unit(6, "pt"))
  save_wilke(p_hm, file.path(outdir, "fig_nhis_heatmap"), width = 7.5, height = 4.4)
}

# ═════════════════════════════════════════════════════════════════════
#  fig_nhis_overlay_10seasons — the ten seasons other than 2018-19, in one figure
#    Each season occupies one row (all-age cases and R_j medians, as in the 2018-19 overlay);
#    two blocks of five rows sit side by side, giving 4 columns x 5 rows = 20 sub-panels.
#    Both axes are shared across seasons so epidemic sizes can be compared directly; the price is
#      that the two weakest seasons (2010-11, 2012-13) are compressed.
#    x axis and the week-53 merge follow fig_nhis_heatmap.
#    Sized to fit one A4 page with room for the caption: 7.5 x 8.4 in.
# ═════════════════════════════════════════════════════════════════════
S10 <- c(2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2022)
if (all(file.exists(file.path(here, "results", sprintf("nhis_rt_%d.csv", S10))))) {
  mw53 <- function(d) {                        # average week 52 and 53 by day of week into week 52, then drop week 53
    d <- d[order(d$day), ]; i53 <- which(d$wk == 53)
    if (length(i53) > 0) {
      i52 <- which(d$wk == 52)
      d$obs_pc[i52] <- (d$obs_pc[i52] + d$obs_pc[i53]) / 2
      d$Rj_med[i52] <- (d$Rj_med[i52] + d$Rj_med[i53]) / 2
      d <- d[-i53, ]
    }
    d$day <- seq_len(nrow(d)); d
  }
  s10 <- bind_rows(lapply(S10, function(sy)
      read_csv(file.path(here, "results", sprintf("nhis_rt_%d.csv", sy)), show_col_types = FALSE) %>%
        mutate(start_year = sy, obs_pc = obs * 1e5 / POP[age_group]))) %>%
    group_by(start_year, age_group) %>% group_modify(~ mw53(.x)) %>% ungroup() %>%
    mutate(age_group = factor(age_group, levels = AGE),
           season = factor(sprintf("%02d–%02d", start_year %% 100, (start_year + 1) %% 100),
                           levels = sprintf("%02d–%02d", S10 %% 100, (S10 + 1) %% 100)))
  s10_lw <- c("0-5"=0.5, "6-11"=0.5, "12-17"=0.5, "18-44"=0.34, "45-64"=0.34, "65+"=0.5)  # thinner lines for the dense layout
  s10_scale <- list(scale_colour_manual(values = age_cols, breaks = AGE),
                    scale_linewidth_manual(values = s10_lw, guide = "none"),
                    guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.1))))
  s10_ref <- distinct(filter(s10, start_year == 2009), day, wk) %>% arrange(day)
  s10_bd <- s10_ref$day[seq(1, nrow(s10_ref), by = 112)]                      # a tick every 16 weeks in this dense layout
  s10_bl <- s10_ref$wk[seq(1, nrow(s10_ref), by = 112)]
  s10_case_lim <- c(0, max(s10$obs_pc) * 1.05); s10_rj_lim <- c(0, max(s10$Rj_med) * 1.05)
  # tick labels are set smaller than the shared rule here: 20 sub-panels on one page leave no room
  s10_theme <- function() pub_theme(7) +
    theme(legend.position = "none", plot.margin = margin(3, 4, 3, 1),   # top/bottom: air between rows
          axis.text = element_text(size = 9))
  use_k(1.39)                                 # Figure 6
  s10_pan <- function(sea, kind, show_x) {
    d <- filter(s10, season == sea)
    p <- ggplot(d, aes(day, if (kind == "cases") obs_pc else Rj_med,
                       colour = age_group, linewidth = age_group))
    if (kind == "rj") p <- p + geom_hline(yintercept = 1, linetype = "dashed",
                                          colour = REF_COL, linewidth = 0.3)
    p <- p + geom_line() + s10_scale +
      scale_x_continuous(breaks = s10_bd, labels = s10_bl, expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.06)),
                         labels = scales::label_comma()) +
      coord_cartesian(ylim = if (kind == "cases") s10_case_lim else s10_rj_lim) +
      labs(x = NULL, y = NULL) +          # the quantity is named once per column, in the header
      s10_theme()
    if (!show_x) p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p
  }
  # A rotated y title repeated on every row would be taller than the row itself, so each quantity
  # is named once, horizontally, above its column.
  # Drawn on a canvas of its own rather than as a plot title, so that the vertical position is set
  # directly: the two headers are at different point sizes, and a title is placed from the TOP of
  # its cell, which left the taller one sitting a baseline lower. Centring both boxes puts them on
  # one line (measured: the two baselines then differ by 1 px at 150 dpi); dy is left for tuning.
  s10_hdr <- function(lab, mult = 1, dy = 0) ggdraw() +
    draw_label(lab, fontfamily = FONT, size = 9 * FIG_K * mult, y = 0.55 + dy, vjust = 0.5)
  # column widths of one block, shared by the header and the data rows: the two spacers are the
  # gap between the season label and the tick labels, and the gap between the two panels
  RW <- c(0.10, 0.05, 1, 0.055, 1)
  s10_block <- function(seas) {                 # header + five rows of [season label | cases | R_j]
    hdr <- plot_grid(NULL, NULL, s10_hdr("Cases (per 100,000)"),
                     NULL, s10_hdr(LAB_RJT_H(), R_LAB_MULT),
                     nrow = 1, rel_widths = RW)
    rows <- lapply(seq_along(seas), function(k) {
      last <- k == length(seas)
      plot_grid(ggdraw() + draw_label(seas[k], fontfamily = FONT, fontface = "bold",
                                      size = 10 * FIG_K, angle = 90),
                NULL, s10_pan(seas[k], "cases", last), NULL, s10_pan(seas[k], "rj", last),
                nrow = 1, rel_widths = RW, align = "h", axis = "tb")   # wider label column
    })
    hs <- rep(1, length(seas)); hs[length(seas)] <- 1.16   # extra height for the x axis labels on the last row
    plot_grid(hdr, plot_grid(plotlist = rows, ncol = 1, rel_heights = hs),
              ncol = 1, rel_heights = c(0.055, 1))
  }
  s10_leg <- get_legend(ggplot(s10, aes(day, obs_pc, colour = age_group, linewidth = age_group)) +
                          geom_line() + s10_scale + pub_theme(8) +
                          theme(legend.position = "top", legend.justification = "center"))
  SL10 <- levels(s10$season)
  fig10 <- plot_grid(s10_leg,
                     plot_grid(s10_block(SL10[1:5]), NULL, s10_block(SL10[6:10]),
                               ncol = 3, rel_widths = c(1, 0.05, 1)),
                     ggdraw() + draw_label("Week", fontfamily = FONT, size = 15),
                     ncol = 1, rel_heights = c(0.05, 1, 0.028)) +
    theme(plot.background = element_rect(fill = "white", colour = NA))   # the assembled plot needs an explicit white background
  save_wilke(fig10, file.path(outdir, "fig_nhis_overlay_10seasons"), width = 7.5, height = 8.4)
}
