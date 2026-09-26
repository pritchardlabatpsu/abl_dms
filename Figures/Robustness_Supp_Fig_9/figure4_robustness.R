# =============================================================================
# Figure 4 (revised): Robustness of variant resistance classification to
# simulated inter-patient PK variability, stratified by clinical dose class.
#
# Inputs : ABL_Supplemental_tables.xlsx  (sheets Table_S10 and Table_S14)
# Outputs: figure4_robustness_panels.pdf / .png
#
# Reads S14 (resistance_probability, class_robustness per variant) and joins it
# to S10 (Highest_Resistant_Dose = Sensitive / 444nM / 760nM / 916nM). The
# scientific claim it tests: variants that are resistant at *higher* doses are
# *more robustly* resistant to PK variability, while 444nM-only variants are the
# PK-sensitive/borderline population.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)   # optional; only used to compose panels
})

xlsx_path <- "ABL_Supplemental_tables.xlsx"   # adjust path as needed

# --- Load. Both sheets carry a title in row 1 and the header in row 3, so skip 2.
s10 <- read_excel(xlsx_path, sheet = "Table_S10", skip = 2)
s14 <- read_excel(xlsx_path, sheet = "Table_S14", skip = 2)

mk_key <- function(ref, pos, alt) paste0(ref, as.integer(pos), alt)

s10 <- s10 %>%
  mutate(variant = mk_key(ref_aa, protein_start, alt_aa)) %>%
  select(variant, Highest_Resistant_Dose)

s14 <- s14 %>%
  mutate(variant = mk_key(ref_aa, protein_start, alt_aa))

dat <- s14 %>%
  inner_join(s10, by = "variant") %>%
  filter(Highest_Resistant_Dose %in% c("Sensitive", "444nM", "760nM", "916nM")) %>%
  mutate(
    dose_class = factor(Highest_Resistant_Dose,
                        levels = c("Sensitive", "444nM", "760nM", "916nM"),
                        labels = c("Sensitive", "444 nM\n(400 mg QD)",
                                   "760 nM\n(400 mg BID)", "916 nM\n(500 mg BID)")),
    class_robustness = factor(class_robustness,
                        levels = c("Robust sensitive",
                                   "PK-sensitive borderline",
                                   "Robust resistant"))
  )

dose_cols <- c("Sensitive"           = "#7fb0d3",
               "444 nM\n(400 mg QD)" = "#2166ac",
               "760 nM\n(400 mg BID)"= "#f4a582",
               "916 nM\n(500 mg BID)"= "#b2182b")

thr <- 0.5   # a variant is called resistant in a patient if resistance_probability
             # crosses; the 0.1 / 0.9 robustness cuts are drawn as guides below.

# ---- Panel A: simulated patient exposure distribution (eff.Cave) at 400 mg QD.
#      Draws the population of serum-adjusted exposures used for the per-variant
#      resistance-probability calculation; the red line is the 444 nM threshold.
mean_eff <- 444; cv_main <- 0.157
set.seed(42)
expo <- data.frame(eff_cave = rnorm(10000, mean_eff, cv_main * mean_eff))
# install.packages("ggplot2")
pA <- ggplot(expo, aes(eff_cave)) +
  geom_density(fill = "#c98b8b", alpha = 0.55, colour = "black") +
  geom_vline(xintercept = mean_eff, colour = "red", linewidth = 1) +
  labs(x = "eff.Cave (nM)", y = "Density",
       title = "A  Simulated patient exposure (400 mg QD)") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# ---- Panel B: violin/box of resistance probability by dose class
pB <- ggplot(dat, aes(dose_class, resistance_probability, fill = dose_class)) +
  geom_violin(scale = "width", alpha = 0.5, colour = NA) +
  geom_boxplot(width = 0.16, outlier.size = 0.4, alpha = 0.9) +
  scale_fill_manual(values = dose_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Resistance probability",
       title = "B  More resistant classes are more robustly resistant") +
  theme_bw(base_size = 11)

# ---- Panel C: stacked robustness composition within each dose class
comp <- dat %>%
  count(dose_class, class_robustness) %>%
  group_by(dose_class) %>%
  mutate(frac = n / sum(n))

rob_cols <- c("Robust sensitive"        = "#7fb0d3",
              "PK-sensitive borderline" = "#fdb863",
              "Robust resistant"        = "#b2182b")

pC <- ggplot(comp, aes(dose_class, frac, fill = class_robustness)) +
  geom_col(width = 0.7, colour = "white") +
  scale_fill_manual(values = rob_cols, name = "Classification\nrobustness") +
  scale_y_continuous(labels = percent) +
  labs(x = NULL, y = "Fraction of variants",
       title = "C  Robustness composition within each dose class") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

# ---- Panel D: ranked median relative viability +/- simulated PK spread.
#      (This is the "improved Figure 4" style plot: one point per variant,
#       ranked, with the interpatient exposure band shown as a vertical range.)
#      Relative viability at the three serum-adjusted doses stands in for the
#      simulated exposure spread; the dashed line is the resistance threshold.
rank_dat <- dat %>%
  rowwise() %>%
  mutate(
    v_lo  = min(c_across(c(rel_viab_444, rel_viab_760, rel_viab_916)), na.rm = TRUE),
    v_hi  = max(c_across(c(rel_viab_444, rel_viab_760, rel_viab_916)), na.rm = TRUE),
    v_med = rel_viab_444
  ) %>%
  ungroup() %>%
  arrange(v_med) %>%
  mutate(rank = row_number())

res_threshold <- 0.420968716

pD <- ggplot(rank_dat, aes(rank, v_med, colour = dose_class)) +
  geom_linerange(aes(ymin = v_lo, ymax = v_hi), alpha = 0.25, linewidth = 0.2) +
  geom_point(size = 0.35) +
  geom_hline(yintercept = res_threshold, linetype = 2) +
  scale_colour_manual(values = dose_cols, name = "Highest\nresistant dose") +
  labs(x = "Variants ranked by relative viability at 444 nM",
       y = "Relative viability",
       title = "D  Ranked variants with dose-exposure spread") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

# ---- Panel E: CV sensitivity sweep.
#      Reproduces each variant's resistance probability at a grid of assumed
#      inter-patient CVs (5-60%) from the fitted curve parameters, then reports
#      the % robustly resistant within each resistant dose class and the overall
#      % of variants receiving a definitive (robust) call.
res_threshold <- 0.420968716
mean_eff      <- 444
set.seed(42)

relviab <- function(conc, ic50, b) 1 / (1 + (conc / ic50)^abs(b))

sweep_one <- function(cv, n = 40000) {
  samp <- rnorm(n, mean_eff, cv * mean_eff); samp <- samp[samp > 0]
  # variants x patients matrix of resistance calls
  rv <- outer(s14$IC50, samp, function(ic, cc) NA)  # placeholder for clarity
  prob <- vapply(seq_len(nrow(s14)), function(i) {
    mean(relviab(samp, s14$IC50[i], s14$Hill.slope[i]) > res_threshold)
  }, numeric(1))
  cls <- dplyr_left <- s10$Highest_Resistant_Dose[match(s14$variant, s10$variant)]
  data.frame(
    cv = cv * 100,
    overall_robust = mean(prob < 0.1 | prob > 0.9) * 100,
    pct_rr_444 = mean(prob[cls == "444nM"] > 0.9) * 100,
    pct_rr_760 = mean(prob[cls == "760nM"] > 0.9) * 100,
    pct_rr_916 = mean(prob[cls == "916nM"] > 0.9) * 100
  )
}
# NOTE: s10 above must retain the Highest_Resistant_Dose column; if you trimmed
# it earlier, re-read Table_S10 keeping that column before running the sweep.
cv_grid <- sort(unique(c(seq(0.05, 0.60, by = 0.025), 0.157)))
sweep   <- do.call(rbind, lapply(cv_grid, sweep_one))

sweep_long <- sweep %>%
  pivot_longer(c(pct_rr_916, pct_rr_760, pct_rr_444, overall_robust),
               names_to = "series", values_to = "pct") %>%
  mutate(series = recode(series,
           pct_rr_916 = "916 nM (500 mg BID)",
           pct_rr_760 = "760 nM (400 mg BID)",
           pct_rr_444 = "444 nM (400 mg QD)",
           overall_robust = "All variants: definitive call"))

series_cols <- c("916 nM (500 mg BID)" = "#b2182b",
                 "760 nM (400 mg BID)" = "#f4a582",
                 "444 nM (400 mg QD)"  = "#2166ac",
                 "All variants: definitive call" = "grey50")

pE <- ggplot(sweep_long, aes(cv, pct, colour = series,
                             linetype = series == "All variants: definitive call")) +
  annotate("rect", xmin = 40, xmax = 60, ymin = -Inf, ymax = Inf,
           fill = "grey70", alpha = 0.10) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
  geom_vline(xintercept = 15.7, linetype = 3) +
  scale_colour_manual(values = series_cols, name = NULL) +
  scale_linetype_manual(values = c(`TRUE` = 2, `FALSE` = 1), guide = "none") +
  labs(x = "Assumed inter-patient PK variability (CV, %)",
       y = "% robustly resistant\n(or % definitive, dashed)",
       title = "E  Dose-escalation calls are stable across PK-variability assumptions") +
  theme_bw(base_size = 11)

# ---- Compose & save
#      Layout: top row = A (exposure) | B (violin) | C (composition);
#              then D (ranked, full width); then E (CV sweep, full width).
#      The redundant per-class histogram panel has been dropped in favor of the
#      violin/box (B), which shows the same per-variant probability distribution.
top    <- pA + pB + pC + plot_layout(widths = c(1, 1, 1))
figure <- top / pD / pE + plot_layout(heights = c(1, 1, 0.95))

ggsave("figure_combined_robustness.pdf", figure, width = 14, height = 12.5)
ggsave("figure_combined_robustness.png", figure, width = 14, height = 12.5, dpi = 200)

# ---- Focused single-panel version (directly answers "are more resistant
#      mutants more robustly resistant"): violin + box + jittered points.
set.seed(0)
pMain <- ggplot(dat, aes(dose_class, resistance_probability)) +
  geom_jitter(aes(colour = dose_class), width = 0.28, size = 0.5, alpha = 0.28) +
  geom_violin(aes(fill = dose_class), scale = "width", alpha = 0.30, colour = NA) +
  geom_boxplot(aes(colour = dose_class), width = 0.14, fill = "white",
               outlier.shape = NA, linewidth = 0.7) +
  geom_hline(yintercept = c(0.1, 0.9), linetype = 3, colour = "grey50") +
  scale_fill_manual(values = dose_cols, guide = "none") +
  scale_colour_manual(values = dose_cols, guide = "none") +
  scale_y_continuous(limits = c(-0.05, 1.05)) +
  labs(x = "Highest clinically resistant dose",
       y = "Resistance probability across simulated patients",
       title = "Higher-dose resistant variants are more robustly resistant to PK variability") +
  theme_bw(base_size = 11)
ggsave("figure_robustness_by_dose_single.pdf", pMain, width = 8, height = 5.5)
ggsave("figure_robustness_by_dose_single.png", pMain, width = 8, height = 5.5, dpi = 200)

# ---- Summary table printed to console (useful for the response text)
summ <- dat %>%
  group_by(dose_class) %>%
  summarise(n = n(),
            mean_prob   = mean(resistance_probability),
            median_prob = median(resistance_probability),
            pct_robust_resistant = mean(resistance_probability > 0.9) * 100,
            pct_robust_sensitive = mean(resistance_probability < 0.1) * 100,
            .groups = "drop")
print(summ)

# Spearman: does IC50 track resistance probability?
message(sprintf("Spearman(log10 IC50, resistance prob) = %.3f",
        suppressWarnings(cor(log10(s14$IC50), s14$resistance_probability,
                             method = "spearman", use = "complete.obs"))))
