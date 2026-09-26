## ---------------------------------------------------------------------------
##  Reviewer Figure R1 -- asciminib pilot DMS (16 residues, 304 variants)
##  Response to Reviewer #2, asciminib comment.
##
##  Inputs:
##    K5.ABL.v101_260813_dr4pl_filtered.csv                       -- 4PL fits
##    K5.ABL.v101_260813_Master_Table_Analysis_Wide_filtered.csv  -- rel. viability x conc
##
##  Output: Figure_R1_asciminib.pdf / .png  (6.5 x 4.6 in)
##
##  Colour scheme inherited from the manuscript figures:
##    heat map        blue - white - red   (as Fig. 1D, Supp. Fig. 5)
##    WT residues     yellow               (as Supp. Fig. 5)
##    V468 curves     #D7191C              (Fig. 3D ">1200 nM" bin)
##    F359 curves     #2C7FB8              (Fig. 3D "(300,600]" bin)
##
##  Asciminib data only -- no imatinib arm is plotted or compared.
## ---------------------------------------------------------------------------
# setwd("OneDrive - The Pennsylvania State University/RProjects/abl_dms/Figures/Asciminib_analysis/")
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)
# read.csv("K5.ABL.v1-01_260813_Master_Table_Analysis_Long_filtered.csv")
f_dr4pl <- "K5.ABL.v1-01_260813_dr4pl_filtered.csv"
f_wide  <- "K5.ABL.v1-01_260813_Master_Table_Analysis_Wide_filtered.csv"

AA_ORDER <- strsplit("HKRDECMNQSTAILVFWYGP", "")[[1]]   # manuscript row order
CAP      <- 3                                            # log2 colour cap (= 8-fold)

COL_RESIST <- "#D7191C"
COL_BLUE   <- "#2C7FB8"
COL_WT     <- "#FFFF00"

base_thm <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text       = element_text(size = 6.4, face = "bold"),
        axis.text        = element_text(size = 5.4),
        axis.title       = element_text(size = 6.4),
        legend.title     = element_text(size = 5.6),
        legend.text      = element_text(size = 5.2),
        legend.key.size  = unit(0.22, "cm"),
        plot.tag         = element_text(size = 11, face = "bold"))

## ---- 1. asciminib dose-response fits --------------------------------------

dr <- read.csv(f_dr4pl, stringsAsFactors = FALSE) %>%
  filter(drug == "Asciminib") %>%
  mutate(pos = as.integer(gsub("\\D", "", mut)),
         ref = substr(mut, 1, 1),
         alt = substr(mut, nchar(mut), nchar(mut))) %>%
  group_by(subpool) %>%
  mutate(bg   = median(ic50, na.rm = TRUE),   # background IC50 of the library region
         fc   = ic50 / bg,
         l2fc = log2(fc)) %>%
  ungroup() %>%
  ## a fit is "low-confidence" if it failed, fits poorly, or has an implausibly
  ## steep slope (the signature of an IC50 pinned to the bottom of the range)
  mutate(flag = fit_status != "ok" | r2 < 0.8 | hill > 4)

## ---- 2. observed relative viability (wide -> long) ------------------------
asc <- read.csv(f_wide, stringsAsFactors = FALSE) %>%
  filter(drug == "Asciminib") %>%
  pivot_longer(starts_with("rel_via_"), names_to = "col",
               values_to = "rel_via", values_drop_na = TRUE) %>%
  mutate(conc = as.numeric(gsub("nM$", "", sub("^rel_via_", "", col))))

## the 2,500 nM arm was only run for 5 variants of one region -- drop it
keep <- asc %>% count(subpool, conc) %>% filter(n >= 100) %>% select(subpool, conc)
asc  <- semi_join(asc, keep, by = c("subpool", "conc"))

## dr4pl 4PL with the asymptotes fixed at 0 and 1, as in the paper:
##   f(x) = 1 / (1 + (x / IC50)^hill)
pred4pl <- function(x, ic50, hill) 1 / (1 + (x / ic50)^hill)

## =========================================================== PANEL A =======
reg_lab <- c(S15 = "ABL 354-361", S29 = "ABL 466-473")

res_levels <- dr %>% distinct(pos, ref) %>% arrange(pos) %>%
  transmute(l = paste0(ref, pos)) %>% pull(l)

hm <- dr %>%
  mutate(alt     = factor(alt, levels = AA_ORDER),
         res_lab = factor(paste0(ref, pos), levels = res_levels),
         reg     = factor(reg_lab[subpool], levels = reg_lab))

wt <- hm %>% distinct(reg, res_lab, ref) %>% mutate(alt = factor(ref, levels = AA_ORDER))

## outline F359 and V468: x index *within* each free_x facet
box <- hm %>% filter(pos %in% c(359, 468)) %>% distinct(reg, subpool, pos) %>%
  left_join(hm %>% distinct(subpool, pos) %>% arrange(subpool, pos) %>%
              group_by(subpool) %>% mutate(xi = row_number()) %>% ungroup(),
            by = c("subpool", "pos")) %>%
  mutate(xmin = xi - 0.5, xmax = xi + 0.5, ymin = 0.5, ymax = length(AA_ORDER) + 0.5)

pA <- ggplot(hm, aes(res_lab, alt)) +
  geom_tile(aes(fill = l2fc)) +
  geom_tile(data = wt, fill = COL_WT) +
  geom_point(data = filter(hm, flag), shape = 1, size = 0.4, stroke = 0.25, colour = "black") +
  geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 0.45) +
  facet_wrap(~reg, nrow = 1, scales = "free_x") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-CAP, CAP), oob = scales::squish,
                       breaks = c(-3, -1, 0, 1, 3),
                       labels = c("≤0.13", "0.5", "1", "2", "≥8"),
                       name = expression(atop("IC"[50], "fold-change"))) +
  scale_y_discrete(limits = rev(AA_ORDER)) +
  labs(x = "Residue on BCR-ABL kinase", y = "Substitution") +
  base_thm +
  theme(axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.4),
        axis.text.y   = element_text(size = 4.8),
        panel.grid    = element_blank(),
        panel.spacing = unit(0.6, "lines"))

## =========================================================== PANEL B =======
foc <- tibble::tribble(
  ~mut,    ~subpool, ~col,
  "V468F", "S29",    COL_RESIST,
  "V468K", "S29",    COL_RESIST,
  "V468D", "S29",    COL_RESIST,
  "V468A", "S29",    COL_RESIST,
  "F359V", "S15",    COL_BLUE,
  "F359C", "S15",    COL_BLUE,
  "F359I", "S15",    COL_BLUE,
  "F359D", "S15",    COL_BLUE) %>%
  left_join(select(dr, mut, ic50, hill, fc), by = "mut") %>%
  mutate(facet = factor(sprintf("%s   %.1f×", mut, fc),
                        levels = sprintf("%s   %.1f×", mut, fc)))

band <- asc %>%
  group_by(subpool, conc) %>%
  summarise(lo = quantile(rel_via, .10), md = median(rel_via),
            hi = quantile(rel_via, .90), .groups = "drop")

band_f <- foc %>% select(facet, subpool, col) %>% left_join(band, by = "subpool")
pts_f  <- foc %>% select(facet, mut, col) %>%
  left_join(select(asc, mut, conc, rel_via), by = "mut")
rug_f  <- foc %>% select(facet, subpool) %>%
  left_join(distinct(asc, subpool, conc), by = "subpool")

grid_x <- 10^seq(log10(3), log10(3000), length.out = 200)
fit_f  <- do.call(rbind, lapply(seq_len(nrow(foc)), function(i)
  data.frame(facet = foc$facet[i], col = foc$col[i],
             x = grid_x, y = pred4pl(grid_x, foc$ic50[i], foc$hill[i]))))

pB <- ggplot() +
  geom_ribbon(data = band_f, aes(conc, ymin = lo, ymax = hi), fill = "grey50", alpha = 0.22) +
  geom_line(data = band_f, aes(conc, md), colour = "grey40", linetype = "dashed", linewidth = 0.3) +
  geom_rug(data = rug_f, aes(x = conc), sides = "b", length = unit(0.04, "npc"), linewidth = 0.2) +
  geom_line(data = fit_f, aes(x, y, colour = col), linewidth = 0.45) +
  geom_point(data = pts_f, aes(conc, rel_via, fill = col), shape = 21, size = 0.9,
             stroke = 0.15, colour = "white") +
  facet_wrap(~facet, nrow = 2) +
  scale_colour_identity() + scale_fill_identity() +
  scale_x_log10(limits = c(3, 3000), breaks = c(10, 100, 1000),
                labels = c("10", "100", "1000")) +
  coord_cartesian(ylim = c(-0.05, 1.35)) +
  labs(x = "Asciminib (nM)", y = "Relative viability") +
  base_thm + theme(panel.spacing = unit(0.35, "lines"))

## =========================================================== assemble ======
fig <- pA / pB +
  plot_layout(heights = c(2.20, 2.35)) +
  plot_annotation(tag_levels = "A")

ggsave("Figure_R1_asciminib.pdf", fig, width = 6.5, height = 4.6, units = "in")
ggsave("Figure_R1_asciminib.png", fig, width = 6.5, height = 4.6, units = "in", dpi = 400)

## ---------------------------------------------------------------------------
## Numbers quoted in the response letter
## ---------------------------------------------------------------------------
cat("\n-- background (regional median) asciminib IC50, nM --\n")
print(distinct(dr, subpool, bg))

for (p in c(468, 359)) {
  sp <- ifelse(p == 468, "S29", "S15")
  g  <- filter(dr, subpool == sp)
  a  <- filter(g, pos == p)$fc; b <- filter(g, pos != p)$fc
  cat(sprintf("\npos %d: %d/19 > 3x, %d/19 > 10x, max %.1fx; median %.2f vs %.2f (rest); Wilcoxon p = %.2g\n",
              p, sum(a > 3), sum(a > 10), max(a), median(a), median(b),
              wilcox.test(a, b, alternative = "greater")$p.value))
  print(filter(g, pos == p) %>% arrange(desc(fc)) %>%
          select(mut, ic50, fc, hill, r2, fit_status))
}
