# ============================================================
# rd_design.R
# RD design and main models for DWI project (cutoff = 0.08)
# ============================================================

# --- Packages ---
library(tidyverse)
library(here)
library(rdrobust)
library(rddensity)
library(fixest)

# --- Cutoff ---
c_cut <- 0.08

# --- Load data ---
dwi_raw <- readr::read_csv(here("data", "DWI_Data.csv"))

# --- Construct running variable and treatment ---
dwi <- dwi_raw %>%
  mutate(
    bac_min = pmin(bac1, bac2, na.rm = TRUE),
    treat   = as.integer(bac_min >= c_cut),
    x       = bac_min - c_cut
  )

# ------------------------------------------------------------
# RD sample (exclude higher 0.15 cutoff region)
# ------------------------------------------------------------
dwi08 <- dwi %>%
  filter(!is.na(bac_min), !is.na(recidivism)) %>%
  filter(bac_min < 0.15)

plot_window <- dwi08 %>%
  filter(bac_min >= 0.03, bac_min <= 0.13)

# ------------------------------------------------------------
# Running variable distribution
# ------------------------------------------------------------
p_hist <- ggplot(plot_window, aes(x = bac_min)) +
  geom_histogram(binwidth = 0.001, color = "white") +
  geom_vline(xintercept = c_cut) +
  labs(
    title = "BAC distribution near 0.08",
    x = "BAC (bac_min)",
    y = "Count"
  )

print(p_hist)

# ------------------------------------------------------------
# Density test at cutoff
# ------------------------------------------------------------
dens_08 <- rddensity(X = dwi08$bac_min, c = c_cut)
print(summary(dens_08))

rdplotdensity(dens_08, X = dwi08$bac_min)

# ------------------------------------------------------------
# RD plot
# ------------------------------------------------------------
rdplot(
  y = dwi08$recidivism,
  x = dwi08$bac_min,
  c = c_cut,
  x.label = "BAC (bac_min)",
  y.label = "Recidivism",
  title   = "Recidivism vs BAC (RD at 0.08)"
)

# ------------------------------------------------------------
# Main RD estimate
# ------------------------------------------------------------
rd_main <- rdrobust(
  y = dwi08$recidivism,
  x = dwi08$bac_min,
  c = c_cut
)

print(summary(rd_main))

# ------------------------------------------------------------
# Bandwidth sensitivity
# ------------------------------------------------------------
rd_small <- rdrobust(
  y = dwi08$recidivism,
  x = dwi08$bac_min,
  c = c_cut,
  h = 0.02
)

rd_big <- rdrobust(
  y = dwi08$recidivism,
  x = dwi08$bac_min,
  c = c_cut,
  h = 0.05
)

cat("\n--- Bandwidth sensitivity ---\n")
cat("\nSMALL h = 0.02\n"); print(summary(rd_small))
cat("\nBIG   h = 0.05\n"); print(summary(rd_big))

# ------------------------------------------------------------
# Covariate balance checks
# ------------------------------------------------------------
covars <- c("male", "white", "aged", "acc")

balance_tests <- purrr::map_dfr(covars, function(v) {
  yv <- dwi08[[v]]
  keep <- !is.na(yv) & !is.na(dwi08$bac_min)
  
  out <- rdrobust(
    y = yv[keep],
    x = dwi08$bac_min[keep],
    c = c_cut
  )
  
  tibble(
    var = v,
    tau = out$Estimate[1],
    se  = out$se[1],
    p   = out$pv[1]
  )
})

print(balance_tests)

# ------------------------------------------------------------
# Parametric local linear check
# ------------------------------------------------------------
h_feols <- 0.05

dwi_h <- dwi08 %>%
  filter(abs(x) <= h_feols)

m_ll <- feols(
  recidivism ~ treat + x + treat:x,
  data = dwi_h,
  vcov = "hetero"
)

etable(m_ll)