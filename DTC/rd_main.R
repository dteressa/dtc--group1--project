# ============================================================
# rd_main.R
# Main RD analysis for DWI project at BAC cutoff = 0.08
# ============================================================

# ----------------------------
# Packages
# ----------------------------
library(tidyverse)
library(here)
library(rdrobust)
library(rddensity)
library(fixest)
library(sandwich)
library(lmtest)

# ----------------------------
# Cutoff
# ----------------------------
c_cut <- 0.08

# ----------------------------
# Load data
# Use a reproducible path so the script runs for everyone
# ----------------------------
dwi_raw <- readr::read_csv(here("data", "DWI_Data.csv"))

# ----------------------------
# Construct RD variables
# Use bac_min as the legally relevant BAC measure
# Set up treatment indicator and centered running variable
# ----------------------------
dwi <- dwi_raw %>%
  mutate(
    bac_min = pmin(bac1, bac2, na.rm = TRUE),
    treat   = as.integer(bac_min >= c_cut),
    x       = bac_min - c_cut
  )

# ----------------------------
# Sample restriction
# Keep observations below 0.15 so the higher cutoff does not interfere
# Also define a local window used for visuals and an added robustness check
# ----------------------------
dwi08 <- dwi %>%
  filter(!is.na(bac_min), !is.na(recidivism)) %>%
  filter(bac_min < 0.15)

plot_window <- dwi08 %>%
  filter(bac_min >= 0.03, bac_min <= 0.13)

# ----------------------------
# Running variable distribution
# Histogram near the cutoff using the local plotting window
# ----------------------------
ggplot(plot_window, aes(x = bac_min)) +
  geom_histogram(binwidth = 0.001, color = "white") +
  geom_vline(xintercept = c_cut) +
  labs(
    title = "BAC distribution near 0.08",
    x = "BAC (bac_min)",
    y = "Count"
  )

# ----------------------------
# Density test at cutoff
# Formal density test and visualization
# Used as part of the RD validity checks
# ----------------------------
dens_08 <- rddensity(X = dwi08$bac_min, c = c_cut)
summary(dens_08)
rdplotdensity(dens_08, X = dwi08$bac_min)

# ----------------------------
# RD plot
# Main RD visualization using the local plotting window
# ----------------------------
rdplot(
  y = plot_window$recidivism,
  x = plot_window$bac_min,
  c = c_cut,
  x.label = "BAC (bac_min)",
  y.label = "Recidivism",
  title   = "Recidivism vs BAC (RD at 0.08)"
)

# ----------------------------
# Main RD estimate
# rdrobust used as the primary specification
# ----------------------------
rd_main <- rdrobust(
  y = dwi08$recidivism,
  x = dwi08$bac_min,
  c = c_cut
)

summary(rd_main)

# ----------------------------
# Local-window RD estimate
# Additional robustness check using the narrower window
# ----------------------------
rd_main_window <- rdrobust(
  y = plot_window$recidivism,
  x = plot_window$bac_min,
  c = c_cut
)

summary(rd_main_window)

# ----------------------------
# Bandwidth sensitivity
# Manual bandwidth checks around the cutoff
# ----------------------------
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

cat("\n--- rdrobust bandwidth sensitivity ---\n")
cat("\nMain estimate\n")
print(summary(rd_main))

cat("\nLocal-window estimate (0.03 to 0.13)\n")
print(summary(rd_main_window))

cat("\nSmall bandwidth h = 0.02\n")
print(summary(rd_small))

cat("\nBig bandwidth h = 0.05\n")
print(summary(rd_big))

# ----------------------------
# Covariate balance checks
# Formal balance tests for pre-treatment covariates
# ----------------------------
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

balance_tests

# ----------------------------
# Parametric local linear checks
# Local linear regression models used as supporting checks
# ----------------------------

# Baseline local window: 0.05 to 0.11
m1 <- lm(
  recidivism ~ treat + x + treat:x,
  data = dwi %>% filter(bac_min >= 0.05, bac_min <= 0.11)
)

summary(m1)

# Narrower window: 0.06 to 0.10
m2 <- lm(
  recidivism ~ treat + x + treat:x,
  data = dwi %>% filter(bac_min >= 0.06, bac_min <= 0.10)
)

summary(m2)

# Wider window: 0.04 to 0.12
m3 <- lm(
  recidivism ~ treat + x + treat:x,
  data = dwi %>% filter(bac_min >= 0.04, bac_min <= 0.12)
)

summary(m3)

# Model with covariates
m_cov <- lm(
  recidivism ~ treat + x + treat:x + male + white + aged + factor(year),
  data = dwi %>% filter(bac_min >= 0.05, bac_min <= 0.11)
)

summary(m_cov)

# Quadratic specification
m_quad <- lm(
  recidivism ~ treat + x + I(x^2) + treat:x + treat:I(x^2),
  data = dwi %>% filter(bac_min >= 0.05, bac_min <= 0.11)
)

summary(m_quad)

# Robust standard errors
coeftest(m1, vcov = vcovHC(m1, type = "HC1"))
coeftest(m_cov, vcov = vcovHC(m_cov, type = "HC1"))

# ----------------------------
# Parametric local linear check in fixest
# Heteroskedasticity-robust linear check in a fixed bandwidth
# ----------------------------
h_feols <- 0.05

dwi_h <- dwi08 %>%
  filter(abs(x) <= h_feols)

m_ll <- feols(
  recidivism ~ treat + x + treat:x,
  data = dwi_h,
  vcov = "hetero"
)

etable(m_ll)

# ----------------------------
# Binned RD visual with fitted local linear line
# Visual representation of the local linear model
# ----------------------------
bins <- aggregate(
  recidivism ~ round(bac_min, 3),
  data = dwi %>% filter(bac_min >= 0.05, bac_min <= 0.11),
  mean
)

grid <- data.frame(bac_min = seq(0.05, 0.11, by = 0.001))
grid$x <- grid$bac_min - c_cut
grid$treat <- ifelse(grid$bac_min >= c_cut, 1, 0)
grid$yhat <- predict(m1, newdata = grid)

ggplot(bins, aes(x = `round(bac_min, 3)`, y = recidivism)) +
  geom_point() +
  geom_vline(xintercept = c_cut, linetype = "dashed") +
  geom_line(data = grid, aes(x = bac_min, y = yhat)) +
  labs(
    x = "BAC (binned)",
    y = "Mean recidivism",
    title = "RD plot with fitted local linear line"
  )