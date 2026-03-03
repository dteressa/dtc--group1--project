library(dplyr)
library(ggplot2)
library(fixest)   # keeping since you already load it; not used below
library(dplyr)

# Load the data
load("DWI_Data.rdata")

dwi <- dwi %>%
  mutate(
    bac_min = pmin(bac1, bac2, na.rm = TRUE),
    running = bac_min - 0.08,
    DUI = if_else(bac_min >= 0.08, 1, 0)
  )

# Summary statistics and descriptive checks
summary(dwi$bac_min)
mean(dwi$bac_min, na.rm = TRUE)
median(dwi$bac_min, na.rm = TRUE)
sum(dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11)

# Recidivism just below vs above cutoff
mean(dwi$recidivism[dwi$bac_min >= 0.05 & dwi$bac_min < 0.08], na.rm = TRUE)
mean(dwi$recidivism[dwi$bac_min >= 0.08 & dwi$bac_min <= 0.11], na.rm = TRUE)

# Diagnostic plots
ggplot(dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ], aes(x = bac_min)) +
  geom_histogram(binwidth = 0.002, color = "black", fill = "lightblue") +
  geom_vline(xintercept = 0.08, color = "red", linetype = "dashed") +
  labs(x = "BAC (min of two tests)", y = "Count",
       title = "Distribution of BAC around 0.08 cutoff")

# Covariate smoothness checks
ggplot(dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ], aes(x = bac_min, y = aged)) +
  geom_smooth() + geom_vline(xintercept = 0.08, linetype = "dashed") +
  labs(y = "Age")

ggplot(dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ], aes(x = bac_min, y = male)) +
  geom_smooth() + geom_vline(xintercept = 0.08, linetype = "dashed") +
  labs(y = "Male Indicator")

ggplot(dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ], aes(x = bac_min, y = acc)) +
  geom_smooth() + geom_vline(xintercept = 0.08, linetype = "dashed") +
  labs(y = "Accident Indicator")

dwi_local <- dwi %>%
  filter(bac_min >= 0.05 & bac_min <= 0.11)

# raw difference around cutoff (descriptive)
mean_below <- mean(dwi_local$recidivism[dwi_local$bac_min < 0.08], na.rm = TRUE)
mean_above <- mean(dwi_local$recidivism[dwi_local$bac_min >= 0.08], na.rm = TRUE)

cat("Average recidivism below cutoff (BAC < 0.08):", round(mean_below,4), "\n")
cat("Average recidivism above cutoff (BAC >= 0.08):", round(mean_above,4), "\n")
cat("Raw difference (above - below):", round(mean_above - mean_below,4), "\n")

# descriptive smooth plot
ggplot(dwi_local, aes(x = bac_min, y = recidivism)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "blue") +
  geom_vline(xintercept = 0.08, linetype = "dashed", color = "red")


# ============================================================
# REPLACEMENT FOR SHARP RDD: Logistic regression approaches
# ============================================================

# -------------------------
# 1) Logistic regression models (binary outcome)
# -------------------------
# logit1: running only (descriptive)
logit1 <- glm(
  recidivism ~ running,
  data = dwi_local,
  family = binomial(link = "logit")
)

# logit2: running + treatment indicator (jump in log-odds)
logit2 <- glm(
  recidivism ~ running + DUI,
  data = dwi_local,
  family = binomial(link = "logit")
)

# logit3: running * DUI (allows different slopes in log-odds)
logit3 <- glm(
  recidivism ~ running * DUI,
  data = dwi_local,
  family = binomial(link = "logit")
)

summary(logit1)
summary(logit2)
summary(logit3)

# Optional: compare models with AIC (lower is better fit)
AIC(logit1, logit2, logit3)


# -------------------------
# 2) Covariate-adjusted logistic regression
# -------------------------
logit_cov <- glm(
  recidivism ~ running * DUI + male + white + aged + factor(year),
  data = dwi_local,
  family = binomial(link = "logit")
)
summary(logit_cov)


# -------------------------
# 3) Bandwidth sensitivity using logistic regression
# -------------------------
windows <- list(c(0.06, 0.10), c(0.05, 0.11), c(0.04, 0.12))
names(windows) <- c("0.06-0.10", "0.05-0.11", "0.04-0.12")

logit_bw <- lapply(windows, function(w){
  dat <- dwi %>% filter(bac_min >= w[1], bac_min <= w[2])
  glm(recidivism ~ running * DUI, data = dat, family = binomial(link = "logit"))
})

# Print just the DUI and interaction terms across bandwidths
for(nm in names(logit_bw)){
  cat("\n====================\nLogit bandwidth:", nm, "\n====================\n")
  print(coef(summary(logit_bw[[nm]]))[c("DUI", "running:DUI"), , drop = FALSE])
}


# -------------------------
# 4) Predicted probability plot around cutoff from logistic model
# -------------------------
# Use logit3 as the main logistic specification
grid <- data.frame(bac_min = seq(min(dwi_local$bac_min), max(dwi_local$bac_min), by = 0.001))
grid$running <- grid$bac_min - 0.08
grid$DUI <- if_else(grid$bac_min >= 0.08, 1, 0)

# Predicted probability of recidivism
grid$phat <- predict(logit3, newdata = grid, type = "response")

# Binned means for visualization (same as your earlier idea)
bins <- aggregate(recidivism ~ round(bac_min, 3), data = dwi_local, mean)

ggplot(bins, aes(x = `round(bac_min, 3)`, y = recidivism)) +
  geom_point() +
  geom_vline(xintercept = 0.08, linetype = "dashed") +
  geom_line(data = grid, aes(x = bac_min, y = phat)) +
  labs(
    x = "BAC (binned to 0.001)",
    y = "Mean recidivism / Predicted probability",
    title = "Logistic fit around cutoff (predicted probability)"
  )