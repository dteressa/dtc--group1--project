library(dplyr)
library(ggplot2)
setwd("/Users/alekhyak/Documents/New project1")
load("DWI_Data.rdata")
ls()
View(dwi)
dwi <- dwi %>% mutate(bac_min = pmin(bac1, bac2), running = bac_min - 0.08, DUI = ifelse(bac_min >= 0.08, 1, 0))
summary(dwi$bac_min)
dwi <- dwi %>% mutate(DUI = ifelse(bac_min >= 0.08, 1, 0))
table(dwi$DUI)
table(dwi$DUI[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11])
#no bins yet
ggplot(dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ],
       aes(x = bac_min, y = recidivism)) +
  geom_smooth() +
  geom_vline(xintercept = 0.08, linetype = "dashed")
#with bins
ggplot(aggregate(recidivism ~ round(bac_min, 3),
                 data = dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ],
                 mean),
       aes(x = `round(bac_min, 3)`, y = recidivism)) +
  geom_point() +
  geom_vline(xintercept = 0.08, linetype = "dashed")
#Running Variable
dwi <- dwi %>% mutate(running = bac_min - 0.08)
#local liner RD
m1 <- lm(recidivism ~ DUI + running + DUI:running,
         data = dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ])

summary(m1)
#m2 - Narrow local window
m2 <- lm(recidivism ~ DUI + running + DUI:running,
         data = dwi[dwi$bac_min >= 0.06 & dwi$bac_min <= 0.10, ])

summary(m2)
#m3 - Wider local window 
m3 <- lm(recidivism ~ DUI + running + DUI:running,
         data = dwi[dwi$bac_min >= 0.04 & dwi$bac_min <= 0.12, ])

summary(m3)
#adding covariates
m_cov <- lm(recidivism ~ DUI + running + DUI:running + male + white + aged + factor(year),
            data = dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ])
summary(m_cov)
#Trying Quadratic 
m_quad <- lm(recidivism ~ DUI + running + I(running^2) + DUI:running + DUI:I(running^2),
             data = dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ])
summary(m_quad)
#robust standard errors
library(sandwich)
library(lmtest)
coeftest(m1, vcov = vcovHC(m1, type = "HC1"))
coeftest(m_cov, vcov = vcovHC(m_cov, type = "HC1"))
#Visual RD for m1(binned means + fitted RD)
bins <- aggregate(recidivism ~ round(bac_min, 3),
                  data = dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ],
                  mean)

grid <- data.frame(bac_min = seq(0.05, 0.11, by = 0.001))
grid$running <- grid$bac_min - 0.08
grid$DUI <- ifelse(grid$bac_min >= 0.08, 1, 0)
grid$yhat <- predict(m1, newdata = grid)

ggplot(bins, aes(x = `round(bac_min, 3)`, y = recidivism)) +
  geom_point() +
  geom_vline(xintercept = 0.08, linetype = "dashed") +
  geom_line(data = grid, aes(x = bac_min, y = yhat)) +
  labs(x = "BAC (binned)", y = "Mean recidivism", title = "RD Plot with Fitted Lines (m1)")
#Visual RD for m_cov
bins <- aggregate(recidivism ~ round(bac_min, 3),
                  data = dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ],
                  mean)

grid2 <- data.frame(bac_min = seq(0.05, 0.11, by = 0.001))
grid2$running <- grid2$bac_min - 0.08
grid2$DUI <- ifelse(grid2$bac_min >= 0.08, 1, 0)
grid2$male <- mean(dwi$male, na.rm = TRUE)
grid2$white <- mean(dwi$white, na.rm = TRUE)
grid2$aged <- mean(dwi$aged, na.rm = TRUE)
grid2$year <- factor(1999, levels = levels(factor(dwi$year)))  # baseline year
grid2$yhat <- predict(m_cov, newdata = grid2)

ggplot(bins, aes(x = `round(bac_min, 3)`, y = recidivism)) +
  geom_point() +
  geom_vline(xintercept = 0.08, linetype = "dashed") +
  geom_line(data = grid2, aes(x = bac_min, y = yhat)) +
  labs(x = "BAC (binned)", y = "Mean recidivism", title = "RD Plot with Fitted Lines (m_cov)")