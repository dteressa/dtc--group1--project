setwd("/Users/alekhyak/Documents/New project1")
load("DWI_Data.rdata")
ls()
View(dwi)
library(dplyr)

dwi <- dwi %>%
  mutate(
    bac_min = pmin(bac1, bac2, na.rm = TRUE),
    running = bac_min - 0.08,
    DUI = if_else(bac_min >= 0.08, 1, 0)
  )
summary(dwi$bac_min)
mean(dwi$bac_min, na.rm = TRUE)
median(dwi$bac_min, na.rm = TRUE)
sum(dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11)
mean(dwi$recidivism[dwi$bac_min >= 0.05 & 
                      dwi$bac_min < 0.08], na.rm = TRUE)

mean(dwi$recidivism[dwi$bac_min >= 0.08 & 
                      dwi$bac_min <= 0.11], na.rm = TRUE)
ggplot(dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ],
       aes(x = bac_min)) +
  geom_histogram(binwidth = 0.002, color = "black", fill = "lightblue") +
  geom_vline(xintercept = 0.08, color = "red", linetype = "dashed")
ggplot(dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ],
       aes(x = bac_min, y = aged)) +
  geom_smooth() +
  geom_vline(xintercept = 0.08, linetype = "dashed")
ggplot(dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ],
       aes(x = bac_min, y = male)) +
  geom_smooth() +
  geom_vline(xintercept = 0.08, linetype = "dashed")
ggplot(dwi[dwi$bac_min >= 0.05 & dwi$bac_min <= 0.11, ],
       aes(x = bac_min, y = acc)) +
  geom_smooth() +
  geom_vline(xintercept = 0.08, linetype = "dashed")
