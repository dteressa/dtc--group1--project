
  
library(tidyverse)
library(fixest)
load("DWI_Data.rdata")


# Get minimum bac of 2 variables and create a variaBle for above 0.08 bac
dwi <- dwi %>%
  mutate(
    bac = pmin(bac1, bac2),
    bac_centered = bac - 0.08,
    above08 = as.numeric(bac >= 0.08)
  )

# Plot for bac between 0.03 and 0.13, filterfor n > 50
plot_data <- dwi %>%
  filter(bac >= 0.03 & bac <= 0.13) %>%
  mutate(bac_bin = round(bac, 3)) %>%
  group_by(bac_bin) %>%
  summarise(
    mean_recid = mean(recidivism),
    n = n()
  ) %>%
  filter(n >= 50)

# Scatter plot forf mean recidivism by bac bins 
ggplot(plot_data, aes(x = bac_bin, y = mean_recid)) +
  geom_point(color = "steelblue", size = 1.5) +
  geom_vline(xintercept = 0.08, linetype = "dashed", 
             color = "red") +
  geom_smooth(data = filter(plot_data, bac_bin < 0.08),
              method = "lm", se = FALSE, color = "thistle") +
  geom_smooth(data = filter(plot_data, bac_bin >= 0.08),
              method = "lm", se = FALSE, color = "thistle")

# Filter to bandwidth of 0.05 around cutoff (0.03 to 0.13)
rdd_data <- dwi %>%
  filter(bac >= 0.03 & bac <= 0.13)

# Baseline model
model1 <- feols(recidivism ~ above08 + bac_centered, data = rdd_data)
summary(model1)


# Adding demographic controls
model2 <- feols(recidivism ~ above08 + bac_centered + male + white + aged + acc, 
                data = rdd_data)
summary(model2)


# Narrow RDD bandwidth 0.055 to 0.105
rdd_narrow <- dwi %>% filter(bac >= 0.055 & bac <= 0.105)
model3 <- feols(recidivism ~ above08 + bac_centered, data = rdd_narrow)
summary(model3)

# Wide RDD bandwidth 0 - 0.18
rdd_wide <- dwi %>% filter(bac >= 0.0 & bac <= 0.18)
model4 <- feols(recidivism ~ above08 + bac_centered, data = rdd_wide)
summary(model4)


# Combined results table
etable(model1, model2, model3, model4,
       headers = c("Baseline", "Controls", "Narrow BW", "Wide BW"),
       dict = c(above08 = "Above 0.08 Cutoff",
                bac_centered = "BAC (centered)",
                "I(bac_centered^2)" = "BAC² (centered)"),
       title = "Effect of Harsher DUI Punishment on Recidivism")

# Our baseline RDD estimate shows that exceeding the 0.08 BAC threshold reduces the 
# probability of future drunk-driving offenses by 2.4%. 
# This almost similar when we add demographics as well
# When we narrow the band to only consider, close to cutoff folks 
# the probability does drop to 1.8% but is not too far from 2.4%. 
# And when we widen the band to get more data we see 2.3% which is again very close to 2.4
# Given a baseline recidivism rate of approximately 12.1% 
# just below the cutoff, this represents 2.4 / 12.1 = 0.198 i.e. a roughly ~20% relative reduction.
