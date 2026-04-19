rm(list=ls())

# Load required packages
library(mgcv)
library(ggplot2)

# Read the data file
data = read.csv('your_parameter_path.csv')
# Convert the 'Group' column to a factor
data$Group <- as.factor(data$Group)

# Model: control for both Group and mean_amp_FC simultaneously
# Note: mean_amp_FC is a numeric variable; s() is not needed unless a smooth term is desired
devMod <- gam(Indicator ~ s(Age, k=40) + Group + mean_amp_FC, data = data, method = "REML", gamma = 2)

# View model summary
summary(devMod)

# Extract Group effect coefficients for subsequent adjustment
group_effects <- coef(devMod)[grep("Group", names(coef(devMod)))]

# Add baseline group effect as 0
baseline_effect <- 0
group_effects <- c(baseline_effect, group_effects)
names(group_effects) <- levels(data$Group)

# Extract coefficient for mean_amp_FC
FC_coef <- coef(devMod)["mean_amp_FC"]

# Compute adjusted indicator values (removing effects of both Group and mean_amp_FC)
# Note: The effect of Age is retained, as we are interested in its relationship with the indicator
data$adjusted_indicator <- data$Indicator - group_effects[data$Group] - FC_coef * data$mean_amp_FC

# Predicted values for visualization (with Group and FC effects removed)
data$predicted_indicator <- predict(devMod) - group_effects[data$Group] - FC_coef * data$mean_amp_FC

# Visualization
ggplot(data, aes(x = Age, y = adjusted_indicator, color = Group)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = predicted_indicator), size = 1) +
  labs(title = "Age-related trends of indicator after regressing out Group and FC effects",
       x = "Age",
       y = "Adjusted indicator value") +
  theme_minimal()

# Save results
write.table(data,file ="your_save_path.csv", sep =",", row.names = FALSE)
