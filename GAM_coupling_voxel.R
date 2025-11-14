rm(list=ls())

# Load required packages
library(mgcv)
library(ggplot2)

# Read the data file
data = read.csv('your_parameter_path.csv')
# Convert the 'Group' column to a factor
data$Group <- as.factor(data$Group)

# Fit a GAM model with smoothing spline for Age and fixed effect for Group
# Increased gamma to further penalize complexity and avoid overfitting
devMod<-gam(Indicator ~ s(Age, k=40) + Group, data = data,  method = "REML", gamma = 2)
data$predicted_indicator <- predict(devMod)
group_effects <- coef(devMod)[grep("Group", names(coef(devMod)))]

# Set the baseline group effect to 0 (for the reference category)
baseline_effect <- 0

# Combine baseline effect with other group effects
group_effects <- c(baseline_effect, group_effects)

# Assign group names to the effect vector for clarity
names(group_effects) <- levels(data$Group) 

# Calculate adjusted indicator values by removing the group-specific effect
data$adjusted_indicator <- data$Indicator - group_effects[data$Group]# 计算调整后的值

# Adjust the predicted values similarly for plotting
data$predicted_indicator <- data$predicted_indicator- group_effects[data$Group]

# Create the visualization
ggplot(data, aes(x = Age, y = adjusted_indicator, color = Group)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(y = predicted_indicator), size = 1) +
  labs(title = "Indicator vs. Age After Systematic Bias Adjustment",
       x = "Age",
       y = "Indicator Value") +
  theme_minimal()
write.table(data,file ="your_save_path.csv", sep =",", row.names = FALSE)
