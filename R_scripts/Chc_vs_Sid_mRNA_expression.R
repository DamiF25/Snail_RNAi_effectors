# Install necessary packages

install.packages("ggplot2")
install.packages("dplyr")
install.packages("reshape2")
install.packages("tidyr") #For data transformation to long data

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggpubr) # for stat_compare_means

#set working directory

# Load the data
data <- read.csv("Chc_vs_Sid.csv")

# Perform Shapiro-Wilk test for each group
shapiro_results <- data %>%
  group_by(Tissue) %>%
  summarise(
    shapiro_Chc1_p = shapiro.test(Chc1)$p.value,
    shapiro_Sid1_p = shapiro.test(Sid1)$p.value
  )

# Print results
print(shapiro_results)

# Reshape data to long format
data_long <- melt(data, id.vars = "Tissue", variable.name = "Gene", value.name = "Expression")

# Perform Wilcoxon rank-sum test for each tissue group
wilcoxon_results <- data_long %>% 
  group_by(Tissue) %>% 
  summarise(p_value = wilcox.test(Expression ~ Gene)$p.value)

# Print Wilcoxon test results
print(wilcoxon_results)

View(wilcoxon_results)

# Convert Tissue to a factor for proper ordering
data_long$Tissue <- factor(data_long$Tissue, levels = c("HF", "OV", "AG", "HP/DG", "T"))

# Create ggplot with bivariate geom points with error bars
ggplot(data_long, aes(x = Tissue, y = Expression, fill = Gene, color = Gene)) +
  # Add horizontal lines to show the mean
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.4), 
               width = 1.5) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.4), 
               width = 0.9, size = 0.5) +  # Error bars showing standard error
      geom_jitter(aes(color = Gene), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.4), 
              size = 3.0, alpha = 0.7) +
  stat_compare_means(aes(x = Tissue, y = Expression, group = Gene), 
                     method = "wilcox.test", label = "p.signif", 
                     size = 6, inherit.aes = FALSE) +
  scale_fill_manual(values = c("Chc1" = "black", "Sid1" = "#66b8fd")) +  # Custom colors for genes
  scale_color_manual(values = c("Chc1" = "black", "Sid1" = "#66b8fd")) +  # Ensuring consistent color mapping
  labs(title = "Chc1 vs Sid1L", x = "Tissue", 
       y = "Relative expression") +
  theme_test() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),  # Title size & bold
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20, face = "bold"),
        panel.border = element_rect(size = 2.0, fill = NA))

#NOTE: Ignore the warning message or reduce bar width

# Save the final figure
ggsave("TS_Chc1_vs_SID1 mRNA expression levels2.png", width = 7.25, height = 6.5, dpi = 300)
