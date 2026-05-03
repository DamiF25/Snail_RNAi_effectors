# Install necessary packages

install.packages("ggplot2")
install.packages("dplyr")
install.packages("ggpubr")
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
data <- read.csv("Chc_vs_SidL.csv")

# Perform Shapiro-Wilk test for each group
shapiro_results <- data %>%
  group_by(Tissue) %>%
  summarise(
    shapiro_Chc1_p = shapiro.test(Chc1)$p.value,
    shapiro_Sid1_p = shapiro.test(Sid1L)$p.value
  )

# Print results
print(shapiro_results)

# Reshape data to long format
data_long <- melt(data, id.vars = "Tissue", variable.name = "Gene", value.name = "Expression")

# Define y-axis expansion height
y_axis_upper <- max(data_long$Expression)

# Perform Wilcoxon rank-sum test for each tissue group
wilcoxon_results <- data_long %>% 
  group_by(Tissue) %>% 
  summarise(p_value = wilcox.test(Expression ~ Gene)$p.value)

# Print Wilcoxon test results
print(wilcoxon_results)

View(wilcoxon_results)

# Convert Tissue to a factor for proper ordering
data_long$Tissue <- factor(data_long$Tissue, levels = c("AG", "HF", "HP/DG", "OV", "T"))

# Define y-axis expansion height
y_axis_upper <- max(data_long$Expression)


# Create ggplot with bivariate geom points with error bars
p <- ggplot(data_long, aes(x = Tissue, y = Expression, fill = Gene, color = Gene)) +
  # Add horizontal lines to show the mean
  stat_summary(fun = mean, geom = "bar", position = position_dodge(width = 0.4), 
               width = 1.2) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.4), 
               width = 0.9, size = 0.5) +  # Error bars showing standard error
      geom_jitter(aes(color = Gene), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.4), 
              size = 3.0, alpha = 0.7) +
  stat_compare_means(aes(x = Tissue, y = Expression, group = Gene), 
                     method = "wilcox.test", label = "p.signif",
                     symnum.args = list(
                       cutpoints = c(0, 0.001, 0.01, 0.05),
                       symbols = c("***", "**", "*")),
                     size = 10, inherit.aes = FALSE,
                     label.y = y_axis_upper * 1.17) + # Asterisk position
  geom_segment(
    aes(x = as.numeric(Tissue) - 0.2,
        xend = as.numeric(Tissue) + 0.2,
        y = y_axis_upper * 1.17,       # Significance bar position
        yend = y_axis_upper * 1.17),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  # Left vertical tip
  geom_segment(
    aes(
      x = as.numeric(Tissue) - 0.2,
      xend = as.numeric(Tissue) - 0.2,
      y = y_axis_upper * 1.17,
      yend = y_axis_upper * 1.15    # tip drops down by 0.02
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  
  # Right vertical tip
  geom_segment(
    aes(
      x = as.numeric(Tissue) + 0.2,
      xend = as.numeric(Tissue) + 0.2,
      y = y_axis_upper * 1.17,
      yend = y_axis_upper * 1.1    # tip drops down by 0.07
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
# Manual color by gene
  scale_fill_manual(values = c("Chc1" = "black", "Sid1L" = "#66b8fd")) +  # Custom colors for genes
  scale_color_manual(values = c("Chc1" = "black", "Sid1L" = "#66b8fd")) +  # Ensuring consistent color mapping
  scale_y_continuous(
    expand = c(0, 0),          
    limits = c(0, y_axis_upper * 1.3)   
  ) +
  labs(title = "Chc1 vs Sid1L", x = "Tissue", 
       y = "Relative expression") +
  theme_test() +
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),  # Title size & bold
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 22, color = "black"),
        legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20, face = "bold"),
        legend.key.height = unit(0.7, "cm"),
        legend.key.spacing.y = unit(0.2, "cm"),
        panel.border = element_rect(linewidth = 1.0, fill = NA))

# Print plot
print(p)

#NOTE: Ignore the warning message or reduce bar width

# Save figure
ggsave(filename = "TS_Chc1_vs_SID1_expression_levels.tif",
       plot = p,
       width = 7,
       height = 6,
       dpi = 1000,
       units = "in",
       compression = "lzw")
