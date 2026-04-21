# Install necessary packages
install.packages("ggplot2")
install.packages("ggpubr") # for comparison significance on graph
install.packages("dplyr") 
install.packages("car")
install.packages("multcompView") #To extract Tukey's test result
install.packages("ggstatsplot")
install.packages("FSA") 

# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(car)
library(multcompView)
library(ggstatsplot)
library(FSA)

# Set working directory

# Load the data
data <- read.csv("Eri1_mRNA_expression.csv")

data$Size <- as.factor(data$Size)  # Ensure Size is a factor

# Shapiro-Wilk test for normality of eggs per snail
shapiro_test <- shapiro.test(data$Delta_Ct)

# Check results (if p<0.05, data is NOT normally distributed)
print(shapiro_test)

# Perform a Levene's test (if p<0.05, variances are NOT equal)
leveneTest(Delta_Ct ~ Size, data = data)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(Delta_Ct ~ Size, data = data)
print(kruskal_result)

# Perform Dunn’s post-hoc test with Bonferroni correction
dunn_results <- dunnTest(Delta_Ct ~ Size, data = data, method = "bonferroni")
print(dunn_results)

# Extract significant comparisons
signif_comparisons <- dunn_results$res %>%
  filter(P.adj < 0.05) %>%
  select(Comparison, P.adj)

# Convert significant comparisons into list format for plotting
comparisons_list <- strsplit(signif_comparisons$Comparison, " - ")

# Convert Size to a factor for proper ordering
data$Size <- factor(data$Size, levels = c("J1", "J2", "A1", "A2"))

#Create bar plot with significance annotation 
p <- ggplot(data, aes(x = Size, y = Delta_Ct, fill = Size)) +
  stat_summary(fun = mean, geom = "bar", width = 0.4, position = position_dodge(), 
               color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, 
               position = position_dodge(0.9), color = "red") + # Red error bars
  geom_signif(comparisons = comparisons_list, 
              map_signif_level = TRUE, 
              textsize = 5) +
  geom_jitter(color = "red", position = position_jitter(width = 0.15), size = 2.5) +
  scale_fill_manual(values = c(                             
    "J1" = "#FFFFFF",
    "J2" = "grey",
    "A1" = "#494949",
    "A2" = "#000000"
  )) +
  labs(title = "Eri1", x = "Size group", y = "Relative expression") +
  theme_test() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Ensures zero intercepts x-axis
  theme(legend.position = "none") + # Removes legend
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),  # Title size & bold
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 22, color = "black"),
        axis.line = element_line(linewidth = NULL),         # Thicker axis lines
        panel.border = element_rect(size = 1.0, fill = NA)) # Thicker panel border

# Print plot
print(p)

# Save figure
ggsave(filename = "Eri1 mRNA expression levels.tif",
       plot = p,
       width = 5,
       height = 5,
       dpi = 1000,
       units = "in",
       compression = "lzw")
