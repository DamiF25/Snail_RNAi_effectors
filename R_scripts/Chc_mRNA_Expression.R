# Install necessary packages

install.packages("ggplot2")
install.packages("ggpubr") # for comparison significance on graph
install.packages("dplyr") #for %>% function
install.packages("car")
install.packages("FSA") 
install.packages("rstatix")
install.packages("ggsignif")

# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(car)
library(FSA)
library(rstatix)
library(ggsignif)


#set working directory

# Load data
data <- read.csv("TS_Chc.csv")

# Convert Group to factor
data$Tissue <- as.factor(data$Tissue)

# Shapiro-Wilk test for normality of EM size per diet
shapiro_test <- shapiro.test(data$Delta_Ct)

# Check results (if p<0.05, data is NOT normally distributed)
print(shapiro_test)

# Perform a Levene's test (if p<0.05, variances are NOT equal)
leveneTest(Delta_Ct ~ Tissue, data = data)

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(Delta_Ct ~ Tissue, data = data)
print(kruskal_result)

# Perform Dunn’s post-hoc test with Bonferroni correction
dunn_results <- dunnTest(Delta_Ct ~ Tissue, data = data, method = "bonferroni")
print(dunn_results)

# Extract significant comparisons
signif_comparisons <- dunn_results$res %>%
  filter(P.adj < 0.05) %>%
  select(Comparison, P.adj)

# Convert significant comparisons into list format for plotting
comparisons_list <- strsplit(signif_comparisons$Comparison, " - ")

# Convert p-values to asterisks
p_to_stars <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

annotations_stars <- sapply(signif_comparisons$P.adj, p_to_stars)

# Define colours for borders and points
tissue_colors <- c(
  "HF" = "#00c4b3",
  "OV" = "#9b9b9b",
  "AG" = "#e33b3b",
  "HP/DG" = "#66b8fd",
  "T" = "#d49d00"
)

# Y-axis expansion
y_axis_upper <- max(data$Delta_Ct) * 1.1

# Base plot
p <- ggplot(data, aes(x = Tissue, y = Delta_Ct, color = Tissue, fill = Tissue)) +
  geom_boxplot(size = 1.2, alpha = 0.5) + 
  geom_jitter(position = position_jitter(width = 0.15), size = 3.5) + 
  scale_color_manual(values = tissue_colors) +
  scale_fill_manual(values = tissue_colors) + 
  scale_y_continuous(limits = c(NA, y_axis_upper)) +
  labs(title = "Chc1", x = "Tissue", y = "Relative expression") +
  theme_test() +
  theme(legend.position = "none",
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 24, face = "bold"),
        axis.text = element_text(size = 22, color = "black"),
        panel.border = element_rect(linewidth = 1.0, fill = NA, color = "black")  )

# Add significance only if there are significant comparisons
if (nrow(signif_comparisons) > 0) {
  p <- p + geom_signif(
    comparisons = comparisons_list,
    annotations = annotations_stars,
    textsize = 10,
    color = "black",
    tip_length = 0.02,
    y_position = seq(
      from = max(data$Delta_Ct) * 1.0,
      length.out = length(comparisons_list),
      by = max(data$Delta_Ct) * 0.08))
}

# Print plot
print(p)

# Save figure
ggsave(filename = "TS_Chc.tif",
       plot = p,
       width = 6,
       height = 6,
       dpi = 1000,
       units = "in",
       compression = "lzw")