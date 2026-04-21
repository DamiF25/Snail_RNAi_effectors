# Install necessary packages
install.packages("ggplot2")
install.packages("ggpubr") # for comparison significance on graph
install.packages("dplyr") 
install.packages("car")
install.packages("multcompView") #To extract Tukey's test result
install.packages("ggstatsplot")


# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(car)
library(multcompView)
library(ggstatsplot)

# Set working directory

# Load the data
data <- read.csv("Ago2_mRNA_expression.csv")

data$Size <- as.factor(data$Size)  # Ensure Size is a factor

# Shapiro-Wilk test for normality of eggs per snail
shapiro_test <- shapiro.test(data$Delta_Ct)

# Check results (if p<0.05, data is NOT normally distributed)
print(shapiro_test)

# Perform a Levene's test (if p<0.05, variances are NOT equal)
leveneTest(Delta_Ct ~ Size, data = data)


# If data were normally distributed, and variances are equal, we can use ANOVA
anova_result <- aov(Delta_Ct ~ Size, data = data)

anova_result

# Perform Tukey's post hoc test (all groups vs each other)
tukey_result <- TukeyHSD(anova_result, conf.level=.95)

# Print the results
print(tukey_result)

# Convert Tukey's results to a data frame
tukey_df <- as.data.frame(tukey_result$Size)

# Add group1 and group2 columns by splitting the row names
tukey_annotations <- data.frame(
  group1 = sapply(strsplit(rownames(tukey_df), "-"), `[`, 1),
  group2 = sapply(strsplit(rownames(tukey_df), "-"), `[`, 2),
  p.adj = tukey_df$`p adj`)

colnames(tukey_annotations)

# Filter significant comparisons
significant_comparisons <- tukey_annotations %>%
  filter(p.adj < 0.05)

print(significant_comparisons)

# Round p-values to 3 decimal places
significant_comparisons$p.adj <- round(significant_comparisons$p.adj, 3)

print(significant_comparisons)

# Order bars based on size group
data$Size <- factor(data$Size, levels = c("J1", "J2", "A1", "A2"))

# Create bar plot (with no significance yet)
p <- ggplot(data, aes(x = Size, y = Delta_Ct, fill = Size)) +
  stat_summary(fun = mean, geom = "bar", width = 0.4, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.2, color = "red") +
  geom_jitter(color = "red", width = 0.15, size = 2.5) +
  scale_fill_manual(values = c(                             
    "J1" = "#FFFFFF",
    "J2" = "grey",
    "A1" = "#494949",
    "A2" = "#000000"
  )) +
  labs(title = "Ago2", x = "Size group", y = "Relative expression") +
  theme_test() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(legend.position = "none",
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 22, color = "black"),
    panel.border = element_rect(linewidth = 1.0, fill = NA))

# Add significance only if comparisons exist
if (nrow(significant_comparisons) > 0) {
  
  comparisons_list <- apply(
    significant_comparisons[, c("group1", "group2")],
    1,
    function(x) as.character(unlist(x)))
  
  comparisons_list <- lapply(comparisons_list, identity)
  
  p <- p + geom_signif(
    comparisons = comparisons_list,
    annotations = significant_comparisons$p.adj,
    textsize = 5,
    tip_length = 0.02)}

# Print final plot
print(p)

# Save figure
ggsave(filename = "Ago2_mRNA_expression_levels.tif",
  plot = p,
  width = 5,
  height = 5,
  dpi = 1000,
  units = "in",
  compression = "lzw")
