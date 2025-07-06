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
data <- read.csv("TS_Sid.csv")

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

# Convert Tissue to a factor for proper ordering
data$Tissue <- factor(data$Tissue, levels = c("HF", "OV", "AG", "HP/DG", "T"))

# Define colors
size_colors <- c("HF" = "#00c4b3", "OV" = "#9b9b9b", "AG" = "#e33b3b", "HP/DG" = "#66b8fd",
                 "T" = "#d49d00")

# Extend y-axis upper limit just a bit higher (if there is any sig. diff. in the comparisons)
y_axis_upper <- max(data$Delta_Ct) * 1.15

# Create a box plot with significance annotations
ggplot(data, aes(x = Tissue, y = Delta_Ct)) +
  geom_boxplot(fill = size_colors) +
  geom_jitter(color = "black", position = position_jitter(width = 0.3), size = 2.5) +
  geom_signif(comparisons = comparisons_list, 
              map_signif_level = FALSE, # FALSE shows actual p-values instead of asterisks
              textsize = 6,
              y_position = c(max(data$Delta_Ct) * 1.0,  
                             max(data$Delta_Ct) * 1.09),  # Adjust step heights
              annotations = c("p = 0.0399", "p = 0.0009")) +  # Manually set p-values
  labs(title = "Sid1L", x = "Tissue", y = "Relative expression") +
  scale_y_continuous(limits = c(NA, y_axis_upper)) +  # Expand y-axis manually
  theme_test() +
  theme(legend.position = "none") +  # Removes legend
  theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),  # Title size & bold
      axis.title = element_text(size = 24, face = "bold"),
      axis.text = element_text(size = 22, face = "bold"),
      panel.border = element_rect(size = 2.0, fill = NA))

# Save the final figure
ggsave("TS_Sid2.png", width = 7.25, height = 6.5, dpi = 300)
