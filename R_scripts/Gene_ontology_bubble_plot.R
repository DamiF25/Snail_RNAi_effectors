# Install required packages, if required
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidytext")
install.packages("forcats")

# Load libraries
library(ggplot2)
library(dplyr)
library(tidytext)  # Needed for reorder_within()
library(forcats)   # For better factor handling

# Set working directory

# Load data
data <- read.csv("GO_data.csv")

# Extract GO enrichment results with appropriate columns
go_cc <- data.frame(ID = c("GO:0016442", "GO:0000932", "GO:0036464", "GO:1902555",
                           "GO:1990904", "GO:0030136", "GO:0032991", "GO:0070578",
                           "GO:0099080", "GO:0030139"), 
      Description = c("RISC complex", "P-body", "Cytoplasmic ribonucleoprotein granule",
                      "Endoribonuclease complex", "Ribonucleoprotein complex",
                      "Clathrin-coated vesicle", "Protein-containing complex",
                      "RISC-loading complex", "Supramolecular complex", "Endocytic vesicle"),
                    FDR = c(0.000000000099, 0.000000000851, 0.000000000261, 0.0000703, 
                            0.00000981, 0.00058, 0.0000000000213, 0.0019, 0.0000362, 0.0011),
                    Count = c(6, 8, 10, 4, 12, 4, 38, 2, 12, 4),
                    Signal = c(3.14, 2.53, 2.42, 1.26, 0.98, 0.95, 0.93, 0.87, 0.86, 0.86),
                    Category = "CC")

go_bp <- data.frame(ID = c("GO:0031047", "GO:0010629", "GO:0035196", "GO:0031054", "GO:0006401",
                           "GO:0006402", "GO:0016441", "GO:0035195", "GO:0017148", "GO:0034655"), 
      Description = c("Gene silencing by RNA", "Negative regulation of gene expression",
                      "miRNA processing", "pre-miRNA processing", "RNA catabolic process",
                      "mRNA catabolic process", "Post-transcriptional gene silencing", 
                      "miRNA-mediated gene silencing", "Negative regulation of translation",
                      "Nucleobase-containing compound catabolic process"),
                    FDR = c(0.0000000000000000983, 0.0000000000000000509, 0.000000000364, 
                            0.0000000695, 0.00000000254, 0.000000012, 0.000000256, 0.000000688, 
                            0.0000000665, 0.0000000265),
                    Count = c(13, 20, 7, 5, 10, 9, 5, 4, 8, 11),
                    Signal = c(4.29, 3.06, 2.88, 2.24, 2.18, 2.09, 2.03, 1.96, 1.93, 1.75),
                    Category = "BP")

go_mf <- data.frame(ID = c("GO:0003725", "GO:0003723", "GO:0004540", "GO:0004518", "GO:0140098",
                           "GO:0016896", "GO:0070883", "GO:0004534", "GO:0004527", "GO:0042626"), 
      Description = c("Double-stranded RNA binding", "RNA binding", "Ribonuclease activity",
                      "Nuclease activity", "Catalytic activity, acting on RNA",
                      "Exoribonuclease activity, producing 5-phosphomonoesters",
                      "pre-miRNA binding", "5-3 exoribonuclease activity", "Exonuclease activity",
                      "ATPase-coupled transmembrane transporter activity"),
                    FDR = c(0.00000000000355, 0.0000000000000000856, 0.000000131, 0.0000000295,
                            0.00000000638, 0.0000132, 0.000039, 0.0000574, 0.0000148, 0.00000265),
                    Count = c(9, 28, 8, 10, 14, 5, 3, 3, 6, 9),
                    Signal = c(3.29, 2.01, 1.82, 1.81, 1.6, 1.44, 1.41, 1.35, 1.33, 1.32),
                    Category = "MF")


# Combine all results
go_data <- bind_rows(go_cc, go_bp, go_mf)

# Verify data after merging
str(go_data)


# If not, convert FDR column to numeric in all data frames
go_cc$FDR <- as.numeric(go_cc$FDR)
go_bp$FDR <- as.numeric(go_bp$FDR)
go_mf$FDR <- as.numeric(go_mf$FDR)


# Convert FDR values to -log10 scale for better visualization
go_data$logFDR <- -log10(go_data$FDR)

# Convert 'Category' to a factor before plotting
go_data$Category <- factor(go_data$Category, levels = c("CC", "BP", "MF"))

# Bubble plot with correct FDR legend
ggplot(go_data, aes(x = Signal, y = reorder(Description, as.numeric(Category)),
                    size = Count, color = logFDR, shape = Category)) +
  geom_point(alpha = NA) +
  scale_size_continuous(range = c(3, 20), breaks = c(5, 10, 20, 30), name = "Gene Count") +  # Adjust bubble size + Custom breaks
  scale_color_gradientn(colors = c("#2ca8f4", "black", "orange"), 
                        name = "-log10(FDR)") +  # Better color scaling
  scale_shape_manual(values = c(16, 17, 15),  # Different shapes for CC, BP, MF
                     guide = guide_legend(override.aes = list(size = 6))) +  # Increase legend shape size
  labs(x = "Signal", 
       y = "GO Terms", 
       size = "Gene Count", 
       color = "logFDR", 
       shape = "GO Category") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 24), 
        axis.text.x = element_text(size = 24, face = "bold"),
        axis.title = element_text(size = 26),
        legend.text = element_text(size = 22), 
        legend.title = element_text(size = 24),
        legend.key.height = unit(1, "cm"),  # Adjust legend height
        legend.spacing.y = unit(1, "cm"))  # Adjust vertical spacing in legend

# Save the final figure
ggsave("GO_enrich_bubble_plot.png", width = 18, height = 16, dpi = 300)