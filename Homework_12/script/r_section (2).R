setwd("/home/group/bioe-591-genomics/students/osterhoudt/Homework_12")

library(easypackages)

libraries("pcadapt", "vcfR", "adegenet", 
          "ggplot2", "dplyr", "readr", "patchwork", "viridis")

#loading data
meta <- read_csv("metadata.csv")
x <- read.pcadapt("salvelinus.vcf", type="vcf")

#scree plotting; to choose K, where is elbow visible along x axis?
pca <- pcadapt(x, K=20)
plot(pca, option = "screeplot")
#elbow from Meet et al 2025 is K=4.
pca_4 <- pcadapt(x, K=4)
plot(pca_4, option = "screeplot")

#identify outliers using pvalues; need to correct for miltiple comparison
#to redice rate of false positive, Bonferroni correction (STATS 511)
alpha <- 0.05
padj <- p.adjust(pca_4$pvalues, method = "bonferroni")
outlier_idx <- which(padj < alpha)
neutral_idx <- which(padj >= alpha)
cat("Outlier loci:", length(outlier_idx), "\n")
#306 outlier loci
cat("Neutral loci:", length(neutral_idx), "\n")
#5223 neutral loci

#plotting Manhattan plot to see outlier distribution
plot(pca_4)

#explore relationship of outliers to inferred patterns of variation across subgroups

vcf <- read.vcfR("salvelinus.vcf", verbose = FALSE)
vcf_neutral <- vcf[neutral_idx, ]
vcf_outlier <- vcf[outlier_idx, ]

#neutral and outlier SNP need to conver to genind objects
dna_neutral <- vcfR2DNAbin(vcf_neutral, unphased_as_NA = FALSE, consensus = TRUE, extract.haps = FALSE)
#5223 variants remain
gi_neutral <- DNAbin2genind(dna_neutral)
dna_outlier <- vcfR2DNAbin(vcf_outlier, unphased_as_NA = FALSE, consensus = TRUE, extract.haps = FALSE)
#306 variants remain
gi_outlier <- DNAbin2genind(dna_outlier)

#scaling genotypes to prior running PCA with prcomp()
gi_neutral_scaled <- scaleGen(gi_neutral, NA.method = "mean", scale = FALSE)
pca_neutral <- prcomp(gi_neutral_scaled, center = FALSE, scale. = FALSE)
gi_outlier_scaled <- scaleGen(gi_outlier, NA.method = "mean", scale = FALSE)
pca_outlier <- prcomp(gi_outlier_scaled, center = FALSE, scale. = FALSE)

#transform PCA objects to dataframe for plotting
# turn lists into dataframes
neutral_df <- data.frame(pca_neutral$x[, 1:2]) %>%
  mutate(sample = rownames(.))
outlier_df <- data.frame(pca_outlier$x[, 1:2]) %>%
  mutate(sample = rownames(.))

# join metadata
neutral_df <- left_join(neutral_df, meta, by = c("sample" = "Sample.ID"))
outlier_df <- left_join(outlier_df, meta, by = c("sample" = "Sample.ID"))

# shared axis limits
xlims <- range(c(outlier_df$PC1, neutral_df$PC1), na.rm = TRUE)
ylims <- range(c(outlier_df$PC2, neutral_df$PC2), na.rm = TRUE)

# percent variance explained
neutral_var <- 100 * summary(pca_neutral)$importance[2, ]
outlier_var <- 100 * summary(pca_outlier)$importance[2, ]

#compare PCA plots between two subsets
p1 <- ggplot(neutral_df, aes(x = PC1, y = PC2, color=Site.Latitude)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis() +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(
    title = "PCA: neutral loci",
    x = paste0("PC1 (", round(neutral_var[1], 1), "%)"),
    y = paste0("PC2 (", round(neutral_var[2], 1), "%)")
  ) 

p2 <- ggplot(outlier_df, aes(x = PC1, y = PC2, color=Site.Latitude)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis() +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(
    title = "PCA: outlier loci only",
    x = paste0("PC1 (", round(outlier_var[1], 1), "%)"),
    y = paste0("PC2 (", round(outlier_var[2], 1), "%)")
  )

p1 + p2

#**Interpretation**
#The comparitive PCA plots for outlier and neutral loci plotted with 
#latitude shows clustering in both loci groups. The neutral loca PCA 
#is displaying a gradient from almost completely homozygous PC2 and then
#rapid declines to 0 around 2 on the PC1 axis.There is a relationnship
#with the site latitiude. As latitude decreases, the association with 
#PC2 decreases and PC1 increases. There is also about two clusters outside
#of this distinct gradient. There are some neutral loci found at latitudes 
#between 35 - 40 that cluster in the top right corner, suggesting potential
#for heterozygosity. The outlier loci has 3 clusters. The first impression is that
#individuals around 35-40 site latitude is somewhat heterozygotic based on placement
#within the PCA plot. Also, individuals alongside 45-50 site latitude
#also are heterozygous with placements around 0 at both axis. This is interesting


#plotting against longitude
p1_long <- ggplot(neutral_df, aes(x = PC1, y = PC2, color=Site.Longitude)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis() +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(
    title = "PCA: neutral loci",
    x = paste0("PC1 (", round(neutral_var[1], 1), "%)"),
    y = paste0("PC2 (", round(neutral_var[2], 1), "%)")
  ) 

p2_long <- ggplot(outlier_df, aes(x = PC1, y = PC2, color=Site.Longitude)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis() +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(
    title = "PCA: outlier loci only",
    x = paste0("PC1 (", round(outlier_var[1], 1), "%)"),
    y = paste0("PC2 (", round(outlier_var[2], 1), "%)")
  )

p1_long + p2_long

#**Interpretation**
#*My interpretation for my PCA plot with longitude shows a more distinguished
#*difference between neutral and outlier loci. The neutral loci also exhibits
#*the gradient. As longitude decreases from -60 to -90, there is a
#*fast decline in association with PCR. At -70 and -60 longitude
#*there are more homozygous neutral loci than those from -80. There are also
#*clusters in similar areas like in the latitude PCA plot. At site longitudes
#*between -80 and -70, there is a stronger heterozygous pattern than compared to the 
#*larger site longitudes. The outlier loci seems to only be found between -75
#* and -85 site longitude. Also, they show heterozygous clustering 
#* between PC1 and PC2.

