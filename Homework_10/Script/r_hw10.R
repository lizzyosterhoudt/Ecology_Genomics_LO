setwd("~/evo_genomics/students/osterhoudt/Homework_10")
library(tidyverse)
library(adegenet)
library(vcfR)

#reading suite of files
# read sample names and extract
fam <- read_table("species.int.fam", 
                  col_names = FALSE,
                  show_col_types = FALSE
)
samples <- fam$X2   # individual IDs are usually column 2

# choose your K value (will be 2 for this demo)
K <- 2

# read Q matrix
q2 <- read_table("species.int.2.Q",
                col_names = FALSE,
                show_col_types = FALSE
)

# name the ancestry columns
colnames(q2) <- paste0("Cluster", 1:K)

# combine with sample names
q_df2 <- q2 %>%
  mutate(sample = samples) %>%
  relocate(sample)

# convert to long format for ggplot
q_long2 <- q_df2 %>%
  pivot_longer(
    cols = starts_with("Cluster"),
    names_to = "cluster",
    values_to = "ancestry"
  ) %>%
  mutate(sample = factor(sample, levels = samples))

#plotting
ggplot(q_long2, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color="white") +
  theme_bw() +
  labs(x = "Individual", y = "Ancestry proportion") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#extracting cross-validation values from log
cv_df2 <- tibble(file = "species.2.log") %>%
  mutate(
    text = map_chr(file, read_file),
    cv_line = str_extract(text, "CV error \\(K=\\d+\\):\\s*[-0-9.eE]+"),
    K  = str_match(cv_line, "CV error \\(K=(\\d+)\\)")[,2] |> as.integer(),
    CV = str_match(cv_line, ":\\s*([-0-9.eE]+)")[,2] |> as.numeric()
  ) %>%
  select(file, K, CV) %>%
  arrange(K)
cv_df2

#PCA and DAPC with adegenet
vcf <- read.vcfR(file = "squirrels.vcf", verbose = TRUE)

dna <- vcfR2DNAbin(vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
species_genind <- DNAbin2genind(dna)
species_genind

#make pca behave; treating missing data to making homo alternate equivalent
species_genind_scaled <- scaleGen(species_genind,NA.method="mean",scale=F)
species_pca <- prcomp(species_genind_scaled, center=F,scale=F)

#display load of each pc indicating whether majority of var can be on 1 axis
screeplot(species_pca)
pc2 <- data.frame(species_pca$x[,1:3])
pc2$sample <- rownames(pc)

ggplot(data=pc,aes(x=PC1,y=PC2))+
  geom_text(aes(label=sample))

#before applying DAPC, needs priori
grp2 <- find.clusters(species_genind, n.pca = 50, n.clust = 2)
grp2

#visualize clusters
pc2$cluster <- grp2$grp
ggplot(pc, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2) +
  geom_text(aes(label = sample), vjust = -0.5, size = 3) 

#DAPC determine which pc contribute and ancestry uncertainty
dapc2 <- dapc(species_genind, pop = grp2$grp, n.pca = 50, n.da = 2)
dapc2

#posterior prob of ancestry extraction
q <- as.data.frame(dapc2$posterior)
q$sample <- rownames(q)
q_long <- q |>
  pivot_longer(
    cols = -sample,
    names_to = "cluster",
    values_to = "ancestry"
  )
q_long


#plotting
ggplot(q_long, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color = "white") +
  theme_bw() +
  labs(x = "Individual", y = "Assignment probability") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


#compare results across diff k

# identify clusters
grp_auto2 <- find.clusters(
  species_genind,
  n.pca = 50,
  choose.n.clust = FALSE, 
  max.n.clust = 10,
  stat = "BIC"
)

# print BIC values
grp_auto2$Kstat


# plot!
plot(
  1:length(grp_auto$Kstat),
  grp_auto$Kstat,
  type = "b",
  xlab = "K",
  ylab = "BIC"
)

#to find k where is elbow?




#**Homework run**
  
fam <- read_table("species.int.fam", 
                  col_names = FALSE,
                  show_col_types = FALSE
)
samples <- fam$X2   # individual IDs are usually column 2

#kval
K <- 3

# read Q matrix
q3 <- read_table("species.int.3.Q",
                col_names = FALSE,
                show_col_types = FALSE
)

# name the ancestry columns
colnames(q3) <- paste0("Cluster", 1:K)

# combine with sample names
q_df3 <- q3 %>%
  mutate(sample = samples) %>%
  relocate(sample)

# convert to long format for ggplot
q_long3 <- q_df3 %>%
  pivot_longer(
    cols = starts_with("Cluster"),
    names_to = "cluster",
    values_to = "ancestry"
  ) %>%
  mutate(sample = factor(sample, levels = samples))

#plotting
ggplot(q_long, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color="white") +
  theme_bw() +
  labs(x = "Individual", y = "Ancestry proportion") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#extracting cross-validation values from log
cv_df3 <- tibble(file = "species.2.log") %>%
  mutate(
    text = map_chr(file, read_file),
    cv_line = str_extract(text, "CV error \\(K=\\d+\\):\\s*[-0-9.eE]+"),
    K  = str_match(cv_line, "CV error \\(K=(\\d+)\\)")[,2] |> as.integer(),
    CV = str_match(cv_line, ":\\s*([-0-9.eE]+)")[,2] |> as.numeric()
  ) %>%
  select(file, K, CV) %>%
  arrange(K)
cv_df3

#PCA and DAPC with adegenet
vcf <- read.vcfR(file = "squirrels.vcf", verbose = TRUE)

dna <- vcfR2DNAbin(vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
species_genind <- DNAbin2genind(dna)
species_genind

#make pca behave; treating missing data to making homo alternate equivalent
species_genind_scaled <- scaleGen(species_genind,NA.method="mean",scale=F)
species_pca <- prcomp(species_genind_scaled, center=F,scale=F)

#display load of each pc indicating whether majority of var can be on 1 axis
screeplot(species_pca)
pc3 <- data.frame(species_pca$x[,1:3])
pc3$sample <- rownames(pc)

ggplot(data=pc,aes(x=PC1,y=PC2))+
  geom_text(aes(label=sample))

#before applying DAPC, needs priori
grp3 <- find.clusters(species_genind, n.pca = 50, n.clust = 3)
grp3

#visualize clusters
pc3$cluster <- grp3$grp
ggplot(pc, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2) +
  geom_text(aes(label = sample), vjust = -0.5, size = 3) 

#DAPC determine which pc contribute and ancestry uncertainty
dapc3 <- dapc(species_genind, pop = grp3$grp, n.pca = 50, n.da = 3)
dapc3

#posterior prob of ancestry extraction
q <- as.data.frame(dapc1$posterior)
q$sample <- rownames(q)
q_long <- q |>
  pivot_longer(
    cols = -sample,
    names_to = "cluster",
    values_to = "ancestry"
  )
q_long


#plotting
ggplot(q_long, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color = "white") +
  theme_bw() +
  labs(x = "Individual", y = "Assignment probability") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


#compare results across diff k

# identify clusters
grp_auto <- find.clusters(
  species_genind,
  n.pca = 50,
  choose.n.clust = FALSE, 
  max.n.clust = 10,
  stat = "BIC"
)

# print BIC values
grp_auto$Kstat


# plot!
plot(
  1:length(grp_auto$Kstat),
  grp_auto$Kstat,
  type = "b",
  xlab = "K",
  ylab = "BIC"
)

#to find k where is elbow?

#**K = 4**
fam <- read_table("species.int.fam", 
                  col_names = FALSE,
                  show_col_types = FALSE
)
samples <- fam$X2   # individual IDs are usually column 2

# choose your K value (will be 2 for this demo)
K <- 4

# read Q matrix
q <- read_table("species.int.4.Q",
                col_names = FALSE,
                show_col_types = FALSE
)

# name the ancestry columns
colnames(q) <- paste0("Cluster", 1:K)

# combine with sample names
q_df <- q %>%
  mutate(sample = samples) %>%
  relocate(sample)

# convert to long format for ggplot
q_long <- q_df %>%
  pivot_longer(
    cols = starts_with("Cluster"),
    names_to = "cluster",
    values_to = "ancestry"
  ) %>%
  mutate(sample = factor(sample, levels = samples))

#plotting
ggplot(q_long, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color="white") +
  theme_bw() +
  labs(x = "Individual", y = "Ancestry proportion") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#extracting cross-validation values from log
cv_df <- tibble(file = "species.2.log") %>%
  mutate(
    text = map_chr(file, read_file),
    cv_line = str_extract(text, "CV error \\(K=\\d+\\):\\s*[-0-9.eE]+"),
    K  = str_match(cv_line, "CV error \\(K=(\\d+)\\)")[,2] |> as.integer(),
    CV = str_match(cv_line, ":\\s*([-0-9.eE]+)")[,2] |> as.numeric()
  ) %>%
  select(file, K, CV) %>%
  arrange(K)
cv_df

#PCA and DAPC with adegenet
vcf <- read.vcfR(file = "squirrels.vcf", verbose = TRUE)

dna <- vcfR2DNAbin(vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
species_genind <- DNAbin2genind(dna)
species_genind

#make pca behave; treating missing data to making homo alternate equivalent
species_genind_scaled <- scaleGen(species_genind,NA.method="mean",scale=F)
species_pca <- prcomp(species_genind_scaled, center=F,scale=F)

#display load of each pc indicating whether majority of var can be on 1 axis
screeplot(species_pca)
pc <- data.frame(species_pca$x[,1:4])
pc$sample <- rownames(pc)

ggplot(data=pc,aes(x=PC1,y=PC2))+
  geom_text(aes(label=sample))

#before applying DAPC, needs priori
grp <- find.clusters(species_genind, n.pca = 50, n.clust = 4)
grp

#visualize clusters
pc$cluster <- grp$grp
ggplot(pc, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2) +
  geom_text(aes(label = sample), vjust = -0.5, size = 3) 

#DAPC determine which pc contribute and ancestry uncertainty
dapc1 <- dapc(species_genind, pop = grp$grp, n.pca = 50, n.da = 4)
dapc1

#posterior prob of ancestry extraction
q <- as.data.frame(dapc1$posterior)
q$sample <- rownames(q)
q_long <- q |>
  pivot_longer(
    cols = -sample,
    names_to = "cluster",
    values_to = "ancestry"
  )
q_long


#plotting
ggplot(q_long, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color = "white") +
  theme_bw() +
  labs(x = "Individual", y = "Assignment probability") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


#compare results across diff k

# identify clusters
grp_auto <- find.clusters(
  species_genind,
  n.pca = 50,
  choose.n.clust = FALSE, 
  max.n.clust = 10,
  stat = "BIC"
)

# print BIC values
grp_auto$Kstat


# plot!
plot(
  1:length(grp_auto$Kstat),
  grp_auto$Kstat,
  type = "b",
  xlab = "K",
  ylab = "BIC"
)




#**k=5**

fam <- read_table("species.int.fam", 
                  col_names = FALSE,
                  show_col_types = FALSE
)
samples <- fam$X2 


K <- 5

# read Q matrix
q <- read_table("species.int.5.Q",
                col_names = FALSE,
                show_col_types = FALSE
)

# name the ancestry columns
colnames(q) <- paste0("Cluster", 1:K)

# combine with sample names
q_df <- q %>%
  mutate(sample = samples) %>%
  relocate(sample)

# convert to long format for ggplot
q_long <- q_df %>%
  pivot_longer(
    cols = starts_with("Cluster"),
    names_to = "cluster",
    values_to = "ancestry"
  ) %>%
  mutate(sample = factor(sample, levels = samples))

#plotting
ggplot(q_long, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color="white") +
  theme_bw() +
  labs(x = "Individual", y = "Ancestry proportion") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


#extracting cross-validation values from log
cv_df <- tibble(file = "species.2.log") %>%
  mutate(
    text = map_chr(file, read_file),
    cv_line = str_extract(text, "CV error \\(K=\\d+\\):\\s*[-0-9.eE]+"),
    K  = str_match(cv_line, "CV error \\(K=(\\d+)\\)")[,2] |> as.integer(),
    CV = str_match(cv_line, ":\\s*([-0-9.eE]+)")[,2] |> as.numeric()
  ) %>%
  select(file, K, CV) %>%
  arrange(K)
cv_df

#PCA and DAPC with adegenet
vcf <- read.vcfR(file = "squirrels.vcf", verbose = TRUE)

dna <- vcfR2DNAbin(vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
species_genind <- DNAbin2genind(dna)
species_genind

#make pca behave; treating missing data to making homo alternate equivalent
species_genind_scaled <- scaleGen(species_genind,NA.method="mean",scale=F)
species_pca <- prcomp(species_genind_scaled, center=F,scale=F)

#display load of each pc indicating whether majority of var can be on 1 axis
screeplot(species_pca)
pc <- data.frame(species_pca$x[,1:5])
pc$sample <- rownames(pc)

ggplot(data=pc,aes(x=PC1,y=PC2))+
  geom_text(aes(label=sample))

#before applying DAPC, needs priori
grp <- find.clusters(species_genind, n.pca = 50, n.clust = 5)
grp

#visualize clusters
pc$cluster <- grp$grp
ggplot(pc, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2) +
  geom_text(aes(label = sample), vjust = -0.5, size = 3) 

#DAPC determine which pc contribute and ancestry uncertainty
dapc1 <- dapc(species_genind, pop = grp$grp, n.pca = 50, n.da = 5)
dapc1

#posterior prob of ancestry extraction
q <- as.data.frame(dapc1$posterior)
q$sample <- rownames(q)
q_long <- q |>
  pivot_longer(
    cols = -sample,
    names_to = "cluster",
    values_to = "ancestry"
  )
q_long


#plotting
ggplot(q_long, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color = "white") +
  theme_bw() +
  labs(x = "Individual", y = "Assignment probability") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


#compare results across diff k

# identify clusters
grp_auto <- find.clusters(
  species_genind,
  n.pca = 50,
  choose.n.clust = FALSE, 
  max.n.clust = 10,
  stat = "BIC"
)

# print BIC values
grp_auto$Kstat


# plot!
plot(
  1:length(grp_auto$Kstat),
  grp_auto$Kstat,
  type = "b",
  xlab = "K",
  ylab = "BIC"
)




#**Comparing BIC across diff k-values**


library(mosaic)
favstats(grp_auto$Kstat)


library(ggplot2)

#make data plotable
bic_data <- data.frame(
  K = factor(1:10),  
  BIC = grp_auto$Kstat,
  cluster_assignment = 1:10 == which.min(grp_auto$Kstat)
)

ggplot(bic_data, aes(x = K, y = BIC, fill = cluster_assignment)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c("steelblue", "red"), 
                    labels = c(" ", paste("Best fit =", which.min(grp_auto$Kstat)))) +
  labs(title = "Comparison of BIC values across cluster assignments \n in northern and southern Idaho ground squirrels",
       x = "Number of Clusters (K)",
       y = "BIC Value") +
  theme_minimal() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = round(BIC, 2)), vjust = -0.5, size = 3.5)
