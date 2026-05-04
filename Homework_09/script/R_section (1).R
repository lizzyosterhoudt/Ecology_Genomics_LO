
#loading necessary libraries 
library(adegenet)
library(vcfR)
library(pegas)
library(tidyr)
library(reshape2)
library(ggplot2)
library(related)
library(tidyr)
library(mosaic)

#uploading dataset
vcf <- read.vcfR(file= "Koala_MaxMissing10.recode.vcf", verbose = TRUE)

#viewing vcf file
View(vcf)
vcf
vcf@meta
vcf@fix
vcf@gt

#converting data to fit formats for adegenet and pegas packages
genind_obj <- vcfR2genind(vcf)
View(genind_obj)
summary(genind_obj@loc.n.all)
summary(genind_obj@loc.fac)
summary(genind_obj@type)

#calculating He and Ho for locus (SNP site) 
af_summary <- adegenet::summary(genind_obj)
af_summary
H_o <- af_summary$Hobs
H_e <- af_summary$Hexp

#taking first look at He and Ho and creating dataframe 
head(H_o)
head(H_e)
het_df <- data.frame(locus = names(H_o), H_o = H_o, H_e = H_e)
head(het_df)


#How to interpret locus name: 
#A:B:(+/-)
#A = scaffold, B = position, +/- = forward or reverse strand 

#calculating Fis, removing the NA values can calc mean Fis for all loci
Fis_per_locus <- 1 - (H_o / H_e)
Fis_per_locus
mean(Fis_per_locus, na.rm = TRUE)
#mean of Fis per locus is -0.0154, how do i interpret this?

#checking to see if loci are in HWE
loci_obj <- genind2loci(genind_obj)
hwe_results <- pegas::hw.test(loci_obj, B = 100)
hwe_results

#isolating HWE deviations by p-value 
library(tibble)
head(hwe_results)
class(hwe_results)
#hwe_results is in a matrix, changing to a dataframe 
hwe_results <- as.data.frame(hwe_results)
write.csv(hwe_results, "hwe_results.csv")
hwe_results <- read.csv("hwe_results.csv")
#this line of code was not working, repeatedly got object Pr.exact not found no matter how i used tibble
#hwe_results %>% as.tibble() %>% filter(Pr.exact<0.05)

#manually filtering by pvalue for hwe
hwe_results[hwe_results[, "Pr.exact"] < 0.05, ]

#transforming data to fit related package
gt_filtered <- vcfR::extract.gt(vcf, element = "GT")

#
# sample ids
sample_ids <- colnames(gt_filtered)

gt_to_alleles <- function(gt_vector) {
  # split "0/1" or "0|1" into two integer alleles, returning a 2-column matrix (samples x 2 alleles)
  allele1 <- integer(length(gt_vector))
  allele2 <- integer(length(gt_vector))
  
  for (i in seq_along(gt_vector)) {
    g <- gt_vector[i]
    if (is.na(g) || g %in% c("./.", ".", "./", "/.")) {
      allele1[i] <- 0
      allele2[i] <- 0
    } else {
      parts <- as.integer(strsplit(g, "[/|]")[[1]])
      allele1[i] <- parts[1] + 1L    # shift: 0->1 (ref), 1->2 (alt)
      allele2[i] <- parts[2] + 1L
    }
  }
  cbind(allele1, allele2)
}

allele_list <- vector("list", nrow(gt_filtered))

for (v in seq_len(nrow(gt_filtered))) {
  allele_list[[v]] <- gt_to_alleles(gt_filtered[v, ])
}

# combine: each element is (n_samples x 2); bind column-wise
allele_matrix <- do.call(cbind, allele_list)

# add individual IDs as the first column
coancestry_input <- data.frame(IndID = sample_ids, allele_matrix,
                               stringsAsFactors = FALSE)

# column names: IndID, L1_a, L1_b, L2_a, L2_b, ...
locus_names <- paste0(rep(paste0("L", seq_len(nrow(gt_filtered))),
                          each = 2),
                      rep(c("_a", "_b"), nrow(gt_filtered)))
colnames(coancestry_input) <- c("IndID", locus_names)



#slicing down to manageble
coancestry_input[1:5, 1:7]


#identifyinh coancestry

kin_results <- related::coancestry(
  genotype.data = coancestry_input,
  wang          = 1,      # 1 = compute; 0 = skip
  dyadml        = 1,
  quellergt     = 1)




# View overall structure
str(kin_results)

# Summary of what's in the results
names(kin_results)

#saving data
kin_results <- readRDS("kin_results.RDS")
relatedness_df <- kin_results$relatedness

#taking vars
favstats(wang~group, data=relatedness_df)
favstats(dyadml~group, data=relatedness_df)
favstats(quellergt~group, data=relatedness_df)


# group stats 
group_stats <- relatedness_df %>%
  group_by(group) %>%
  summarise(
    n_pairs = n(),
    mean_wang = mean(wang),
    sd_wang = sd(wang),
    mean_dyadml = mean(dyadml),
    mean_quellergt = mean(quellergt),
    n_positive_wang = sum(wang > 0),
    n_close_r = sum(wang > 0.25)
  ) %>%
  arrange(desc(mean_wang))

print(group_stats)



#whats diff between diff estimators
consistency_data <- relatedness_df %>%
  rowwise() %>%
  mutate(
    mean_estimate = mean(c(wang, dyadml, quellergt), na.rm = TRUE),
    sd_estimate = sd(c(wang, dyadml, quellergt), na.rm = TRUE),
    cv_estimate = sd_estimate / abs(mean_estimate) * 100  # Coefficient of variation
  ) %>%
  ungroup()

# Summary stats for each estimator
summary_stats <- relatedness_df %>%
  summarise(
    Estimator = c("Wang", "Dyadic ML", "Queller & Goodnight"),
    Mean = c(mean(wang, na.rm = TRUE), mean(dyadml, na.rm = TRUE), mean(quellergt, na.rm = TRUE)),
    Median = c(median(wang, na.rm = TRUE), median(dyadml, na.rm = TRUE), median(quellergt, na.rm = TRUE)),
    SD = c(sd(wang, na.rm = TRUE), sd(dyadml, na.rm = TRUE), sd(quellergt, na.rm = TRUE)),
    Min = c(min(wang, na.rm = TRUE), min(dyadml, na.rm = TRUE), min(quellergt, na.rm = TRUE)),
    Max = c(max(wang, na.rm = TRUE), max(dyadml, na.rm = TRUE), max(quellergt, na.rm = TRUE)),
    Prop_Negative = c(mean(wang < 0, na.rm = TRUE), mean(dyadml < 0, na.rm = TRUE), mean(quellergt < 0, na.rm = TRUE)),
    Prop_Positive = c(mean(wang > 0, na.rm = TRUE), mean(dyadml > 0, na.rm = TRUE), mean(quellergt > 0, na.rm = TRUE))
  )

print(summary_stats)


#wang relatedness
relatedness_df <- relatedness_df %>%
  mutate(relationship = case_when(
    wang < 0 ~ "Unrelated",
    wang >= 0 & wang < 0.125 ~ "Distant",
    wang >= 0.125 & wang < 0.25 ~ "Half-sibling",
    wang >= 0.25 & wang < 0.375 ~ "3/4-sibling",
    wang >= 0.375 ~ "Full-sibling"
  ))

#mean diffs of relatedness across coeffecients 
comparison_data <- relatedness_df %>%
  group_by(relationship) %>%
  summarise(
    Wang = mean(wang, na.rm = TRUE),
    DyadicML = mean(dyadml, na.rm = TRUE),
    QuellerGT = mean(quellergt, na.rm = TRUE),
    n_pairs = n(),
    .groups = 'drop'
  ) %>%
  mutate(relationship = factor(relationship, 
                               levels = c("Unrelated", "Distant", "Half-sibling", 
                                          "3/4-sibling", "Full-sibling")))

#plotting
comparison_long <- comparison_data %>%
  pivot_longer(cols = c(Wang, DyadicML, QuellerGT),
               names_to = "estimator",
               values_to = "mean_value")

]

#distribution
hist(group_stats$mean_wang, breaks = 20, 
     main = "Distribution of Mean Relatedness (Wang)",
     xlab = "Relatedness (r)", col = "skyblue")
abline(v = 0, col = "red", lwd = 2)
abline(v = 0.25, col = "green", lwd = 2, lty = 2)


#anyone most related??
top_pairs <- relatedness_df %>%
  arrange(desc(wang)) %>%
  select(pair.no, ind1.id, ind2.id, group, wang, dyadml, quellergt) %>%
  head(10)

print(top_pairs)

kin_wang <- relatedness_df %>%
  mutate(relationship = case_when(
    wang < 0 ~ "Unrelated",
    wang >= 0 & wang < 0.125 ~ "Distant/Unclear",
    wang >= 0.125 & wang < 0.25 ~ "Half-sibling range",
    wang >= 0.25 & wang < 0.375 ~ "Half-sibling to 3/4-sibling",
    wang >= 0.375 & wang < 0.5 ~ "Full-sibling range",
    wang >= 0.5 ~ "Duplicate/Inbred"
  )) %>%
  count(relationship)
head(kin_wang)



# dyadml relatedness hist
ggplot(relatedness_df, aes(x = dyadml)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, 
                 fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Dyadml Estimates",
       x = "Relatedness (r)", y = "Density") +
  theme_minimal()

#comp between wang and dyadml
plot(relatedness_df$wang, relatedness_df$dyadml,
     xlab = "Wang estimator", ylab = "DyadML",
     main = "Wang vs DyadML Comparison")
abline(0, 1, col = "red", lwd = 2)
abline(h = 0, v = 0, lty = 2)

# working w dyadml
kin_dyadml <- relatedness_df %>%
  mutate(dyadml_relationship = case_when(
    dyadml < 0 ~ "Unrelated/Negative (sampling error)",
    dyadml >= 0 & dyadml < 0.125 ~ "Unrelated to Distant",
    dyadml >= 0.125 & dyadml < 0.25 ~ "Distant Relatives (e.g., cousins)",
    dyadml >= 0.25 & dyadml < 0.375 ~ "Half-sibling range",
    dyadml >= 0.375 & dyadml < 0.5 ~ "Full-sibling range",
    dyadml >= 0.5 ~ "Very close (twins/duplicates?)"
  ))


table(kin_dyadml$dyadml_relationship)


table(kin_wang > 0.25, kin_dyadml > 0.25,
      dnn = c("Wang > 0.25", "DyadML > 0.25"))


  top_pairs <- kin_combined %>%
    arrange(desc(dyadml)) %>%
    head(10)
  
  ggplot(kin_dyadml, aes(x = reorder(pair.no, dyadml), y = dyadml)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.25, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(title = "Top 10 DyadML Estimates with 95% CI",
         x = "Pair Number", y = "DyadML Relatedness (r)") +
    theme_minimal()
}

print(top_pairs)

hist(relatedness_df$dyadml, breaks = 20, 
     main = "Distribution of Relatedness (dyadml)",
     xlab = "Relatedness (r)", col = "skyblue")


hist(relatedness_df$wang, breaks = 20, 
     main = "Distribution of Relatedness (dyadml)",
     xlab = "Relatedness (r)", col = "pink")


#group summary
group_summary <- relatedness_df %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_wang = mean(wang),
    sd_wang = sd(wang),
    se_wang = sd_wang / sqrt(n)
  ) %>%
  arrange(desc(mean_wang))

print(group_summary)


ggplot(relatedness_df, aes(x = reorder(group, wang, FUN = median), y = wang)) +
  geom_boxplot(fill = "lightblue") +
  coord_flip() +
  labs(title = "Wang Relatedness by Group",
       x = "Group", 
       y = "Wang Relatedness (r)") +
  theme_minimal()

#summary of all 
relationship_summary <- relatedness_df %>%
  group_by(relationship) %>%
  summarise(
    count = n(),
    mean_wang = mean(wang),
    sd_wang = sd(wang),
    min_wang = min(wang),
    max_wang = max(wang)
  ) %>%
  arrange(desc(mean_wang))

print(relationship_summary)

 