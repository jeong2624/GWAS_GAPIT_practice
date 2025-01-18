# Load necessary libraries
packages <- c("data.table", "tidyverse", "GAPIT", "GWAS.utils", 
              "CMplot", "foreach", "doParallel")
sapply(packages, require, character.only = TRUE)

# Set working directory
setwd("~/Desktop/GWAS_test")

# Register parallel backend (using one less core to avoid overloading)
nCores <- detectCores() - 1  
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Load genotype and phenotype data
df_geno <- fread("HapMap_3_r3_1.hmp.txt", header = FALSE, data.table = FALSE)
df_pheno <- fread("Phenotype.csv", data.table = FALSE) %>%
  filter(Phenotype != -9) %>% # Remove missing phenotype
  mutate(Phenotype = ifelse(Phenotype == 1, 0, 1)) # For logistic regression, replace the values.

# Load covariate and kinship data
df_pop <- fread("pops_HapMap_3_r3", data.table = FALSE)

# Load eigenvector data and merge with population data
df_covariate <- fread("HapMap_3_r3_1.eigenvec", data.table = FALSE) %>%
  mutate(Group = df_pop[df_pop$IID %in% .$SampleName, "population"],
         Sex = df_pop[df_pop$IID %in% .$SampleName, "sex"]) %>%
  select(SampleName, Group, Cluster, Sex, starts_with("PC")) %>%
  rename(IND = SampleName, Pop = Group) %>%
  filter(IND %in% df_pheno$IND)

df_kinship <- fread("HapMap_3_r3_1.Normalized_IBS.Kinship", header = FALSE) %>%
  rename(IND = V1)

# PCA analysis and plot
pca_plot <- df_covariate %>%
  ggplot(aes(x = PC1, y = PC2, color = Cluster, shape = Pop)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2") +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Save PCA plot
ggsave("PCA_plot.jpeg", plot = pca_plot)

# Set model list for GWAS
model_names <- c("GLM")

# Start time measurement for GWAS
start_time <- proc.time()

# Perform cross-validation for each model in parallel
results <- foreach(i = model_names, .packages = 'GAPIT') %dopar% {
  myGAPIT <- GAPIT(
    Y = df_pheno,
    G = df_geno,
    CV = df_covariate[, c(1, 4:7)],
    KI = df_kinship,
    SNP.MAF = 0.05,
    file.output = FALSE,
    model = i
  )
  list(GWAS = myGAPIT$GWAS)  # Extract GWAS results
}

# Stop the cluster after computation
stopCluster(cl)

# End time measurement
end_time <- proc.time()
execution_time <- end_time - start_time
print(paste("Total time:", execution_time['elapsed'], "seconds"))

# Extract GWAS results for each model
gwas_results_list <- setNames(
  lapply(results, function(res) res[["GWAS"]]),
  model_names
)

# Prepare GWAS data for Manhattan plot
GWAS_table <- gwas_results_list$GLM
colnames(GWAS_table)[3] = "Position"
rownames(GWAS_table) <- NULL
GWAS_table <- GWAS_table %>%
  select(c("SNP", "Chromosome", "Position", "P.value"))

# Clean up unnecessary data
rm(df_pop, execution_time, end_time, start_time, cl, nCores, model_names, pca_plot, results,
   df_geno, df_pheno, df_kinship, df_covariate)

# Calculate genomic inflation factor (lambda)
lambda <- round(genomic_inflation(P = GWAS_table$P.value), 2)

# Define significant threshold
threshold <- 0.05 / nrow(GWAS_table)
Marker_SNPs <- GWAS_table[GWAS_table$P.value < threshold, "SNP"]

# Create a Manhattan plot using CMplot
CMplot(GWAS_table,
       type = "p", 
       plot.type = "m",
       LOG10 = TRUE,
       threshold = threshold,
       file.output = TRUE,
       chr.labels.angle = 45,
       highlight = Marker_SNPs,
       signal.col = "red",
       highlight.text = Marker_SNPs,
       main = "Manhattan plot")

# Create a Q-Q plot using CMplot
CMplot(GWAS_table,
       plot.type = "q",
       conf.int = TRUE,
       conf.int.col = NULL,
       threshold.col = "red",
       threshold.lty = 2,
       file.output = TRUE,
       main = paste0("Q-Q plot (λ = ", lambda, ")"))

# Create a SNP-density plot using CMplot
CMplot(GWAS_table,
       plot.type = "d",
       chr.den.col = c("darkgreen", "yellow", "red"),
       file.output = TRUE)
