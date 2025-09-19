setwd("~/OneDrive - University of California, San Diego Health/Projects/Autism_Pierce")

library(tidyverse)
library(mixOmics)
library(ggpubr)
library(vegan)
library(caret)
library(limma)
library(patchwork)
library(rstatix)
library(readxl)
library(Spectra)
library(MsBackendMgf)

# Read in data
feature_table <- read_csv("mzmine/autism_gnps_quant.csv")
metadata_metabolomics <- read_csv("metadata_metabolomics.csv") %>%
  dplyr::mutate(race = replace_na(as.character(race), "unknown")) %>%
  dplyr::mutate(household_income_code = replace_na(as.character(household_income_code), "unknown"))
annotations <- read.delim("fbmn_annotations.tsv")
annotations$X.Scan. <- as.character(annotations$X.Scan.)

info_feature <- feature_table %>% dplyr::select(1:3,7)
colnames(info_feature) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature$Feature <- as.character(info_feature$Feature)

info_feature_complete <- info_feature %>% 
  left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  dplyr::select(1:4,18,24)

# Data table
data <- feature_table %>%
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE) %>%
  dplyr::mutate(SampleID = gsub(".mzML Peak area", "", SampleID)) %>%
  dplyr::mutate(SampleID = gsub(".*_(36)", "\\1", SampleID))

# Investigate total peak area
data_TIC <- data.frame(TIC = rowSums(data %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

data_TIC %>%
  ggscatter("Order", "TIC", add = "reg.line") +
  stat_cor()

metadata_metabolomics <- metadata_metabolomics %>% 
  left_join(data_TIC %>% dplyr::select(SampleID, TIC))

# Check sample type
sample_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "3644", SampleID)) %>% summarise(mean(TIC))
pool_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = regex("Pool", ignore_case = TRUE), SampleID)) %>% summarise(mean(TIC))
six_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "mix", SampleID)) %>% summarise(mean(TIC))
blank_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Blank", SampleID)) %>% summarise(mean(TIC))

ratio_tic_pb <- pool_tic/sample_tic
ratio_tic_sb <- sample_tic/blank_tic
ratio_tic_6p <- pool_tic/six_tic

# Check TIC overtime for QCpool, QCmix, Blank, SRM
data_TIC %>% dplyr::filter(str_detect(pattern = regex("Pool", ignore_case = TRUE), SampleID)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 6e8) + stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "mix", SampleID)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 6e8) + stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "Blank", SampleID)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 2.5e8) + stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = regex("3644", ignore_case = TRUE), SampleID)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 2.5e9) + stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = regex("3643", ignore_case = TRUE), SampleID)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 2.5e8) + stat_cor()

# Check internal standard (IS)
fbmn_IS <- annotations %>% dplyr::filter(str_detect(Compound_Name, regex("sulf", ignore_case = TRUE))) %>% 
  distinct(X.Scan., .keep_all = TRUE) %>% dplyr::filter(Organism != "BILELIB19")

# Extract IS 
table_IS <- data %>% column_to_rownames("SampleID") %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% 
  dplyr::filter(ID %in% fbmn_IS$X.Scan.) %>% column_to_rownames("ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("SampleID") %>% dplyr::filter(!(str_detect(SampleID, "mix|Blank|Pool|pool|wash|3643"))) %>%
  dplyr::select(SampleID, `7963`) %>%
  left_join(metadata_metabolomics)

colnames(table_IS)[2] <- "Sulfamethazine"

table_IS %>% ggscatter(x = "Order", y = "Sulfamethazine", add = c("reg.line")) + ylim(0, 5e6) + stat_cor()

table_IS %>% ggbarplot(x = "Order", y = "Sulfamethazine", xlab = "Run Order", 
                       ylab = "Peak Area Sulfamethazine", title = "Internal Standard Acquisition") +
  geom_hline(yintercept = mean(table_IS$Sulfamethazine, na.rm = TRUE), linetype = "dashed", color = "blue")

cv_is <- sd(table_IS$Sulfamethazine)/mean(table_IS$Sulfamethazine) 

# Check features per sample type
data_blank <- data %>% dplyr::filter(str_detect(pattern = "Blank", SampleID))
data_pool <- data %>% dplyr::filter(str_detect(pattern = "Pool|pool", SampleID))
data_sixmix <- data %>% dplyr::filter(str_detect(pattern = "mix", SampleID))

# Blank
blanks_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                  Mean_blank = data_blank %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_blank =  data_blank %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_blank, SD_blank, CV_blank) %>% 
  dplyr::filter(Mean_blank > 0) %>% arrange(desc(Mean_blank))

# Six mix
sixmix_feature_info <- data.frame(Feature = colnames(data_sixmix)[-1],
                                  Mean_sixmix = data_sixmix %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_sixmix = data_sixmix %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sixmix = SD_sixmix/Mean_sixmix) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_sixmix, SD_sixmix, CV_sixmix) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% arrange(desc(Mean_sixmix))

# Pool
pools_feature_info <- data.frame(Feature = colnames(data_pool)[-1],
                                 Mean_pool = data_pool %>% column_to_rownames("SampleID") %>% colMeans(), 
                                 SD_pool =  data_pool %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_pool = SD_pool/Mean_pool) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_pool, SD_pool, CV_pool) %>% 
  dplyr::filter(Mean_pool > 0) %>% arrange(desc(Mean_pool))

# Features to be removed Pools/Blank < 5
feature_to_remove <- blanks_feature_info %>% left_join(pools_feature_info) %>%
  dplyr::filter(Mean_blank > 0) %>% 
  dplyr::mutate(Pool_Blank = Mean_pool/Mean_blank) %>% 
  dplyr::filter(Pool_Blank < 5 | is.na(Pool_Blank))

# Data with blank removal
data_clean <- data %>% dplyr::select(-c(feature_to_remove$Feature)) %>%
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")

# Features to be removed Pool/QCmix < 5
feature_to_remove_mix <- sixmix_feature_info %>% left_join(pools_feature_info) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% 
  dplyr::mutate(Pool_Mix = Mean_pool/Mean_sixmix) %>% 
  dplyr::filter(Pool_Mix < 5 | is.na(Pool_Mix)) %>% 
  dplyr::filter(!(Feature %in% feature_to_remove$Feature))

# Data with QCmix removal
data_clean2 <- data_clean %>% dplyr::select(-c(feature_to_remove_mix$Feature))

# PCA raw data
PCA_raw <- mixOmics::pca(data_clean2 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("SampleID") %>%
  dplyr::mutate(ID = case_when(str_detect(pattern = "3643", SampleID) ~ "Blank",
                               str_detect(pattern = "3644|3641", SampleID) ~ "Sample",
                               str_detect(pattern = "Blank", SampleID) ~ "Blank_run",
                               str_detect(pattern = "pool|Pool", SampleID) ~ "Pool",
                               str_detect(pattern = "mix", SampleID) ~ "Mix",
                               TRUE ~ "wash"))

i <- "ID"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA Raw -", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Keep only samples
data_sample <- data_clean2 %>% 
  dplyr::filter(!(str_detect(pattern = "mix|Pool|Blank|pool|wash|364309", SampleID)))

# PCA sample raw
PCA_raw <- mixOmics::pca(data_sample %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_metabolomics)

i <- "dx"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA -", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# RCLR transformation
data_sample_clr <- decostand(data_sample %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics) %>% 
  dplyr::mutate_at("host_age", as.numeric)

i <- "dx"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = "PCA - Fecal Metabolome", palette = c("#05629c", "#fecd06"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic(), legend = "none") +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) + coord_fixed()

#ggsave(plot = PCA_plot, filename = "pca_metabolome.svg", device = "svg", dpi = "retina", height = 2.5, width = 2.5)

# PERMANOVA
dist_metabolites <- vegdist(data_sample_clr, method = "euclidean")
disper_study <- betadisper(dist_metabolites, PCA_whole_scores$dx)
anova(disper_study)
permanova_age <- adonis2(dist_metabolites ~ dx + host_age, PCA_whole_scores, na.action = na.omit, by = "term")
permanova_cov <- adonis2(dist_metabolites ~ dx + host_age + race + household_income_code, PCA_whole_scores, na.action = na.omit, by = "term")

# PLSDA
PLSDA_diagnosis <- mixOmics::plsda(data_sample_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                PCA_whole_scores$dx, ncomp = 3, scale = TRUE)
PLSDA_diagnosis_scores <- data.frame(PLSDA_diagnosis$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_metabolomics)

PLSDA_diagnosis_plot <- PLSDA_diagnosis_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "dx", alpha = 0.6, title = "PLSDA - Fecal Metabolome",
            xlab = paste("Component 1 (", round(PLSDA_diagnosis$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_diagnosis$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend = "none", ggtheme = theme_classic(), palette = c("#05629c", "#fecd06")) +
  geom_point(data = PLSDA_diagnosis_scores %>% group_by(dx) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = dx), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(plot = PLSDA_diagnosis_plot, filename = "plsda_metabolome.svg", device = "svg", dpi = "retina", height = 2.5, width = 2.5)

Loadings_diagnosis <- plotLoadings(PLSDA_diagnosis, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_diagnosis <- perf(PLSDA_diagnosis, validation = "Mfold", folds = 4, nrepeat = 999, progressBar = TRUE, auc = TRUE) 
#plot(perf_plsda_diagnosis, legend = FALSE)

VIPs_diagnosis <- as.data.frame(mixOmics::vip(PLSDA_diagnosis))
VIPs_diagnosis_filter <- dplyr::filter(VIPs_diagnosis, VIPs_diagnosis$comp1 > 1)
VIPs_diagnosis_filter$ID <- rownames(VIPs_diagnosis_filter)
VIPs_diagnosis_select <- VIPs_diagnosis_filter %>% dplyr::select(ID, comp1)
VIPs_diagnosis_Load <- VIPs_diagnosis_select %>% 
  left_join(Loadings_diagnosis, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1))

#######################
# Univariate analysis #
#######################

trypt_plot <- data_sample_clr %>%
  dplyr::select(`5043`) %>%
  rownames_to_column("SampleID") %>%
  dplyr::rename(tryptophan = `5043`) %>% 
  left_join(metadata_metabolomics) %>%
  ggboxplot(x = "dx", y = "tryptophan", add = "jitter") + 
  stat_compare_means()

ky_plot <- data_sample_clr %>%
  dplyr::select(`5980`) %>%
  rownames_to_column("SampleID") %>%
  dplyr::rename(kynurenic_acid = `5980`) %>% 
  dplyr::filter(kynurenic_acid != 0) %>%
  left_join(metadata_metabolomics) %>%
  ggboxplot(x = "dx", y = "kynurenic_acid", add = "jitter") + 
  stat_compare_means()

ky_trp_plot <- data_sample %>%
  dplyr::select(SampleID, `5043`, `5980`) %>%
  dplyr::mutate(Ratio = `5980` / `5043`) %>%
  left_join(metadata_metabolomics) %>% 
  dplyr::filter(Ratio < 10) %>%
  dplyr::filter(Ratio > 0) %>%
  ggboxplot(x = "dx", y = "Ratio", add = "jitter", add.params = list(color = "dx", alpha = 0.6), 
            palette = c("#fecd06", "#05629c"), legend = "none", title = "KYNA/TRP") + 
  stat_compare_means() +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(plot = ky_trp_plot, filename = "ratio_kyna_trp.svg", device = "svg", dpi = "retina", height = 2.5, width = 1.5)

model_linear_simple <- data_sample %>%
  dplyr::select(SampleID, `5043`, `5980`) %>%
  dplyr::mutate(Ratio = `5980` / `5043`) %>%
  left_join(metadata_metabolomics) %>% 
  dplyr::filter(Ratio < 10) %>% 
  dplyr::filter(Ratio > 0) %>%
  dplyr::mutate(Ratio = log(Ratio)) %>%
  lm(formula = Ratio ~ dx + host_age)

summary(model_linear_simple)

model_linear_complex <- data_sample %>%
  dplyr::select(SampleID, `5043`, `5980`) %>%
  dplyr::mutate(Ratio = `5980` / `5043`) %>%
  left_join(metadata_metabolomics) %>% 
  dplyr::filter(Ratio < 10) %>% 
  dplyr::filter(Ratio > 0) %>%
  dplyr::mutate(Ratio = log(Ratio)) %>%
  lm(formula = Ratio ~ dx + host_age + race + household_income_code)

anova(model_linear_complex)

model_logistic_simple <- data_sample %>%
  dplyr::select(SampleID, `5043`, `5980`) %>%
  dplyr::mutate(Ratio = `5980` / `5043`) %>%
  left_join(metadata_metabolomics) %>% 
  dplyr::filter(Ratio < 10) %>% 
  dplyr::filter(Ratio > 0) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  dplyr::mutate(Ratio = log(Ratio)) %>%
  glm(formula = dx ~ Ratio + host_age,family = binomial(link = "logit"))

summary(model_logistic_simple)

model_logistic_complex <- data_sample %>%
  dplyr::select(SampleID, `5043`, `5980`) %>%
  dplyr::mutate(Ratio = `5980` / `5043`) %>%
  left_join(metadata_metabolomics) %>% 
  dplyr::filter(Ratio < 10) %>% 
  dplyr::filter(Ratio > 0) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  dplyr::mutate(Ratio = log(Ratio)) %>%
  glm(formula = dx ~ Ratio + host_age + race + household_income_code,family = binomial(link = "logit"))

summary(model_logistic_complex)

######################
# Check all features #
######################

univariate_nopar <- data_sample_clr %>%
  dplyr:: select(where(~ mean(.x == 0, na.rm = TRUE) <= 0.5)) %>%               # keep only feature found in at least 50% of the samples
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%                    # remove features with near zero variance
  rename_with(~ paste0("Metabolite_", .x)) %>%
  rownames_to_column("SampleID") %>%
  left_join(metadata_metabolomics %>% select(SampleID, dx), by = "SampleID") %>%
  pivot_longer(cols = starts_with("Metabolite"), names_to = "Metabolite", values_to = "Abundance") %>%
  group_by(Metabolite) %>%
  summarise(cliff_delta = effsize::cliff.delta(Abundance, dx)$estimate,         # calculate effect size
            p_value = wilcox.test(Abundance ~ dx)$p.value,                      # check p values using Wilcoxon
            .groups = "drop") %>%
  dplyr::mutate(Metabolite = gsub("Metabolite_", "", Metabolite)) %>%
  dplyr::filter(abs(cliff_delta) > 0.3) %>%                                     # keep only feature with effect size > 0.3
  dplyr::mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  dplyr::arrange(cliff_delta) %>%
  left_join(info_feature_complete, by = c("Metabolite" = "Feature")) %>%
  group_by(Corr_ID) %>% dplyr::arrange(is.na(Compound_Name), Metabolite) %>%
  slice(if (is.na(Corr_ID[1])) 1:n() else 1) %>% ungroup() %>%
  dplyr::arrange(cliff_delta)

#write_csv(x = univariate_nopar, file = "univariate_nopar.csv")

univariate_par <- data_sample_clr %>%
  dplyr:: select(where(~ mean(.x == 0, na.rm = TRUE) <= 0.5)) %>%               # keep only feature found in at least 50% of the samples
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%                    # remove features with near zero variance
  rename_with(~ paste0("Metabolite_", .x)) %>%
  rownames_to_column("SampleID") %>%
  left_join(metadata_metabolomics %>% select(SampleID, dx), by = "SampleID") %>%
  pivot_longer(cols = starts_with("Metabolite"), names_to = "Metabolite", values_to = "Abundance") %>%
  group_by(Metabolite) %>%
  summarise(cohen = effsize::cohen.d(Abundance, dx)$estimate,                   # calculate effect size
            p_value = t.test(Abundance ~ dx)$p.value,                           # check p values using t test
            .groups = "drop") %>%
  dplyr::mutate(Metabolite = gsub("Metabolite_", "", Metabolite)) %>%
  dplyr::filter(abs(cohen) > 0.3) %>%                                           # keep only feature with effect size > 0.3
  dplyr::mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::arrange(cohen) %>%
  left_join(info_feature_complete, by = c("Metabolite" = "Feature")) %>%
  group_by(Corr_ID) %>% dplyr::arrange(is.na(Compound_Name), Metabolite) %>%
  slice(if (is.na(Corr_ID[1])) 1:n() else 1) %>% ungroup() %>%
  dplyr::arrange(cohen)

#write_csv(x = univariate_par, file = "univariate_par.csv")

univariate_common <- univariate_nopar %>% inner_join(univariate_par, by = "Metabolite")

# Plot
interest_univariate_boxplot <- data_sample_clr %>%
  dplyr::select(`24977`, `26310`, `28209`, `9466`, `5974`, `23013`) %>%
  rownames_to_column("SampleID") %>%
  rename_with(~paste0("Metabolite_", .x), .cols = (2:7)) %>%
  left_join(metadata_metabolomics %>% dplyr::select(SampleID, dx)) %>%
  pivot_longer(cols = (2:7), names_to = "Metabolite", values_to = "Abundance") %>%
  ggboxplot(x = "dx", y = "Abundance", facet.by = "Metabolite", add = "jitter", scale = "free_y",
            add.params = list(color = "dx", alpha = 0.6), 
            palette = c("#fecd06", "#05629c"), ylab = "RCLR(Abundance)") + 
  stat_compare_means(label.y.npc = "middle") +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 6),
        axis.text = element_text(size = 4))

#ggsave(plot = interest_univariate_boxplot, filename = "univariate_interest.svg", device = "svg", dpi = "retina", width = 4, height = 4)


# Check covariates effect
model_24977_linear_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `24977`) %>%
  dplyr::rename(Metabolite = `24977`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age)

summary(model_24977_linear_simple)

model_24977_linear_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `24977`) %>%
  dplyr::rename(Metabolite = `24977`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age + race + household_income_code)

anova(model_24977_linear_complex)

model_24977_logistic_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `24977`) %>%
  dplyr::rename(Metabolite = `24977`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age,family = binomial(link = "logit"))

summary(model_24977_logistic_simple)

model_24977_logistic_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `24977`) %>%
  dplyr::rename(Metabolite = `24977`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age + race + household_income_code, family = binomial(link = "logit"))

summary(model_24977_logistic_complex)


model_26310_linear_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `26310`) %>%
  dplyr::rename(Metabolite = `26310`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age)

summary(model_26310_linear_simple)

model_26310_linear_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `26310`) %>%
  dplyr::rename(Metabolite = `26310`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age + race + household_income_code)

anova(model_26310_linear_complex)

model_26310_logistic_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `26310`) %>%
  dplyr::rename(Metabolite = `26310`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age,family = binomial(link = "logit"))

summary(model_26310_logistic_simple)

model_26310_logistic_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `26310`) %>%
  dplyr::rename(Metabolite = `26310`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age + race + household_income_code, family = binomial(link = "logit"))

summary(model_26310_logistic_complex)


model_28209_linear_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `28209`) %>%
  dplyr::rename(Metabolite = `28209`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age)

summary(model_28209_linear_simple)

model_28209_linear_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `28209`) %>%
  dplyr::rename(Metabolite = `28209`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age + race + household_income_code)

anova(model_28209_linear_complex)

model_28209_logistic_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `28209`) %>%
  dplyr::rename(Metabolite = `28209`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age,family = binomial(link = "logit"))

summary(model_28209_logistic_simple)

model_28209_logistic_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `28209`) %>%
  dplyr::rename(Metabolite = `28209`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age + race + household_income_code, family = binomial(link = "logit"))

summary(model_28209_logistic_complex)


model_9466_linear_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `9466`) %>%
  dplyr::rename(Metabolite = `9466`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age)

summary(model_9466_linear_simple)

model_9466_linear_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `9466`) %>%
  dplyr::rename(Metabolite = `9466`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age + race + household_income_code)

anova(model_9466_linear_complex)

model_9466_logistic_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `9466`) %>%
  dplyr::rename(Metabolite = `9466`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age,family = binomial(link = "logit"))

summary(model_9466_logistic_simple)

model_9466_logistic_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `9466`) %>%
  dplyr::rename(Metabolite = `9466`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age + race + household_income_code, family = binomial(link = "logit"))

summary(model_9466_logistic_complex)


model_5974_linear_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `5974`) %>%
  dplyr::rename(Metabolite = `5974`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age)

summary(model_5974_linear_simple)

model_5974_linear_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `5974`) %>%
  dplyr::rename(Metabolite = `5974`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age + race + household_income_code)

anova(model_5974_linear_complex)

model_5974_logistic_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `5974`) %>%
  dplyr::rename(Metabolite = `5974`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age,family = binomial(link = "logit"))

summary(model_5974_logistic_simple)

model_5974_logistic_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `5974`) %>%
  dplyr::rename(Metabolite = `5974`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age + race + household_income_code, family = binomial(link = "logit"))

summary(model_5974_logistic_complex)


model_23013_linear_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `23013`) %>%
  dplyr::rename(Metabolite = `23013`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age)

summary(model_23013_linear_simple)

model_23013_linear_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `23013`) %>%
  dplyr::rename(Metabolite = `23013`) %>%
  left_join(metadata_metabolomics) %>%
  lm(formula = Metabolite ~ dx + host_age + race + household_income_code)

anova(model_23013_linear_complex)

model_23013_logistic_simple <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `23013`) %>%
  dplyr::rename(Metabolite = `23013`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age,family = binomial(link = "logit"))

summary(model_23013_logistic_simple)

model_23013_logistic_complex <- data_sample_clr %>%
  rownames_to_column("SampleID") %>%
  dplyr::select(SampleID, `23013`) %>%
  dplyr::rename(Metabolite = `23013`) %>%
  left_join(metadata_metabolomics) %>%
  dplyr::mutate(dx = as.factor(dx)) %>%
  glm(formula = dx ~ Metabolite + host_age + race + household_income_code, family = binomial(link = "logit"))

summary(model_23013_logistic_complex)
