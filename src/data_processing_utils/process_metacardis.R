library(readr)
library(tidyr)
library(dplyr)
library(logger)

source('src/ml_pipeline/clustering.R')
source('src/ml_pipeline/preprocessing.R')
log_threshold(DEBUG)

# ------------------------
# Metacardis study groups
# ------------------------
# Group 1	(266)	
# Subjects with MetS according to IDF-criteria without T2D or CAD	
# Group 2a (293)	
# Subjects with severe obesity (BMI > 35 kg/m^2) without diabetes or CAD	
# Group 2b (193)	
# Subjects eligible for bariatric surgery 	
# Group 3	(644)	
# Subjects with T2D without CAD	
# Group 4	(159)
# Subjects with a first acute CAD event	
# Group 5	(173)	
# Subjects suffering from chronic CAD without heart failure 
# Group 6	(114)	
# Subjects suffering from chronic CAD with heart failure 
# Group 7 (36)	
# Subjects with non-CAD related heart failure 	
# Group 8 (295)
# Control group, subjects exempt of CAD, T2D, MS and with BMI<25

# -----------------------
# Read metacardis files
# -----------------------
demographics <- read.table("~/MICROBIOME_METABOLOME/METACARDIS_original_processing/metadata/demographic_20201210.r", sep= "\t")
antibiotics <- read.table("~/MICROBIOME_METABOLOME/METACARDIS_original_processing/metadata/antibiotics_20201210.r", sep= "\t")
drugs <- read.table("~/MICROBIOME_METABOLOME/METACARDIS_original_processing/metadata/cmd_drugs_20201210.r", sep = '\t')
species_raw <- read.table("~/MICROBIOME_METABOLOME/METACARDIS_original_processing/input_features/hub.cellcount.motu.Species.v2.data.frame.r", sep= "\t", header = TRUE)
modules_raw <- read.table("~/MICROBIOME_METABOLOME/METACARDIS_original_processing/input_features/hub.adjusted_KEGG_module.down.10000000.v3.data.frame.r", sep= "\t", header = TRUE)
serum_mtb_ms_raw <- read.table("~/MICROBIOME_METABOLOME/METACARDIS_original_processing/input_features/hub.serum_absolute.v1.data.frame.r", sep= "\t", header = TRUE)
serum_mtb_nmr_raw <- read.table("~/MICROBIOME_METABOLOME/METACARDIS_original_processing/input_features/hub.serum_nmr_absolute.v1.data.frame.r", sep= "\t", header = TRUE)
serum_tmao_raw <- read.table("~/MICROBIOME_METABOLOME/METACARDIS_original_processing/input_features/hub.serum_tmao.v1.data.frame.r", sep= "\t", header = TRUE)

# Paths for output files
clusters_1_8_file <- 'data/ml_input/metacardis_1_8/clusters.tsv'
clusters_3_8_file <- 'data/ml_input/metacardis_3_8/clusters.tsv'
final_1_8_file <- 'data/ml_input/metacardis_1_8/all_data.tsv'
final_3_8_file <- 'data/ml_input/metacardis_3_8/all_data.tsv'

cluster_init_clusters_table(clusters_1_8_file)
cluster_init_clusters_table(clusters_3_8_file)

# -----------------------------------------------
# Create separate metadata files 
# (for each case-control setting we want to use)
# -----------------------------------------------

# Similar to the exclusion criteria in the original study [Fromentin et al.], 
#  We excluded subjects with recent antibiotics (ANTIBIOTICS_TOTAL>0) and those 
#  treated by multiple drugs (DRUGTOTAL>=3). 

metadata <- demographics %>%
  left_join(antibiotics %>% select(1:2), by = 'SampleID') %>%
  left_join(drugs %>% select(SampleID, DRUGTOTAL), by = 'SampleID') %>%
  rename(sample_id__ = SampleID) %>%
  # No antibiotics/drugs
  filter(ANTIBIOTICS_TOTAL == 0) %>%
  filter(DRUGTOTAL < 3)

## TODO: consider matching cases and controls by BMI (will require removing missing BMI samples... leaves us with very small groups)
# library(MatchIt)
# m.out0 <- matchit(PATGROUPFINAL_C ~ BMI_C, data = metadata %>% filter(PATGROUPFINAL_C %in% c(1,8)) %>% mutate(PATGROUPFINAL_C = (PATGROUPFINAL_C==1)) %>% filter(!is.na(BMI_C)), method = NULL, distance = "glm")
# summary(m.out0) # See stats without matching
# m.out1 <- matchit(PATGROUPFINAL_C ~ BMI_C, data = metadata, method = "nearest", distance = "glm")
# summary(m.out1, un = FALSE) # See stats with matching
# plot(m.out1, type = "jitter", interactive = FALSE)

# 1 vs. 8 (MetS_without_T2D_CAD)
metadata_1_8 <- metadata %>%
  filter(PATGROUPFINAL_C %in% c(1,8)) %>%
  mutate(DiseaseState = ifelse(PATGROUPFINAL_C == 1, 'disease', 'healthy'))

# 3 vs. 8 (T2D_without_CAD)
metadata_3_8 <- metadata %>%
  filter(PATGROUPFINAL_C %in% c(3,8)) %>%
  mutate(DiseaseState = ifelse(PATGROUPFINAL_C == 3, 'disease', 'healthy'))

# metadata %>% group_by(PATGROUPFINAL_C) %>% summarise(N = n(), avg_age = mean(AGE, na.rm = T), avg_bmi = mean(BMI_C, na.rm=T)) %>% View()

# -----------------------------------------------
# Organize species table (add T__ prefix)
# -----------------------------------------------

species <- species_raw %>%
  select(SampleID, FeatureDisplayName, FeatureValue) %>%
  rename(sample_id__ = SampleID) %>%
  mutate(FeatureDisplayName = make.names(FeatureDisplayName)) %>%
  pivot_wider(id_cols = sample_id__, names_from = FeatureDisplayName, values_from = FeatureValue)

## Add prefix
names(species) <- c('sample_id__', paste0('T__', names(species)[-1]))

## Convert to relative abundances
species <- prep_convert_to_relative_abundances(species)

## Split 
species_1_8 <- species %>% filter(sample_id__ %in% metadata_1_8$sample_id__)
species_3_8 <- species %>% filter(sample_id__ %in% metadata_3_8$sample_id__)

## Remove rare features + constant features
species_1_8 <- prep_sanitize_dataset(species_1_8, 'T', rare_feature_cutoff = 0.15, mean_abundance_cutoff = 0.00005)
species_3_8 <- prep_sanitize_dataset(species_3_8, 'T', rare_feature_cutoff = 0.15, mean_abundance_cutoff = 0.00005)

# -----------------------------------------------
# Organize pathways table (add P__ prefix)
# -----------------------------------------------

pwy <- modules_raw %>%
  select(SampleID, FeatureDisplayName, FeatureValue) %>%
  rename(sample_id__ = SampleID) %>%
  mutate(FeatureDisplayName = make.names(FeatureDisplayName)) %>%
  pivot_wider(id_cols = sample_id__, names_from = FeatureDisplayName, values_from = FeatureValue)

## Add prefix
names(pwy) <- c('sample_id__', paste0('P__', names(pwy)[-1]))

## Convert to relative abundances
pwy <- prep_convert_to_relative_abundances(pwy)

## Split 
pwy_1_8 <- pwy %>% filter(sample_id__ %in% metadata_1_8$sample_id__)
pwy_3_8 <- pwy %>% filter(sample_id__ %in% metadata_3_8$sample_id__)

## Remove rare features + constant features
pwy_1_8 <- prep_sanitize_dataset(pwy_1_8, 'P', rare_feature_cutoff = 0.15, mean_abundance_cutoff = 0.00005)
pwy_3_8 <- prep_sanitize_dataset(pwy_3_8, 'P', rare_feature_cutoff = 0.15, mean_abundance_cutoff = 0.00005)

  
# -----------------------------------------------
# Organize serum metabolomics (add S__ prefix)
# -----------------------------------------------

mtb <- bind_rows(serum_mtb_ms_raw, serum_mtb_nmr_raw, serum_tmao_raw) %>%
  select(SampleID, Feature, FeatureValue) %>%
  rename(sample_id__ = SampleID) %>%
  mutate(Feature = make.names(Feature)) %>%
  mutate(Feature = chartr("???¹²³??????????????????", "0123456789", Feature)) %>%
  pivot_wider(id_cols = sample_id__, names_from = Feature, values_from = FeatureValue) %>%
  na.omit()

## Add prefix
names(mtb) <- c('sample_id__', paste0('S__', names(mtb)[-1]))

## Log transform
mtb <- mtb %>% mutate_if(is.numeric, log)

## Split 
mtb_1_8 <- mtb %>% filter(sample_id__ %in% metadata_1_8$sample_id__)
mtb_3_8 <- mtb %>% filter(sample_id__ %in% metadata_3_8$sample_id__)

## Remove rare features + constant features
mtb_1_8 <- prep_sanitize_dataset(mtb_1_8, 'S', rare_feature_cutoff = 0.15)
mtb_3_8 <- prep_sanitize_dataset(mtb_3_8, 'S', rare_feature_cutoff = 0.15)

# -----------------------------------------------
# Cluster highly correlated features
# -----------------------------------------------

all_data_1_8 <- metadata_1_8 %>%
  select(sample_id__, DiseaseState) %>%
  inner_join(species_1_8, by = 'sample_id__') %>%
  inner_join(pwy_1_8, by = 'sample_id__') %>%
  inner_join(mtb_1_8, by = 'sample_id__') 
  
all_data_1_8_clustered <- cluster_cluster_features(
  dataset = all_data_1_8 %>% select(-sample_id__),
  feature_set_type = 'T+P+S',
  cluster_type = 'clustering99',
  clusters_output = clusters_1_8_file
)
all_data_1_8 <- bind_cols(all_data_1_8 %>% select(sample_id__), all_data_1_8_clustered)

all_data_3_8 <- metadata_3_8 %>%
  select(sample_id__, DiseaseState) %>%
  inner_join(species_3_8, by = 'sample_id__') %>%
  inner_join(pwy_3_8, by = 'sample_id__') %>%
  inner_join(mtb_3_8, by = 'sample_id__') 

all_data_3_8_clustered <- cluster_cluster_features(
  dataset = all_data_3_8 %>% select(-sample_id__),
  feature_set_type = 'T+P+S',
  cluster_type = 'clustering99',
  clusters_output = clusters_1_8_file
)
all_data_3_8 <- bind_cols(all_data_3_8 %>% select(sample_id__), all_data_3_8_clustered)

# ---------------------------------------------------
# Save to files
# ---------------------------------------------------

# Save for intermediate integration (one large table)
write_tsv(all_data_1_8, final_1_8_file)
write_tsv(all_data_3_8, final_3_8_file)

# Print summaries (numbers of cases/controls)
table(all_data_1_8$DiseaseState)
table(all_data_3_8$DiseaseState)

