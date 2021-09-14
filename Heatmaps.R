
# for unified and processed RNA-seq data
library(recount3)
# to normalize the RNA-seq data 
library(DESeq2) 
# for access to TCGA data
library(TCGAbiolinks)
# to look at the data
library(tidyverse)
# to visualize the mutation data
library(maftools)
# to create heatmaps
library(ComplexHeatmap)

scale2 <- function(mat, ...) {
  t(scale(t(mat), ...))
}


rse_gene <- 
  create_rse(
    subset(
      available_projects(),
      project == "GBM" & project_type == "data_sources"
    )
  )

#Now, let's explore what are we looking at?
assayNames(rse_gene)

#We need to scale the reads to be able to use them in `DESeq2` processing.
assay(rse_gene, "counts") <- 
  transform_counts(rse_gene)

#The attached colData contains a lot of information that we can use on top of the expression data. 
sample_sheet <-
  colData(rse_gene) %>%
  data.frame() %>%
  rownames_to_column("sample_id")

sample_sheet %>%
  head(n = 3)


#For plotting, we will use the variance stabilizing transformation (**vst**) normalized counts
#(as it is quicker than **rlog**). To read more about normalization, please read
#[Data transformations and visualization]
#(https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
  #data-transformations-and-visualization) section of RNA-seq data anlysis guide from Bioconductor.

#сводит дисперсии к однородному виду, нормализует распределениеы

normalized_counts <- 
  vst(assays(rse_gene)$counts)


## First visualization
#Let's see the expression of top variable genes.
#rowVArs - выбирает топ самых эксперсисвных генов

row_var <-
  rowVars(normalized_counts)


ht <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          col = viridis::viridis(100))

ht

#make sample type
sample_sheet %>% 
  select(tcga.tcga_barcode, tcga.cgc_sample_sample_type) %>% 
  mutate(patient_id = str_extract(tcga.tcga_barcode, 
                                  "[^-]{4}-[^-]{2}-[^-]{4}")) %>% 
  group_by(patient_id) %>% 
  summarise(count = n(), 
            sample_type = paste(unique(sort(tcga.cgc_sample_sample_type)),
                                collapse = ", ")) %>% 
  filter(count > 1)

#make gender
sample_sheet %>% 
  select(tcga.tcga_barcode,tcga.gdc_cases.demographic.gender ) %>% 
  mutate(patient_id = str_extract(tcga.tcga_barcode, 
                                  "[^-]{4}-[^-]{2}-[^-]{4}")) %>% 
  group_by(patient_id) %>% 
  summarise(count = n(), 
            gender = paste(unique(sort(tcga.gdc_cases.demographic.gender)),
                                collapse = ", ")) 

sample_sheet %>%
  select(starts_with("tcga.gdc_cases.demographic.gender"))

#heatmap annotation
ha <-
  sample_sheet %>% 
  select(sample_id, `Sample type` = tcga.cgc_sample_sample_type) %>%
  column_to_rownames("sample_id") %>%
  HeatmapAnnotation(df = .)

ht2 <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          col = viridis::viridis(100), top_annotation = ha)

ht2

###########gen
ha_gen <-
  sample_sheet %>% 
  #select(sample_id, `Sample type` = tcga.cgc_sample_sample_type) %>%
  select(sample_id, `gender` = tcga.gdc_cases.demographic.gender) %>%
  column_to_rownames("sample_id") %>%
  HeatmapAnnotation(df = .)

htg <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          top_annotation = ha_gen)

htg

ha_race <-
  sample_sheet %>% 
  #select(sample_id, `Sample type` = tcga.cgc_sample_sample_type) %>%
  select(sample_id, `race` = tcga.gdc_cases.demographic.race) %>%
  column_to_rownames("sample_id") %>%
  HeatmapAnnotation(df = .)

ha_ethno <-
  sample_sheet %>% 
  #select(sample_id, `Sample type` = tcga.cgc_sample_sample_type) %>%
  select(sample_id, `race` = tcga.gdc_cases.demographic.ethnicity) %>%
  column_to_rownames("sample_id") %>%
  HeatmapAnnotation(df = .)

############################


tryCatch(maf <- tcgaLoad(study = "GBM"), 
         error = function(e) {
           print(paste(rep("#", 50), collapse = ""))
           print(paste0("# ERROR! Read the message below!", 
                        paste(rep(" ", 17), collapse = ""),
                        "#"))
           print(paste(rep("#", 50), collapse = ""))
           print(e)
           print(paste("If you're seeing this message you probably don't have",
                       "maftools package loaded, or have an older version.", 
                       "This function is available with v2.8.",
                       "Install the new version of maftools package with",
                       "BiocManager::install('PoisonAlien/maftools')", 
                       "and try again!"))
         })
plotmafSummary(maf = maf, rmOutlier = TRUE, 
               addStat = 'median', dashboard = TRUE, log_scale = TRUE)

tcga_mutation_burden <- 
  tcgaCompare(maf = maf, cohortName = "GBM")
#Now, we can look at the top mutated genes
oncoplot(maf = maf, top = 10)
oncoplot(maf = maf, top = 10,removeNonMutated = TRUE)



# classifies Single Nucleotide Variants into Transitions and Transversions
titv = titv(maf = maf, 
            plot = FALSE, 
            useSyn = TRUE)
plotTiTv(res = titv)
getSampleSummary(maf) %>% head

lollipopPlot(maf, "PTEN", labelPos = 'all')
lollipopPlot(maf, "TP53", labelPos = 'all')
lollipopPlot(maf, "TTN", labelPos = 'all')

somaticInteractions(maf, top = 10, pvalue = c(0.01, 0.05))
#TTN - PTEN
#ARTX - TP53


#clinics
sample_sheet %>%
  select(starts_with("tcga.gdc_cases.demographic.gender"))

tcga_subtype_data <-
  TCGAquery_subtype(tumor = "gbm")

tcga_subtype_data %>%
  select(ends_with("subtype"))


#############
sample_sheet_with_patient_id <-
  sample_sheet %>% 
  mutate(patient_id = str_extract(tcga.tcga_barcode, 
                                  "[^-]{4}-[^-]{2}-[^-]{4}")) %>%
  select(sample_id, patient_id)

pten_mutation <-
  maf@data %>%
  mutate(patient_id = str_extract(Tumor_Sample_Barcode, 
                                  "[^-]{4}-[^-]{2}-[^-]{4}")) %>%
  select(Hugo_Symbol, patient_id) %>%
  filter(Hugo_Symbol == "PTEN")

tp53_mutation <-
  maf@data %>%
  mutate(patient_id = str_extract(Tumor_Sample_Barcode, 
                                  "[^-]{4}-[^-]{2}-[^-]{4}")) %>%
  select(Hugo_Symbol, patient_id) %>%
  filter(Hugo_Symbol == "TP53")

ttn_mutation <-
  maf@data %>%
  mutate(patient_id = str_extract(Tumor_Sample_Barcode, 
                                  "[^-]{4}-[^-]{2}-[^-]{4}")) %>%
  select(Hugo_Symbol, patient_id) %>%
  filter(Hugo_Symbol == "TTN")

egfr_mutation <-
  maf@data %>%
  mutate(patient_id = str_extract(Tumor_Sample_Barcode, 
                                  "[^-]{4}-[^-]{2}-[^-]{4}")) %>%
  select(Hugo_Symbol, patient_id) %>%
  filter(Hugo_Symbol == "EGFR")
  
####
ha_bottom <-
  sample_sheet_with_patient_id %>%
  left_join(pten_mutation) %>%
  mutate(present = ifelse(is.na(Hugo_Symbol), FALSE, TRUE)) %>%
  group_by(sample_id) %>%
  summarise(PTEN = any(present)) %>%
  column_to_rownames("sample_id") %>%
  HeatmapAnnotation(df = ., col = list("PTEN" = c(`TRUE` = "black", 
                                                  `FALSE` = "white")))
tp53_ha_bottom <-
  sample_sheet_with_patient_id %>%
  left_join(tp53_mutation) %>%
  mutate(present = ifelse(is.na(Hugo_Symbol), FALSE, TRUE)) %>%
  group_by(sample_id) %>%
  summarise(TP53 = any(present)) %>%
  column_to_rownames("sample_id") %>%
  HeatmapAnnotation(df = ., col = list("TP53" = c(`TRUE` = "black", 
                                                  `FALSE` = "white")))

ttn_ha_bottom <-
  sample_sheet_with_patient_id %>%
  left_join(ttn_mutation) %>%
  mutate(present = ifelse(is.na(Hugo_Symbol), FALSE, TRUE)) %>%
  group_by(sample_id) %>%
  summarise(TTN = any(present)) %>%
  column_to_rownames("sample_id") %>%
  HeatmapAnnotation(df = ., col = list("TTN" = c(`TRUE` = "black", 
                                                  `FALSE` = "white")))

egfr_ha_bottom <-
  sample_sheet_with_patient_id %>%
  left_join(egfr_mutation) %>%
  mutate(present = ifelse(is.na(Hugo_Symbol), FALSE, TRUE)) %>%
  group_by(sample_id) %>%
  summarise(EGFR = any(present)) %>%
  column_to_rownames("sample_id") %>%
  HeatmapAnnotation(df = ., col = list("EGFR" = c(`TRUE` = "black", 
                                                 `FALSE` = "white")))




#######
ht3 <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          col = viridis::viridis(100), 
          top_annotation = ha_gen, bottom_annotation = ha_bottom)

ht3


ht_tp53 <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
                  show_row_names = FALSE, show_column_names = FALSE,
                  clustering_distance_rows = "pearson", name = "gene expression",
                  col = viridis::viridis(100), 
                  top_annotation = ha_gen, bottom_annotation = tp53_ha_bottom)
ht_tp53

ht_ttn <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          top_annotation = ha_race, bottom_annotation = ttn_ha_bottom)
ht_ttn

#попробовать на топ 10
ht_egfr <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          col = viridis::viridis(100), 
          top_annotation = ha, bottom_annotation = egfr_ha_bottom)
ht_egfr



#############################

plotVaf(maf = maf, vafCol = 'i_TumorVAF_WU')

sample_sheet %>%
  select(starts_with("tcga.gdc_cases.demographic."))






#######survival
# install.packages("survminer")
#https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/analysis.html
library('survminer')
clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
TCGAanalyze_survival(clin.gbm,
                     "gender",
                     main = "TCGA Set\n GBM",height = 10, width=10)


#####PCA
