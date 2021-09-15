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
library(plotly)

#demogr, tcga_subtype_data, maf_demogr, maf, sample_sheet
load("vars.RData")
#############
############
#############

#TOP MUTATED GENES
top15_g <- head(maf@gene.summary$Hugo_Symbol,15)

######################################################################
####################TOP15 BY RACE###################################
######################################################################
maf_demogr %>%
  filter(Hugo_Symbol %in% top15_g) %>%
  #filter(tcga.gdc_cases.demographic.gender == "male") %>%
  group_by(., tcga.gdc_cases.demographic.race) %>%
  group_by(., tcga.gdc_cases.demographic.gender) %>%
  summarise(count = n(),
            Hugo_Symbol=Hugo_Symbol,
            tcga.gdc_cases.demographic.race=tcga.gdc_cases.demographic.race,
            tcga.gdc_cases.demographic.gender=tcga.gdc_cases.demographic.gender)%>%
  plot_ly(x = ~Hugo_Symbol, y = ~count, color=~tcga.gdc_cases.demographic.race, type = 'bar') %>% 
  #add_trace(y = ~tcga.gdc_cases.demographic.gender) %>% 
  layout(yaxis = list(title = 'Count'), barmode = 'stack')


###PLOT1 TOP15 BY GENDER & RACE#####################################################################################
fem <- maf_demogr %>%
  filter(Hugo_Symbol %in% top15_g) %>%
  filter(tcga.gdc_cases.demographic.gender == "female") %>%
  group_by(., tcga.gdc_cases.demographic.race) %>%
  group_by(., tcga.gdc_cases.demographic.gender) %>%
  summarise(count = n(),
            Hugo_Symbol=Hugo_Symbol,
            tcga.gdc_cases.demographic.race=tcga.gdc_cases.demographic.race,
            tcga.gdc_cases.demographic.gender=tcga.gdc_cases.demographic.gender)%>%
  plot_ly(x = ~Hugo_Symbol, y = ~count, color=~tcga.gdc_cases.demographic.race, type = 'bar') %>% 
  layout(yaxis = list(title = 'Count'), barmode = 'stack')

male <- maf_demogr %>%
  filter(Hugo_Symbol %in% top15_g) %>%
  filter(tcga.gdc_cases.demographic.gender == "male") %>%
  group_by(., tcga.gdc_cases.demographic.race) %>%
  group_by(., tcga.gdc_cases.demographic.gender) %>%
  summarise(count = n(),
            Hugo_Symbol=Hugo_Symbol,
            tcga.gdc_cases.demographic.race=tcga.gdc_cases.demographic.race,
            tcga.gdc_cases.demographic.gender=tcga.gdc_cases.demographic.gender)%>%
  plot_ly(x = ~Hugo_Symbol, y = ~count, color=~tcga.gdc_cases.demographic.race, type = 'bar') %>% 
  layout(yaxis = list(title = 'Count'), barmode = 'stack')

annotations = list( 
  list( 
    x = 0.1,  
    y = 1.0,  
    text = "Female",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ),  
  list( 
    x = 0.1,  
    y = 0.45,  
    text = "Male",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ))

subplot(fem, male, nrows = 2, shareX = TRUE) %>% layout(annotations = annotations) %>% 
  layout(title = 'Mutations by race and gender')

###############################################################################################################
#PLOT2 GENES BY GENDER 
maf_demogr %>%
  filter(Hugo_Symbol %in% top15_g) %>%
  #filter(tcga.gdc_cases.demographic.gender == "female") %>%
  group_by(., tcga.gdc_cases.demographic.race) %>%
  group_by(., tcga.gdc_cases.demographic.gender) %>%
  summarise(count = n(),
            Hugo_Symbol=Hugo_Symbol,
            tcga.gdc_cases.demographic.race=tcga.gdc_cases.demographic.race,
            tcga.gdc_cases.demographic.gender=tcga.gdc_cases.demographic.gender)%>%
  plot_ly(x = ~Hugo_Symbol, y = ~count, color=~tcga.gdc_cases.demographic.gender, type = 'bar') %>% 
  #add_trace(y = ~tcga.gdc_cases.demographic.gender) %>% 
  layout(yaxis = list(title = 'Count'), barmode = 'stack') 
#+ facet_grid(~ tcga.gdc_cases.demographic.race)
##################################################################################3  
#PLOT GENES BY SUBTYPE
########################################################################################
#subt %>%
#  mutate(patient_id = str_extract(tcga.tcga_barcode, "[^-]{4}-[^-]{2}-[^-]{4}")) %>%
#  {. ->> subt }
#maf_demogr<-inner_join(maf_demogr, subt, by=c("patient_id" = "patient"))

maf_demogr %>%
  filter(Hugo_Symbol %in% top15_g) %>%
  group_by(., tcga.gdc_cases.demographic.race) %>%
  group_by(., tcga.gdc_cases.demographic.gender) %>%
  group_by(., Original.Subtype) %>%
  summarise(count = n(),
            Hugo_Symbol=Hugo_Symbol,
            tcga.gdc_cases.demographic.race=tcga.gdc_cases.demographic.race,
            tcga.gdc_cases.demographic.gender=tcga.gdc_cases.demographic.gender,
            Original.Subtype=Original.Subtype)%>%
  plot_ly(x = ~Hugo_Symbol, y = ~count, color=~Original.Subtype, type = 'bar') %>% 
  layout(yaxis = list(title = 'Count'), barmode = 'stack') 

##same but FREQ################################################################################
#ORIGINAL
maf_demogr %>%
  filter(complete.cases(.$Original.Subtype)) %>%
  filter(Hugo_Symbol %in% top15_g) %>%
  group_by(., Hugo_Symbol,Original.Subtype) %>%
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n))%>%
  plot_ly(x = ~Hugo_Symbol, y = ~freq, color=~Original.Subtype, type = 'bar') %>% 
  layout(yaxis = list(title = 'freq'), barmode = 'stack') 
#TRANSCR
maf_demogr %>%
  filter(complete.cases(.$Transcriptome.Subtype)) %>%
  filter(Hugo_Symbol %in% top15_g) %>%
  group_by(., Hugo_Symbol,Transcriptome.Subtype) %>%
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n))%>%
  plot_ly(x = ~Hugo_Symbol, y = ~freq, color=~Transcriptome.Subtype, type = 'bar') %>% 
  layout(yaxis = list(title = 'freq'), barmode = 'stack')

##################################################################################################3
####SUBTYPES IN MALES VS FEMALES
####################################
maf_demogr %>%
  filter(complete.cases(.$Original.Subtype)) %>%
  filter(Hugo_Symbol %in% top15_g) %>%
  group_by(., tcga.gdc_cases.demographic.gender, Original.Subtype) %>%
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n))%>%
  plot_ly(x = ~tcga.gdc_cases.demographic.gender, y = ~freq, color=~Original.Subtype, type = 'bar') %>% 
  layout(yaxis = list(title = 'Frequency'), barmode = 'stack') 

table(maf_demogr$tcga.gdc_cases.demographic.gender)


#################################
############DEATH
################################

maf_death<-inner_join(maf@data, tcga_subtype_data, by=c("patient_id" = "patient"))

tcga_subtype_data %>%
  colnames()
#Age..years.at.diagnosis.
#Survival..months.

ggplot(tcga_subtype_data, aes(Age..years.at.diagnosis., fill = Original.Subtype))+
  geom_histogram()

ggplot(tcga_subtype_data, aes(Survival..months., fill = Transcriptome.Subtype))+
  geom_histogram()


#demogr$tcga.gdc_cases.demographic.year_of_death
ggplot(demogr, aes(demogr$tcga.gdc_cases.demographic.year_of_death, fill = tcga.gdc_cases.demographic.gender))+
  geom_bar()
ggplotly()

#save(demogr, tcga_subtype_data, maf_demogr, maf, sample_sheet, file="C:/Users/limon/OneDrive/Documents/NGSprint_data_viz-main/vars.RData")
###################
##AGE BY SUBTYPE
##################
tcga_subtype_data %>%
  filter(complete.cases(.$Gender)) %>%
  ggplot(aes(Age..years.at.diagnosis., fill = Original.Subtype)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~Gender, scales = "free", ncol = 1) +
  scale_fill_jco() +
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  labs(y = "Density")
ggplotly()

tcga_subtype_data %>%
  filter(complete.cases(.$Gender)) %>%
  ggplot(aes(Survival..months., fill = Original.Subtype)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~Gender, scales = "free", ncol = 1) +
  scale_fill_jco() +
  theme_linedraw() +
  theme(axis.title.x = element_blank()) +
  labs(y = "Density")
ggplotly()


#############
############
#############

#maf_d,tSNE_df,maf_d_demogr
load("vars_tsne.RData")
#############
############
#############

#FILTER AND PLOT
########################UPSETR#######################
maf_d %>%
  filter(., Hugo_Symbol == "TP53") %>%
  group_by(., Tumor_Sample_Barcode) %>%
  summarize(n_Frame_Shift_Del = max(Frame_Shift_Del),
            n_Frame_Shift_Ins = max(Frame_Shift_Ins),
            n_In_Frame_Del = max(In_Frame_Del), 
            n_In_Frame_Ins = max(In_Frame_Ins), 
            n_Missense_Mutation = max(Missense_Mutation), 
            n_Nonsense_Mutation = max(Nonsense_Mutation), 
            n_Nonstop_Mutation = max(Nonstop_Mutation), 
            n_Splice_Site = max(Splice_Site), 
            n_Translation_Start_Site = max(Translation_Start_Site),
            sample_type_description=sample_type_description) %>%
  {. ->> b }

unique(maf_d$sample_type_description)
upset(as.data.frame(b), queries = list(list(query = elements, params = list("sample_type_description", "Primary Solid Tumor"), color = "blue", active = T,query.name = "Primary Solid Tumor")), query.legend = "bottom")

maf_d %>%
  filter(., Hugo_Symbol == "TTN") %>%
  group_by(., Tumor_Sample_Barcode) %>%
  summarize(n_Frame_Shift_Del = max(Frame_Shift_Del),
            n_Frame_Shift_Ins = max(Frame_Shift_Ins),
            n_In_Frame_Del = max(In_Frame_Del), 
            n_In_Frame_Ins = max(In_Frame_Ins), 
            n_Missense_Mutation = max(Missense_Mutation), 
            n_Nonsense_Mutation = max(Nonsense_Mutation), 
            n_Nonstop_Mutation = max(Nonstop_Mutation), 
            n_Splice_Site = max(Splice_Site), 
            n_Translation_Start_Site = max(Translation_Start_Site),
            sample_type_description=sample_type_description) %>%
  {. ->> b }

unique(maf_d$sample_type_description)
upset(as.data.frame(b), queries = list(list(query = elements, params = list("sample_type_description", "Primary Solid Tumor"), color = "blue", active = T,query.name = "Primary Solid Tumor")), query.legend = "bottom")

##################ADDING DEMOGRAPHICS
#####################################
maf_d_demogr %>%
  filter(., Hugo_Symbol == "TTN") %>%
  group_by(., Tumor_Sample_Barcode) %>%
  summarize(n_Frame_Shift_Del = max(Frame_Shift_Del),
            n_Frame_Shift_Ins = max(Frame_Shift_Ins),
            n_In_Frame_Del = max(In_Frame_Del), 
            n_In_Frame_Ins = max(In_Frame_Ins), 
            n_Missense_Mutation = max(Missense_Mutation), 
            n_Nonsense_Mutation = max(Nonsense_Mutation), 
            n_Nonstop_Mutation = max(Nonstop_Mutation), 
            n_Splice_Site = max(Splice_Site), 
            n_Translation_Start_Site = max(Translation_Start_Site),
            sample_type_description=sample_type_description,
            tcga.gdc_cases.demographic.ethnicity=tcga.gdc_cases.demographic.ethnicity) %>%
  {. ->> b }
upset(as.data.frame(b), queries = list(list(query = elements, params = list("tcga.gdc_cases.demographic.ethnicity", "hispanic or latino"), color = "blue", active = T,query.name = "hispanic or latino")), query.legend = "bottom")
#################################################3
#################################################
###sample_sheet_with_patient_id$sample_id
###tcga_subtype_data$patient

tSNE_df %>%
  ggplot(aes(x = V1, 
             y = V2,
             color = Original.Subtype,
             shape = Transcriptome.Subtype, size=4))+geom_point()+theme(legend.position="bottom") %>%
  {. ->> fig }

ggplotly()



tSNE_df %>%
  plot_ly(x = ~V1, y = ~V2, color=~Original.Subtype, type = 'scatter',mode = 'markers', 
          symbol = ~Transcriptome.Subtype, marker = list(size = 10)) 


tSNE_df %>%
  plot_ly(x = ~V1, y = ~V2, color=~Transcriptome.Subtype, type = 'scatter',mode = 'markers', marker = list(size = 10)) 



#############
############
#############
###Aryuna's part - Heatmaps

###Upload data

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


normalized_counts <- 
  vst(assays(rse_gene)$counts)


## First visualization
#Let's see the expression of top variable genes.
#rowVArs - ???????? ??? ????? ????????????? ?????

row_var <-
  rowVars(normalized_counts)


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

ht <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          col = viridis::viridis(100), top_annotation = ha)

ht


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

###########race
ha_race <-
  sample_sheet %>% 
  #select(sample_id, `Sample type` = tcga.cgc_sample_sample_type) %>%
  select(sample_id, `race` = tcga.gdc_cases.demographic.race) %>%
  column_to_rownames("sample_id") %>%
  HeatmapAnnotation(df = .)


htrace <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          top_annotation = ha_race)

htrace

##########ethno

ha_ethno <-
  sample_sheet %>% 
  #select(sample_id, `Sample type` = tcga.cgc_sample_sample_type) %>%
  select(sample_id, `race` = tcga.gdc_cases.demographic.ethnicity) %>%
  column_to_rownames("sample_id") %>%
  HeatmapAnnotation(df = .)

htethno <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          top_annotation = ha_ethno)
htethno

#########



#####mut_burden
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
lollipopPlot(maf, "EGFR", labelPos = 'all')

#good graph!!!!!
somaticInteractions(maf, top = 10, pvalue = c(0.01, 0.05))

############# create mutation
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

############## create mutation annotation
pten_ha_bottom <-
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
#############mutatiuon heatmaps

ht_pten <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha, bottom_annotation = pten_ha_bottom)

ht_pten


ht_tp53 <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha, bottom_annotation = tp53_ha_bottom)
ht_tp53

ht_ttn <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          top_annotation = ha, bottom_annotation = ttn_ha_bottom)
ht_ttn

#??????????? ?? ??? 10
ht_egfr <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha, bottom_annotation = egfr_ha_bottom)
ht_egfr




#cool heatmaps######perv gender
ht_pten_gen <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha_gen, bottom_annotation = pten_ha_bottom)


ht_pten_gen
ht_tp53_gen <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha_gen, bottom_annotation = tp53_ha_bottom)
ht_tp53_gen

ht_ttn_gen <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          top_annotation = ha_gen, bottom_annotation = ttn_ha_bottom)
ht_ttn_gen

ht_egfr_gen <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha_gen, bottom_annotation = egfr_ha_bottom)
ht_egfr_gen



#########################race
ht_pten_race <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha_race, bottom_annotation = pten_ha_bottom)


ht_pten_race
ht_tp53_race <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha_race, bottom_annotation = tp53_ha_bottom)
ht_tp53_race

ht_ttn_race <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          top_annotation = ha_race, bottom_annotation = ttn_ha_bottom)
ht_ttn_race

ht_egfr_race <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha_race, bottom_annotation = egfr_ha_bottom)
ht_egfr_race



############ethno
ht_pten_ethno <-
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha_ethno, bottom_annotation = pten_ha_bottom)


ht_pten_ethno
ht_tp53_ethno <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha_ethno, bottom_annotation = tp53_ha_bottom)
ht_tp53_ethno

ht_ttn_ethno <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          top_annotation = ha_ethno, bottom_annotation = ttn_ha_bottom)
ht_ttn_ethno

ht_egfr_ethno <- 
  Heatmap(scale2(normalized_counts[row_var > quantile(row_var, 0.995),]),
          show_row_names = FALSE, show_column_names = FALSE,
          clustering_distance_rows = "pearson", name = "gene expression",
          #col = viridis::viridis(100), 
          top_annotation = ha_ethno, bottom_annotation = egfr_ha_bottom)
ht_egfr_ethno


#Survival_analysis############################################################################
clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
TCGAanalyze_survival(clin.gbm,
                     "gender",
                     main = "TCGA Set\n GBM",height = 10, width=10)






















