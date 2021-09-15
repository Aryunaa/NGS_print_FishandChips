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
