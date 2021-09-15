#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# Setup

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
library(shiny)
library(tidyverse)
library(plotly)
library(survminer)
setwd("C:/Users/limon/OneDrive/Documents/NGSprint_data_viz-main")
load("vars.RData")
load("vars_tsne.RData")
load("counts.RData")
# Define UI for application that draws a histogram
ui <- fluidPage(
    fluidRow("Survival by subtype",
                 selectInput("data", label = h2("Data"), 
                             choices = list("Survival"="Survival..months." ,"Age of diagnosis"="Age..years.at.diagnosis."), 
                             selected = 1),
                 selectInput("trans", label = h2("Subtype"), 
                             choices = list("Original"="Original.Subtype" ,"Transcriptome"="Transcriptome.Subtype" ), 
                             selected = 1),
                 #selectInput("gene", label = h2("Select box"), 
                 #            choices = list("PTEN","TP53","TTN","EGFR","MUC16","FLG","NF1","RYR2","ATRX","SPTA1","PIK3R1", "PIK3CA", "RB1","SYNE1","OBSCN"), 
                 #            selected = 1),
                     
                 hr(),
                 mainPanel(
                     plotlyOutput("Plot1")
                     
                 )
                 ),
        
    fluidRow("Frequency o gene mutations by subr",
                 # Show a plot of the generated distribution
                 mainPanel(
                     plotlyOutput("Plot2")
                     
                 )),
    
    fluidRow("Frequency o gene mutations by subr",
             # Show a plot of the generated distribution
             mainPanel(
                 plotlyOutput("Plot3")
                 
             )),
    
    fluidRow("Genes", "Name2",
             # Show a plot of the generated distribution
             mainPanel(
                 plotlyOutput("Plot4")
                 
             )),
    fluidRow("Name1", "Name2",
             # Show a plot of the generated distribution
             mainPanel(
                 plotlyOutput("Plot5")
                 
             )),
    )

# Define server logic required to draw a histogram
server <- function(input, output) {
    tcga_subtype_data <-
        TCGAquery_subtype(tumor = "GBM")
    g <-tcga_subtype_data %>%
        filter(complete.cases(.$Gender))

    output$Plot1 <- renderPlotly({
        # generate bins based on input$bins from ui.R
        pc <-input$data
        pc1 <-input$trans
        g %>%
            ggplot(aes_(as.name(pc), fill = as.name(pc1))) +
            geom_density(alpha = 0.4) +
            facet_wrap(~Gender, scales = "free", ncol = 1) +
            scale_fill_jco() +
            theme_linedraw() +
            theme(axis.text.x = element_text(angle = 90, size = 12)) +
            labs(y = "Density")
    })
            
       # ggplotly()
    output$Plot2 <- renderPlotly({
        maf_demogr %>%
            filter(complete.cases(.$Original.Subtype)) %>%
            filter(Hugo_Symbol %in% top15_g) %>%
            group_by(., Hugo_Symbol,Original.Subtype) %>%
            summarise(n = n()) %>% 
            mutate(freq = n/sum(n))%>%
            plot_ly(x = ~Hugo_Symbol, y = ~freq, color=~Original.Subtype, type = 'bar') %>% 
            layout(yaxis = list(title = 'freq'), barmode = 'stack')
        
     })
    output$Plot3 <- renderPlotly({
        tSNE_df %>%
            ggplot(aes(x = V1, 
                       y = V2,
                       color = Original.Subtype,
                       shape = Transcriptome.Subtype, size=4))+geom_point()+theme(legend.position="bottom")
    })
    
    output$Plot4 <- renderPlotly({
    })
    
    output$Plot5 <- renderPlotly({

    })
    
    output$Plot6 <- renderPlotly({
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
