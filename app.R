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

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Survival by subtype"),

    # Sidebar with a slider input for number of bins 
    selectInput("data", label = h2("Select box"), 
                choices = list("Survival"="Survival..months." ,"Age of diagnosis"="Age..years.at.diagnosis."), 
                selected = 1),
    selectInput("trans", label = h2("Subtype"), 
                choices = list("Original"="Original.Subtype" ,"Transcriptome"="Transcriptome.Subtype" ), 
                selected = 1),
    
    hr(),
    fluidRow(column(2, verbatimTextOutput("value"))),

        # Show a plot of the generated distribution
        mainPanel(
           #plotOutput("distPlot")
           plotlyOutput("distPlot")
           
        )
    )


# Define server logic required to draw a histogram
server <- function(input, output) {
    tcga_subtype_data <-
        TCGAquery_subtype(tumor = "GBM")
    g <-tcga_subtype_data %>%
        filter(complete.cases(.$Gender))

    output$distPlot <- renderPlotly({
        # generate bins based on input$bins from ui.R
        pc <-input$data
        pc1 <-input$trans
        g %>%
            ggplot(aes_(as.name(pc), fill = as.name(pc1))) +
            geom_density(alpha = 0.4) +
            facet_wrap(~Gender, scales = "free", ncol = 1) +
            scale_fill_jco() +
            theme_linedraw() +
            theme(axis.title.x = element_blank()) +
            labs(y = "Density")
            
       # ggplotly()
            

        
     })
}

# Run the application 
shinyApp(ui = ui, server = server)
