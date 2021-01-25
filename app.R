#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyWidgets)
library(DT)
library(dplyr)
library(reactable)
allben1 <- read.delim("~/Documents/Andersen Lab /Annotation Concordance /cendr_incorporation /ben_1_table_rev.csv", header=FALSE) %>%  
  parse_BCSQ() 
  

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("BCSQ Anotation Browser"),
  # fluidRow(
  #  column(4,
  #         pickerInput("SAMPLES",
  #                     "SAMPLES:", 
  #                     choices = c(unique(shiny_testdata$SAMPLE)), 
  #                     options = list(`actions-box` = TRUE),
  #                     multiple = T, 
  #                     selected = c(unique(shiny_testdata$SAMPLE)) )
  # )
  # ), 
  fluidRow(
    column(4,
           selectizeInput("GENE", #May need new input widget with new filtering style
                          "Gene Search (e.g. trt-1):", 
                          choices = c(unique(shiny_testdata$GENE)), 
                          selected = shiny_testdata$GENE[1],
                          multiple = TRUE, # allow for multiple inputs
                          options = list(create = FALSE)) # if TRUE, allows newly created inputs
    )
    
  ), 
 
  fluidRow(
    column(4,
           textInput("GENOME_POS",
                     'Enter Genomic Location (example: II:4942126 - 4942168'
                     #chrII:4942126 - 4942168' #Can have a starting input
                     
                     
           )
    )
    
  ),
  

  fluidRow(
    column(4,
           checkboxGroupInput("IMPACT",
                              "Select Variant Impact:", 
                              choices = c(unique(shiny_testdata$VARIANT_IMPACT)),
                              selected = ""
                              
           )
    )
    
  ), 
  fluidRow(
    column(4,
           pickerInput("CONSEQUENCE",
                       "CONSEQUENCE:", 
                       choices = c(unique(shiny_testdata$CONSEQUENCE)), 
                       options = list(`actions-box` = TRUE),
                       selected = NULL,
                       multiple = T )
    )
  ),
  fluidRow(
    column(4,
           sliderTextInput("BLOSSUM",
                       "BLOSSUM:", 
                       choices = c(sort(unique(shiny_testdata$BSCORE))))
          
                      
    )
  ), 
  

  reactableOutput("table_filter_samples")

)





# Server Side
server <- function(input, output) {


  
  
#Function that allows gene to accept blanks - Creates too large of a file ATM
  

  
  
# gene_filtered <- reactive ({
#   if(input$GENE == ""){
#   shiny_testdata()
#   } else {
#     genomic_pos_filtered() %>%
#       filter(shiny_testdata() %in% input$GENE)
#   }})


gene_filtered <- reactive({ #Does not Accept blank at the moment *FIX
   shiny_testdata %>% filter(GENE %in% input$GENE)
})
#Genomic Posistion Filtering
 
 genomic_pos_filtered <- reactive ({ #Currently Supports single chromosome queries
   if(input$GENOME_POS == ""){
     gene_filtered()
   } else {
     chr_split <- strsplit(input$GENOME_POS, ":")
     chr <- chr_split[[1]][1]
     range_split <- strsplit(chr_split[[1]][2], "-") %>% as.data.frame() #Maybe could be more efficent
     gene_filtered() %>%
     filter( CHROM == chr) %>% 
     filter(dplyr::between(POS, range_split[1,1], range_split[2,1]))
     
   }})


variantimpact_filtered <- reactive ({
  if(is.null(input$IMPACT)){
    genomic_pos_filtered()
  } else {
    genomic_pos_filtered() %>%
      filter(VARIANT_IMPACT %in% input$IMPACT)
  }})

consequence_filtered <- reactive ({ 
  if(is.null(input$CONSEQUENCE)){
    variantimpact_filtered()
  } else {
    variantimpact_filtered() %>%  
      filter(CONSEQUENCE %in% input$CONSEQUENCE)
  }})

    output$table_filter_samples <- renderReactable({reactable(consequence_filtered())})
}


# Run the application 
shinyApp(ui = ui, server = server)
