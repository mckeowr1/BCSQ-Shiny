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
library(DBI)
library(RMySQL)
library(glue)



# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("BCSQ Anotation Browser"),
  fluidRow(
   column(4,
          shinyWidgets::pickerInput("SAMPLES",
                      "SAMPLES:",
                      choices = c("AB1", "CB4856", "DL238", "JU2007", "JU258", "QX1211"),
                      selected = "CB4856",
                      multiple = TRUE)

  )
  ),
  fluidRow(
    column(4,
           textInput("GENE",
                     'Enter GENE (ex ben-1)',
                     value="ben-1"



           )
    )

  ),



                    #Complicated to provide a list of options based on samples input

  # fluidRow(
  #   column(4,
  #          selectizeInput("GENE", #May need new input widget with new filtering style
  #                         "Gene Search (e.g. trt-1):",
  #                         choices = "ben-1",
  #                         selected = "ben-1",
  #                         multiple = FALSE, # allow for multiple inputs
  #                         options = list(create = FALSE)) # if TRUE, allows newly created inputs
  #   )
  #
  # ),

  fluidRow(
    column(4,
           textInput("GENOME_POS",
                     'Enter Genomic Location (example: II:4942126 - 4942168',
                     value="IV:119712 - 7094479"
                     #chrII:4942126 - 4942168' #Can have a starting input


           )
    )

  ),


  fluidRow(
    column(4,
           checkboxGroupInput("IMPACT",
                              "Select Variant Impact:",
                              choices = c("HIGH", "MODERATE", "LOW", "MODIFIER", "No Annotation"),
                              selected = c("HIGH", "MODERATE")

           )
    )

  ),
  # fluidRow(
  #   column(4,
  #          pickerInput("CONSEQUENCE",
  #                      "CONSEQUENCE:",
  #                      choices = c(unique(shiny_testdata$CONSEQUENCE)),
  #                      options = list(`actions-box` = TRUE),
  #                      selected = NULL,
  #                      multiple = T )
  #   )
  # ),
  # fluidRow(
  #   column(4,
  #          sliderTextInput("BLOSSUM",
  #                      "BLOSSUM:",
  #                      choices = c(sort(unique(shiny_testdata$BSCORE))))
  #
  #
  #   )
  # ),


  reactableOutput("table_filter_samples")

  # tableOutput("tbl")

)





# Server Side
server <- function(input, output) {




#ATM Initial Query is capable of filtering 5 samples or gene


samples_filtered <- reactive({
  x <- length(input$SAMPLES)
  chr_split <- strsplit(input$GENOME_POS, ":")
  chr <- chr_split[[1]][1]
  range_split <- strsplit(chr_split[[1]][2], "-") %>% as.data.frame() #Maybe could be more efficient
  conn <- dbConnect(
        drv = RMySQL::MySQL(),
        user = 'root',
        password = '<put in password>', #*******PUT INT PW BEFORE RUNNING*******
        dbname = "bcsq_annotation",
        host = 'localhost')
      on.exit(dbDisconnect(conn), add = TRUE)

      #Return Default Gene if genomic coordinates are left blank
  if(input$GENOME_POS == ""){
    variable_gene_query <- if(x == 1){
      dbSendQuery(conn,
                  glue::glue(
                    "SELECT * FROM {input$SAMPLES[1]}
                    WHERE GENE = '{input$GENE}'"
                  ))
    }else if (x ==2 ){
      dbSendQuery(conn,
                  glue::glue(
                    "SELECT * FROM (
                        select * from {input$SAMPLES[1]}
                        union all
                        select * from {input$SAMPLES[2]}
                    ) foo where GENE = '{input$GENE}'
                      "
                    ))


    }else if (x == 3){
      dbSendQuery(conn,
                  glue::glue(
                    "SELECT * FROM (
                        select * from {input$SAMPLES[1]}
                        union all
                        select * from {input$SAMPLES[2]}
                        union all
                        select * from {input$SAMPLES[3]}
                    ) foo where GENE = '{input$GENE}'
                      "
                  ))
    }else if (x == 4){
      dbSendQuery(conn,
                  glue::glue(
                    "SELECT * FROM (
                        select * from {input$SAMPLES[1]}
                        union all
                        select * from {input$SAMPLES[2]}
                        union all
                        select * from {input$SAMPLES[3]}
                        union all
                        select * from {input$SAMPLES[4]}
                    ) foo where GENE = '{input$GENE}'
                      "
                  ))
    }else if (x == 5){
      dbSendQuery(conn,
                  glue::glue(
                    "SELECT * FROM (
                        select * from {input$SAMPLES[1]}
                        union all
                        select * from {input$SAMPLES[2]}
                        union all
                        select * from {input$SAMPLES[3]}
                        union all
                        select * from {input$SAMPLES[4]}
                        union all
                        select * from {input$SAMPLES[5]}
                    ) foo where GENE = '{input$GENE}'
                      "
                  ))
    }
    dbFetch(variable_gene_query, n = Inf)
  } else{
        variable_query <- if(x == 1){
        dbSendQuery(conn,
                    glue::glue(
                      "SELECT * FROM {input$SAMPLES[1]}
                      WHERE CHROM = '{chr}'
                      AND POS >= {range_split[1,1]}
                      AND POS <= {range_split[2,1]}"
                    ))
      }else if (x == 2){
        dbSendQuery(conn,
                    glue::glue(
                      "SELECT * FROM (
                        select * from {input$SAMPLES[1]}
                        union all
                        select * from {input$SAMPLES[2]}
                    ) foo where CHROM = '{chr}'
                      AND POS >= {range_split[1,1]}
                      AND POS <= {range_split[2,1]} "
                    ))
      }else if (x == 3){
        dbSendQuery(conn,
                    glue::glue(
                      "SELECT * FROM (
                        select * from {input$SAMPLES[1]}
                        union all
                        select * from {input$SAMPLES[2]}
                        union all
                        select * from {input$SAMPLES[3]}
                    ) foo where CHROM = '{chr}'
                      AND POS >= {range_split[1,1]}
                      AND POS <= {range_split[2,1]} "
                    ))
      }else if (x == 4){
        dbSendQuery(conn,
                    glue::glue(
                      "SELECT * FROM (
                        select * from {input$SAMPLES[1]}
                        union all
                        select * from {input$SAMPLES[2]}
                        union all
                        select * from {input$SAMPLES[3]}
                        union all
                        select * from {input$SAMPLES[4]}
                    ) foo where CHROM = '{chr}'
                      AND POS >= {range_split[1,1]}
                      AND POS <= {range_split[2,1]} "
                    ))
      }else if (x == 5){
        dbSendQuery(conn,
                    glue::glue(
                      "SELECT * FROM (
                        select * from {input$SAMPLES[1]}
                        union all
                        select * from {input$SAMPLES[2]}
                        union all
                        select * from {input$SAMPLES[3]}
                        union all
                        select * from {input$SAMPLES[4]}
                         union all
                        select * from {input$SAMPLES[5]}
                    ) foo where CHROM = '{chr}'
                      AND POS >= {range_split[1,1]}
                      AND POS <= {range_split[2,1]} "
                    ))
      }



      dbFetch(variable_query, n = Inf)
}
   })



variantimpact_filtered <- reactive ({
  if(is.null(input$IMPACT)){
    samples_filtered()
  } else {
    samples_filtered() %>%
      filter(VARIANT_IMPACT %in% input$IMPACT)
  }})


#Genomic Posistion Filtering

# genomic_pos_filtered <- reactive ({ #Currently Supports single chromosome queries
#    if(input$GENOME_POS == ""){
#      reactable(pool %>% tbl(input$SAMPLES) %>% head(n= 100) %>% as.data.frame())
#    } else {
#      chr_split <- strsplit(input$GENOME_POS, ":")
#      chr <- chr_split[[1]][1]
#      range_split <- strsplit(chr_split[[1]][2], "-") %>% as.data.frame() #Maybe could be more efficent
#      samples_filtered() %>%
#      filter( CHROM == chr) %>%
#      filter(dplyr::between(POS, range_split[1,1], range_split[2,1]))
#
#    }})
#
#

# consequence_filtered <- reactive ({
#   if(is.null(input$CONSEQUENCE)){
#     variantimpact_filtered()
#   } else {
#     variantimpact_filtered() %>%
#       filter(CONSEQUENCE %in% input$CONSEQUENCE)
#   }})

    output$table_filter_samples <- renderReactable({reactable(variantimpact_filtered(),
                                                              columns = list (
                                                              VARIANT_IMPACT = colDef (style =  function(value){
                                                                if (value == "MODIFIER") {
                                                                  color <-  "#777"
                                                                } else if (value == "LOW") {
                                                                  color <- "#15A306"
                                                                } else if (value == "MODERATE") {
                                                                  color <- "#DF9310"
                                                                } else if (value == "HIGH") {
                                                                  color <- "#FF0000"
                                                                }else if (value == "No Annotation") {
                                                                  color <- "#f0eee1"
                                                                }
                                                                list(color = color, fontWeight = "bold")
                                                              })),


                                                              bordered = TRUE




                                                              )})
}


# Run the application
shinyApp(ui = ui, server = server)
