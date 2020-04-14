library(ggplot2);library(shiny);library(plotly)
library(Cairo)   # For nicer ggplot2 output when deployed on Linux
setwd("~/Dropbox/SARS-CoV-2/GISAID_12_04_2020/Homoplasy/")
stat.file <- data.frame(read.csv("SNP_homoplasy_counts_table.csv",as.is=T))

ui <- fluidPage(
  fluidRow(
    column(width = 12,
           plotOutput("SNP_count", height = 350,hover = hoverOpts(id ="plot_hover"))
    )
  ),
  fluidRow(
    column(width = 5,
           verbatimTextOutput("hover_info")
    )
  )
)

server <- function(input, output) {
  
  
  output$SNP_count <- renderPlot({
    
    ggplotly(ggplot(stat.file,aes(bp,SNP_count,color=genome_feature)) + geom_line())
    
  })
  
  output$hover_info <- renderPrint({
    if(!is.null(input$plot_hover)){
      hover=input$plot_hover
      dist=sqrt((hover$x-stat.file$bp)^2+(hover$y-stat.file$SNP_count)^2)
      cat("Nucleotide Pos")
      if(min(dist) < 3)
        stat.file$bp[which.min(dist)]
    }
    
    
  })
}

shinyApp(ui, server)