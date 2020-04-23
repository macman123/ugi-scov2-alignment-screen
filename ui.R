library(ape)
library(ggplot2)
library(shiny)
library(plotly)
library(Cairo)
library(gggenes)
library(shinycssloaders)
library(data.table)
library(grid)
library(gridExtra)
library(scales)
library(ggtree)
library(ape)
library(seqinr)
library(tidytree)
library(plyr)
library(shinydashboard)
library(shinyWidgets)
library(shinycustomloader)

setwd("~/Desktop/__PhD__/ucl-genetics-institute-alignment-screen/")
setkey(stat.file,bp)
stat.file <- fread("data/SNP_homoplasy_counts_20-04-2020.csv")
genes <-  fread("data/Genes_gggenes_20-04-2020.tsv")
setkey(genes,gene)
t <- read.tree("data/tree.tree")
t <- fortify(t)
t_phylo <- as.phylo(t)
t <- t[t$isTip,]
t$loc <- unlist(sapply(strsplit(t$label,"\\/"),"[", 2))
aln <- read.alignment("data/alignment.aln", format = "fasta")

# alignment to matrix
# awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' aln.fasta.aln > aln.fasta
#aln <- read.alignment("data/aln.fasta",format = "fasta",forceToLower = F)
#aln <- as.matrix(aln)
#write.csv(aln,file = "data/aln_matrix.csv",quote = F)


clrs <- c("NC" = '#a1a48c', #beige
          "ORF1ab" = '#e4c826', # yellow
          "S" = '#c82d19', # red
          "ORF3a" = '#28316b', #purple
          "E" = '#cd0164', #pink
          "M" = '#840304', #cherry red
          "ORF6" = '#b2df8a', #light green
          "ORF7a" = '#97abc3', #teal/blue/grey
          "ORF7b" = '#7a8a9e', #teal/blue/grey
          "ORF8" = '#33a02c', #green
          "N" = '#335a93', #blue
          "ORF10" = '#ffff99') #light yellow

ui <- fluidPage(
  
  ## Custom CSS
  tags$head(tags$style(
    HTML('#snpCount, #homPlot, #tree{height:70vh !important;},
          #tree{height:75vh !important;},
         .half-fill { width: 50%; height: 100%; },
         #sidebar {background-color: #b6d1bf;}
         body, label, input, button, select { 
         font-family: "Arial";color: #3e4240}
         #panel {background-color: #386f56}
         .tabbable > .nav > li > a {color:#3e4240}
         .tabbable > .nav > li > a:focus {
         background-color: white}
         .tabbable > .nav > li > a:hover {
         background-color: #f0f0f0}
         ')
  )),
  setBackgroundColor("#386f56"),
  
  titlePanel("SARS-Cov-2 Alignment Screen"),
  
  sidebarLayout(
    sidebarPanel(
      id="sidebar",
      div(style="font-weight:bold","Displaying"),
      verbatimTextOutput("info"),
      selectizeInput('gene', 'Select gene(s)', choices = genes$gene[order(genes$start)], 
                     multiple = TRUE),
      div(style="font-weight:bold","Custom range"),
      div(style="display: inline-block;vertical-align:middle", "From: "),
      div(style="display: inline-block;vertical-align:middle;width:35%",
          textInput("from", label = NULL,value = 1)),
      div(style="display: inline-block;vertical-align:middle", "   To: "),
      div(style="display: inline-block;vertical-align:middle;width:35%",
          textInput("to", label = NULL, value = nrow(stat.file))),
      selectizeInput('loc', 'Prune tree by sampling location', 
                     choices = locs, multiple = TRUE),
      actionButton("reset", "Reset"),
      actionButton("about", "About")
    ),
  
  
    mainPanel(
      id="panel",
      # Gene arrows plot
      fluidRow(
        div(style="align:center",
            plotOutput("geneDirection", height = "75px", width = "auto",
                       click = "plot_click") %>% withSpinner(color="#0dc5c1")
               )
        ),
      
      # PLOT TABSET PANEL
      tabsetPanel(
        id='tabset',
        tabPanel("SNPs", withLoader(plotlyOutput("snpCount"))),
        tabPanel('Homoplasies', withLoader(plotlyOutput("homPlot"))),
        tabPanel('ML Tree', withLoader(plotOutput("tree")))
        )
      )
    )
)


shinyApp(ui, server)

