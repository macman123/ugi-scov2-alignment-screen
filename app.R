library(ggplot2);
library(shiny);
library(plotly)
library(Cairo)
library(gggenes)
library(shinycssloaders)

setwd("~/Repositories/ucl-genetics-institute-alignment-screen/")
stat.file <- data.frame(read.csv("SNP_homoplasy_counts_table_20200414.csv", as.is=T, sep='\t'))
genes <-  read.delim("Genes_gggenes_20200414.tsv", as.is=T)

genomeFeatures <- subset(unique(stat.file[2]), genome_feature!="NC")

getGenomeFeatureRange <- function(genomeFeature){
  if(genomeFeature=="All") {
    return(c(min(stat.file$bp), max(stat.file$bp)))
  }
  else {
    selection = subset(stat.file, genome_feature==genomeFeature)
    return(c(min(selection$bp), max(selection$bp)))
  }
}

ui <- fluidPage(
  titlePanel("SARS-Cov-2 Alignment Screen"),
  tags$p(
    "Visualisation toolproviding the distribution of SNPs and homoplasies across all high coverage, complete, SARS-CoV-2 assemblies, updated in close to real time. Made possible by the open data sharing platform of The GISAID Initiative [1] which makes available many thousands of uploads from the global research community. For a full list of contributors see gisaid.org.",
    tags$br(),    tags$br(),
    "Alignments are generated using the real-time phylodynamic analysis pipeline implemented by Augur [2]. Augur align employs MAFFT [3] to align all assemblies against the Wuhan-Hu-1 reference genome (NCBI RefSeq NC_045512.2).",
    tags$br(),    tags$br(),
    "SNPs are identified in the alignment using the R package adegenet [4].",
    tags$br(),    tags$br(),
    "Homoplasies are identified following maximum parsimony tree inference implemented in mpboot [5] and screening with HomoplasyFinder [6].",
    tags$br(),    tags$br(),
    "Provided plots are generated using the R packages ggplotly and gggenes with annotations as per the coordinates provided for Wuhan-Hu-1."
    ),
  
  fluidRow(
    br(),
    br(),
    column(12, align="center",
      plotOutput("geneDirection", height = 450, width = 1200) %>% withSpinner(color="#0dc5c1")
    )
  ),
  
  fluidRow(
    br(),
    br(),
    column(12, align="center",
      h4("Header"),
      br(),
      radioButtons("genome_feature", "Select Genome Feature for Inspection",
                   c('All', genomeFeatures[,]), inline=TRUE
                   ),
    ),
  ),
    fluidRow(
        plotlyOutput("snpCount", height = 450) %>% withSpinner(color="#0dc5c1")
      ),
    fluidRow(
        plotlyOutput("homPlot", height = 450) %>% withSpinner(color="#0dc5c1")
      ),
  
  fluidRow(
    h4("References"),
    tags$p(
    "[1] The GISAID Initiative. https://www.gisaid.org. [2] Augur. https://github.com/nextstrain/augur. [3] Nakamura, Yamada, Tomii, Katoh (2018). Parallelization of MAFFT for large-scale multiple sequence alignments. Bioinformatics. 34,14:2490–2492.
    [4] Jombart & Ahmed (2011). adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics. 27,21:3070-3071.
    [5] Hoang, Vinh, Flouri, Stamatakis, von Haeseler, Minh (2018). MPBoot: fast phylogenetic maximum parsimony tree inference and bootstrap approximation. BMC Evol. Biol. 18, 11.
    [6] Crispell, Balaz, Vincent Gordon (2019). HomoplasyFinder: a simple tool to identify homoplasies on a phylogeny. Microb Genom. 5(1):e000245."
    )
  )
)


server <- function(input, output) {
  
  output$geneDirection <- renderPlot({
    p <- ggplot(genes, aes(xmin = start, xmax = end, y = molecule, fill = gene)) +
      geom_gene_arrow(arrowhead_height = grid::unit(15, "mm"),arrow_body_height = grid::unit(10, "mm")) + theme_genes()
    p <- p + theme(legend.position = "top",axis.title.y = element_blank(),axis.title.x = element_blank(),legend.title=element_blank(),plot.title = element_text(size = 20, face = "bold"))
    p <- p + ggtitle("SARS-CoV-2 Genome (29903bp)")
    p + facet_wrap(~ molecule, scales = "free_y", ncol = 1) 
  })
  
  output$snpCount <- renderPlotly({
    print(getGenomeFeatureRange(input$genome_feature))
    ggplotly(
      ggplot(
        stat.file, aes(bp,SNP_count,color=genome_feature, "protein"=Protein)
        ) 
      + geom_line() 
      + theme(legend.position="bottom", plot.title= element_text(hjust = 0.5))
      + coord_cartesian(xlim = getGenomeFeatureRange(input$genome_feature))
      + ggtitle(paste("Some plot title ", input$genome_feature))
      ) %>%
      layout(legend = list(
        orientation = "h", x = 0.25, y =-0.2
      )
      )
  })
  
  output$homPlot <- renderPlotly({
    ggplotly(
      ggplot(
        stat.file, aes(bp,Min.No.ChangesonTree,color=genome_feature, "consis_index"=consistency_index, "protein"=Protein)
        ) 
      + geom_line() 
      + theme(plot.title= element_text(hjust = 0.5))
      + coord_cartesian(xlim = getGenomeFeatureRange(input$genome_feature))
      + ggtitle(paste("Some plot title ", input$genome_feature))
    )  %>%
      layout(legend = list(orientation = "h", x = 0.25, y =-0.2))
  })
  
}

shinyApp(ui, server)