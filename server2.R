server <- function(input, output, session) {
  
  # plotting the genes
  # jsonlite outputs a warning message because of this plot
  # In the future this will have to be plotted using ggplotly
  # Unfortunately, ggplotly does not support gggenes plots at the moment
  # Once this is supported we should switch to ggplotly
  output$geneDirection <- renderPlot({
    ggplot(genes, aes(xmin = start, xmax = end, y = molecule, 
                      fill = gene)) +
      geom_segment(x = 1,y=1,xend=29903,yend=1,size=2.5) +
      geom_gene_arrow(arrowhead_height = unit(0.9, "npc"),
                      arrow_body_height = unit(0.9, "npc"),
                      color="#3e4240",
                      show.legend = FALSE,size=1,
                      arrowhead_width = unit(0.035, "npc")) + 
      geom_gene_label(aes(label=gene),grow = T,
                      color="#3e4240",
                      height = unit(0.9, "npc"),
                      padding.y = unit(0, "lines"),
                      padding.x = unit(0.01, "npc")) +
      scale_fill_manual(values = hlight$color) +
      scale_y_discrete() +
      scale_x_continuous(breaks = NULL) +
      theme_genes() +
      theme(plot.background = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid = element_blank()) +
      scale_color_manual(values=hlight$color)
  }, height=50, width = "auto", bg="transparent" )
  
  ### REACTIVE ###
  #selection variable
  info <- reactiveValues(gene=NULL,SNP=NULL,country="All")
  output$info <- renderPrint(info)
  
  # recording the click on the gene & higlighting it
  hlight <- reactiveValues(color= clrs)
  observeEvent(input$plot_click,{
    gene <- genes$gene[genes$start <= input$plot_click$x &
                         genes$end > input$plot_click$x]
    updateSelectizeInput(session,'gene',selected = gene)
  })
  observeEvent(input$gene,{
    if(is.null(input$gene)){
      hlight$color <- clrs
      updateTextInput(session, "from", value = 1)
      updateTextInput(session, "to", value = nrow(stat.file))
    }
    if(!is.null(input$gene)){
      hlight$color <- clrs
      p <- match(input$gene,names(hlight$color))
      p <- min(p):max(p)
      hlight$color[-p] <- "#5b615e"
      updateTextInput(session, "from", value = min(genes[input$gene,start]))
      updateTextInput(session, "to", value = max(genes[input$gene,end]))
    }
  },ignoreNULL = F)

  # reset plots
  observeEvent(input$reset, {
    hlight$color <- rep("black",length(genes$gene))
    updateTextInput(session, "from", value = 1)
    updateTextInput(session, "to", value = nrow(stat.file))
    updateSelectizeInput(session,'gene',selected = '')
    updateSelectizeInput(session,'loc',selected = '')
  })
  
  # subsetting stat.file with gene clicked or with input in text
  stat.file.sub <- reactiveValues(data = stat.file)
  observeEvent(c(input$from,input$to), {
    if(!is.null(input$to) & !is.null(input$from))
       stat.file.sub$data <- stat.file[input$from:input$to]
  })
  
  
  ### PLOTTING
  
  # plot snpCont
  output$snpCount <- renderPlotly({
    ggplotly(
      ggplot(data = stat.file.sub$data, 
             aes(x=bp,y=SNP_count,color=genome_feature, "protein"=Protein)) + 
        geom_line() +
        scale_x_continuous(breaks = pretty_breaks(n = 8),minor_breaks = pretty_breaks(n=16)) +
        theme_classic() + 
        theme(plot.title= element_text(hjust = 0.5)) +
        ggtitle(paste("Some plot title ", input$genome_feature)) +
        scale_color_manual(values = clrs),
      source="SNP"
    ) %>%
      layout(legend = list(
        orientation = "h", xanchor = "center", yanchor="top", x = 0.5, y = -0.15))
  })
  
  # plot homoplasies
  output$homPlot <- renderPlotly({
    ggplotly(
      ggplot(stat.file.sub$data, 
             aes(x=bp, y=Min.No.ChangesonTree, color=genome_feature, 
                 "consis_index"=consistency_index, "protein"=Protein)) +
        geom_line() +
        scale_x_continuous(breaks = pretty_breaks(n = 8),
                           minor_breaks = pretty_breaks(n=16)) +
        scale_color_manual(values = clrs) +
        theme_classic() + 
        theme(plot.title= element_text(hjust = 0.5)) +
        ggtitle(paste("Some plot title ", input$genome_feature)),
      source="SNP"
    )  %>%
      layout(legend = list(
        orientation = "h", xanchor = "center", yanchor="top", x = 0.5, y = -0.15))
  })
  
  # plot tree
  output$tree <- renderPlot({
    if(!is.null(input$loc)){
      samples <- t$label[!(t$loc %in% input$loc)]
      tree <- drop.tip(t_phylo,samples)
    }
    else
      tree <- t_phylo
    
    SNP <- rep("None selected",aln$nb)
    position <- event_data("plotly_click", source = "SNP")$x
    if (!is.null(position)){
      position <- round(position)
      SNP <- substr(aln$seq,position,position)
    }
      
    names(SNP) <- aln$nam
    SNP <- SNP[tree$tip.label]
    tree <- fortify(tree)
    tree$SNP <- SNP[match(tree$label, names(SNP))]
    ggtree(tree)+ 
      geom_tippoint(aes(fill=SNP),show.legend = T,pch=21) +
      theme(legend.position="right")
    })
}



