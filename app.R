# load libraries
for( package in c('shiny', 'shinyjs', 'shinyBS', 'ggplot2', 'DESeq2', 'rnaseqVis', 'rprojroot', 'DT' ) ){
  library( package, character.only = TRUE )
}
source('functions.R')

# globals
testing <- FALSE
rootPath <- find_root(is_rstudio_project)

# Define UI for application
ui <- fluidPage(
  useShinyjs(), # Include shinyjs
  navbarPage("geneExpr",
             tabPanel("Heatmap",
                      sidebarLayout(
                        sidebarPanel(
                          fileInput('sampleFile', 'Load Sample File'),
                          fileInput('countFile', 'Load Count File'),
                          hr(),
                          fileInput('geneIdsFile', 'Subset by Gene Id'),
                          actionButton("subsetReset", "Reset"),
                          hr(),
                          radioButtons("transform", label = h4("Transform Counts"),
                                       choices = list("Raw" = 1, "Max Scaled" = 2, "log 10" = 3), 
                                       selected = 1),
                          hr(),
                          checkboxGroupInput("clusterCheckGroup", 
                                             label = h4("Clustering"), 
                                             choices = list("By Genes" = 1, "By Samples" = 2 ),
                                             selected = c() ),
                          hr(),
                          h4('Filter Genes'),
                          sliderInput("minMeanCount", "Mean Count Minimum Threshold:",
                                      min=0, max=1000, value=0),
                          sliderInput("maxMeanCount", "Mean Count Maximum Threshold:",
                                      min=100, max=100000, value=100000),
                          hr(),
                          h4('Downloads'),
                          radioButtons("plotFormat", label = h5("Plot File"),
                                       choices = list("pdf" = 'pdf', "png" = 'png'), 
                                       selected = 'pdf'),
                          downloadButton('downloadPlot', 'Download plot'),
                          h5('Count File'),
                          downloadButton('downloadData', 'Download counts (tsv)'),
                          h5('Gene List'),
                          downloadButton('downloadGenes', 'Download genes (tsv)'),
                          hr(),
                          fileInput('dataFile', 'Load .RData File'),
                          width = 3
                        ),
                        mainPanel(
                          bsAlert("geneIdsWarningAlert"),
                          plotOutput("exprHeatmap",
                                     dblclick = "heatmap_dblclick",
                                     brush = brushOpts(
                                       id = "heatmap_brush",
                                       resetOnNew = TRUE
                                     ),
                                     height = "1000px"
                          ),
                          width = 9
                        )
                      )
             ),
             tabPanel("Counts",
                      DT::dataTableOutput("table")
             ),
             tabPanel("Help",includeMarkdown("README.md"))             
  )
)

# server side stuff
server <- function(input, output, session) {
  # increase max upload size to 30 MB
  options(shiny.maxRequestSize=30*1024^2)
  # ranges object for zooming plot
  ranges <- reactiveValues(x = NULL, 
                           y = NULL)
  ids2Names <- reactiveValues() # this will contain genes and samples, both named character vectors
  selected <- reactiveValues() # this will contain genes and samples, both named character vectors
  warnings <- reactiveValues() # this will contain warning messages to display

  DeSeqCounts <- reactive({
    if( testing ){
      testDataFile <- file.path(rootPath, 'data', 'DESeq.shield.testdata.RData')
      load(testDataFile)
    } else{
      dataFileInfo <- input$dataFile
      sampleFileInfo <- input$sampleFile
      countFileInfo <- input$countFile
      if( !is.null(sampleFileInfo) & 
                 !is.null(countFileInfo) ){
        DESeqData <- FilesToDESeqObj( sampleFileInfo$datapath, countFileInfo$datapath )
      } else if( !is.null(dataFileInfo) ){
        load(dataFileInfo$datapath)
      }
      else{
        return( NULL )
      }
    }
    genes <- as.character(rowData(DESeqData)$Gene.name)
    names(genes) <- as.character(rowData(DESeqData)$Gene.ID)
    samples <- rownames(colData(DESeqData))
    names(samples) <- rownames(colData(DESeqData))
    if( testing ){
      print(head(genes))
      print(head(samples))
    }
    selected$genes <- genes
    selected$samples <- samples
    ids2Names$genes <- genes
    ids2Names$samples <- samples
    return( DESeqData )
  })
  
  subsetGeneList <- observe({
    geneIdsFileInfo <- input$geneIdsFile
    if( !is.null(geneIdsFileInfo) ){
      geneInfo <- read.table(geneIdsFileInfo$datapath)
      geneIds <- as.character(geneInfo[[1]])
      if(testing){ print( head(geneIds) ) }
      # check for any ids that don't exist
      nonexistentIds <- vector( 'list', length = length(geneIds) )
      Ids <- vector( 'list', length = length(geneIds) )
      for( i in seq_len( length(geneIds) ) ){
        if( sum( names(ids2Names$genes) == geneIds[i] ) == 0 ){
          nonexistentIds[[i]] <- geneIds[i]
          Ids[[i]] <- NULL
        } else{
          nonexistentIds[i] <- NULL
          Ids[[i]] <- geneIds[i]
        }
      }
      nonexistentIds <- do.call(c, nonexistentIds)
      Ids <- do.call(c, Ids)
      if( testing ){ 
        print( sprintf('Missing Ids: %s', paste(nonexistentIds, collapse=" ") ) )
        print( sprintf('Matched Ids: %s', paste(Ids, collapse = " " ) ) )
      }
      # and warn
      if( length(nonexistentIds) > 0 ){
        missingGenesWarning <- paste0("Some of the gene ids couldn't be matched! Ids: ", paste(nonexistentIds, collapse=", ") )
        createAlert(session, "geneIdsWarningAlert", "geneIdsAlert", title = "Non-matching Ids",
                    content = missingGenesWarning, append = FALSE, style = 'warning' )
      }
      # set selected genes to ids
      selected$genes <- isolate(selected$genes[ Ids ])
    }
  })
  
  # respond to reset button
  observeEvent(input$subsetReset, {
    selected$genes <- ids2Names$genes
    selected$samples <- ids2Names$samples
    warnings$geneSubset <- NULL
    reset("geneIdsFile") 
  })
  
  # retrieve normalised counts
  normalisedCounts <- reactive({
    counts <- DeSeqCounts()
    if( testing ){
      print( sprintf('Num Genes: %d', nrow(counts) ) )
    }
    if( is.null(counts) ){
      return(NULL)
    } else {
      return( counts(DeSeqCounts(), normalized=TRUE) )
    }
  })
  
  # filter counts by expression level
  filteredCounts <- reactive({
    counts <- normalisedCounts()
    if( is.null(counts) ){
      return(NULL)
    } else {
      meanCount <- apply( counts, 1, mean )
      if( testing ){
        print( sprintf('Max mean value: %f', max(meanCount)) )
      }
      counts <- counts[ meanCount >= input$minMeanCount &
                          meanCount <= input$maxMeanCount, ]
      if( testing ){
        print( sprintf('Num Genes Filtered: %d', nrow(counts) ) )
      }
      # reset selected genes to everything in counts
      geneIds <- rownames(counts)
      selected$genes <- ids2Names$genes[ geneIds ]

      return( counts )
    }
  })
  
  # subset to genes/samples
  subsettedCounts <- reactive({
    counts <- filteredCounts()
    if( is.null(counts) ){
      return(NULL)
    } else {
      if( testing ){
        print( sprintf('selected Genes length = %d', length(reactiveValuesToList(selected)$genes) ) )
        print( sprintf('selected Samples length = %d', length(reactiveValuesToList(selected)$samples) ) )
        print( sprintf('selected Gene Ids = %s', head(names(selected$genes)) ) )
        print( sprintf('selected Gene Names = %s', head(selected$genes) ) )
        print( sprintf('selected Sample Ids = %s', head(names(selected$samples)) ))
      }
      numRow <- length(selected$genes)
      numCol <- length(selected$samples)
      # If there is only one row or col the matrix gets simplified into a vector.
      # Needs to force it to stay as a matrix
      if( numRow == 1 & numCol == 1 ){
        count <- matrix( counts[ names(selected$genes), names(selected$samples) ],
                         dimnames = list( names(selected$genes), names(selected$samples) )
                        )
      } else if( numRow == 1 ){
        count <- matrix( counts[ names(selected$genes), names(selected$samples) ],
                         nrow = 1,
                         dimnames = list( names(selected$genes), names(selected$samples) )
                      )
      } else if( numCol == 1 ){
        count <- matrix( counts[ names(selected$genes), names(selected$samples) ],
                         ncol = 1,
                         dimnames = list( names(selected$genes), names(selected$samples) )
                        )
      } else{
        count <- counts[ names(selected$genes), names(selected$samples) ]
      }
      return( count )
    }
  })
  
  # cluster counts by genes/samples
  clusteredCounts <- reactive({
    counts <- subsettedCounts()
    if( is.null(counts) ){
      return(NULL)
    }
    if( testing & length(input$clusterCheckGroup) > 0 ){
      print( sprintf('Cluster checkbox value: %s', input$clusterCheckGroup ) )
    }
    if( any( input$clusterCheckGroup == "1" ) ){
      if( any( input$clusterCheckGroup == "2" ) ){
        counts <- clusterMatrix(counts, byRow = TRUE, byCol = TRUE )
      } else{
        counts <- clusterMatrix(counts, byRow = TRUE, byCol = FALSE )
      }
    } else if( any( input$clusterCheckGroup == "2" ) ){
      counts <- clusterMatrix(counts, byRow = FALSE, byCol = TRUE )
    }
    return(counts)
  })
  
  # transform counts. max scaled/log10
  transformedCounts <- reactive({
    counts <- clusteredCounts()
    if( is.null(counts) ){
      return(NULL)
    } else {
      if( testing ){
        print( sprintf('Transform checkbox value: %s', input$transform ) )
      }
      if( input$transform == 2 ){
        geneMaxCounts <- apply(counts, 1, max)
        # scale operates on the column so need to transpose, scale and then transpose back
        counts <- t( scale( t(counts), scale = geneMaxCounts, center = FALSE ) )
      } else if( input$transform == 3 ){
        counts <- log10( counts + 1 )
      }
      return(counts)
    }
  })
  
  # create plot object
  heatmapObj <- reactive({
    counts <- transformedCounts()
    if( is.null(counts) ){
      return(NULL)
    } else {
      plot <- ggplotExprHeatmap(counts)
      if( nrow(counts) <= 80 ){
        genes <- ids2Names$genes[ rownames(counts) ]
        if( ncol(counts) <= 48 ){
          plot <- plot + 
            scale_y_discrete( labels = genes ) + 
            theme( axis.text.x = element_text(colour="black", angle = 90, vjust = 0.5, hjust = 1, debug = FALSE),
                   axis.text.y = element_text(colour="black", angle = 0, debug = FALSE ) )
        } else{
          plot <- plot + 
            scale_y_discrete( labels = genes ) + 
            theme( axis.text.y = element_text(colour="black", angle = 0, debug = FALSE ) )
        }
      } else if( ncol(counts) <= 48 ){
        plot <- plot + 
          theme( axis.text.x = element_text(colour="black", angle = 90, vjust = 0.5, hjust = 1, debug = FALSE) )
      }
      return( plot )
    }
  })
  
  # # create text for warning message
  # output$warningsText <- renderText({
  #   if( !is.null(warnings) ){
  #     if( !is.null(warnings$geneSubset) ){
  #       return( paste(warnings$geneSubset, missingGenes$genes, sep="\n") )
  #     }
  #   } else{
  #     return( NULL )
  #   }
  # })
  
  # render heatmap plot
  output$exprHeatmap <- renderPlot({
    return( heatmapObj() )
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$heatmap_dblclick, {
    brush <- input$heatmap_brush
    if( testing ){
      print( sprintf('Brush: Xmin: %s Xmax: %s Ymin: %s Ymax: %s', brush$xmin, brush$xmax, brush$ymin, brush$ymax ) )
    }
    if (!is.null(brush)) {
      plotRanges <- list(
        x = floor( c(brush$xmin, brush$xmax) - c( 0.5, -0.5 ) ),
        y = floor( c(brush$ymin, brush$ymax) - c( 0.5, -0.5 ) )
      )
      # make sure values have not gone out of range
      # i.e. less than 0 or greater than ncol/nrow of the data
      plotRanges$x[1] <- ifelse(plotRanges$x[1] < 1, 1, plotRanges$x[1] )
      plotRanges$x[2] <- ifelse(plotRanges$x[2] > ncol(normalisedCounts()), ncol(normalisedCounts()), plotRanges$x[2] )
      plotRanges$y[1] <- ifelse(plotRanges$y[1] < 1, 1, plotRanges$y[1] )
      plotRanges$y[2] <- ifelse(plotRanges$y[2] > nrow(normalisedCounts()), nrow(normalisedCounts()), plotRanges$y[2] )
      
      if( testing ){
        print( sprintf('Plot Ranges X: %f %f', plotRanges$x[1], plotRanges$x[2] ) )
        print( sprintf('Plot Ranges Y: %f %f', plotRanges$y[1], plotRanges$y[2] ) )
      }
      selection <- reactiveValuesToList(selected)
      selectedGeneIds <- rev( rev( names(selection$genes) )[ seq(plotRanges$y[1], plotRanges$y[2]) ] )
      selectedSampleIds <- names(selection$samples)[ seq(plotRanges$x[1], plotRanges$x[2]) ]
      if( testing ){
        print( sprintf('Num Genes = %d, Num Samples = %d', 
                       length(selectedGeneIds),
                       length(selectedSampleIds) ) )        
      }

      selected$genes <- selection$genes[ selectedGeneIds ]
      selected$samples <- selection$samples[ selectedSampleIds ]
    } else {
      selected$genes <- ids2Names$genes
      selected$samples <- ids2Names$samples
    }
  })
  
  # for downloading the plot as a pdf/png
  output$downloadPlot <- downloadHandler(
    filename = function(){ paste('geneExpr', Sys.Date(), input$plotFormat, sep='.') },
    content = function(file) {
      heatmapPlot <- heatmapObj()
        # theme( axis.text = element_text( size = rel(0.5) ) )
      if(input$plotFormat == "pdf"){
        pdf(file, paper="special", height = 16.5, width = 11.7 ) # open the pdf device
      } else if(input$plotFormat == "png"){
        png(file, height = 960, width = 480, res = 100 ) # open the png device
      }
      print( heatmapPlot )
      dev.off()  # close device 
    }
  )
  
  # for downloading the counts file
  output$downloadData <- downloadHandler(
    filename = function(){ paste('geneExpr', Sys.Date(), 'tsv', sep='.') },
    content = function(file) {
      write.table(clusteredCounts(), file, quote=FALSE, col.names = NA, sep="\t")
    },
    contentType = 'text/tsv'
  )
  # for downloading a gene list
  output$downloadGenes <- downloadHandler(
    filename = function(){ paste('geneExpr', Sys.Date(), 'genes', 'tsv', sep='.') },
    content = function(file) {
      write.table(rownames(clusteredCounts()), file, quote=FALSE, 
                  row.names = FALSE, col.names = FALSE, sep="\t")
    },
    contentType = 'text/tsv'
  )
  
  # render the counts matrix as a table
  output$table <- DT::renderDataTable({
    data <- clusteredCounts()
    if( is.null(data) ){
      return(NULL)
    }
    DT::datatable(data)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
