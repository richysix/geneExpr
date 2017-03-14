# load libraries
for( package in c('shiny', 'ggplot2', 'DESeq2', 'rnaseqVis', 'rprojroot', 'DT' ) ){
  library( package, character.only = TRUE )
}

# globals
testing <- FALSE
rootPath <- find_root(is_rstudio_project)


maxScale <- function( counts ){
  geneMaxCounts <- apply(counts, 1, max)
  # scale operates on the column so need to transpose, scale and then transpose back
  return( t( scale( t(counts), scale = geneMaxCounts, center = FALSE ) ) )
}

clipRangeToMinMax <- function( range, minValue, maxValue ){
  range[1] <- ifelse(range[1] < minValue, minValue, range[1] )
  range[2] <- ifelse(range[2] > maxValue, maxValue, range[2] )
  return(range)
}

# Define UI for application
ui <- fluidPage(
  navbarPage("geneExpr",
             tabPanel("Heatmap",
                      sidebarLayout(
                        sidebarPanel(
                          fileInput('dataFile', 'Load Data File'),
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
                                      min=0, max=1000, value=10),
                          sliderInput("maxMeanCount", "Mean Count Maximum Threshold:",
                                      min=100, max=10000, value=1000),
                          hr(),
                          h4('Download Count File'),
                          downloadButton('downloadData', 'Download counts (tsv)'),
                          hr(),
                          h4('Download Gene List'),
                          downloadButton('downloadGenes', 'Download genes (tsv)'),
                          width = 3
                        ),
                        mainPanel(
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
             )
  )
)

# server side stuff
server <- function(input, output) {
  # increase max upload size to 30 MB
  options(shiny.maxRequestSize=30*1024^2)
  # ranges object for zooming plot
  ranges <- reactiveValues(x = NULL, 
                           y = NULL )

  DeSeqCounts <- reactive({
    if( testing ){
      testDataFile <- file.path(rootPath, 'data', 'DESeq.shield.testdata.RData')
      load(testDataFile)
      return( DESeqData )
    } else{
      fileInfo <- input$dataFile
      if (is.null(fileInfo)){
        return(NULL)
      } else{
        load(fileInfo$datapath)
        return( DESeqData )
      }
    }
  })
    
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
      return( counts )
    }
  })
  
  subsettedCounts <- reactive({
    counts <- filteredCounts()
    if( is.null(counts) ){
      return(NULL)
    } else {
      if( testing ){
        print( sprintf('Num Genes: %d', nrow(counts) ) )
      }
      # subset counts based on plot brush
      # set ranges to max extent if they are null
      if( is.null(ranges$x) ){
        ranges$x <- c(1,ncol(counts))
      }
      if( is.null(ranges$y) ){
        ranges$y <- c(1,nrow(counts))
      }
      # set ranges to max extent if they are larger than current ranges
      ranges$x[2] <- ifelse( ranges$x[2] > ncol(counts), ncol(counts), ranges$x[2] )
      ranges$y[2] <- ifelse( ranges$y[2] > nrow(counts), nrow(counts), ranges$y[2] )
      
      # y coordinates from the brush are reversed compared to the count matrix
      # need to reverse the matrix, then subset, then reverse back to get correct genes
      if( testing ){
        print( seq(ranges$y[1], ranges$y[2] ) )
        print( seq(ranges$x[1], ranges$x[2] ) )
      }
      revCounts <- counts[ rev( rownames(counts) ), ]
      count <- revCounts[ seq(ranges$y[1], ranges$y[2] ), 
                          seq(ranges$x[1], ranges$x[2] ) ]
      return( count[ rev( rownames(count) ), ] )
    }
  })
  
  clusteredCounts <- reactive({
    counts <- subsettedCounts()
    if( is.null(counts) ){
      return(NULL)
    }
    if( testing ){
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
  
  transformedCounts <- reactive({
    counts <- clusteredCounts()
    if( is.null(counts) ){
      return(NULL)
    } else {
      if( testing ){
        print( sprintf('Transform checkbox value: %s', input$transform ) )
      }
      if( input$transform == 2 ){
        counts <- maxScale( counts )
      } else if( input$transform == 3 ){
        counts <- log10( counts + 1 )
      }
      return(counts)
    }
  })
  
  geneNames <- reactive({
    DESeqData <- DeSeqCounts()
    genes <- as.character(rowData(DESeqData)$Gene.name)
    names(genes) <- as.character(rowData(DESeqData)$Gene.ID)
    return(genes)
  })
  
  output$exprHeatmap <- renderPlot({
    counts <- transformedCounts()
    if( is.null(counts) ){
      return(NULL)
    } else {
      plot <- ggplotExprHeatmap(counts)
      if( nrow(counts) <= 80 ){
        genes <- geneNames()[ rownames(counts) ]
        if( ncol(counts) <= 96 ){
          plot <- plot + 
            scale_y_discrete( labels = genes ) + 
            theme( axis.text.x = element_text(colour="black", angle = 90, hjust = 1, debug = FALSE),
                   axis.text.y = element_text(colour="black", angle = 0, debug = FALSE ) )
        } else{
          plot <- plot + 
            scale_y_discrete( labels = genes ) + 
            theme( axis.text.y = element_text(colour="black", angle = 0, debug = FALSE ) )
        }
      } else if( ncol(counts) <= 96 ){
        plot <- plot + 
          theme( axis.text.x = element_text(colour="black", angle = 90, hjust = 1, debug = FALSE) )
      }
      return( plot )
    }
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
        x = floor( c(brush$xmin, brush$xmax) - c( 0.5, -0.5 ) ) + c( 0, -1 ),
        y = floor( c(brush$ymin, brush$ymax) - c( 0.5, -0.5 ) ) + c( 0, -1 )
      )
      plotRanges$x <- clipRangeToMinMax( plotRanges$x, 0, ncol(normalisedCounts()))
      plotRanges$y <- clipRangeToMinMax( plotRanges$y, 0, nrow(normalisedCounts()))
      if( testing ){
        print( sprintf('Plot Ranges X: %f %f', plotRanges$x[1], plotRanges$x[2] ) )
        print( sprintf('Plot Ranges Y: %f %f', plotRanges$y[1], plotRanges$y[2] ) )
      }
      
      ranges$x <- c(ranges$x[1], ranges$x[1]) + plotRanges$x
      ranges$y <- c(ranges$y[1], ranges$y[1]) + plotRanges$y
      if( testing ){
        print( sprintf('Ranges X: %f %f', ranges$x[1], ranges$x[2] ) )
        print( sprintf('Ranges Y: %f %f', ranges$y[1], ranges$y[2] ) )
      }
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  output$downloadData <- downloadHandler(
    filename = function(){ paste('geneExpr', Sys.Date(), 'tsv', sep='.') },
    content = function(file) {
      write.table(dataForTable(), file, quote=FALSE, col.names = NA, sep="\t")
    },
    contentType = 'text/tsv'
  )
  output$downloadGenes <- downloadHandler(
    filename = function(){ paste('geneExpr', Sys.Date(), 'genes', 'tsv', sep='.') },
    content = function(file) {
      write.table(rownames(dataForTable()), file, quote=FALSE, 
                  row.names = FALSE, col.names = FALSE, sep="\t")
    },
    contentType = 'text/tsv'
  )
  
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
