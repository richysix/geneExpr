# load libraries
for( package in c('shiny', 'ggplot2', 'DESeq2', 'rnaseqVis' ) ){
  require( package, character.only = TRUE )
}

testing <- FALSE

maxScale <- function( counts ){
  geneMaxCounts <- apply(counts, 1, max)
  # scale operates on the column so need to transpose, scale and then transpose back
  return( t( scale( t(counts), scale = geneMaxCounts, center = FALSE ) ) )
}

# loadData <- function(){
#   
# }

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
                          downloadButton('downloadData', 'Download tsv'),
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
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  DeSeqCounts <- reactive({
    if( testing ){
      testDataFile <- '/Users/rw4/sanger/lustre/scratch117/maz/team31/infection/zmp_ph263/DESeq.Dr.rnaseq.grcz10.RData'
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
      counts <- counts[ meanCount >= input$minMeanCount &
                          meanCount <= input$maxMeanCount, ]
      if( any( input$clusterCheckGroup == "1" ) ){
        if( any( input$clusterCheckGroup == "2" ) ){
          counts <- clusterMatrix(counts, byRow = TRUE, byCol = TRUE )
        } else{
          counts <- clusterMatrix(counts, byRow = TRUE, byCol = FALSE )
        }
      }
      if( any( input$clusterCheckGroup == "2" ) ){
        counts <- clusterMatrix(counts, byRow = FALSE, byCol = TRUE )
      }
      return( counts )
    }
  })
  
  transformedCounts <- reactive({
    if( testing ){
      print( input$transform )
    }
    counts <- filteredCounts()
    if( is.null(counts) ){
      return(NULL)
    } else {
      if( input$transform == 2 ){
        counts <- maxScale( counts )
      } else if( input$transform == 3 ){
        counts <- log10( counts + 1 )
      }
      return(counts)
    }
  })
  
  dataForTable <- reactive({
    # reverse order of rows
    revCounts <- transformedCounts()[ rev( rownames(transformedCounts()) ), ]
    if( is.null(ranges$x) ){
      ranges$x <- c(1,ncol(revCounts))
    }
    if( is.null(ranges$y) ){
      ranges$x <- c(1,nrow(revCounts))
    }
    count <- revCounts[ seq(ranges$y[1]+0.5, ranges$y[2]-0.5 ), 
                        seq(ranges$x[1]+0.5, ranges$x[2]-0.5 ) ]
    return( count[ rev( rownames(count) ), ] )
  })

  output$exprHeatmap <- renderPlot({
    counts <- transformedCounts()
    if( is.null(counts) ){
      return(NULL)
    } else {
      return( ggplotExprHeatmap(counts) + 
              coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand=FALSE) +
              theme( axis.text.x = element_text(colour="black", angle = 90, hjust = 1) )
          )
    }
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$heatmap_dblclick, {
    brush <- input$heatmap_brush
    if( testing ){
      print( brush )
    }
    if (!is.null(brush)) {
      ranges$x <- floor( c(brush$xmin, brush$xmax) - c( 0.5, -0.5 ) ) + c( 0.5, 0.5 )
      ranges$y <- floor( c(brush$ymin, brush$ymax) - c( 0.5, -0.5 ) ) + c( 0.5, 0.5 )
      ranges$x[2] <- ifelse(ranges$x[2] > ncol(normalisedCounts()), 
                            ncol(normalisedCounts()) + 0.5, ranges$x[2] )
      ranges$y[2] <- ifelse(ranges$y[2] > nrow(normalisedCounts()), 
                            nrow(normalisedCounts()) + 0.5, ranges$y[2] )
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
    # output for testing
    if( testing ){
      print( ranges$x )
      print( ranges$y )
    }
  })
  output$downloadData <- downloadHandler(
    filename = function(){ paste('geneExpr', Sys.Date(), 'tsv', sep='.') },
    content = function(file) {
      write.table(dataForTable(), file, quote=FALSE, col.names = NA, sep="\t")
    },
    contentType = 'text/tsv'
  )
  
  output$table <- DT::renderDataTable({
    data <- dataForTable()
    print( class(data) )
    print( dim(data) )
    DT::datatable(data)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
