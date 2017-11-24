# load libraries
for (package in c('shiny',
                  'shinyjs',
                  'shinyBS',
                  'ggplot2',
                  'DESeq2',
                  'rnaseqVis',
                  'rprojroot',
                  'DT',
                  'genefilter')) {
  library(package, character.only = TRUE)
}
source('R/functions.R')

# globals
testing <- FALSE
debug <- TRUE
rootPath <- find_root(is_rstudio_project)

# Define UI for application
ui <- fluidPage(useShinyjs(),
                # Include shinyjs
                navbarPage(
                  "geneExpr",
                  tabPanel("Heatmap",
                           sidebarLayout(
                             sidebarPanel(
                               radioButtons(
                                 "dataType",
                                 label = h5("Data Type"),
                                 choices = list("RNA-seq" = 'rnaseq', "DeTCT" = 'detct'),
                                 selected = 'rnaseq'
                               ),
                               fileInput('sampleFile', 'Load Sample File'),
                               fileInput('countFile', 'Load Count File'),
                               hr(),
                               fileInput('geneIdsFile', 'Subset by Gene Id'),
                               actionButton("subsetReset", "Reset"),
                               hr(),
                               radioButtons(
                                 "transform",
                                 label = h4("Transform Counts"),
                                 choices = list(
                                   "Raw" = 1,
                                   "Max Scaled" = 2,
                                   "log 10" = 3,
                                   "Mean Centred and Scaled" = 4
                                 ),
                                 selected = 1
                               ),
                               hr(),
                               checkboxGroupInput(
                                 "clusterCheckGroup",
                                 label = h4("Clustering"),
                                 choices = list("By Genes" = 1, "By Samples" = 2),
                                 selected = c()
                               ),
                               hr(),
                               h4('Filter Genes'),
                               sliderInput(
                                 "minMeanCount",
                                 "Mean Count Minimum Threshold:",
                                 min = 0,
                                 max = 1000,
                                 value = 0
                               ),
                               sliderInput(
                                 "maxMeanCount",
                                 "Mean Count Maximum Threshold:",
                                 min = 100,
                                 max = 1000000,
                                 value = 1000000
                               ),
                               hr(),
                               h4('Downloads'),
                               radioButtons(
                                 "plotFormat",
                                 label = h5("Plot File"),
                                 choices = list("pdf" = 'pdf', "png" = 'png'),
                                 selected = 'pdf'
                               ),
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
                               bsAlert("HeatmapAlert"),
                               plotOutput(
                                 "exprHeatmap",
                                 dblclick = "heatmap_dblclick",
                                 brush = brushOpts(id = "heatmap_brush",
                                                   resetOnNew = TRUE),
                                 height = "1000px"
                               ),
                               width = 9
                             )
                           )),
                  tabPanel("Counts",
                           DT::dataTableOutput("table")),
                  tabPanel("Help", includeMarkdown("README.md"))
                ))

# server side stuff
server <- function(input, output, session) {
  # increase max upload size to 100 MB
  options(shiny.maxRequestSize = 100 * 1024 ^ 2)
  # ranges object for zooming plot
  ranges <- reactiveValues(x = NULL,
                           y = NULL)
  ids2Names <-
    reactiveValues() # this will contain genes and samples, both named character vectors
  selected <-
    reactiveValues() # this will contain genes, samples, all named character vectors
  clustered <- 
    reactiveValues() # this will contain genes, samples, all named character vectors
  
  normalisedCounts <- reactive({
    if (debug) {
      print('Function: normalisedCounts')
    }
    if (testing) {
      # testDataFile <- file.path(rootPath, 'data', 'DESeq.shield.testdata.RData')
      # load(testDataFile)
      if (input$dataType == 'rnaseq') {
        sampleFile <-
          file.path(rootPath, 'data', 'zfs-rnaseq-sampleInfo.tsv')
        countFile <-
          file.path(rootPath, 'data', 'counts.shield-subset.tsv')
      } else if (input$dataType == 'detct') {
        sampleFile <-
          file.path(rootPath, 'data', 'zfs-detct-sampleInfo.tsv')
        countFile <-
          file.path(rootPath, 'data', 'zfs-detct.subset.tsv')
      }
      exptData <-
        load_data(sampleFile, countFile, input$dataType, session)
    } else{
      dataFileInfo <- input$dataFile
      sampleFileInfo <- input$sampleFile
      countFileInfo <- input$countFile
      if (!is.null(sampleFileInfo) &
          !is.null(countFileInfo)) {
        exptData <-
          load_data(sampleFileInfo$datapath,
                          countFileInfo$datapath,
                          input$dataType,
                          session)
      } else if (!is.null(dataFileInfo)) {
        load(dataFileInfo$datapath)
      }
      else{
        return(NULL)
      }
    }
    
    # calculate the max mean value and set input slider
    maxMean <-
      ceiling(max(apply(assays(exptData)$norm_counts, 1, mean)))
    if (debug) {
      print(sprintf('Max mean value: %f', maxMean))
    }
    updateSliderInput(session, "maxMeanCount", value = maxMean, max = maxMean)
    
    if (input$dataType == 'rnaseq') {
      genes <- rowData(exptData)$Gene.name
      names(genes) <- rowData(exptData)$Gene.ID
    } else{
      genes <- rowData(exptData)$Gene.name
      names(genes) <- rownames(exptData)
    }
    samples <- rownames(colData(exptData))
    names(samples) <- rownames(colData(exptData))
    if (debug) {
      cat("Genes and Samples\n")
      print(sprintf('Initial genes: %s', paste0(head(genes), collapse = " ")))
      print(sprintf('Initial samples: %s', paste0(head(samples), collapse = " ")))
    }
    selected$genes <- genes
    selected$samples <- samples
    ids2Names$genes <- genes
    ids2Names$samples <- samples
    return(assays(exptData)$norm_counts)
  })
  
  subsetGeneList <- observe({
    geneIdsFileInfo <- input$geneIdsFile
    if (!is.null(geneIdsFileInfo)) {
      geneInfo <- read.table(geneIdsFileInfo$datapath)
      geneIds <- as.character(geneInfo[[1]])
      if (debug) {
        print('Function: subsetGeneList')
        print(sprintf(
          "Genes ids for subsetting: %s",
          paste0(head(geneIds), collapse = ", ")
        ))
      }
      ids <- check_for_missing_ids( geneIds, names(isolate(ids2Names$genes)) )
      Ids <- ids$ids
      nonexistentIds <- ids$missing_ids
      if (debug) {
        print(sprintf('Missing Ids: %s', paste(nonexistentIds, collapse = " ")))
        print(sprintf('Matched Ids: %s', paste(head(Ids), collapse = " ")))
      }
      # and warn
      if (length(nonexistentIds) > 0) {
        missing_genes_warning(nonexistentIds, session)
      }
      # set selected genes to ids
      selected$genes <- isolate(selected$genes[Ids])
      selected$currentSubset <- isolate(selected$genes[Ids])
    }
  })
  
  # respond to reset button
  observeEvent(input$subsetReset, {
    if (debug) {
      print('Function: subsetReset')
    }
    selected$genes <- ids2Names$genes
    selected$samples <- ids2Names$samples
    selected$currentSubset <- NULL
    reset("geneIdsFile")
  })
  
  # filter counts by expression level
  filteredCounts <- reactive({
    counts <- normalisedCounts()
    if (is.null(counts)) {
      return(NULL)
    } else {
      meanCount <- apply(counts, 1, mean)
      counts <- counts[meanCount >= input$minMeanCount &
                         meanCount <= input$maxMeanCount,]
      if (debug) {
        print('Function: filteredCounts')
        print(sprintf('Num Genes Filtered: %d', nrow(counts)))
      }
      # reset selected genes to everything in counts
      geneIds <- rownames(counts)
      selected$genes <- ids2Names$genes[geneIds]
      
      return(counts)
    }
  })
  
  # subset to genes/samples
  subsettedCounts <- reactive({
    counts <- filteredCounts()
    if (is.null(counts)) {
      return(NULL)
    } else {
      if (debug) {
        print('Function: subsettedCounts')
        print(sprintf(
          'selected Genes length = %d',
          length(isolate(reactiveValuesToList(selected)$genes))
        ))
        print(sprintf(
          'selected Samples length = %d',
          length(isolate(reactiveValuesToList(selected)$samples))
        ))
        print(sprintf('selected Gene Ids = %s', head(names(
          isolate(selected$genes)
        ))))
        print(sprintf('selected Gene Names = %s', head(isolate(selected$genes))))
        print(sprintf('selected Sample Ids = %s', head(names(
          isolate(selected$samples)
        ))))
        print(sprintf(
          'Filtered Counts dimensions: %d, %d',
          dim(counts)[1],
          dim(counts)[2]
        ))
      }
      numRow <- length(isolate(selected$genes))
      numCol <- length(isolate(selected$samples))
      # If the filtered set is too small go back to all genes/samples
      if (numRow > nrow(counts) | numCol > ncol(counts)) {
        counts <- normalisedCounts()
      }
      # check for missing ids
      ids <- check_for_missing_ids( names(isolate(selected$genes)), names(isolate(ids2Names$genes)) )
      Ids <- ids$ids
      nonexistentIds <- ids$missing_ids
      if (debug) {
        print(sprintf('Missing Ids: %s', paste(nonexistentIds, collapse = " ")))
        print(sprintf('Matched Ids: %s', paste(head(Ids), collapse = " ")))
      }
      # and warn
      if (length(nonexistentIds) > 0) {
        missing_genes_warning(nonexistentIds, session)
      }
      # set selected genes to ids
      selected$genes <- isolate(selected$genes[Ids])
      selected$currentSubset <- isolate(selected$genes[Ids])
      
      # If there is only one row or col the matrix gets simplified into a vector.
      # Needs to force it to stay as a matrix using drop = FALSE
      count <-
        tryCatch(
          counts[names(selected$genes), names(selected$samples), drop = FALSE],
          error = function(err) {
            if (debug) {
              print("Subsetting Error")
              print(err)
              print(sprintf('Selected Genes: %s', paste(
                names(selected$genes), collapse = " "
              )))
              print(sprintf('Selected Samples: %s', paste(
                names(selected$samples), collapse = " "
              )))
            }
            if (err$message == 'subscript out of bounds') {
              subsetErrorMsg <-
                'There was an error subsetting the counts matrix. Try resetting or changing the data type.'
              createAlert(
                session,
                "HeatmapAlert",
                "subsetErrorAlert",
                title = "Subsetting error",
                content = subsetErrorMsg,
                append = FALSE,
                style = 'danger'
              )
            } else{
              stop(err)
            }
          }
        )
      # if (numRow == 1 & numCol == 1) {
      #   count <- matrix(count,
      #                   dimnames = list(names(selected$genes), names(selected$samples)))
      # } else if (numRow == 1) {
      #   count <- matrix(count,
      #                   nrow = 1,
      #                   dimnames = list(names(selected$genes), names(selected$samples)))
      # } else if (numCol == 1) {
      #   count <- matrix(count,
      #                   ncol = 1,
      #                   dimnames = list(names(selected$genes), names(selected$samples)))
      # }
      if (debug) {
        print(sprintf(
          'Subset Counts dimensions: %d, %d',
          dim(count)[1],
          dim(count)[2]
        ))
      }
      return(count)
    }
  })
  
  # cluster counts by genes/samples
  clusteredCounts <- reactive({
    counts <- subsettedCounts()
    if (is.null(counts)) {
      return(NULL)
    }
    if (debug & length(input$clusterCheckGroup) > 0) {
      print('Function: clusteredCounts')
      print(sprintf('Cluster checkbox value: %s', paste0(isolate(input$clusterCheckGroup), collapse=" ") ) )
    }
    # get gene ids
    selection <- isolate(reactiveValuesToList(selected))
    err <- FALSE
    if (any(input$clusterCheckGroup == "1")) {
      if (any(input$clusterCheckGroup == "2")) {
        # check genes and samples for ones with zero variance
        zeroVar <- RemoveZeroVariance( counts, rows = TRUE, cols = TRUE )
        counts <- clusterMatrix(zeroVar$matrix, byRow = TRUE, byCol = TRUE)
        if( !is.null(zeroVar$rowsRemoved) | !is.null(zeroVar$colsRemoved) ){
          err <- TRUE
        }
      } else{
        # check genes for ones with zero variance
        zeroVar <- RemoveZeroVariance( counts, rows = TRUE, cols = FALSE )
        counts <- clusterMatrix(zeroVar$matrix, byRow = TRUE, byCol = FALSE)
        if( !is.null(zeroVar$rowsRemoved) ){
          err <- TRUE
        }
      }
      clustered$genes <- selection$genes[ rownames(counts) ]
      clustered$samples <- selection$samples[ colnames(counts) ]
    } else if (any(input$clusterCheckGroup == "2")) {
      # check samples for ones with zero variance
      zeroVar <- RemoveZeroVariance( counts, rows = FALSE, cols = TRUE )
      counts <- clusterMatrix(zeroVar$matrix, byRow = FALSE, byCol = TRUE)
      clustered$genes <- selection$genes[ rownames(counts) ]
      clustered$samples <- selection$samples[ colnames(counts) ]
      if( !is.null(zeroVar$colsRemoved) ){
        err <- TRUE
      }
    } else{
      clustered$samples <- NULL
      clustered$genes <- NULL
    }
    
    if(err){
      # first close any open alerts
      closeAlert(session, "zeroVarErrorAlert")
      
      # show a warning saying some genes/columns have been removed
      clusterErrorMsg <- ''
      if( !is.null(zeroVar$rowsRemoved) ){
        clusterErrorMsg <-
          paste0(
            'Some of the genes that you are trying to cluster have zero variance across the selected samples and have been removed: ',
            paste0(
              zeroVar$rowsRemoved,
              collapse = ", "
            ),
            '<br>'
          )
      }
      if( !is.null(zeroVar$colsRemoved) ){
        clusterErrorMsg <-
          paste0(
            clusterErrorMsg,
            'Some of the samples that you are trying to cluster have zero variance across the selected genes and have been removed: ',
            paste0(
              zeroVar$colsRemoved,
              collapse = ", "
            )
          )
      }
      
      createAlert(
        session,
        "HeatmapAlert",
        "zeroVarErrorAlert",
        title = "Clustering error",
        content = clusterErrorMsg,
        append = FALSE,
        style = 'danger'
      )
      if (debug) {
        print(sprintf(
          'Rows removed due to zero variance = %d, Num Genes left = %d',
          length(zeroVar$rowsRemoved),
          length(zeroVar$rowsKept) 
        ))
        print( sprintf(
          'Cols removed due to zero variance = %d, Num Samples left = %d',
          length(zeroVar$colsRemoved),
          length(zeroVar$colsKept)
        ))
      }
    }
    
    return(counts)
  })
  
  # transform counts. max scaled/log10/centred and scaled
  transformedCounts <- reactive({
    counts <- clusteredCounts()
    if (is.null(counts)) {
      return(NULL)
    } else {
      if (debug) {
        print('Function: transformedCounts')
        print(sprintf('Transform checkbox value: %s', input$transform))
      }
      if (input$transform == 2) {
        geneMaxCounts <- apply(counts, 1, max)
        # scale operates on the column so need to transpose, scale and then transpose back
        counts <-
          t(scale(t(counts), scale = geneMaxCounts, center = FALSE))
        # genes with a max of zero get converted to NAs
        # reset to zeros
        counts[ geneMaxCounts == 0, ] <- matrix( rep(0, sum(geneMaxCounts == 0)*ncol(counts) ), ncol = ncol(counts) )
      } else if (input$transform == 3) {
        counts <- log10(counts + 1)
      } else if (input$transform == 4) {
        # scale operates on the column so need to transpose, scale and then transpose back
        counts <- t(scale(t(counts)))
      }
      return(counts)
    }
  })
  
  # create plot object
  heatmapObj <- reactive({
    counts <- transformedCounts()
    print('Function: heatmapObj')
    if (is.null(counts)) {
      return(NULL)
    } else {
      plot <- ggplotExprHeatmap(counts)
      if (nrow(counts) <= 80) {
        genes <- ids2Names$genes[rownames(counts)]
        if (ncol(counts) <= 48) {
          plot <- plot +
            scale_y_discrete(labels = genes) +
            theme(
              axis.text.x = element_text(
                colour = "black",
                angle = 90,
                vjust = 0.5,
                hjust = 1,
                debug = FALSE
              ),
              axis.text.y = element_text(
                colour = "black",
                angle = 0,
                debug = FALSE
              )
            )
        } else{
          plot <- plot +
            scale_y_discrete(labels = genes) +
            theme(axis.text.y = element_text(
              colour = "black",
              angle = 0,
              debug = FALSE
            ))
        }
      } else if (ncol(counts) <= 48) {
        plot <- plot +
          theme(axis.text.x = element_text(
            colour = "black",
            angle = 90,
            vjust = 0.5,
            hjust = 1,
            debug = FALSE
          ))
      }
      if (input$transform == 4) {
        plot <- plot +
          scale_fill_distiller(type= 'div', palette = "RdBu")
      }
      return(plot)
    }
  })
  
  # render heatmap plot
  output$exprHeatmap <- renderPlot({
    return(heatmapObj())
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$heatmap_dblclick, {
    if (debug) {
      print('Function: heatmap dbl-click')
    }
    brush <- input$heatmap_brush
    if (debug) {
      print(
        sprintf(
          'Brush: Xmin: %s Xmax: %s Ymin: %s Ymax: %s',
          brush$xmin,
          brush$xmax,
          brush$ymin,
          brush$ymax
        )
      )
    }
    if (!is.null(brush)) {
      plotRanges <- list(x = floor(c(brush$xmin, brush$xmax) - c(0.5,-0.5)) + c(1,0),
                         y = floor(c(brush$ymin, brush$ymax) - c(0.5,-0.5)) + c(1,0))
      
      selection <- isolate(reactiveValuesToList(selected))
      clusteredSelection <- isolate(reactiveValuesToList(clustered))
      subset <- list( genes = NULL, samples = NULL )
      if( is.null(clusteredSelection$genes) | is.null(clusteredSelection$samples) ){
        subset$genes <- selection$genes
        subset$samples <- selection$samples
      } else{
        subset$genes <- clusteredSelection$genes
        subset$samples <- clusteredSelection$samples
      }
      
      if (debug) {
        print(sprintf('Plot Ranges X: %f %f', plotRanges$x[1], plotRanges$x[2]))
        print(sprintf('Plot Ranges Y: %f %f', plotRanges$y[1], plotRanges$y[2]))
        print(sprintf('Selected genes: %s', paste0(head(selection$genes), collapse=", ") ) )
        print(sprintf('Selected samples: %s', paste0(head(selection$samples), collapse=", ") ) )
        print(sprintf('Clustered genes: %s', paste0(head(clusteredSelection$genes), collapse=", ") ) )
        print(sprintf('Clustered samples: %s', paste0(head(clusteredSelection$samples), collapse=", ") ) )
        print(sprintf('Subset genes: %s', paste0(head(subset$genes), collapse=", ") ) )
        print(sprintf('Subset samples: %s', paste0(head(subset$samples), collapse=", ") ) )
      }

      # make sure values have not gone out of range
      # i.e. less than 1 or greater than num genes/samples in the current subset
      plotRanges$x[1] <-
        ifelse(plotRanges$x[1] < 1, 1, plotRanges$x[1])
      plotRanges$x[2] <-
        ifelse(
          plotRanges$x[2] > length(subset$samples),
          length(subset$samples),
          plotRanges$x[2]
        )
      plotRanges$y[1] <-
        ifelse(plotRanges$y[1] < 1, 1, plotRanges$y[1])
      plotRanges$y[2] <-
        ifelse(
          plotRanges$y[2] > length(subset$genes),
          length(subset$genes),
          plotRanges$y[2]
        )
      
      selectedGeneIds <-
        rev(rev(names(subset$genes))[seq(plotRanges$y[1], plotRanges$y[2])])
      selectedSampleIds <-
        names(subset$samples)[seq(plotRanges$x[1], plotRanges$x[2])]
      if (debug) {
        print(sprintf(
          'Num Genes = %d, Num Samples = %d',
          length(selectedGeneIds),
          length(selectedSampleIds)
        ))
      }
      
      # reset order to original order
      genePosition <- integer( length = length(selectedGeneIds) )
      if( is.null(selected$currentSubset) ){
        genes <- isolate(ids2Names$genes)
      } else{
        genes <- isolate(selected$currentSubset)
      }
      for( i in seq_len(length(selectedGeneIds)) ){
        geneId <- selectedGeneIds[i]
        genePosition[i] <- which( names(genes) == geneId )
      }
      selected$genes <- selection$genes[ selectedGeneIds[ order(genePosition) ] ]
      
      samplePosition <- integer( length = length(selectedSampleIds) )
      samples <- isolate(ids2Names$samples)
      for( i in seq_len(length(selectedSampleIds)) ){
        sampleId <- selectedSampleIds[i]
        samplePosition[i] <- which( names(samples) == sampleId )
      }
      selected$samples <- selection$samples[ selectedSampleIds[ order(samplePosition) ] ]
      
    } else {
      if( is.null(selected$currentSubset) ){
        selected$genes <- ids2Names$genes
      } else{
        selected$genes <- selected$currentSubset
      }
      selected$samples <- ids2Names$samples
    }
  })
  
  # for downloading the plot as a pdf/png
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste('geneExpr', Sys.Date(), input$plotFormat, sep = '.')
    },
    content = function(file) {
      heatmapPlot <- heatmapObj()
      # theme( axis.text = element_text( size = rel(0.5) ) )
      if (input$plotFormat == "pdf") {
        pdf(file,
            paper = "special",
            height = 16.5,
            width = 11.7) # open the pdf device
      } else if (input$plotFormat == "png") {
        png(file,
            height = 960,
            width = 480,
            res = 100) # open the png device
      }
      print(heatmapPlot)
      dev.off()  # close device
    }
  )
  
  # for downloading the counts file
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('geneExpr', Sys.Date(), 'tsv', sep = '.')
    },
    content = function(file) {
      write.table(
        clusteredCounts(),
        file,
        quote = FALSE,
        col.names = NA,
        sep = "\t"
      )
    },
    contentType = 'text/tsv'
  )
  # for downloading a gene list
  output$downloadGenes <- downloadHandler(
    filename = function() {
      paste('geneExpr', Sys.Date(), 'genes', 'tsv', sep = '.')
    },
    content = function(file) {
      write.table(
        rownames(clusteredCounts()),
        file,
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t"
      )
    },
    contentType = 'text/tsv'
  )
  
  # render the counts matrix as a table
  output$table <- DT::renderDataTable({
    data <- clusteredCounts()
    if (is.null(data)) {
      return(NULL)
    }
    DT::datatable(data)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
