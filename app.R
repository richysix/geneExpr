# load libraries
for (package in c('shiny',
                  'shinyjs',
                  'shinyBS',
                  'ggplot2',
                  'DESeq2',
                  'rnaseqVis',
                  'rprojroot',
                  'DT')) {
  library(package, character.only = TRUE)
}
source('functions.R')

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
                                   "log 10" = 3
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
    reactiveValues() # this will contain genes and samples, both named character vectors
  
  DeSeqCounts <- reactive({
    if (debug) {
      print('Function: DeSeqCounts')
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
      DESeqData <-
        FilesToDESeqObj(sampleFile, countFile, input$dataType)
    } else{
      dataFileInfo <- input$dataFile
      sampleFileInfo <- input$sampleFile
      countFileInfo <- input$countFile
      if (!is.null(sampleFileInfo) &
          !is.null(countFileInfo)) {
        DESeqData <-
          FilesToDESeqObj(sampleFileInfo$datapath,
                          countFileInfo$datapath,
                          input$dataType)
      } else if (!is.null(dataFileInfo)) {
        load(dataFileInfo$datapath)
      }
      else{
        return(NULL)
      }
    }
    
    # calculate the max mean value and set input slider
    maxMean <-
      ceiling(max(apply(
        counts(DESeqData, normalized = TRUE), 1, mean
      )))
    if (debug) {
      print(sprintf('Max mean value: %f', maxMean))
    }
    updateSliderInput(session, "maxMeanCount", value = maxMean, max = maxMean)
    
    if (input$dataType == 'rnaseq') {
      genes <- as.character(rowData(DESeqData)$Gene.name)
      names(genes) <- as.character(rowData(DESeqData)$Gene.ID)
    } else{
      genes <- as.character(rowData(DESeqData)$Gene.name)
      names(genes) <- rownames(DESeqData)
    }
    samples <- rownames(colData(DESeqData))
    names(samples) <- rownames(colData(DESeqData))
    if (debug) {
      cat("Genes and Samples\n")
      print(sprintf('Initial genes: %s', paste0(head(genes), collapse = " ")))
      print(sprintf('Initial samples: %s', paste0(head(samples), collapse = " ")))
    }
    selected$genes <- genes
    selected$samples <- samples
    ids2Names$genes <- genes
    ids2Names$samples <- samples
    return(DESeqData)
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
      # check for any ids that don't exist
      nonexistentIds <- vector('list', length = length(geneIds))
      Ids <- vector('list', length = length(geneIds))
      for (i in seq_len(length(geneIds))) {
        if (sum(names(ids2Names$genes) == geneIds[i]) == 0) {
          nonexistentIds[[i]] <- geneIds[i]
          Ids[[i]] <- NULL
        } else{
          nonexistentIds[i] <- NULL
          Ids[[i]] <- geneIds[i]
        }
      }
      nonexistentIds <- do.call(c, nonexistentIds)
      Ids <- do.call(c, Ids)
      if (debug) {
        print(sprintf('Missing Ids: %s', paste(nonexistentIds, collapse = " ")))
        print(sprintf('Matched Ids: %s', paste(head(Ids), collapse = " ")))
      }
      # and warn
      if (length(nonexistentIds) > 0) {
        missingGenesWarning <-
          paste0(
            "Some of the gene ids couldn't be matched! Ids: ",
            paste(nonexistentIds, collapse = ", ")
          )
        createAlert(
          session,
          "HeatmapAlert",
          "geneIdsAlert",
          title = "Non-matching Ids",
          content = missingGenesWarning,
          append = FALSE,
          style = 'warning'
        )
      }
      # set selected genes to ids
      selected$genes <- isolate(selected$genes[Ids])
    }
  })
  
  # respond to reset button
  observeEvent(input$subsetReset, {
    if (debug) {
      print('Function: subsetReset')
    }
    selected$genes <- ids2Names$genes
    selected$samples <- ids2Names$samples
    reset("geneIdsFile")
  })
  
  # retrieve normalised counts
  normalisedCounts <- reactive({
    counts <- DeSeqCounts()
    if (debug) {
      print('Function: normalisedCounts')
      print(sprintf('Num Genes: %d', nrow(counts)))
    }
    if (is.null(counts)) {
      return(NULL)
    } else {
      return(counts(DeSeqCounts(), normalized = TRUE))
    }
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
          length(reactiveValuesToList(selected)$genes)
        ))
        print(sprintf(
          'selected Samples length = %d',
          length(reactiveValuesToList(selected)$samples)
        ))
        print(sprintf('selected Gene Ids = %s', head(names(
          selected$genes
        ))))
        print(sprintf('selected Gene Names = %s', head(selected$genes)))
        print(sprintf('selected Sample Ids = %s', head(names(
          selected$samples
        ))))
        print(sprintf(
          'Filtered Counts dimensions: %d, %d',
          dim(counts)[1],
          dim(counts)[2]
        ))
      }
      numRow <- length(selected$genes)
      numCol <- length(selected$samples)
      # If the filtered set is too small go back to all genes/samples
      if (numRow > nrow(counts) | numCol > ncol(counts)) {
        counts <- normalisedCounts()
      }
      # If there is only one row or col the matrix gets simplified into a vector.
      # Needs to force it to stay as a matrix
      count <-
        tryCatch(
          counts[names(selected$genes), names(selected$samples)],
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
      if (numRow == 1 & numCol == 1) {
        count <- matrix(count,
                        dimnames = list(names(selected$genes), names(selected$samples)))
      } else if (numRow == 1) {
        count <- matrix(count,
                        nrow = 1,
                        dimnames = list(names(selected$genes), names(selected$samples)))
      } else if (numCol == 1) {
        count <- matrix(count,
                        ncol = 1,
                        dimnames = list(names(selected$genes), names(selected$samples)))
      }
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
      print(sprintf('Cluster checkbox value: %s', input$clusterCheckGroup))
    }
    if (any(input$clusterCheckGroup == "1")) {
      # check the genes for ones were sd is zero
      zeroVarRows <- rowSds(counts) == 0
      if (sum(zeroVarRows) > 0) {
        # remove the rows that have zero variance
        selection <- reactiveValuesToList(selected)
        selectedGeneIds <- names(selection$genes)[!zeroVarRows]
        # show a warning saying some rows have been removed
        clusterErrorMsg <-
          paste0(
            'Some of the genes that you are trying to cluster have zero variance across the selected samples and have been removed: ',
            paste(
              names(selection$genes)[zeroVarRows],
              sep = "",
              collapse = ", "
            )
          )
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
            sum(zeroVarRows),
            sum(!zeroVarRows)
          ))
        }
        selected$genes <- selection$genes[selectedGeneIds]
      }
      if (any(input$clusterCheckGroup == "2")) {
        counts <- clusterMatrix(counts, byRow = TRUE, byCol = TRUE)
      } else{
        counts <- clusterMatrix(counts, byRow = TRUE, byCol = FALSE)
      }
    } else if (any(input$clusterCheckGroup == "2")) {
      counts <- clusterMatrix(counts, byRow = FALSE, byCol = TRUE)
    }
    return(counts)
  })
  
  # transform counts. max scaled/log10
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
      } else if (input$transform == 3) {
        counts <- log10(counts + 1)
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
      plotRanges <- list(x = floor(c(brush$xmin, brush$xmax) - c(0.5,-0.5)),
                         y = floor(c(brush$ymin, brush$ymax) - c(0.5,-0.5)))
      if (debug) {
        print(sprintf('Plot Ranges X: %f %f', plotRanges$x[1], plotRanges$x[2]))
        print(sprintf('Plot Ranges Y: %f %f', plotRanges$y[1], plotRanges$y[2]))
      }
      selection <- reactiveValuesToList(selected)
      # make sure values have not gone out of range
      # i.e. less than 1 or greater than num genes/samples in the current subset
      plotRanges$x[1] <-
        ifelse(plotRanges$x[1] < 1, 1, plotRanges$x[1])
      plotRanges$x[2] <-
        ifelse(
          plotRanges$x[2] > length(selection$samples),
          length(selection$samples),
          plotRanges$x[2]
        )
      plotRanges$y[1] <-
        ifelse(plotRanges$y[1] < 1, 1, plotRanges$y[1])
      plotRanges$y[2] <-
        ifelse(
          plotRanges$y[2] > length(selection$genes),
          length(selection$genes),
          plotRanges$y[2]
        )
      
      selectedGeneIds <-
        rev(rev(names(selection$genes))[seq(plotRanges$y[1], plotRanges$y[2])])
      selectedSampleIds <-
        names(selection$samples)[seq(plotRanges$x[1], plotRanges$x[2])]
      if (debug) {
        print(sprintf(
          'Num Genes = %d, Num Samples = %d',
          length(selectedGeneIds),
          length(selectedSampleIds)
        ))
      }
      
      selected$genes <- selection$genes[selectedGeneIds]
      selected$samples <- selection$samples[selectedSampleIds]
    } else {
      selected$genes <- ids2Names$genes
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
