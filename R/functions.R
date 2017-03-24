#' FilesToDESeqObj
#'
#' \code{FilesToDESeqObj} returns a DESeqDataSet object produced from the supplied sample and count files 
#'
#'    The first column of the sample file must match the sample names in the column names of the count file
#'    The sample file must also have a column labelled "condition" for the DESeq2 design formula.
#'    If the samples file contains a column labelled "sampleName" the samples will be renamed.
#'    
#' @param sampleFile character - name of the sample file
#' @param countFile character - name of the count file
#'
#' @return a DESeqDataSet object
#'
#' @examples
#' FilesToDESeqObj( 'samples.txt', 'all.tsv' )
#'
#' @export
#'
FilesToDESeqObj <- function( sampleFile, countFile, dataType, session ){
  # read in sample info
  sampleInfo <- read.table(sampleFile, sep="\t", header=TRUE, row.names = 1)
  # Read in count data
  data <- read.delim(countFile, header=TRUE, check.names=FALSE)
  
  # Support different column names
  names(data)[names(data) == 'chr']     <- 'Chr'
  names(data)[names(data) == '#Chr']     <- 'Chr'
  names(data)[names(data) == 'start']   <- 'Start'
  names(data)[names(data) == 'end']     <- 'End'
  names(data)[names(data) == 'Gene ID']      <- 'Gene.ID'
  names(data)[names(data) == 'ID']      <- 'Gene.ID'
  names(data)[ grepl("e[0-9][0-9] Ensembl Gene ID", names(data), perl=TRUE ) ]    <- 'Gene.ID'
  names(data)[names(data) == 'Name']      <- 'Gene.name'
  names(data)[names(data) == 'Gene name']      <- 'Gene.name'
  names(data)[names(data) == 'adjpval'] <- 'adjp'
  
  # Get count data
  countData <- data[,grepl(" count$", names(data)) &
                      !grepl(" normalised count$", names(data)) &
                      !grepl(" end read count$", names(data)) ]
  names(countData) <- gsub(" count$", "", names(countData))
  if( dataType == 'rnaseq' ){
    rownames(countData) <- data[ , "Gene.ID" ]
  } else {
    rownames(countData) <- paste( data[['Chr']], data[['Region start']], data[['Region end']], data[["3' end strand"]], sep=":")
  }
  
  # check column names of countData match rownames of sample info
  # check for samples that don't exist
  missingSamples <- vector('list', length = length(rownames(sampleInfo)))
  Samples <- vector('list', length = length(rownames(sampleInfo)))
  for (i in seq_len(length(rownames(sampleInfo)))) {
    if (sum(colnames(countData) == rownames(sampleInfo)[i]) == 0) {
      missingSamples[[i]] <- rownames(sampleInfo)[i]
      Samples[[i]] <- NULL
    } else{
      missingSamples[i] <- NULL
      Samples[[i]] <- rownames(sampleInfo)[i]
    }
  }
  missingSamples <- do.call(c, missingSamples)
  Samples <- do.call(c, Samples)
  # and warn
  if (length(missingSamples) > 0) {
    missingSamplesWarning <-
      paste0(
        "Some of the sample Ids couldn't be matched in the count file header! <br>Ids: ",
        paste(missingSamples, collapse = ", ")
      )
    shinyBS::createAlert(
      session,
      "HeatmapAlert",
      "sampleIdsAlert",
      title = "Non-matching Sample Ids",
      content = missingSamplesWarning,
      append = FALSE,
      style = 'warning'
    )
  }
  # reorder columns of countData based on sample names
  countData <- countData[ , Samples ]
  # and subset sample data to Samples
  sampleInfo <- sampleInfo[ Samples, ]

  # change sample names if sample names present in sample file
  if( any( colnames( sampleInfo ) == 'sampleName' ) ){
    colnames(countData) <- sampleInfo$sampleName
    rownames(sampleInfo) <- sampleInfo$sampleName
    if( sum(!grepl("sampleName", colnames(sampleInfo) ) ) == 1 ){
      samples <- data.frame( condition = sampleInfo[ , !grepl("sampleName", colnames(sampleInfo) ) ],
                             row.names = rownames(sampleInfo) )
    } else {
      samples <- sampleInfo[ , !grepl("sampleName", colnames(sampleInfo) ) ]
    }
  } else {
    samples <- sampleInfo
  }
  
  # make DESeq2 object
  DESeqData <- DESeq2::DESeqDataSetFromMatrix(countData, samples, design = ~ condition)
  DESeqData <- DESeq2::estimateSizeFactors(DESeqData)
  # add metadata
  featureData <- data[ , !grepl("count", colnames(data) ) ]
  mcols(DESeqData) <- DataFrame(mcols(DESeqData), featureData)
  
  return( DESeqData )
}

#' RemoveZeroVariance
#'
#' \code{RemoveZeroVariance} removes rows/columns from a matrix that have zero variance
#'
#'    
#'    
#' @param x - matrix
#' @param rows - logical, default TRUE
#' @param cols - logical, default FALSE
#'
#' @return a list containing the following elements
#' 
#'    matrix - the supplied matrix with the rows/columns removed
#'    
#'    rowsKept - a character vector with the names of the rows that were not removed
#'      NULL if rows = FALSE
#'    
#'    rowsRemoved - a character vector with the names of the rows that were removed
#'      NULL if rows = FALSE
#'      
#'    colsKept - a character vector with the names of the columns that were not removed
#'      NULL if cols = FALSE
#'      
#'    colsRemoved - a character vector with the names of the columns that were removed
#'      NULL if cols = FALSE
#'
#' @examples
#' countMatrix <- matrix( sample(1:100, 100), ncol = 10)
#' countMatrix[ 3, ] <- rep(0,10)
#' countMatrix[ 6, ] <- rep(0,10)
#' RemoveZeroVariance( countMatrix, rows = TRUE, cols = TRUE )
#'
#' @export
#'
RemoveZeroVariance <- function( x, rows = TRUE, cols = FALSE ){
  # check matrix has row and col names
  if( rows & is.null(rownames(x)) ){
    stop("The supplied matrix does not have row names!")
  } else if( cols & is.null(colnames(x)) ){
    stop("The supplied matrix does not have column names!")
  }
  rowsKept <- NULL
  rowsRemoved <- NULL
  colsKept <- NULL
  colsRemoved <- NULL
  if( rows ){
    zeroVarRows <- genefilter::rowSds(x) == 0
    if (sum(zeroVarRows) > 0) {
      rowsKept <- rownames(x)[ !zeroVarRows ]
      rowsRemoved <- rownames(x)[ zeroVarRows ]
      x <- x[ !zeroVarRows, ]
    }
  } else{
    rowsKept <- rownames(x)
  }
  if( cols ){
    zeroVarCols <- genefilter::rowSds(t(x)) == 0
    if (sum(zeroVarCols) > 0) {
      colsKept <- colnames(x)[ !zeroVarCols ]
      colsRemoved <- colnames(x)[ zeroVarCols ]
      x <- x[ , !zeroVarCols ]
    }
  } else{
    colsKept <- colnames(x)
  }
  return(
    list(
      matrix = x,
      rowsKept = rowsKept,
      rowsRemoved = rowsRemoved,
      colsKept = colsKept,
      colsRemoved = colsRemoved
    )
  )
}

