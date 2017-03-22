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
FilesToDESeqObj <- function( sampleFile, countFile, dataType ){
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
  if( length(colnames(countData)) != length(rownames(sampleInfo)) ){
      stop("Number of columns in counts file does not match number of rows in samples file")
  } else if( !identical(colnames(countData), rownames(sampleInfo)) ){
    # order count data in the same order as the samples
    sampleOrder <- sapply( rownames(sampleInfo), 
                           function(x){ 
                             Idx <- which( colnames(countData) == x )
                             if( length(Idx) == 0 ){
                               stop( sprintf('Could not find sample, %s, in count data columns!', x) )
                             }
                             return(Idx)
                           } )
    countData <- countData[ , sampleOrder ]
  }
  
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
  DESeqData <- DESeqDataSetFromMatrix(countData, samples, design = ~ condition)
  DESeqData <- estimateSizeFactors(DESeqData)
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
    zeroVarRows <- rowSds(x) == 0
    if (sum(zeroVarRows) > 0) {
      rowsKept <- rownames(x)[ !zeroVarRows ]
      rowsRemoved <- rownames(x)[ zeroVarRows ]
      x <- x[ !zeroVarRows, ]
    }
  }
  if( cols ){
    zeroVarCols <- rowSds(t(x)) == 0
    if (sum(zeroVarCols) > 0) {
      colsKept <- colnames(x)[ !zeroVarCols ]
      colsRemoved <- colnames(x)[ zeroVarCols ]
      x <- x[ , !zeroVarCols ]
    }
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

