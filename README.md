# geneExpr

A shiny app for visualising gene expression data.

## Instructions

The app is designed to display gene expression data produced by RNASeq.

### Loading data

Data can be loaded using the buttons at the top of the side bar control panel. Users can supply a count file along with a sample file.

#### Count File
The count file should be a tab-separated file containing genes in rows and samples in columns.

e.g.
<table>
  <tr>
    <th>Gene ID</th>
    <th>Gene Name</th>
    <th>012536 count</th>
    <th>012537 count</th>
    <th>012538 count</th>
    <th>012525 count</th>
    <th>012526 count</th>
    <th>012527 count</th>
  </tr>
  <tr>
    <td>ENSDARG00000000212</td>
    <td>krt97</td>
    <td>10</td>
    <td>12</td>
    <td>15</td>
    <td>35</td>
    <td>42</td>
    <td>23</td>
  </tr>
  <tr>
    <td>ENSDARG00000000567</td>
    <td>znf281a</td>
    <td>345</td>
    <td>333</td>
    <td>325</td>
    <td>357</td>
    <td>365</td>
    <td>354</td>
  </tr>
</table>

#### Sample File

The sample file is a tab-separated file containg at least two columns. The required columns are the sample names that match the sample column names in the count file. This must be the first column. The other required column should be labelled "condition" and details how the samples are divided into groups (e.g. Control and Treated). If the file contains a column labelled "sampleName" these names will be displayed on the heatmap.

e.g.
<table>
  <tr>
    <th></th>
    <th>condition</th>
    <th>sampleName</th>
    <th>batch</th>
  </tr>
  <tr>
    <td>012536</td>
    <td>Control</td>
    <td>Ctrl1</td>
    <td>A</td>
  </tr>
  <tr>
    <td>012537</td>
    <td>Control</td>
    <td>Ctrl2</td>
    <td>B</td>
  </tr>
  <tr>
    <td>012538</td>
    <td>Control</td>
    <td>Ctrl3</td>
    <td>C</td>
  </tr>
  <tr>
    <td>012525</td>
    <td>Treated</td>
    <td>Trt1</td>
    <td>A</td>
  </tr>
  <tr>
    <td>012526</td>
    <td>Treated</td>
    <td>Trt2</td>
    <td>B</td>
  </tr>
  <tr>
    <td>012527</td>
    <td>Treated</td>
    <td>Trt3</td>
    <td>C</td>
  </tr>
</table>

#### Subset by Gene id

A text file of Ensembl gene ids to subset the heatmap to can be uploaded using this button. Any ids that can't be matched will be listed in a warning alert above the heatmap.

### Heatmap

The default plot is a heatmap showing all the genes and all the samples in the count data. Users can zoom in to a section of the plot by selecting an area of the plot and double-clicking on it. Double-clicking without highlighting resets the view back to all genes and all samples.

#### Transform counts

The counts can be transformed as follows

* Raw - The untransformed normalised counts (Default)
* Max Scaled - Each row is scaled to the maximum value for that row.
* log10

#### Clustering

Both the rows (genes) and columns (samples) of the heatmap can be clustered. Currently, the clustering is done by hierarchical clustering of the Pearson correlation coefficients between genes/samples. Only the genes/samples currently displayed in the heatmap are used for the clustering

#### Filter Genes

The genes can be filtered according to the mean counts across all samples. Both minimum and maximum thresholds are available to set with the sliders.

#### Downloads

The **Download plot** button downloads the current plot as either a pdf or png.

The **Download Count File** button downloads the count data for the genes and samples displayed in the heatmap.

The **Download Gene List** button downloads a file containing just the gene ids for the genes and samples displayed in the heatmap.

#### Load RData file

As an alternative to sample/count files, a pre-computed R data file containing a DESeq2 DESeqDataSet object (named DESeqData) can be uploaded using the **Load .RData File** button 

## Prerequisites

* shiny
* shinyjs
* shinyBS
* ggplot2
* DESeq2
* DT
* [rnaseqVis](https://github.com/richysix/rnaseqVis)
* rprojroot
* devtools
