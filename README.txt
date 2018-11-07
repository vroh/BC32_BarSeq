This is a package for BC32 barcode analysis. See details in Roh et al. , Cell Reports 2018, 25:1-15.
Below is a brief description of the analysis pipeline.

1. Prepare a table for your samples
===================================

In addition to the sequencing files, it is recommended to provide an information textfile in your working directory.
This tab-delimited table (without header) named sampname.txt must contain 3 columns in the following order:
name (name of your sample), file (corresponding FASTQ file), index (multiplexing index sequence)
 
2. Compute barcodes summary for all samples
===========================================

To compute a BC32 barcode summary for each of your samples, follow the instructions in the Vignette.
This step will also generate a merged summary table.

3. Explore barcodes distributions
=================================

Use the function "exploreBC()" to start a shinyapp and explore the barcodes distribution and most common clones in
your samples.
This app uses the tables generated in Step 2 as input.

Quick guide:
------------
Select which samples you want to explore.
Switch between barcode distribution and topclones graphs using the tab above the plot.
Distribution mode:
    Adjust the minimum number of reads required to include a barcode in the graph.
    Use the checkbox to switch between pooled or faceted graph.
TopClones mode:
    Adjust the rank down to which barcodes are displayed (default: top 5).
Click 'Save Plot' to save the current plot to the working directory.

See source code for further details.