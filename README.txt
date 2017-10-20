This is a package for BC32 barcode analysis. See details in Roh et al. (2017).

1. Prepare a table for your samples
===================================

In addition to the sequencing file, it is suggested to provide an information textfile in your working directory.
This tab-delimited table (without header) named sampname.txt must contain 3 columns in the following order:
name (name of your sample), file (corresponding FASTQ file), index (multiplexing index sequence)
 
2. Compute barcodes summary for all samples
===========================================

To compute a BC32 barcode summary for each of your samples, follow the instructions in the Vignette.
You will also generate a merged summary table.

3. Explore barcodes distributions
=================================

The function exploreBC() starts a shinyapp for exploration of barcodes distribution and most common clones in your samples.
This app uses the tables generated in Step 2 as input.

Guide:
Select which samples you want to explore
Switch between barcode distribution and topclones graphs using the tab above the plot
Distribution mode:
    Adjust the minimum number of reads required to include a barcode in the graph
    Use the checkbox to switch between pooled or faceted graph
TopClones mode:
    Adjust the rank down to which barcodes are displayed (default: top 5)
Click 'Save Plot' to save the current plot to the working directory

See source code for details.