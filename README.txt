This is the final protocol for BC32 barcode analysis. See details in Roh et al. (2017).
Files in toy_data folder can be used to test the scripts.

1. Prepare a table for your samples
===================================

This tab-delimited table (without header) must be named sampname.txt and contain 3 columns in the following order:
name (name of your sample), file (corresponding FASTQ file), index (multiplexing index sequence)
 
2. Compute barcodes summary for all samples
===========================================

To compute a BC32 barcode summary for each of your samples, run the barcode_processing.R script.
This will also generate a merged summary table.

3. Explore barcodes distributions
=================================

barcode_explorer.R is a shinyapp for exploration of barcodes distribution and most common clones in your samples.
This app uses the output from Step 2 (barcode_processing.R) as input.

Guide:
Select which samples you want to explore
Switch between barcode distribution and topclones graphs using the tab above the plot
Distribution mode:
    Adjust the minimum number of reads required to include a bacorde in the graph
    Use the checkbox to switch between pooled or faceted graph
TopClones mode:
    Adjust the rank down to which barcodes are displayed (default: top 5)

See barcode_explorer.R for details.