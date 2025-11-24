# Winter Lab Bulk RNA-seq app
This is an R Shiny app created by the lab of Dr. Deborah Winter at Northwestern University: a database of bulk RNA-seq data. I was tasked with debugging it. The code was highly un-annotated, making the fixes difficult. Here are the changes I made:

✅ Input files no longer case sensitive

✅ Reorganized heatmap customization options:

      - Made "No Clustering" the default option in the Clustering Method dropdown
      
      - Group Or Samples dropdown moved above Clustering Method dropdown
      
✅ Made heatmap smaller

✅ Implemented the ability to open heatmap tab before gene selection

✅ Moved all 3 overarching tabs ("Introduction", "Data Visualization", "Gene Expression Heatmaps") next to sidebar

✅ Implemented deletion of heatmap upon deselection of gene

✅ Sped up loading time for gene selection

✅ Eliminated the freezing bug when a second gene is selected before the first gene finishes loading

✅ Repaired the glitch for box plots where a sub-tab isn't always created when a new dataset is selected

✅ Eliminated the crash that happens when a dataset is selected, then deselected


# Winter Lab Single Cell RNA-seq app
In this app, you upload an expression matrix and an accompanying metadata file. You select one or more genes, and can view one of 5 different plot types: feature plot (UMAP), bubble plot, violin plot, heat map, box plot. You can choose metadata variables to group or split by.

The most major improvement I made to this app was to convert the code, which was written in ggplot, to plotly. Other improvements include:

- speeding up slow runtime

- fixing multiple plot types that initially failed to run

- the addition of a second heatmap to visualize pseudobulked data.

Once the plots were all functional, I made more edits:

- min-max scaling the expression levels represented in the heatmaps

- fixing the split feature plots (bottom panel on the Feature Plot tab) so they don't become squashed or stretched

- fixing the shrinking violins bug by only allowing one violin plot per row instead of two per row

- swapping the axes of the bubble plots (selected genes on y-axis, group variable on x-axis), and creating one plot per split variable instead of one plot per selected gene