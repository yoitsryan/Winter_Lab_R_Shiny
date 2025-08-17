# Winter Lab R Shiny app
Here is an R Shiny app created by the lab of Dr. Deborah Winter at Northwestern University: a database of bulk RNA-seq data. I was tasked with debugging it. The code was highly un-annotated, making the fixes difficult. Here are the changes I made:

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
