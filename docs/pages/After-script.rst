
============================
What to do with your results
============================

Once the script is finished running you can use the results to perform Gene Set Enrichment Analysis (GSEA). 
An example could be to run a R script that takes the .csv file and performs the analysis like

.. code::

    if (!requireNamespace("BiocManager", quietly=TRUE)) 
    install.packages("BiocManager")

    BiocManager::install("GeneOverlap")

    library("GeneOverlap")

    x <- read.csv("~/Desktop/yourgenes.csv", sep=",") - assuming it contains one single column with the genes you want to test enrichment for (for example your Down's syndrome genes)

    PLS <- read.csv("~/Desktop/PLS1.csv", sep=",") - assuming you have two columns one for positively weighted (Pos) genes and another for negatively weighted (Neg) genes



    go.obj_pos <- newGeneOverlap(PLS$Pos,x, genome.size = 15633)

    PLS1_pos <- testGeneOverlap(go.obj_pos)



    go.obj_neg <- newGeneOverlap(PLS$Neg,x, genome.size = 15633)
    PLS1_neg <- testGeneOverlap(go.obj_neg)
