install_ez_packages <- function(answer = T) {
  if (answer == T) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkgs = c("tidyverse", "R.matlab", "tidyr", "rlang", "flowCore", "ks", "flowVS", "flowViz", "RColorBrewer", "gtools", "gplots","ggpubr","ggrastr", 
                           "ggplot2", "openxlsx", "samr", "lattice", "flowStats", "gdata", "Rtsne", "umap",
                           "FlowSOM", "dplyr", "plyr", "pryr", "doBy", "scales", "mixOmics", "reshape2", 
                           "plotly", "Rmisc", "Hmisc", 
                           "data.table", "EBImage", "magick", "phonTools", "heatmaply", "superheat", "esquisse", "robustbase","corrplot","scales","pheatmap","SCORPIUS"))
  }
  print("Done")
}

load_ez_packages <- function (answer = T) {
  if (answer == T) {
    library("tidyverse")
    library("R.matlab")
    library('tidyr')
    library("rlang")
    library("flowCore")
    library("ks")
    library("flowVS") #
    library("flowViz")
    library("RColorBrewer")
    library("gtools")
    library("gplots")
    library("ggplot2")
    library("ggpubr")
    library("ggrastr")
    library("openxlsx") #
    library("samr") #
    library("lattice")
    library("flowStats") #
    library("gdata")
    library("Rtsne")
    library("umap")
    library("FlowSOM") #
    library("dplyr")
    library('plyr')
    library("pryr")
    library("doBy") #
    library("scales")
    library("mixOmics") #
    library("reshape2")
    library("plotly") #
    library("Rmisc")
    library("Hmisc") #
    library("data.table")
    library("EBImage")
    library("magick")
    library("phonTools")
    library("heatmaply")
    library("superheat")
    library("esquisse")
    library("robustbase")
    library("corrplot")
    library("scales")
    library("pheatmap")
    library("SCORPIUS")
    
  
    #https://support.bioconductor.org/p/109128/ --> explains why use Biobase::exprs
    exprs = Biobase::exprs
    
  }
}
