ReadXenium2 <- function (data.dir, outs = c("matrix", "microns"), type = "centroids", 
                         mols.qv.threshold = 20) 
{
  library(progressr)
  
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"), 
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"), 
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      pmtx <- progressor()
      pmtx(message = "Reading counts matrix", class = "sticky", 
           amount = 0)
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir, 
                                                              "cell_feature_matrix/")))
      pmtx(type = "finish")
      matrix
    }, centroids = {
      pcents <- progressor()
      pcents(message = "Loading cell centroids", class = "sticky", 
             amount = 0)
      if (has_dt) {
        cell_info <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                               "cells.csv.gz")))
      } else {
        cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
      }
      cell_centroid_df <- data.frame(x = cell_info$x_centroid, 
                                     y = cell_info$y_centroid, cell = cell_info$cell_id, 
                                     stringsAsFactors = FALSE)
      pcents(type = "finish")
      cell_centroid_df
    }, segmentations = {
      psegs <- progressor()
      psegs(message = "Loading cell segmentations", class = "sticky", 
            amount = 0)
      if (has_dt) {
        cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                                        "cell_boundaries.csv.gz")))
      } else {
        cell_boundaries_df <- read.csv(file.path(data.dir, 
                                                 "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
      }
      names(cell_boundaries_df) <- c("cell", "x", "y")
      psegs(type = "finish")
      cell_boundaries_df
    }, microns = {
      pmicrons <- progressor()
      pmicrons(message = "Loading molecule coordinates", 
               class = "sticky", amount = 0)
      
      transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
      transcripts <- subset(transcripts, qv >= mols.qv.threshold)
      
      df <- data.frame(x = transcripts$x_location, y = transcripts$y_location, 
                       gene = transcripts$feature_name, stringsAsFactors = FALSE)
      pmicrons(type = "finish")
      df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)
  return(data)
}
