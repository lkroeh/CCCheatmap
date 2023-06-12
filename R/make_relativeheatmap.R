#' Preps df for relative heatmap
#'
#' @param df.netPx1outgoing df
#'
#' @return df for heatmap
#' @import reshape stats
#'
#' @examples
preprelativehm <- function(df.netPx1outgoing) {
  dfoutgoing <- reshape::cast(df.netPx1outgoing, pathway_name~source, value = 'prob', fun.aggregate = 'sum')
  dfincoming<- reshape::cast(df.netPx1outgoing, pathway_name~target, value = 'prob', fun.aggregate = 'sum')

  rownames(dfoutgoing) <- dfoutgoing$pathway_name
  rownames(dfincoming) <- dfincoming$pathway_name
  dfall <- as.data.frame(merge(dfoutgoing, dfincoming))
  rownames(dfall) <- dfall$pathway_name
  dfall1 <- dfall[,-1]

  #exclude low probability interactions
  dfall2 <- stats::na.omit(dfall1)
  zeros = apply(dfall2, 1, function(row) all(row ==0 ))
  dfall3 <- dfall2[!zeros,]
  #pt1 = which(rowSums(dfall2)>=0.1)
  #dfall3 <- dfall2[pt1,]

  #order by decreasing y axis
  order <- order(rowSums(dfall3), decreasing = T)
  dfall5 <- dfall3[order(-rowSums(dfall3)),]

  mat <- dfall5
  mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
  dfall4 <- mat
  test1 <- list()
  test1[[1]] <- dfall4
  test1[[2]] <- dfall5
  return(test1)
}


#' This function makes a heatmap with relative pathway singnaling strength
#'
#' @param cellchatobj a cellchat object
#' @param col_fun color function
#' @param sources list of cell types that are sender cells
#' @param targets list of cell types that are receiver cells
#' @param fontsize numeric
#' @param hmtitle string as title
#'
#' @return a heatmap in which row (pathway) values are scaled between 0 and 1.
#' @export
#' @import ComplexHeatmap reshape2 circlize
#' @examples
makerelativehm <- function(cellchatobj,
                           col_fun,
                           sources,
                           targets,
                           fontsize,
                           hmtitle) {
  hmobj <- cellchatobj

  if(!is.na(col_fun)) {
    col_fun = cols
  } else {
    col_fun <- colorRamp2(c(0, 0.005, 1), c("white", "darkseagreen2", "darkgreen"))
  }

  size <- as.numeric(fontsize)
  title <- as.character(hmtitle)

  #add if not centr, run compute centrality
  df.netPx <- reshape2::melt(hmobj@netP$centr, value.name = "outdeg")
  colnames(df.netPx)[1:3] <- c("source","target","pathway_name")

  #choose cells
  df.netPx1outgoing <- df.netPx[(df.netPx$source %in% sources) & (df.netPx$target %in% targets),]

  dfall4 <- as.matrix(prephm(df.netPx1outgoing))

  ht_allctrl = ComplexHeatmap::Heatmap(dfall4[[1]], name = hmtitle,
                                       top_annotation = HeatmapAnnotation(Ctrl = anno_barplot(as.numeric(colSums(dfall4[[2]])))),
                                       show_column_dend = FALSE,
                                       cluster_rows = FALSE,
                                       cluster_columns = FALSE,
                                       cluster_column_slices = FALSE,
                                       show_row_dend = FALSE,
                                       row_title = "Pathways",
                                       row_names_gp = gpar(fontsize = size),
                                       column_split = paste0(c(rep("outgoing,", length(sources)), rep("incoming,", length(targets)-1), rep("incoming", 1))),
                                       column_names_gp = grid::gpar(fontsize = size),
                                       column_title_side="bottom",
                                       col = col_fun) +
    rowAnnotation(pathway = anno_barplot(rowSums(dfall4[[2]]))) + rowAnnotation(rn = anno_text(rownames(dfall4[[2]]), gp = gpar(fontsize = size)))

  return(ht_allctrl)
}
