# This function preps the df
#' Title
#'
#' @param df.netPx1outgoing df
#'
#' @return a df for input to the heatmap
#'
#' @import reshape stats reshape2
#'
#' @examples
#' df <- prephm(df)
prephm <- function(df.netPx1outgoing) {
  dfoutgoing <- reshape::cast(df.netPx1outgoing, pathway_name~source, value = 'prob', fun.aggregate = 'sum')
  dfincoming<- reshape::cast(df.netPx1outgoing, pathway_name~target, value = 'prob', fun.aggregate = 'sum')

  rownames(dfoutgoing) <- dfoutgoing$pathway_name
  rownames(dfincoming) <- dfincoming$pathway_name
  dfall <- as.data.frame(merge(dfoutgoing, dfincoming, by = 'pathway_name'))
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
  return(dfall5)
}


# This function creates a heatmap
#' Title
#'
#' @param cellchatobj a cellchat object
#' @param col_fun color function
#' @param sources a list of cell types to be used as sender cells
#' @param targets a list of cell types to be used as receiver cells
#' @param fontsize a number
#' @param hmtitle a string for a title
#'
#' @return a list of heatmaps
#' @export
#'
#' @import ComplexHeatmap reshape2 reshape grid circlize
#' @examples
#'
#' hm <- makehm(cellchatobj = cellchatobj, col_fun = cols)
#' col_fun = colorRamp2(c(0, 0.005, 1), c("white", "darkseagreen2", "darkgreen"))
#' hm[[1]]
makehm <- function(cellchatobj,
                   col_fun,
                   sources,
                   targets,
                   fontsize,
                   hmtitle) {
  hmobj <- cellchatobj
  col_fun = col_fun
  size <- as.numeric(fontsize)
  title <- as.character(hmtitle)

  df.netPx <- reshape2::melt(hmobj@netP$prob, value.name = "prob")
  colnames(df.netPx)[1:3] <- c("source","target","pathway_name")

  #choose cells
  df.netPx1outgoing <- df.netPx[(df.netPx$source %in% sources) & (df.netPx$target %in% targets),]

  dfall4 <- as.matrix(prephm(df.netPx1outgoing))

  ht_allctrl = ComplexHeatmap::Heatmap(dfall4, name = hmtitle,
                                       top_annotation = HeatmapAnnotation(Ctrl = anno_barplot(as.numeric(colSums(dfall4)))),
                                       show_column_dend = FALSE,
                                       cluster_rows = FALSE,
                                       cluster_columns = FALSE,
                                       cluster_column_slices = FALSE,
                                       show_row_dend = FALSE,
                                       row_title = "Pathways",
                                       row_names_gp = grid::gpar(fontsize = size),
                                       column_split = paste0(c(rep("outgoing", length(sources)), rep("incoming", length(targets)-1))),
                                       column_names_gp = grid::gpar(fontsize = size),
                                       column_title_side="bottom",
                                       col = col_fun) +
    rowAnnotation(pathway = anno_barplot(rowSums(dfall4))) + rowAnnotation(rn = anno_text(rownames(dfall4), gp = grid::gpar(fontsize = size)))

  return(ht_allctrl)
}



