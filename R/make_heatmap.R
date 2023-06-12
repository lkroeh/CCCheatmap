# This function preps the df
#' Title
#'
#' @param df.netPx1outgoing df
#'
#' @return a df for input to the heatmap
#'
#' @import reshape2 stats
#'
#' @examples
#' df <- prephm(df)
prephm <- function(df.netPx1outgoing) {
  dfoutgoing <- reshape2::cast(df.netPx1outgoing, pathway_name~source, value = 'prob', fun.aggregate = 'sum')
  dfincoming<- reshape2::cast(df.netPx1outgoing, pathway_name~target, value = 'prob', fun.aggregate = 'sum')

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
#' @import ComplexHeatmap reshape2
#' @examples
#'
#' hm <- makehm(cellchatobj = cellchatobj, col_fun = cols)
#' hm[[1]]
makehm <- function(cellchatobj = cellchatobj,
                   col_fun = cols,
                   sources = c(),
                   targets = c(),
                   fontsize = 4,
                   hmtitle = "title") {
  hmobj <- cellchatobj
  col_fun = cols
  df.netPx <- reshape2::melt(hmobj@netP$prob, value.name = "prob")
  colnames(df.netPx)[1:3] <- c("source","target","pathway_name")

  #choose cells
  lensources <- length(sources)
  df.netPx1outgoing <- df.netPx[(df.netPx$source==c('iCAF') | df.netPx$source==c('Lymphocyte') | df.netPx$source==c('Endothelial')) & (df.netPx$target=='CTL Lum ER-' | df.netPx$target=='CTL Lum ER+' | df.netPx$target=='CTL Basal'),]

  dfall4 <- prephm(df.netPx1outgoing)

  ht_allctrl = ComplexHeatmap::Heatmap(dfall4, name = title,
                       top_annotation = HeatmapAnnotation(Ctrl = anno_barplot(as.numeric(colSums(dfall4)))),
                       show_column_dend = FALSE,
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       cluster_column_slices = FALSE,
                       show_row_dend = FALSE,
                       row_title = "Pathways",
                       row_names_gp = gpar(fontsize = fontsize),
                       column_split = paste0(c(rep("outgoing,", 3), rep("incoming,", length(targets-1)), rep("incoming", 1))),
                       column_names_gp = grid::gpar(fontsize = 8),
                       column_title_side="bottom",
                       col = col_fun) +
    rowAnnotation(pathway = anno_barplot(rowSums(dfall4))) + rowAnnotation(rn = anno_text(rownames(dfall4), gp = gpar(fontsize = fontsize)))

  return(ht_allctrl)
}


#if want other relative hm
#add if not centr, run compute centrality
#df.netPx1 <- reshape2::melt(hmobj@netP$centr, value.name = "outdeg")

#mat <- dfall5
#mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
#dfall4 <- mat
#test1 <- list()
#test1[[1]] <- dfall4
#test1[[2]] <- dfall5
