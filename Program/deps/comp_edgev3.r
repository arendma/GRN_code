comp_edge <- function (e_lista, e_listb) {
  #takes to edgelists chr(node1) -> chr(node2) and returns the logical vector indicating if edge of lista is contained in listb
  ce_lista <- paste(e_lista[,1], e_lista[, 2], sep=' ')
  ce_listb <- paste(e_listb[,1], e_listb[, 2], sep=' ')
  
  return(ce_lista %in% ce_listb)
}