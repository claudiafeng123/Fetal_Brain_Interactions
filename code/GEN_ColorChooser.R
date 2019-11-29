
ColorChooser <- function(chrs = cellTypeList,
                         cellTypeList = c("End", "ExM", "ExDp1", "ExDp2", "ExM", "ExM-U", "ExN",
                                          "InCGE", "InMGE", "IP", 
                                          "Mic", 
                                          "OPC", "oRG", 
                                          "Per", "PgG2M", "PgS", 
                                          "vRG"),
                         colors = hue_pal()(length(cellTypeList))){
  rtn <- colors[match(chrs, cellTypeList)]
  return(rtn)
}
