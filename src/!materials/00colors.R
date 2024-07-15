## Transforms a character or numeric vector into colors
.colorF <- function(namVcn) {
  ## 16 color palette without 'gray'
  palVc <- c("red",  "blue", "yellow","green", "magenta", "#FF7F00", "#6A3D9A", "#B15928", "aquamarine4", "yellow4", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#FFFF99")
  
  if(is.null(namVcn) || all(is.na(namVcn))) {
    
    if(!is.null(namVcn)) {
      
      dev.new()
      
      palNamVc <- paste0(1:length(palVc),
                         "_",
                         palVc)
      
      pie(rep(1, length(palVc)),
          col = palVc,
          labels = palNamVc)
      
      print(matrix(palNamVc, ncol = 1))
      
    }
    
    return(palVc)
    
  } else {
    
    if(is.character(namVcn)) {
      
      namFcn <- factor(namVcn)
      
      if(length(levels(namFcn)) <= length(palVc)) {
        scaVc <- palVc[1:length(levels(namFcn))]
      } else
        scaVc <- c(palVc,
                   rep("gray",
                       length(levels(namFcn)) - length(palVc)))
      
      names(scaVc) <- levels(namFcn)
      
      colVc <- scaVc[unlist(sapply(namVcn,
                                   function(scaleC) {
                                     if(is.na(scaleC))
                                       return(NA)
                                     else
                                       which(levels(namFcn) == scaleC)
                                   }))]
      
    } else if(is.numeric(namVcn)) {
      
      scaVc <- rev(rainbow(100, end = 4/6))
      if(length(namVcn) > 1) {
#### initial string        colVc <- scaVc[round((namVcn - min(namVcn, na.rm = TRUE)) / diff(range(namVcn, na.rm = TRUE)) * 99) + 1]
        colVc <- palVc[namVcn + 1]
      } else
        colVc <- rep("black", length(namVcn))
      
    } else
      stop("'namVcn' argument must be a vector of either character or numeric mode", call. = FALSE)
    
    colVc[is.na(colVc)] <- "grey"
    names(colVc) <- namVcn
    
  }
  
  return(list(colVc = colVc,
              scaVc = scaVc))
  
}  ## end of .colorF()


