# Scores Plot Drawing -----------------------------------------------------


opls.plot <- function(opls.model, put_labs = F, draw_ellipse = F, fig_name = "tst.tiff"){

  ############### Plot preparation
  opl <- opls.model
  typeVc = "x-score"
  parCompVi = c(1, 2)
  if (grepl("^OPLS", opl@typeC)) {
    tCompMN <- cbind(opl@scoreMN[, 1], opl@orthoScoreMN[, 1])
  } else{
    tCompMN <- cbind(opl@scoreMN[, 1], opl@scoreMN[, 2])
  }
  
  hotFisN <- (nrow(tCompMN) - 1) * 2 * (nrow(tCompMN)^2 - 1) / (nrow(tCompMN) * nrow(tCompMN) * (nrow(tCompMN) - 2)) * qf(0.95, 2, nrow(tCompMN) - 2)
  radVn <- seq(0, 2 * pi, length.out = 100)
  obsColVc <- .colorF(c(opl@suppLs$data_groups))[["colVc"]]
  
  #obsColVc <- .colorF(c(opl@suppLs[["yMCN"]]))[["colVc"]]
  #obsColVc <- .colorF(endoTable[,1])[["colVc"]]

  #obsLegVc <- as.vector(parAsColFcVn)
  ploColVc <- obsColVc
  obsLabVc <- rownames(opl@suppLs[["yMCN"]])
  parCexN = 0.8
  
  layL <- F
  cexRqcN <- ifelse(layL, 0.7, 1)
  parTitleL = TRUE
  
#### VC added labeles shift ####
  if(opl@suppLs$scaleC == "standard"){
    labeles_shift <- 0.5
  } else{
    labeles_shift <- 0.05
  }
  
 
#####
  
  
  
  ###############
  
  maiC <- paste0("Scores (", opl@typeC, ")")
  
  xLabC <- paste0("t",
                  parCompVi[1],
                  " (",
                  round(opl@modelDF[parCompVi[1], "R2X"] * 100),
                  "%)")
  
  yLabC <- paste0("t",
                  parCompVi[2],
                  " (",
                  round(opl@modelDF[parCompVi[2], "R2X"] * 100),
                  "%)")
  
  ploMN <- tCompMN
  
  if(grepl("^OPLS", opl@typeC))
    yLabC <- paste0("to", parCompVi[2] - 1)
  
  xLimVn <- c(-1, 1) * max(sqrt(var(ploMN[, 1]) * hotFisN), max(abs(ploMN[, 1])))
  yLimVn <- c(-1, 1) *max(sqrt(var(ploMN[, 2]) * hotFisN), max(abs(ploMN[, 2])))
  
  ploColVc <- obsColVc
  
  
  
  #### Alisa's part ####

if(draw_ellipse){  
    param<-unique(opl@suppLs[["yMCN"]])
    index<-which(param[1]==opl@suppLs$yMCN)
    circle1<-ploMN[index,]
    index<-which(param[2]==opl@suppLs$yMCN)
    circle2<-ploMN[index,]
    cov1<-cov(circle1)
    cov2<-cov(circle2)
    eig1<-eigen(cov1)
    eig2<-eigen(cov2)
    angle<-rbind(cos(radVn),sin(radVn))
    coord_opls1<-matrix(data=c(mean(circle1[,1]),mean(circle1[,2])),
                        nrow=2,ncol=100)+(eig1$vectors%*%sqrt(diag(eig1$values)*qchisq(0.95,2)))%*%angle
    coord_opls2<-matrix(data=c(mean(circle2[,1]),mean(circle2[,2])),
                        nrow=2,ncol=100)+(eig2$vectors%*%sqrt(diag(eig2$values)*qchisq(0.95,2)))%*%angle
  }
    
  
   ############### Plot drawing
  layRowN <- ceiling(sqrt(length(typeVc)))

  tiff(paste("./figures/", fig_name),
      w = 1100,
      h = 1000,
      pointsize = 20) 
  
  
  layout(matrix(1:layRowN^2, byrow = TRUE, nrow = layRowN))
  marVn <- c(5.1, 4.1, 4.1, 2.1)
  par(font=2, font.axis=2, font.lab=2, lwd=2,
      mar=marVn,
      pch=18)
  ##################################
  
  if(is.null(xLimVn))
    xLimVn <- range(ploMN[, 1])
  if(is.null(yLimVn))
    yLimVn <- range(ploMN[, 2])
 

   
  plot(ploMN,
       main=ifelse(parTitleL, maiC, ""),
       type = "p",
       pch = 19,
       cex = 2,
       xlab = xLabC,
       ylab = yLabC,
       xlim = xLimVn,
       ylim = yLimVn,
       col = ploColVc)
  
  
  abline(v = axTicks(1),
         col = "grey")
  
  abline(h = axTicks(2),
         col = "grey")
  
  abline(v = 0)
  abline(h = 0)
  
  ###############
  
  lines(sqrt(var(ploMN[, 1]) * hotFisN) * cos(radVn),
        sqrt(var(ploMN[, 2]) * hotFisN) * sin(radVn))
  
  if(draw_ellipse){ 
      #### Alisa's part ####
    lines(coord_opls2[1,],coord_opls2[2,],col="black")
    lines(coord_opls1[1,],coord_opls1[2,],col="black")
  }
  
  ## Tenenhaus98, p87
  
  # if(!is.null(obsLegVc))
  #   .legendF(obsLegVc,
  #            ploMN)
  
  
if(put_labs){  
  text(ploMN + labeles_shift,
       cex = parCexN,
       col = ploColVc,
       #     labels = obsLabVc)
       labels = opl@suppLs$data_labels)
}

  dev.off()  
return()  
}



#### Colors Definitions ####

## Transforms a character or numeric vector into colors
.colorF <- function(namVcn) {
  ## 16 color palette without 'gray'
  palVc <- c("red", "green",  "blue", "yellow","magenta", "#FF7F00", "#6A3D9A", "#B15928", "aquamarine4", "yellow4", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#FFFF99")
  
  if(is.null(namVcn) || all(is.na(namVcn))) {
    
    if(!is.null(namVcn)) {
      
      dev.new()
      
      palNamVc <- paste0(1:length(palVc),
                         "_",
                         palVc)
      
      pie(rep(1, length(palVc)),
          col = palVc,
          labels = palNamVc)
      
#      print(matrix(palNamVc, ncol = 1))
      
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
#####



