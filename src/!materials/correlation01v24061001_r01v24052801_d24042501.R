# Include Libraries -------------------------------------------------------
library(ggplot2)
library(FactoMineR)
library(corrplot)
library(xlsx)


Sys.setlocale(category = "LC_ALL", locale = "Russian, Russia")

# Read in Data ------------------------------------------------------------


source("./src/read01v24052801_d24042501.R")

rm(list = setdiff(ls(), c("clin_data_cleaned_01",
                          "init_aa_names",
                          "t1",
                          "t2",
                          "t3")))

# Functions ---------------------------------------------------------------

# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = 'spearman', ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

chi_square_fun <- function(chi_square_data, v1, v2){
  result <- chisq.test(unlist(chi_square_data[v1]), 
                       unlist(chi_square_data[v2]), 
                       simulate.p.value = T, correct = T)
  return(result$p.value)
  
}

cv.test = function(chi_square_data, v1, v2) {
  x <- unlist(chi_square_data[v1])
  y <- unlist(chi_square_data[v2])
  CV = sqrt(chisq.test(x, y, correct=T)$statistic /
              (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  return(as.numeric(CV))
}

numerical_correlation <- function(numerical_data_filled_na, cleaned_names, 
                                  file_name_pref,
                                  file_name_temp,
                                  file_name_suff,
                                  plot_tit){
  logic_num_data_clean <- sapply(numerical_data_filled_na, function(x) sd(x, na.rm = T) != 0)
  
  num_data_clean <- numerical_data_filled_na[, logic_num_data_clean]
  
  # all ----------------------------------------------------------------
  data_cor_res <- round(cor(num_data_clean, method = 'spearman'), 2)
  rownames(data_cor_res) <-  as.character(factor(rownames(data_cor_res),
                                                 levels = cleaned_names$name2,
                                                 labels = cleaned_names$name1))
  colnames(data_cor_res) <-  factor(colnames(data_cor_res),
                                    levels = cleaned_names$name2,
                                    labels = cleaned_names$name1)
  write.xlsx(data_cor_res, 
             paste0("./writeup/", file_name_pref, "_numeric_spear_", file_name_suff, ".xlsx"))
  
  p.mat_res <- round(cor.mtest(num_data_clean), 4)
  
  rownames(p.mat_res) <-  as.character(factor(rownames(p.mat_res),
                                              levels = cleaned_names$name2,
                                              labels = cleaned_names$name1))
  colnames(p.mat_res) <-  as.character(factor(colnames(p.mat_res),
                                              levels = cleaned_names$name2,
                                              labels = cleaned_names$name1))
  write.xlsx(p.mat_res, 
             paste0("./writeup/", file_name_pref, "_p-value_numeric_spear_", file_name_suff, ".xlsx"))
  
  
  ## Make Figures ##
  png(paste0("./figures/", file_name_pref, "_numeric_", file_name_suff, ".png"), w=1100, h=800)
  corrplot(data_cor_res, method = 'circle', type = 'upper', title = plot_tit,
           order = 'original', p.mat = p.mat_res, sig.level = 0.05,
           tl.col="black", tl.srt=45, tl.cex = 1.3, #Text label color and rotation
           cl.cex = 1, number.cex = 1, pch.cex = 2,
           mar = c(1, 1, 1, 1))
  dev.off()
  
  
  
}

categorical_correlation <- function(categorical_data_filled_na, cleaned_names, 
                                    file_name_pref,
                                    file_name_temp,
                                    file_name_suff,
                                    plot_tit){
  
  logic_cat_data_clean <- sapply(categorical_data_filled_na, function(x) sd(x, na.rm = T) != 0)
  cat_data_clean <- categorical_data_filled_na[, logic_cat_data_clean]
  
  # all ----------------------------------------------------------------
  cramers_v_result <- matrix(ncol = ncol(cat_data_clean),  
                             nrow = ncol(cat_data_clean))
  chi_square_result <- matrix(ncol = ncol(cat_data_clean),  
                              nrow = ncol(cat_data_clean))
  
  for(i in seq(ncol(cat_data_clean))){
    for(j in seq(ncol(cat_data_clean))){
      cramers_v_result[i, j] <- cv.test(cat_data_clean, i, j)
      chi_square_result[i, j] <- chi_square_fun(cat_data_clean, i, j)
      
    }
  }
  
  rownames(cramers_v_result) <-  as.character(factor(colnames(cat_data_clean),
                                                     levels = cleaned_names$name2,
                                                     labels = cleaned_names$name1))
  colnames(cramers_v_result) <-  as.character(factor(colnames(cat_data_clean),
                                                     levels = cleaned_names$name2,
                                                     labels = cleaned_names$name1))
  
  rownames(chi_square_result) <-  as.character(factor(colnames(cat_data_clean),
                                                      levels = cleaned_names$name2,
                                                      labels = cleaned_names$name1))
  colnames(chi_square_result) <-  as.character(factor(colnames(cat_data_clean),
                                                      levels = cleaned_names$name2,
                                                      labels = cleaned_names$name1))
  
  
  write.xlsx(cramers_v_result, 
             paste0("./writeup/", file_name_pref, "_cramers_v_categ_", file_name_suff, ".xlsx"))
  
  write.xlsx(chi_square_result, 
             paste0("./writeup/", file_name_pref, "_chi_square_correlation_", file_name_suff, ".xlsx"))
  
  
  #  logic_filter_chi_sq <- colSums(chi_square_result < 0.05) >= 1
  logic_filter_chi_sq <- colSums(chi_square_result < 0.05) >= 0
  filt_chi_sq <- chi_square_result[logic_filter_chi_sq, logic_filter_chi_sq]
  filt_cramers_v <- cramers_v_result[logic_filter_chi_sq, logic_filter_chi_sq]
  
  #### Make Figure
  png(paste0("./figures/", file_name_pref, "_categorical_", file_name_suff, ".png"), w=4000, h=4000, pointsize=30)
  corrplot(filt_cramers_v, method = 'circle', type = 'upper', title = plot_tit, 
           order = 'original', p.mat = filt_chi_sq, sig.level = 0.05,
           tl.col="black", tl.srt=45, tl.cex = 1.5, #Text label color and rotation
           number.cex = 3, pch.cex = 2,
           cl.lim=c(0, 1), cl.cex = 1,
           mar = c(1, 1, 1, 1))
  dev.off()
  
  
  
}


############ Fill NA ######################
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



filling_na_wrapper <- function(data, groupped = T, fill_by = "mean"){
  columns_to_clean <- colnames(data)
  
  if(groupped == T){
    filled_na <-  data.frame(lapply(columns_to_clean, function(x)
      fill_na_grouped(data, x, fill_by)))
  } else {
    filled_na <-  data.frame(lapply(columns_to_clean, function(x)
      fill_na_ungrouped(data, x, fill_by)))
  }
  colnames(filled_na) <- columns_to_clean
  return(filled_na)
  
}

fill_na_grouped <- function(data, outer_column, fill_by = "mean") {
  column <- sym(outer_column)
  if (fill_by == "mean"){
    grouped_data <- data %>%
      group_by(group) %>%
      dplyr::mutate(!!column :=  ifelse(is.na(!!column), mean(!!column, na.rm = T), !!column))
  } else if (fill_by == "mode"){
    grouped_data <- data %>%
      group_by(group) %>%
      dplyr::mutate(!!column :=  ifelse(is.na(!!column), Mode(!!column), !!column))
  }
  else if (fill_by == "median"){
    grouped_data <- data %>%
      group_by(group) %>%
      dplyr::mutate(!!column :=  ifelse(is.na(!!column), median(!!column, na.rm = T), !!column))
  }
  return(grouped_data[[column]])
}

fill_na_ungrouped <- function(data, outer_column, fill_by = "mean") {
  column <- sym(outer_column)
  if (fill_by == "mean"){
    grouped_data <- data %>%
      dplyr::mutate(!!column :=  ifelse(is.na(!!column), mean(!!column, na.rm = T), !!column))
  } else if (fill_by == "mode"){
    grouped_data <- data %>%
      dplyr::mutate(!!column :=  ifelse(is.na(!!column), Mode(!!column), !!column))
  }
  else if (fill_by == "median"){
    grouped_data <- data %>%
      dplyr::mutate(!!column :=  ifelse(is.na(!!column), median(!!column, na.rm = T), !!column))
  }
  return(grouped_data[[column]])
}


# Clean Data --------------------------------------------------------------
init_names <- init_aa_names

## v258, v262, v263 - 0/1
t1_cor <- t1 %>%
  select(v003, v004, v005, v014, v017, v019, v040, v043, 
         v044, v053, 
         v056, v057, v063, 
         v079, v157, v160,  
         v163, v166, v169,  
         v172, v241, v244, v247, v250, v253)

# Correlation, by groups -------------------------------------
t1g1_cor <- t1_cor %>%
  filter(v003 == 1) %>%
  select(-v003)

t1g1_cor_filled_na <- filling_na_wrapper(t1g1_cor, groupped = F, fill_by = "mode")

numerical_correlation(t1g1_cor_filled_na, init_names, 
                      file_name_pref = "correlation01v24061001_t1g1",
                      file_name_temp = "t1g1",
                      file_name_suff = "r01v24052801_d24042501",
                      plot_tit = "")

t1g2_cor <- t1_cor %>%
  filter(v003 == 2) %>%
  select(-v003)

t1g2_cor_filled_na <- filling_na_wrapper(t1g2_cor, groupped = F, fill_by = "mode")

numerical_correlation(t1g2_cor_filled_na, init_names, 
                      file_name_pref = "correlation01v24061001_t1g2",
                      file_name_temp = "t1g2",
                      file_name_suff = "r01v24052801_d24042501",
                      plot_tit = "")


# Correlation, all groups -------------------------------------
t1_cor_filled_na <- rbind(t1g1_cor_filled_na,
                          t1g2_cor_filled_na)

numerical_correlation(t1_cor_filled_na, init_names, 
                      file_name_pref = "correlation01v24061001_t1",
                      file_name_temp = "t1",
                      file_name_suff = "r01v24052801_d24042501",
                      plot_tit = "")


# testing area ------------------------------------------------------------







