# Import Libraries --------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(dplyr)
library(plyr)
library(stringr)
library(RColorBrewer)
library(readtext)
library(tidyr)
library(pROC)
library(purrr)

Sys.setlocale(category = "LC_ALL", locale = "Russian, Russia")

source("./src/read01v24052801_d24042501.R")

rm(list = setdiff(ls(), c("clin_data_cleaned_01",
                          "init_aa_names",
                          "t1",
                          "t2",
                          "t3")))

# Functions --------------------------------------------------------------- 

logit <- function(data_logit){
  model <- glm(group ~ ., family = binomial(link = 'logit'), data = data_logit)
  return(model)
  
}


logit_modeling <- function(data,
                           group_a, 
                           group_b,
                           file_name_pref,
                           file_name_temp,
                           file_name_suff,
                           plot_tit){
  
  # Log Regression ----------------------------------------------------------
  
  group1 <- group_a
  group2 <- group_b
  
  transferred_all <- data %>%
    filter(group == group1 | group == group2) %>%
    mutate(group = ifelse(group == group1, 0,1))

  
  file_name_prefix <- file_name_pref
  file_name_template <- file_name_temp
  file_name_suffix <- file_name_suff
  plot_title <- plot_tit
  
  auc_all_comb <- data.frame(stringsAsFactors = F)
  roc_plot = list()
  model_parameters <- data.frame()
  model_coef <- list()
  model_parameters_summary <- matrix(nrow = 0, ncol = 10)
  n <- 0
  
  for (i in seq(1, 4)) {
    combinations = combn(names(transferred_all[,-1]), m = i)
    for (j in seq(ncol(combinations))) {
      n <- n + 1
      data_selected <- transferred_all %>%
        select(group, combinations[, j])
      
      model_result <- logit(data_selected)
      model_prediction <- predict(model_result, type = "response")
      area <- round(pROC::auc(data_selected$group, model_prediction), 4)
      roc.plot <- pROC::roc(data_selected$group, model_prediction)
      
      model_parameters <- rbind(model_parameters, data.frame(model_id = n,
                                                             vars_list = paste(combinations[, j], collapse = ", "),
                                                             AUC = area,
                                                             stringsAsFactors = F))
      
      roc_plot[n] <- list(roc.plot)
    }
  }
  
  # Select 4 models with the biggest AUC ------------------------------------
  
  top_4_data <- model_parameters[order(-model_parameters$AUC),][seq(1, 4),]
  top_4_models <- roc_plot[top_4_data$model_id]
  
  # Additional statistics for a model ---------------------------------------
  
  model_params <- data.frame()
  model_coefficients_summary <- matrix(nrow = 0, ncol = 11)
  n <- 0
  model_tmp <- data.frame()
  
  for(i in top_4_data$vars_list) {
    n = n + 1  
    selected_vars <- unlist(strsplit(i, ', '))
    data_selected <- transferred_all %>%
      select(group, all_of(selected_vars))
    
    model_result <- logit(data_selected)
    model_prediction <- predict(model_result, type = "response")
    area <- round(pROC::auc(data_selected$group, model_prediction), 4)
    roc.plot <- pROC::roc(data_selected$group, model_prediction)
    
    model_threshold <-  round(pROC::coords(roc.plot, "best", ret = c("threshold"),
                                           transpose = T
    ), 4)
    
    
    parameters_list_formatted <- paste(as.character(factor(str_trim(str_split(i, ",", simplify = T)),
                                                           levels = init_names$name2,
                                                           labels = init_names$name1)), collapse = ", ")
    
    
    
    add_params <- ci.coords(roc.plot, x="best", best.policy = "random", input = "threshold", ret=c("sensitivity", "specificity", "ppv", "tp"))
    add_params1 <- coords(roc.plot, x="best", ret = c("fn", "fp", "tn", "tp", "tpr", "fpr", "tnr", "fnr"))
    transformed_add_params <- add_params$sensitivity[1]
    model_params <-
      rbind(
        model_params,
        data.frame(
          parameters_list = parameters_list_formatted,
          AUC = round(area, 2),
          threshold = round(model_threshold, 2),
          sensitivity = paste0(round(add_params$sensitivity[2], 2), ' (',
                               round(add_params$sensitivity[1], 2), '; ',
                               round(add_params$sensitivity[3], 2), ')'),
          sensitivity = paste0(round(add_params$specificity[2], 2), ' (',
                               round(add_params$specificity[1], 2), '; ',
                               round(add_params$specificity[3], 2), ')'),
          ppv = paste0(round(add_params$ppv[2], 2), ' (',
                       round(add_params$ppv[1], 2), '; ',
                       round(add_params$ppv[3], 2), ')'),
          tp = paste0(round(add_params$tp[2], 2), ' (',
                      round(add_params$tp[1], 2), '; ',
                      round(add_params$tp[3], 2), ')'),
          
          FN = add_params1[1, 1],
          FP = add_params1[1, 2],
          TN = add_params1[1, 3],
          TP = add_params1[1, 4],
          TPR = round(add_params1[1, 5], 2),
          FPR = round(add_params1[1, 6], 2),
          TNR = round(add_params1[1, 7], 2),
          FNR = round(add_params1[1, 8], 2),
          lr_plus = round(add_params1[1, 5] / add_params1[1, 6], 2),
          lr_minus = round(add_params1[1, 8] / add_params1[1, 7], 2),
          stringsAsFactors = F
        )
      )
    
    glm_coefficients <- round(coef(summary(model_result)), 4)
    glm_coefficients <- cbind(row.names(glm_coefficients),
                              glm_coefficients, confint(model_result), 
                              or = exp(coef(model_result)), 
                              exp(confint(model_result)), n)
    colnames(glm_coefficients) <- c("Coefficients",
                                    "Estimate",
                                    "Error",
                                    "Wald",
                                    "p_value",
                                    "ci_2p5",
                                    "ci_97p5",
                                    "or",
                                    "or_2p5",
                                    "or_97p5",
                                    "#")
    model_coefficients_summary <- rbind(model_coefficients_summary, glm_coefficients)
    row.names(model_coefficients_summary) <- str_replace_all(row.names(model_coefficients_summary), '\\.', '-')
    row.names(model_coefficients_summary) <- str_replace_all(row.names(model_coefficients_summary), "^X", "")
    
    
    
    
    add_params1 <- coords(roc.plot, x="best", ret = c("fn", "fp", "tn", "tp", "tpr", "fpr"), transpose = FALSE)
    #    coords(roc.plot, x="best", best.policy = "random", input = "threshold", ret=c("sensitivity", "specificity", "ppv", "tp"))
    
    a <- data.frame(group = data_selected$group, prognosis = model_prediction) %>%
      mutate(prediction_0_1 = ifelse(prognosis>= model_threshold[1], 1, 0)) %>%
      mutate(group = factor(group, levels = c(group1, group2),
                            labels = c(0, 1))) %>%
      mutate(group = as.numeric(as.character(group)))
    
    
    positives_amount <- sum(a$group)
    negatives_amount <- nrow(a) - sum(a$group)
    false_negative <- sum(ifelse((a == 1) & (a$prediction_0_1 == 0), 1, 0))
    false_positive <- sum(ifelse((a$group == 0) & (a$prediction_0_1 == 1), 1, 0))
    true_negative <- sum(ifelse((a$group == 0) & (a$prediction_0_1 == 0), 1, 0))
    true_positive <- sum(ifelse((a$group == 1) & (a$prediction_0_1 == 1), 1, 0))
    TPR <- true_positive/positives_amount
    FPR <- false_positive/negatives_amount
    TNR <-true_negative/negatives_amount
    FNR <- false_negative/positives_amount
    lr_plus <- TPR / FPR
    lr_minus <- FNR / TNR
    
    model_tmp <- rbind(
      model_tmp,
      data.frame(
        TPR = TPR,
        FPR = FPR,
        TNR = TNR,
        FNR = FNR,
        lr_plus = lr_plus,
        lr_minus = lr_minus,
        stringsAsFactors = F
      ))
    
    
    
  }
  
  model_coefficients_summary <-  as.data.frame(model_coefficients_summary, stringsAsFactors = F) %>%
    mutate(Coefficients = factor(Coefficients,
                                 levels = c("(Intercept)", init_names$name2),
                                 labels = c("(Intercept)", init_names$name1)))
  
  write.csv(model_params, 
            paste0('./writeup/', file_name_prefix, file_name_template, '_model_params_', file_name_suffix, '.csv'), 
            row.names = F)
  
  write.csv(model_coefficients_summary, 
            paste0('./writeup/', file_name_prefix, file_name_template, '_model_coeff_', file_name_suffix, '.csv'), 
            row.names = F)
  
  # Make ROC figure ---------------------------------------------------------
  
  for(i in seq(length(top_4_data$vars_list))){
    top_4_data$vars_list[i] <- paste(as.character(factor(str_trim(str_split(top_4_data$vars_list[i], ",", simplify = T)),
                                                         levels = init_names$name2,
                                                         labels = init_names$name1)), collapse = ", ")
    
  }
  
  png(paste0('./figures/', file_name_prefix, file_name_template, '_rocr_', file_name_suffix, '.png'),
      w = 1100,
      h = 1000,
      pointsize = 20)
  
  l = 3
  plot(top_4_models[[1]] #smooth(top_4_models[[1]], method="binormal")
       , col = 1, lty = 1, lwd = l, smooth = F, xlab = '1 - Specificity', ylab = 'Sensitivity')
  legend_text <- c(paste0(top_4_data[1, 2], '; AUC = ', round(top_4_data$AUC[1], 3)))
  
  
  for(i in seq(2, 4)){
    plot(top_4_models[[i]]  #smooth(top_4_models[[i]], method="binormal")
         , col = i, lty = i, lwd = l, add = T, smooth = F)
    legend_text <- c(legend_text,
                     paste0(top_4_data[i, 2], '; AUC = ', round(top_4_data$AUC[i], 3)))
  }
  
  legend("bottomright", legend = legend_text, col = 1:9, lty = 1:9, lwd = 2, merge = F)
  
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

# End Functions -----------------------------------------------------------

# Modeling ----------------------------------------------------------------
init_names <- init_aa_names

# t1 ------------------------------------------------

## v258, v262, v263 - 0/1
t1_stat <- t1 %>%
  select(v003, v040, v043, 
         v044, v053, 
         v056, v057, v063, 
         v079)


t1 <- t1_stat %>%
  dplyr::rename(group=v003) %>%
  mutate(group=factor(group, levels = c(1, 2),
                      labels = c("1",
                                 "2")))
## Fill NA
t1 <- cbind(t1[1],
                 filling_na_wrapper(t1[-1], groupped = F, fill_by = "mean"))



logit_modeling(
  t1,
  group_a = "1",
  group_b = "2",
  file_name_pref = "log_regr01v24061101_",
  file_name_temp = "g1_vs_g2",
  file_name_suff = "r01v24052801_d24042501",
  plot_tit = ""
)



# testing area ------------------------------------------------------------


