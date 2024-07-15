
# Import Libraries --------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(dplyr)
library(plyr)
library(stringr)
library(RColorBrewer)
library(readtext)
library(tidyr)
library(ggsignif)
library(reshape2)

library(rstatix)

Sys.setlocale(category = "LC_ALL", locale = "Russian, Russia")

source("./src/read01v24052801_d24042501.R")

rm(list = setdiff(ls(), c("clin_data_cleaned_01",
                          "init_aa_names",
                          "t1",
                          "t2",
                          "t3")))

# Functions ---------------------------------------------------------------

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))


wilcox_test_multiple <- function(data){
  mir_list <- unique(data$metabolite)
  p_value <- data.frame()
  for(i in mir_list){
    mir_data <- data %>%
      filter(metabolite == i)
    mir_pvalue <- pairwise.wilcox.test(mir_data$level, 
                                       mir_data$group, 
                                       p.adjust.method = "none")$p.value
   
    pvalue_result <- data.frame(metabolite = i, 
                    group = row.names(mir_pvalue), 
                    round(mir_pvalue, 4))
    p_value <- rbind(p_value, pvalue_result)
    
  }
  return(p_value)
}

t_test_multiple <- function(data){
  mir_list <- unique(data$metabolite)
  p_value <- data.frame()
  for(i in mir_list){
    mir_data <- data %>%
      filter(metabolite == i)
    mir_pvalue <- pairwise.t.test(mir_data$level, 
                                       mir_data$group, 
                                       p.adjust.method = "none")$p.value
    
    pvalue_result <- data.frame(metabolite = i, 
                                group = row.names(mir_pvalue), 
                                round(mir_pvalue, 4))
    p_value <- rbind(p_value, pvalue_result)
    
  }
  return(p_value)
}

chisq_test_multiple <- function(data){
  vars_list <- unique(data$vars)
  p_value <- data.frame()
  for(i in vars_list){
    vars_data <- data %>%
      filter(vars == i)
    vars_pvalue <- chisq.test(vars_data$values, 
                              vars_data$group)$p.value
    
    pvalue_result <- data.frame(vars = i, 
                                round(vars_pvalue, 4))
    p_value <- rbind(p_value, pvalue_result)
    
  }
  return(p_value)
}

data_stats <- function(data,
                       data_wide,
                       signif_y_manual_calc = T,
                       signif_y_correction = 0,
                       custom_ylim,
                       file_name_pref,
                       file_name_temp,
                       file_name_suff,
                       x_label,
                       y_label,
                       fill_label,
                       plot_tit){
  
  custom_ylim <- as.numeric(custom_ylim)
  plot_data <- data
  
  file_name_prefix <- file_name_pref
  file_name_template <- file_name_temp
  file_name_suffix <- file_name_suff
  plot_title <- plot_tit
  
  # Student Stats -----------------------------------------------------------
  
  t_stats <- t_test_multiple(plot_data)
  
  mean_sd_result <- plot_data %>%
    group_by(metabolite, group) %>%
    dplyr::summarise(conc_mean = mean(level),
                     conc_sd = sd(level)) %>%
    ungroup() %>%
    mutate(mean_sd = paste(round(conc_mean, 2), round(conc_sd, 2), sep = "\u00B1")) %>%
    dcast(formula = metabolite ~ group, value.var = "mean_sd") 
  
  write.csv(mean_sd_result, paste0('./writeup/', file_name_prefix, file_name_template, '_mean_sd_', file_name_suffix, '.csv'))
  write.csv(t_stats, paste0('./writeup/', file_name_prefix, file_name_template, '_t_test_p_value_', file_name_suffix, '.csv'))
  
  # Wilcox Stats -----------------------------------------------------------
  
  w_stats <- wilcox_test_multiple(plot_data)
  
  
  median_q_result <- plot_data %>%
    group_by(metabolite, group) %>%
    dplyr::summarise(conc_median = median(level),
                     q1 = quantile(level, 0.25),
                     q3 = quantile(level, 0.75)) %>%
    ungroup() %>%
    mutate(median_q = paste0(round(conc_median, 2), " (", round(q1, 2), "; ", round(q3, 2), ")")) %>%
    dcast(formula = metabolite ~ group, value.var = "median_q") 
  
  write.csv(median_q_result, paste0('./writeup/', file_name_prefix, file_name_template, '_median_q_', file_name_suffix, '.csv'))
  write.csv(w_stats, paste0('./writeup/', file_name_prefix, file_name_template, '_w_test_p_value_', file_name_suffix, '.csv'))
  
  
  # Make Plot ---------------------------------------------------------------
  plot_data_for_signif <- data_wide %>%
    gather(metabolite, level, -group)
  
  stat.tst <- plot_data_for_signif %>%
    group_by(metabolite) %>%
    wilcox_test(level ~ group) %>%
    adjust_pvalue(method = "none") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "metabolite", 
                    fun = "max", 
                    #                    dodge = 0.8,
                    step.increase = 0.05)
  
  if(signif_y_manual_calc){
    stat.tst$y.position <- calc_y_position(plot_data_for_signif)
    
  } else {
    stat.tst$y.position <- calc_y_position_corrected(plot_data_for_signif)
  }
  if(signif_y_correction != 0){
    stat.tst$y.position <- stat.tst$y.position * (1 - signif_y_correction) 
    
  }
  
  
  stat.tst1 <- stat.tst %>%
    filter(p <= 0.05)
  #  stat.tst1 <- stat.tst
  
  
  plot.theme = theme(
    text = element_text(size = 20),
    title = element_text(hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5),
    #  legend.position = c(0.9, 0.95),
    axis.text.x = element_text(
      angle = 75,
      hjust = 1,
      vjust = 1,
      size = 20
    )
  )
  
  
  common_labs = labs(
    x = x_label,
    y = y_label,
    fill = fill_label,
    title = plot_title
  )
  
  # p.adj.signif  
  metabolite.plot <-
    ggplot(plot_data, aes(x = metabolite, y = level)) +
    geom_boxplot(aes(fill = group), outlier.shape=NA) +
    theme_classic() +
    ylim(custom_ylim[1], custom_ylim[2]) +
    plot.theme + common_labs +
    stat_pvalue_manual(stat.tst1,  label = "p.adj.signif", tip.length = 0.005, label.size = 6, bracket.size = 0.5) +
    scale_fill_manual(values = c('#999999', '#ef8a62', '#8dd3c7', '#ffffb3', 
                                 '#fb8072', '#80b1d3', '#386cb0', '#7fc97f')) 
  
  ggsave(
    paste0("figures/", file_name_prefix, file_name_template, "_box_", file_name_suffix, '.jpg'),
    plot = metabolite.plot,
    device = "jpeg",
    scale = 1,
    width = 50,
    height = 30,
    units = "cm",
    dpi = 800,
    limitsize = TRUE
  ) 
  return()
  
}

calc_y_position <- function(data){
  mean_sd_result <- data %>%
    group_by(metabolite, group) %>%
    dplyr::summarise(q3 = quantile(level, 0.75),
                     iqr = IQR(level)) %>%
    mutate(max_value = q3 + 1.5 * iqr) %>%
    ungroup()
  
  outlier_filtered <- data.frame(group = character(),
                                 metabolite = character(),
                                 level = numeric())
  for(i in unique(data$group)) {
    for (j in unique(data$metabolite)) {
      outlier_threshold <- mean_sd_result %>%
        dplyr::filter(group == i) %>%
        dplyr::filter(metabolite == j)
      
      outlier_filtered_tmp <- data %>%
        filter(group == i) %>%
        filter(metabolite == j) %>%
        filter(level <= outlier_threshold$max_value)
      outlier_filtered <- rbind(outlier_filtered,
                                outlier_filtered_tmp)
    }
  }
  
  data_for_y_position <- outlier_filtered %>%
    group_by(metabolite) %>%
    wilcox_test(level ~ group) %>%
    add_xy_position(x = "metabolite",
                    fun = "max",
                    dodge = 0.8,
                    step.increase = 0.07)
  
  return(data_for_y_position$y.position)
}

calc_y_position_corrected <- function(data){
  mean_sd_result <- data %>%
    group_by(metabolite, group) %>%
    dplyr::summarise(q3 = quantile(level, 0.75),
                     iqr = IQR(level)) %>%
    mutate(max_value = q3 + 1.5 * iqr) %>%
    ungroup()
  
  outlier_filtered <- data.frame(group = character(),
                                 metabolite = character(),
                                 level = numeric())
  for(i in unique(data$group)) {
    for (j in unique(data$metabolite)) {
      outlier_threshold <- mean_sd_result %>%
        dplyr::filter(group == i) %>%
        dplyr::filter(metabolite == j)
      
      outlier_filtered_tmp <- data %>%
        filter(group == i) %>%
        filter(metabolite == j) %>%
        filter(level <= outlier_threshold$max_value)
      outlier_filtered <- rbind(outlier_filtered,
                                outlier_filtered_tmp)
    }
  }
  
  data_for_y_position <- outlier_filtered %>%
    group_by(metabolite) %>%
    dunn_test(level ~ group) %>%
    add_xy_position(x = "metabolite",
                    fun = "max",
                    dodge = 0.8,
                    step.increase = 0.07)
  
  return(data_for_y_position$y.position)
}


process_categorical_data <- function(data, 
                                     file_name_pref,
                                     file_name_temp,
                                     file_name_suff){
  
  
  plot_data <- data
  
  file_name_prefix <- file_name_pref
  file_name_template <- file_name_temp
  file_name_suffix <- file_name_suff
  
  # Chi-Square Stats -----------------------------------------------------------
  
  chi_stats <- chisq_test_multiple(plot_data)
  
  
  # median_q_result <- plot_data %>%
  #   group_by(vars, group) %>%
  #   dplyr::summarise(values_median = median(values, na.rm = T),
  #                    q1 = quantile(values, 0.25, na.rm = T),
  #                    q3 = quantile(values, 0.75, na.rm = T)) %>%
  #   ungroup() %>%
  #   mutate(median_q = paste0(round(values_median, 2), " (", round(q1, 2), "; ", round(q3, 2), ")")) %>%
  #   dcast(formula = vars ~ group, value.var = "median_q") 
  
  #  write.csv(median_q_result, paste0('./writeup/', file_name_prefix, file_name_template, '_median_q_', file_name_suffix, '.csv'))
  
  colnames(chi_stats) <- c("Variable", "p-value")
  write.csv(chi_stats, paste0('./writeup/', file_name_prefix, file_name_template, '_chisq_test_p_value_', file_name_suffix, '.csv'))
  
  
  return()  
  
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




# Read in Data ------------------------------------------------------------
init_names <- init_aa_names

# t1 ------------------------------------------------

## v258, v262, v263 - 0/1
t1_stat <- t1 %>%
  select(v160,  
         v169,  
         v250,
         v003, v040, v043, 
         v044, v053, 
         v056, v057, v063, 
         v079) %>%
  mutate(v160 = ifelse(v160 == 0, 0, 1)) %>%
  mutate(v169 = ifelse(v169 == 0, 0, 1)) %>%
  mutate(v250 = ifelse(v250 == 0, 0, 1))


# ПНЯ Вазомоторные симптомы -----------------------------------------------

#### g1g2
t1_g1g2_v160 <- t1_stat %>%
  select(-v169, -v250, -v003) %>%
  dplyr::rename(group=v160) %>%
  mutate(group=factor(group, levels = c(0, 1),
                      labels = c("Нет симпт.",
                                 "Есть симпт.")))


## Standard data
t1_g1g2_v160_stdz <- cbind(group = t1_g1g2_v160[, 1],
                               as.data.frame(scale(t1_g1g2_v160[, -1])))

## Fill NA
t1_g1g2_v160_stdz <- cbind(t1_g1g2_v160_stdz[1],
                           filling_na_wrapper(t1_g1g2_v160_stdz[-1], groupped = F, fill_by = "mean"))

t1_g1g2_v160_stdz_wide <- t1_g1g2_v160_stdz

t1_g1g2_v160_stdz_long <- t1_g1g2_v160_stdz_wide %>%
  gather(metabolite, level, -group) %>%
  mutate(metabolite = factor(metabolite, levels = init_names$name2,
                             labels = init_names$name1))

data_stats(t1_g1g2_v160_stdz_long, 
           t1_g1g2_v160_stdz_wide,
           signif_y_manual_calc = F,
           custom_ylim = c(-2, 3),
           file_name_pref = "stats01v24052801_",
           file_name_temp = "t1_g1g2_v160_stdz_box",
           file_name_suff = "r01v24052801_d24042501",
           x_label = "Параметр",
           y_label = "Уровень",
           fill_label = "Группа",
           plot_tit = "")

######################## g1
t1_g1_v160 <- t1_stat %>%
  filter(v003 == 1) %>%
  select(-v169, -v250, -v003) %>%
  dplyr::rename(group=v160) %>%
  mutate(group=factor(group, levels = c(0, 1),
                      labels = c("Нет симпт.",
                                 "Есть симпт.")))


## Standard data
t1_g1_v160_stdz <- cbind(group = t1_g1_v160[, 1],
                           as.data.frame(scale(t1_g1_v160[, -1])))

## Fill NA
t1_g1_v160_stdz <- cbind(t1_g1_v160_stdz[1],
                           filling_na_wrapper(t1_g1_v160_stdz[-1], groupped = F, fill_by = "mean"))

t1_g1_v160_stdz_wide <- t1_g1_v160_stdz

t1_g1_v160_stdz_long <- t1_g1_v160_stdz_wide %>%
  gather(metabolite, level, -group) %>%
  mutate(metabolite = factor(metabolite, levels = init_names$name2,
                             labels = init_names$name1))

data_stats(t1_g1_v160_stdz_long, 
           t1_g1_v160_stdz_wide,
           signif_y_manual_calc = F,
           custom_ylim = c(-2, 3),
           file_name_pref = "stats01v24052801_",
           file_name_temp = "t1_g1_v160_stdz_box",
           file_name_suff = "r01v24052801_d24042501",
           x_label = "Параметр",
           y_label = "Уровень",
           fill_label = "Группа",
           plot_tit = "")


######################## g2
t1_g2_v160 <- t1_stat %>%
  filter(v003 == 2) %>%
  select(-v169, -v250, -v003) %>%
  dplyr::rename(group=v160) %>%
  mutate(group=factor(group, levels = c(0, 1),
                      labels = c("Нет симпт.",
                                 "Есть симпт.")))


## Standard data
t1_g2_v160_stdz <- cbind(group = t1_g2_v160[, 1],
                         as.data.frame(scale(t1_g2_v160[, -1])))

t1_g2_v160_stdz_wide <- t1_g2_v160_stdz

t1_g2_v160_stdz_long <- t1_g2_v160_stdz_wide %>%
  gather(metabolite, level, -group) %>%
  mutate(metabolite = factor(metabolite, levels = init_names$name2,
                             labels = init_names$name1))

data_stats(t1_g2_v160_stdz_long, 
           t1_g2_v160_stdz_wide,
           signif_y_manual_calc = F,
           custom_ylim = c(-2, 3),
           file_name_pref = "stats01v24052801_",
           file_name_temp = "t1_g2_v160_stdz_box",
           file_name_suff = "r01v24052801_d24042501",
           x_label = "Параметр",
           y_label = "Уровень",
           fill_label = "Группа",
           plot_tit = "")


# ПНЯ Симптомы ВВА -----------------------------------------------

#### g1g2
t1_g1g2_v169 <- t1_stat %>%
  select(-v160, -v250, -v003) %>%
  dplyr::rename(group=v169) %>%
  mutate(group=factor(group, levels = c(0, 1),
                      labels = c("Нет симпт.",
                                 "Есть симпт.")))


## Standard data
t1_g1g2_v169_stdz <- cbind(group = t1_g1g2_v169[, 1],
                           as.data.frame(scale(t1_g1g2_v169[, -1])))

## Fill NA
t1_g1g2_v169_stdz <- cbind(t1_g1g2_v169_stdz[1],
                           filling_na_wrapper(t1_g1g2_v169_stdz[-1], groupped = F, fill_by = "mean"))

t1_g1g2_v169_stdz_wide <- t1_g1g2_v169_stdz

t1_g1g2_v169_stdz_long <- t1_g1g2_v169_stdz_wide %>%
  gather(metabolite, level, -group) %>%
  mutate(metabolite = factor(metabolite, levels = init_names$name2,
                             labels = init_names$name1))

data_stats(t1_g1g2_v169_stdz_long, 
           t1_g1g2_v169_stdz_wide,
           signif_y_manual_calc = F,
           custom_ylim = c(-2, 3),
           file_name_pref = "stats01v24052801_",
           file_name_temp = "t1_g1g2_v169_stdz_box",
           file_name_suff = "r01v24052801_d24042501",
           x_label = "Параметр",
           y_label = "Уровень",
           fill_label = "Группа",
           plot_tit = "")

######################## g1
t1_g1_v169 <- t1_stat %>%
  filter(v003 == 1) %>%
  select(-v160, -v250, -v003) %>%
  dplyr::rename(group=v169) %>%
  mutate(group=factor(group, levels = c(0, 1),
                      labels = c("Нет симпт.",
                                 "Есть симпт.")))


## Standard data
t1_g1_v169_stdz <- cbind(group = t1_g1_v169[, 1],
                         as.data.frame(scale(t1_g1_v169[, -1])))

## Fill NA
t1_g1_v169_stdz <- cbind(t1_g1_v169_stdz[1],
                         filling_na_wrapper(t1_g1_v169_stdz[-1], groupped = F, fill_by = "mean"))

t1_g1_v169_stdz_wide <- t1_g1_v169_stdz

t1_g1_v169_stdz_long <- t1_g1_v169_stdz_wide %>%
  gather(metabolite, level, -group) %>%
  mutate(metabolite = factor(metabolite, levels = init_names$name2,
                             labels = init_names$name1))

data_stats(t1_g1_v169_stdz_long, 
           t1_g1_v169_stdz_wide,
           signif_y_manual_calc = F,
           custom_ylim = c(-2, 3),
           file_name_pref = "stats01v24052801_",
           file_name_temp = "t1_g1_v169_stdz_box",
           file_name_suff = "r01v24052801_d24042501",
           x_label = "Параметр",
           y_label = "Уровень",
           fill_label = "Группа",
           plot_tit = "")


######################## g2
t1_g2_v169 <- t1_stat %>%
  filter(v003 == 2) %>%
  select(-v160, -v250, -v003) %>%
  dplyr::rename(group=v169) %>%
  mutate(group=factor(group, levels = c(0, 1),
                      labels = c("Нет симпт.",
                                 "Есть симпт.")))


## Standard data
t1_g2_v169_stdz <- cbind(group = t1_g2_v169[, 1],
                         as.data.frame(scale(t1_g2_v169[, -1])))

t1_g2_v169_stdz_wide <- t1_g2_v169_stdz

t1_g2_v169_stdz_long <- t1_g2_v169_stdz_wide %>%
  gather(metabolite, level, -group) %>%
  mutate(metabolite = factor(metabolite, levels = init_names$name2,
                             labels = init_names$name1))

data_stats(t1_g2_v169_stdz_long, 
           t1_g2_v169_stdz_wide,
           signif_y_manual_calc = F,
           custom_ylim = c(-2, 3),
           file_name_pref = "stats01v24052801_",
           file_name_temp = "t1_g2_v169_stdz_box",
           file_name_suff = "r01v24052801_d24042501",
           x_label = "Параметр",
           y_label = "Уровень",
           fill_label = "Группа",
           plot_tit = "")


# Грина Вазомоторные симптомы -----------------------------------------------

#### g1g2
t1_g1g2_v250 <- t1_stat %>%
  select(-v160, -v169, -v003) %>%
  dplyr::rename(group=v250) %>%
  mutate(group=factor(group, levels = c(0, 1),
                      labels = c("Нет симпт.",
                                 "Есть симпт.")))


## Standard data
t1_g1g2_v250_stdz <- cbind(group = t1_g1g2_v250[, 1],
                           as.data.frame(scale(t1_g1g2_v250[, -1])))

## Fill NA
t1_g1g2_v250_stdz <- cbind(t1_g1g2_v250_stdz[1],
                           filling_na_wrapper(t1_g1g2_v250_stdz[-1], groupped = F, fill_by = "mean"))

t1_g1g2_v250_stdz_wide <- t1_g1g2_v250_stdz

t1_g1g2_v250_stdz_long <- t1_g1g2_v250_stdz_wide %>%
  gather(metabolite, level, -group) %>%
  mutate(metabolite = factor(metabolite, levels = init_names$name2,
                             labels = init_names$name1))

data_stats(t1_g1g2_v250_stdz_long, 
           t1_g1g2_v250_stdz_wide,
           signif_y_manual_calc = F,
           custom_ylim = c(-2, 3),
           file_name_pref = "stats01v24052801_",
           file_name_temp = "t1_g1g2_v250_stdz_box",
           file_name_suff = "r01v24052801_d24042501",
           x_label = "Параметр",
           y_label = "Уровень",
           fill_label = "Группа",
           plot_tit = "")

######################## g1
t1_g1_v250 <- t1_stat %>%
  filter(v003 == 1) %>%
  select(-v160, -v169, -v003) %>%
  dplyr::rename(group=v250) %>%
  mutate(group=factor(group, levels = c(0, 1),
                      labels = c("Нет симпт.",
                                 "Есть симпт.")))


## Standard data
t1_g1_v250_stdz <- cbind(group = t1_g1_v250[, 1],
                         as.data.frame(scale(t1_g1_v250[, -1])))

## Fill NA
t1_g1_v250_stdz <- cbind(t1_g1_v250_stdz[1],
                         filling_na_wrapper(t1_g1_v250_stdz[-1], groupped = F, fill_by = "mean"))

t1_g1_v250_stdz_wide <- t1_g1_v250_stdz

t1_g1_v250_stdz_long <- t1_g1_v250_stdz_wide %>%
  gather(metabolite, level, -group) %>%
  mutate(metabolite = factor(metabolite, levels = init_names$name2,
                             labels = init_names$name1))

data_stats(t1_g1_v250_stdz_long, 
           t1_g1_v250_stdz_wide,
           signif_y_manual_calc = F,
           custom_ylim = c(-2, 3),
           file_name_pref = "stats01v24052801_",
           file_name_temp = "t1_g1_v250_stdz_box",
           file_name_suff = "r01v24052801_d24042501",
           x_label = "Параметр",
           y_label = "Уровень",
           fill_label = "Группа",
           plot_tit = "")


######################## g2
t1_g2_v250 <- t1_stat %>%
  filter(v003 == 2) %>%
  select(-v160, -v169, -v003) %>%
  dplyr::rename(group=v250) %>%
  mutate(group=factor(group, levels = c(0, 1),
                      labels = c("Нет симпт.",
                                 "Есть симпт.")))


## Standard data
t1_g2_v250_stdz <- cbind(group = t1_g2_v250[, 1],
                         as.data.frame(scale(t1_g2_v250[, -1])))

t1_g2_v250_stdz_wide <- t1_g2_v250_stdz

t1_g2_v250_stdz_long <- t1_g2_v250_stdz_wide %>%
  gather(metabolite, level, -group) %>%
  mutate(metabolite = factor(metabolite, levels = init_names$name2,
                             labels = init_names$name1))

data_stats(t1_g2_v250_stdz_long, 
           t1_g2_v250_stdz_wide,
           signif_y_manual_calc = F,
           custom_ylim = c(-2, 3),
           file_name_pref = "stats01v24052801_",
           file_name_temp = "t1_g2_v250_stdz_box",
           file_name_suff = "r01v24052801_d24042501",
           x_label = "Параметр",
           y_label = "Уровень",
           fill_label = "Группа",
           plot_tit = "")










# testing area ------------------------------------------------------------







