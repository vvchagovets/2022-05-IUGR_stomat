
# Import Libraries --------------------------------------------------------
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(reshape2)

rm(list = setdiff(ls(), c()))
Sys.setlocale(category = "LC_ALL", locale = "Russian, Russia")

# Functions ---------------------------------------------------------------
filter_na <- function(data_unfiltered, na_allowed = 0.4, by_group = 0){
  g1 <- data_unfiltered %>%
    filter(group == 1)
  g2 <- data_unfiltered %>%
    filter(group == 2)
  na_number_g1 <- nrow(g1) * na_allowed
  na_number_g2 <- nrow(g2) * na_allowed
  logic_g1 <- colSums(is.na(g1)) < na_number_g1
  logic_g2 <- colSums(is.na(g2)) < na_number_g2
  if(by_group == T){
    data_filtered <- data_unfiltered[, logic_g1 & logic_g2]
  }else{
    data_filtered <- data_unfiltered[, logic_g1 | logic_g2]
  }
  return(data_filtered)
}

filter_na_single_group <- function(data_unfiltered, na_allowed = 0.4){
  na_number <- nrow(data_unfiltered) * na_allowed
  logic_g1 <- colSums(is.na(data_unfiltered)) < na_number
  data_filtered <- data_unfiltered[, logic_g1]
  
  return(data_filtered)
}


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

numeric_group_processing <- function(data, group_name){
  numerical_g <- data %>%
    filter(group == group_name) %>%
    select(-group)
  
  numerical_data_na_filtered <- filter_na_single_group(numerical_g, na_allowed = 0.3)
  
  numerical_data_filled_na <- filling_na_wrapper(numerical_data_na_filtered, F, "mean")
  return(numerical_data_filled_na)
  
}

categorical_group_processing <- function(data, group_name){
  categorical_g <- data %>%
    filter(group == group_name) %>%
    select(-group)
  
  categorical_data_na_filtered <- filter_na_single_group(categorical_g, na_allowed = 0.3)
  
  categorical_data_filled_na <- filling_na_wrapper(categorical_data_na_filtered, F, "mode")
  return(categorical_data_filled_na)
  
}

processing_ms_id <- function(init_id){
  corrected_id <- str_extract(init_id, "[[:alpha:]]{2}_\\d+")
  corrected_id <- str_replace(corrected_id, "_", "")
  return(corrected_id)
  
}

processing_biobank_id <- function(init_id){
  corrected_id <- str_extract(init_id, "[[:alpha:]]{2}\\d+")
  return(corrected_id)
}

normalize_on_total <- function(data){
  
  numeric_data <- data[, -1]
  normalized <- numeric_data / rowSums(numeric_data)
  result <- cbind(sample = data[, 1], normalized)
  return(result)
  
}


# EDA Functions -----------------------------------------------------------
amount_most_frequent <- function(x){
  vector_stats <- as.data.frame(table(x))
  sorted_stats <- vector_stats[order(vector_stats$Freq, decreasing = T),]
  return(sorted_stats[1, 2])  
}

eda_table <- function(data, names, file){
  # NA count
  clean_na_tf <- data.frame(is.na(data))
  na_amount <- sapply(clean_na_tf, sum)
  na_percent <- round(na_amount / nrow(data) * 100, 0)
  
  # Unique count
  unique_amount <- sapply(data, function(x) length(unique(x)))
  unique_percent <- round(unique_amount / nrow(data) * 100, 0)
  
  # single value part
  
  
  most_frequent_amount <- sapply(data, amount_most_frequent)
  most_frequent_percent <- round(most_frequent_amount / nrow(data) * 100, 0)
  
  eda01_result <- data.frame(na_percent = na_percent,
                             unique_percent = unique_percent,
                             most_frequent_percent = most_frequent_percent)
  
  param_alias <- row.names(eda01_result)
  eda01_result <- cbind(param_alias, eda01_result)
  eda01_final <- eda01_result %>%
    mutate(param_name = factor(param_alias, levels=names$name2,
                               labels=names$name1))
  
  write.csv(eda01_final, paste0("./writeup/", file, ".csv"))
  
}

# End Functions -----------------------------------------------------------


# Read in Data ------------------------------------------------------------
init_aa_data <-
  read.table(
    './data/24042501/preproc/init_data_24042501.csv',
    header = T,
    sep = ',', 
    stringsAsFactors = F
  )

init_aa_names <-
  read.table(
    './data/24042501/preproc/init_names_24042501.csv',
    header = T,
    sep = ',', 
    stringsAsFactors = F
  )



# Clean data -----------------------------------------------------

clin_data_cleaned_01 <- init_aa_data %>%
  select(v003, v004, v005, v014, v017, v019, v040, v041, v042, v043, 
         v044, v045, v046, v047, v048, v049, v050, v051, v052, v053, 
         v054, v055, v056, v057, v058, v059, v060, v061, v062, v063, 
         v064, v065, v079, v080, v081, v157, v158, v159, v160, v161, 
         v162, v163, v164, v165, v166, v167, v168, v169, v170, v171, 
         v172, v173, v174, v241, v242, v243, v244, v245, v246, v250, 
         v247, v248, v249, v251, v252, v253, v254, v255, v258, v262, 
         v263)

eda_table(clin_data_cleaned_01, init_aa_names, "read01v24052801_eda_d24042501")



# separate time1, time2, time3 --------------------------------------------

t1 <- clin_data_cleaned_01 %>%
  select(v003, v004, v005, v014, v017, v019, v040, v043, 
         v044, v053, 
         v056, v057, v063, 
         v079, v157, v160,  
         v163, v166, v169,  
         v172, v241, v244, v247, v250, v253, v258,
         v262, v263)

t2 <- init_aa_data %>%
  select(v003, v004, v005, v014, v017, v019, v041, 
         v045, v046,  
         v054, v058, v059, 
         v064, v080, v158, v161, 
         v164, v167, v170,  
         v173, v242, v245, v248, v251,  
         v254, v258, v262, v263)

t3 <- init_aa_data %>%
  select(v003, v004, v005, v014, v017, v019, v042,  
         v047, v048,  
         v055, v060, v061,  
         v065, v081, v159,  
         v162, v165, v168, v171, 
         v174, v243, v246, v249, v252,  
         v255, v258, v262, v263)


# Testing Area ------------------------------------------------------------


