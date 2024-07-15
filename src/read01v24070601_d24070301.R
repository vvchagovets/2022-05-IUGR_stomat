# Positive ions processing

# Import Libraries --------------------------------------------------------
library(dplyr)
library(stringr)

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

cn_dbn_parser <- function(lipid_name) {
  cn_dbn_data <- str_extract_all(lipid_name, "\\d+:\\d+")
  
  carbon_number <- 0
  double_bond_number <- 0
  for (i in cn_dbn_data[[1]]) {
    cn_i <- str_extract(i, "^\\d+")
    dbn_i <- str_extract(i, "\\d+$")
    carbon_number <- carbon_number + as.numeric(cn_i)
    double_bond_number <- double_bond_number + as.numeric(dbn_i)
  }
  cn_dbn_result <- data.frame(id =  lipid_name,
                              cn = carbon_number,
                              dbn = double_bond_number)
  return(cn_dbn_result)
  
}

ms_file_name_parser <- function(ms_data) {
  columns_number  = ncol((ms_data))
  full_name_dot_splitted <- str_split_fixed(colnames(ms_data)[seq(5, columns_number - 6)], "\\.", 3)
  
  full_name_only <- data.frame(str_split_fixed(full_name_dot_splitted[, 1], "_", 8),
                                  stringsAsFactors = F)
  
  for(i in seq(nrow(full_name_only))){

    if (grepl("^[[:alpha:]]{2}", full_name_only[i, 3])) {
      full_name_only[i, 5] <- str_extract(full_name_only[i, 3], "^[[:alpha:]]{2}")
      full_name_only[i, 3] <- str_replace(full_name_only[i, 3], "[[:alpha:]]{2}", '')
      
    } else if(grepl("^[[:digit:]]{2}", full_name_only[i, 3])){
      full_name_only[i, 5] <- "Stomat"
    }
    
    
    if(full_name_only[i, "X2"] == "qc0"){
      full_name_only[i, seq(4, 8)] <- full_name_only[i, seq(3, 7)]
      full_name_only[i, 3] <- full_name_only[i, 2]
      
    }
  }
  ms_file_meta <- full_name_only 

  
  colnames(ms_file_meta) <- c("in_batch_num", "sample_type",
                                 "sample_id", "vial_coord", "s1", 
                                 "ms_id", "s2", "polarity")
  
  ms_file_meta <- ms_file_meta %>%
    mutate(file_id = paste(sample_type, s1, sample_id, sep = "_")) %>%
    mutate(sample_id = as.numeric(sample_id))
  
  return(ms_file_meta)
  
}

# 
# ms_file_name_parser_pl <- function(ms_data) {
#   columns_number  = ncol((ms_data))
#   full_name_dot_splitted <- str_split_fixed(colnames(ms_data)[seq(5, columns_number - 6)], "\\.", 3)
#   
#   full_name_only <- data.frame(str_split_fixed(full_name_dot_splitted[, 1], "_", 10),
#                                stringsAsFactors = F)
#   
#   for(i in seq(nrow(full_name_only))){
#     if(full_name_only[i, "X8"] == "Sample"){
#       full_name_only[i, seq(6, 10)] <- full_name_only[i, seq(5, 9)]
#     }
#     
#     if(grepl("^qc", full_name_only[i, "X2"])){
#       full_name_only[i, seq(5, 10)] <- full_name_only[i, seq(2, 7)]
#       full_name_only[i, 3] <- full_name_only[i, 2]
#       
#     }
#   }
#   
#   ms_file_meta <- full_name_only %>%
#     mutate(X5 = ifelse(X5 == 1, 1, 0))
#   
#   colnames(ms_file_meta) <- c("in_batch_num", "sample", "sample_type", 
#                               "sample_id", "is_repeat", "vial_coord", "s1", 
#                               "ms_id", "s2", "polarity")
#   
#   ms_file_meta <- ms_file_meta %>%
#     mutate(file_id = ms_id) %>%
#     mutate(sample_id = as.numeric(sample_id))
#   
#   return(ms_file_meta)
#   
# }



ms_part_from_init_data <- function(ms_data, meta_data) {
  ms_only <- ms_data[, seq(5, ncol(ms_data) - 6)]
  
  colnames(ms_only) <- c(meta_data$file_id)
  return(ms_only) 

}

lipid_name_parser <- function(ms_data){
  lipid_ids <- sapply(ms_data$ID_Ranked, function(x) strsplit(x, "\\|"))
  selected_lipids <- sapply(lipid_ids, function(x) x[[1]][[1]])
  names(selected_lipids) = ""
  
  id_rank_tmp1 <- sapply(selected_lipids, function(x) str_extract(x, "\\d_"))
  id_rank_tmp2 <- str_replace(id_rank_tmp1, "_", "")
  id_no_rank <- sapply(selected_lipids, function(x) str_replace(x, "\\d_", ""))
  
  id_no_rank_clean <- sapply(id_no_rank, function(x) strsplit(x, ";")[[1]][[1]])
  id_no_rank_clean <- str_replace(id_no_rank_clean, "MSDIAL_", "")
  ion_type <- sapply(id_no_rank_clean, function(x) str_extract(x, "\\+[[:graph:]]+"))
  id_no_rank_no_ion <- sapply(id_no_rank_clean, function(x) str_replace(x, "\\+[[:graph:]]+", ""))
  
  lipids_number <- length(id_no_rank)
  final_digits_number <- str_length(as.character(lipids_number))
  id_alias <- list()
  
  for(i in seq(lipids_number)){
    digits_number <- str_length(as.character(i))
    zeroes <- paste(rep('0', final_digits_number-digits_number), collapse='')
    current_element <- paste0('l', zeroes, i)
    id_alias <- c(id_alias, current_element)
  }
  
  lipids_meta <- data.frame(id_alias = unlist(id_alias),
                            mz = ms_data$row.m.z,
                            rt = ms_data$row.retention.time,
                            id_rank = as.numeric(id_rank_tmp2),
                            id = id_no_rank_no_ion,
                            ion_type = ion_type,
                            stringsAsFactors = F)
  return(lipids_meta)
}

lipids_meta_extended_parser <- function(lipids_meta){
  lipid_classes <- sapply(lipids_meta$id, function(x) str_extract(x, "^[:alpha:]+[-]?[:alpha:]+"))
  
  lipid_cn_dbn <- lapply(lipids_meta$id, function(x) cn_dbn_parser(x))
  cn_dbn_tmp <- do.call(rbind, lipid_cn_dbn) %>%
    select(-id)
  
  lipids_meta_extended <- lipids_meta %>% 
    mutate(lipid_class = lipid_classes) %>%
    cbind(cn_dbn_tmp)
  return(lipids_meta_extended)
  
}

transpose_ms <- function(ms_part, lipids_meta){
  ms_data_clean <- cbind(lipid_id = lipids_meta$id_alias,
                            ms_part)
  
  ms_data_clean_t <- t(ms_data_clean)
  colnames(ms_data_clean_t) <- ms_data_clean_t[1,]
  
  lipids_data <- data.frame(ms_data_clean_t[-1, ])
  lipids_data <- data.frame(sapply(lipids_data, function(x) as.numeric(as.character(x))))
  lipids_data <- cbind(file_id = names(ms_part),
                       lipids_data)
  return(lipids_data)
}

# End Functions -----------------------------------------------------------

# Read in Data ------------------------------------------------------------
init_ms_cm <-
  read.table(
    './data/24070301/preproc/PosIDed.csv',
    header = T,
    sep = ','
  )


init_clin_data <-
  read.table(
    './data/24070301/preproc/clin_data.csv',
    header = T,
    sep = ','
  )

init_clin_names <-
  read.table(
    './data/24070301/preproc/clin_names.csv',
    header = T,
    sep = ','
  )

init_biobank_data <-
  read.table(
    './data/24070301/preproc/biobank_data.csv',
    header = T,
    sep = ','
  )

# Clean FL MS Data -----------------------------------------------------------

## Clean lipid file names; Form MS meta
ms_file_meta_cm <- ms_file_name_parser(init_ms_cm)

ms_part_only_cm <- ms_part_from_init_data(init_ms_cm, ms_file_meta_cm)
  
# Clean lipid names ---
lipids_meta_cm <- lipid_name_parser(init_ms_cm)
lipids_meta_extended_cm <- lipids_meta_extended_parser(lipids_meta_cm)

## transpose ms data
ms_data_transposed_cm <- transpose_ms(ms_part_only_cm, lipids_meta_cm)

## Attach meta info
ms_data_clean_cm <-  ms_file_meta_cm %>%
  left_join(ms_data_transposed_cm, by = "file_id") %>%
  filter(sample_type == "Urine") %>%
  select(-in_batch_num, -vial_coord,
         -s1, -ms_id, -s2, -polarity, -sample_type, -sample_id)



# Join clin and ms-meta ---
init_biobank_data <- init_biobank_data %>%
  dplyr::rename(sample_id = ms_id)
ms_biobank_meta <- ms_file_meta_cm %>%
  filter(s1 == 'Stomat') %>%
  left_join(init_biobank_data, by = 'sample_id')

init_clin_data <- init_clin_data %>%
  dplyr::rename(name = v003)
clin_ms_meta <- ms_biobank_meta %>%
  left_join(init_clin_data, by = 'name') %>%
  filter(!is.na(v001))


# Add IUGR Norm meta
reduced_norm_meta <- ms_file_meta_cm %>%
  filter(s1 != 'Stomat') %>%
  filter(sample_type != 'qc0') %>%
  select(sample_type, sample_id, file_id) %>%
  mutate(group = 'N')
reduced_clin_ms_meta <- clin_ms_meta %>%
  dplyr::rename(group = v001) %>%
  select(sample_type, sample_id, file_id, group)

joint_ms_meta <- rbind(reduced_clin_ms_meta, reduced_norm_meta)

#### Final MS data
joint_ms_meta_cleaned <- joint_ms_meta %>%
  select(file_id, group)

ms_data_cleaned <- joint_ms_meta_cleaned %>%
  left_join(ms_data_clean_cm, by='file_id')
  


# ms_n_clin_cm <- clin_cm %>%
#   right_join(ms_data_clean_pl, by = "sample_id")

## right ovary matching
# clin_cm <- init_clin_data %>%
#   select(v003, v009) %>%
#   dplyr::rename(sample_id = v003) %>%
#   dplyr::rename(group = v009)
# 
# ms_n_clin_cm <- clin_cm %>%
#   right_join(ms_data_clean_pl, by = "sample_id")



# lipids_to_save_cm <- data.frame(mz = round(lipids_meta_cm$mz, 2),
#                              rt = round(lipids_meta_cm$rt, 2),
#                              id = lipids_meta_cm$id)
# write.csv(lipids_to_save_cm, './writeup/!lipids_cm.csv')
# 
# lipids_to_save_cm <- data.frame(mz = round(lipids_meta_cm$mz, 2),
#                                 rt = round(lipids_meta_cm$rt, 2),
#                                 id = lipids_meta_cm$id)
# write.csv(lipids_to_save_cm, './writeup/!lipids_cm.csv')

# preprocessed lipid data to transfer -------------------------------------



# write.csv(lipids_meta_extended_cm, './transfer/preproc_lipids/lipids_meta_cm.csv')
# write.csv(ms_data_clean_cm, './transfer/preproc_lipids/ms_n_clin_cm.csv')
# 
# 
# # Normalize on TIC
# lipids_data_norm_pl <- ms_n_clin_pl %>% 
#   mutate(row_sum = rowSums(select(., 4:ncol(ms_n_clin_pl)))) %>% 
#   mutate_at(4:ncol(ms_n_clin_pl), ~ ./row_sum) %>% 
#   select(-row_sum)
# 
# lipids_data_norm_fl <- ms_n_clin_fl %>% 
#   mutate(row_sum = rowSums(select(., 5:ncol(ms_n_clin_fl)))) %>% 
#   mutate_at(5:ncol(ms_n_clin_fl), ~ ./row_sum) %>% 
#   select(-row_sum)
# 
# write.csv(lipids_data_norm_fl, './transfer/preproc_lipids/ms_n_clin_norm_fl.csv')
# write.csv(lipids_data_norm_pl, './transfer/preproc_lipids/ms_n_clin_norm_pl.csv')

# Testing Area ------------------------------------------------------------






