rm(list=ls())

package_list = c('dplyr', 'here', 'ggplot2', 'tidyr','irr')
for (package in package_list) {
  require(package, character.only = TRUE); 
  library(package, character.only = TRUE) }


filepath_isotypes <- 'N:/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/processed_data/DRAFT_CAVD_G002_Glycan_Microarray_data_processed_2024-10-16.txt'
df_glycan_isotypes <- read.delim(filepath_isotypes, header=T)

usecols <- c('sample_id',
             'isotype', 
             'ptid', 
             'study_week',
             'spot_name', 
             'glycan_m_number', 
             'background_subtraced_mean_signal')

df <- df_glycan_isotypes[,usecols]

centered_mean <- function(x){
  if (length(x) >= 6) {
    return(mean(sort(x)[2:5], na.rm=TRUE))
  } else {
    return(mean(x, na.rm=TRUE))
  }
}

df <- df %>% group_by(sample_id, isotype, glycan_m_number) %>%
  mutate(centered_mean = centered_mean(background_subtraced_mean_signal)) %>%
  ungroup()

# drop sample level results
df <- select(df, -background_subtraced_mean_signal)

# unique rows only
df <- unique(df)

# create table for boxplots
calc_responses <- df
calc_responses$threshold <- 100
calc_responses$response_flag <- calc_responses$centered_mean > 100

calc_responses <- calc_responses %>%
  group_by(isotype, study_week, glycan_m_number) %>%
  mutate(count_of_responses = sum(response_flag, na.rm = TRUE)) %>%
  ungroup()

calc_responses <- calc_responses %>%
  group_by(isotype, study_week, glycan_m_number) %>%
  mutate(prop_of_responses = mean(response_flag, na.rm = TRUE)) %>%
  ungroup()

calc_responses <- calc_responses %>%
  group_by(isotype, glycan_m_number) %>%
  mutate(overall_response_rate = mean(response_flag, na.rm = TRUE)) %>%
  ungroup()

calc_responses <- calc_responses %>%
  arrange(desc(overall_response_rate))

agreement <- calc_responses[
  ,c('isotype','glycan_m_number','ptid','study_week','response_flag')] %>%
  distinct()

agreement <- agreement %>%
  pivot_wider(
    names_from=study_week,
    values_from=response_flag
  ) %>%
  select(isotype, glycan_m_number, ptid, `Wk 0`, `Wk 8`, `Wk 10`)

colnames(agreement) <- gsub(" ", "", tolower(colnames(agreement)))

agreement <- agreement %>% group_by(isotype, glycan_m_number) %>%
  mutate(
    wk0_wk8_kappa = kappa2(na.omit(cbind(wk0, wk8)))$value,
    wk0_wk8_pval = kappa2(na.omit(cbind(wk0, wk8)))$p.value,
    wk0_wk10_kappa = kappa2(na.omit(cbind(wk0, wk10)))$value,
    wk0_wk10_pval = kappa2(na.omit(cbind(wk0, wk10)))$p.value,
    ) %>%
  ungroup()

agreement <- agreement %>% select(-ptid, -wk0, -wk8, -wk10) %>%
  distinct()

agreement <- agreement %>% 
  filter(rowSums(is.na(agreement[c('wk0_wk8_kappa','wk0_wk10_kappa')])) < 2)

agreement <- agreement %>%
  arrange(isotype, desc(wk0_wk8_kappa))

savepath = "N:/vtn/lab/SDMC_labscience/operations/documents/templates/assay/template_testing/agreement.csv"
write.csv(agreement, savepath, row.names=FALSE, na='NaN')
