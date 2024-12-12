rm(list=ls())

#load packages
package_list = c('dplyr', 'here', 'ggplot2', 'tidyr')
for (package in package_list) {
  require(package, character.only = TRUE); 
  library(package, character.only = TRUE) }

# read in data
filepath_isotypes <- 'N:/cavd/VDCs/Schief/Schief_856-G002/SkinReactions/data/Glycan_array_Scripps/processed_data/DRAFT_CAVD_G002_Glycan_Microarray_data_processed_2024-10-16.txt'
df_glycan_isotypes <- read.delim(filepath_isotypes, header=T)

# subset to desired columns
usecols <- c('sample_id',
             'isotype', 
             'ptid', 
             'study_week',
             'spot_name', 
             'glycan_m_number', 
             'background_subtraced_mean_signal')

df <- df_glycan_isotypes[,usecols]

# if there are 6 datapoints per isotype/guspec/glycan, then drop
# outer two and average the middle ones.
# if fewer than 6, just take the average
centered_mean <- function(x){
  if (length(x) >= 6) {
    return(mean(sort(x)[2:5], na.rm=TRUE))
  } else {
    return(mean(x, na.rm=TRUE))
  }
}

df <- df %>% 
  group_by(sample_id, isotype, glycan_m_number) %>%
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

## hi tzu-jung! this isn't the code that generated what we were looking at on that call
## but this should do the same thing.
## calc_responses still has unique rows per each sample/isotype/week/glycan. 
## so i believe you'd get the table you want from:

summary <- calc_responses[,c('isotype', 'study_week', 'glycan_m_number', 'count_of_responses')]
summary <- unique(summary)

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

# save to csv to diff against python version
# savedir <- "N:/vtn/lab/SDMC_labscience/operations/documents/templates/assay/template_testing/test_11_01_2024/"
# write.table(calc_responses, file = file.path(savedir, "calc_responses.txt"), sep = "\t", row.names = FALSE, quote = TRUE)

# Set the isotype and filter for glycan ordering
iso <- "IgG"

# order glycans by overall response rate
glycan_ordering <- calc_responses %>%
  filter(isotype == iso) %>%
  select(glycan_m_number, overall_response_rate) %>%
  distinct() %>%
  pull(glycan_m_number)


dt <- as.data.frame(calc_responses)
dt$study_week <- factor(dt$study_week, levels=c('Wk 0', 'Wk 8', 'Wk 10'), ordered=TRUE)

igg <- dt %>% filter(isotype=="IgG" & glycan_m_number %in% glycan_ordering[1:20] & centered_mean >= 100)

# turn glycan into factor so plots are ordered
igg$glycan_m_number_f = factor(igg$glycan_m_number, levels=glycan_ordering[1:20])

# this plots the first 20 glycans for igg.
ggplot(igg, aes(study_week, centered_mean)) + 
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~glycan_m_number_f)

