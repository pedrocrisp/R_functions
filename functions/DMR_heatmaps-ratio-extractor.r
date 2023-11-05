############ ############ ############ ############ ############ ############ ############ ############ ############ ############ ############ ############ 
# Function #4 extract ratios data for DMR tile list including an out group
# Should probably use the same criteria here as were used for DMR calling eg number of CG sutes and coverage 
ratio_extractor <- function(context, DMR_list, data_folder, sample_list_to_analyse, DMR_list_folder){
############ read in file stack
  # ARGS
  # context= "CG"
  # data_folder = "analysis/tiles_CT_GA_filtered_4C_2x_12pctA"
  # sample_list_to_analyse <-"samples_6.txt"
  # DMR_list_folder <- "analysis/tiles_CT_GA_filtered_4C_2x_12pctA_DSS_DMRs_samples_6"
  # DMR_list = "DMRs_CG_tiles.tsv"

############ setup
sample_list_subset <- read_tsv(sample_list_to_analyse, col_names = c("sample"))

# list the folder with the individual sample outout files output from the make ratios script, should end "_ratio_annotated.csv"
# data_path <- paste0(data_folder, "/")   # path to the data

# DMR list
DMR_tile_list <- read_tsv(file.path(DMR_list_folder, DMR_list))

outFolder_parent = paste0(DMR_list_folder, "_ratios")
dir.create(outFolder_parent, showWarnings = F)

outFolder = paste0(DMR_list_folder, "_ratios", "/Ratios_for_", file_path_sans_ext(DMR_list))
dir.create(outFolder, showWarnings = F)

##############
# subset files
samples_to_filter = pull(sample_list_subset, sample)

ratio_filter_per_sample <- function(target_sample){
  # args
  # target_sample = samples_to_filter[1]

target_sample

ratio_data <- read_tsv(file.path(data_folder, paste0(target_sample, "_BSMAP_out.txt.100.", context,"_filtered.txt")))

sites_col <- colnames(ratio_data %>% select(contains("sites")))

# If I want to add another filter t could be done here
ratio_data_filter <- ratio_data %>%
  mutate(tile = paste0(chr, ":", start, ":", start + 99), ratio = round((C/CT)*100, digits = 1))  %>%
  # filter((!!as.name(sites_col)) >= 2, cov >= 2) %>%
  filter(tile %in% DMR_tile_list$tile) %>%
  select(tile, ratio)

ratio_data_filter

write_tsv(ratio_data_filter, paste0(outFolder, "/",target_sample, "_BSMAP_out.txt.100.", context,"_filtered_DMRs.txt"))

}

walk(samples_to_filter, ratio_filter_per_sample)

}