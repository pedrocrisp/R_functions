############ ############ ############ ############ ############ ############ ############ ############ ############ ############ ############ ############ 
# Function # 2 get list of distinct DMRs per context
common_distinct_tile_list <- function(context, prefix, data_path, sample_list, clustering_group, min_samples_with_DMR){
  
  # context = "CG"
  # prefix = "*CG_DMRs.csv"
  # data_path <- "analysis/tiles_CT_GA_filtered_4C_2x_12pctA_DSS_DMRs"
  # sample_list = "samples_6.txt"

############ chunk options
#folder to put the results
# experiment <- "apple_mC"
clustering_group <- sub("*.txt", "", sample_list)
out_folder <- paste0(data_path, "_", clustering_group, "_DMR_lists_common_", min_samples_with_DMR)
# dir.create(out_folder)

# all samples
# sample_list_to_analyse <-"samples_HC_and_GD_minus3.txt"
# just MN samples?
sample_list_to_analyse <- sample_list

############ setup
dir.create(out_folder, showWarnings = F)
sample_list_subset <- read_tsv(sample_list_to_analyse, col_names = c("sample"))

############ read in file stack
# list the folder with the individual sample outout files output from the make ratios script, should end "_ratio_annotated.csv"
# data_path <- "analysis/tiles_filtered_2_samp_MN_DMR_samples"   # path to the data

# list the files
files_all <- dir(data_path, pattern = prefix) # get file names

# limit to comparisons of interest
files <- tibble(files_all) %>% 
  separate(col = files_all, into = c("sample1", "sample2"), sep = '.vs.', remove = F) %>%
  mutate(sample2 = sub("_C[HG].*.csv", "", sample2)) %>%
  filter(sample1 %in% sample_list_subset$sample, sample2 %in% sample_list_subset$sample) %>% 
  # dplyr::slice(., 1:3) %>%
  pull(files_all)
files

# Read in as a list of dfs
data <- data_frame(filename = files) %>% # create a data frame
                                         # holding the file names
  mutate(file_contents = purrr::map(filename,          # read files into
           ~ read_csv(file.path(data_path, .), col_names = T, cols_only(
  chr = col_character(),
  pos = col_integer()
))) # a new data column
        )  
data

############ make a giant combined table 
# unnest the dfs
big_data <- unnest(data)
big_data

# distinct tile list for DMRs in x number of samples - 
# for filtering to common DMRs - this is a hacky way to try extract common RoyalRed DMRs vs HC mother
DMR_tiles_2 <- big_data %>% 
  group_by(chr, pos) %>%
  mutate(number_of_contrasts = n()/2) %>%
  ungroup() %>%
  filter(number_of_contrasts >= min_samples_with_DMR) %>%
  distinct(chr, pos) %>% mutate(end = pos + 99) %>% mutate(tile = paste0(chr, ":", pos, ":", end)) %>% select(tile)

write_tsv(DMR_tiles_2, paste0(out_folder, "/DMRs_", context, "_tiles_in_", min_samples_with_DMR, "samples", ".tsv"))

}
