################## summary function

get_summary_distribution <- function(context_filter, context_DMR_list_file, data_path, outFolder){
###########
# args

# context_filter = "CG"
# data_path <- "analysis/tiles_CT_GA_filtered_4C_2x_12pctA_DSS_DMRs_samples_6_DMR_lists_ratios/Ratios_for_DMRs_CG_tiles"   # path to the data
# context_DMR_list_file <- "analysis/tiles_CT_GA_filtered_4C_2x_12pctA_DSS_DMRs_samples_6_DMR_lists/DMRs_CG_tiles.tsv"

###########

# list the files
files <- dir(data_path, pattern = "*filtered_DMRs.txt") # get file names

files

# Read in as a list of dfs
data <- data_frame(filename = files) %>% # create a data frame
                                         # holding the file names
  mutate(file_contents = purrr::map(filename,          # read files into
           ~ read_tsv(file.path(data_path, .))) # a new data column
        )  
data

############ make a giant combined table 
# unnest the dfs
big_data <- unnest(data)
big_data

distinct(big_data, tile)
# 372,964 tiles with data

###### context DMR list

context_DMR_list <- read_tsv(context_DMR_list_file)

######## all contexts missing data

big_data_summary <- big_data %>%
  separate(filename, into = c("Sample", "context"), sep = "_BSMAP_out.txt.100.") %>%
  mutate(context = sub("_filtered_DMRs.txt", "", context)) %>%
  filter(context == context_filter, tile %in% context_DMR_list$tile) %>%
  filter(!Sample %in% c("Average_GD", "Average_ST")) %>%
  select(Sample, tile, ratio) %>%
  group_by(tile) %>%
  summarise(number_of_samples = n()) %>%
  group_by(number_of_samples) %>%
  summarise(number_of_tiles = n())

big_data_summary <- big_data_summary %>% 
  mutate(pct = number_of_tiles/sum(number_of_tiles)*100)

big_data_summary

write.csv(big_data_summary, paste0(outFolder, "/", context_filter, "_sample_with_data_distribution.csv"))

}

