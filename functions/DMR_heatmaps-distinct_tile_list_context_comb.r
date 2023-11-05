############ ############ ############ ############ ############ ############ ############ ############ ############ ############ ############ ############
# Function #3 get list of distnct DMR regions all contexts
distinct_tile_list_context_comb <- function(data_path){
############ chunk options
# data_path <- "analysis/tiles_filtered_4C_2x_DSS_DMRs"   # path to the data
###########

# list the files
files <- dir(data_path, pattern = "*tiles.tsv") # get file names

files

# Read in as a list of dfs
data <- data_frame(filename = files) %>% # create a data frame
                                         # holding the file names
  mutate(file_contents = map(filename,          # read files into
           ~ read_tsv(file.path(data_path, .))) # a new data column
        )
data

############ make a giant combined table
# unnest the dfs
big_data <- unnest(data)
big_data

############ distinct tiles

DMR_tiles <- big_data %>% select(tile) %>% distinct(tile)

print(length(DMR_tiles$tile))

write_tsv(DMR_tiles, paste0(data_path, "/DMRs_all_contexts_tiles.tsv"))
}
