
make_heatmaps <- function(context_filter, context_DMR_list_file, data_path, max_rows_to_include, h, w){
###########
# args

# context_filter = "CG"
# data_path <- "analysis/tiles_CT_GA_filtered_4C_2x_12pctA_DSS_DMRs_samples_6_DMR_lists_ratios/Ratios_for_DMRs_CG_tiles"   # path to the data
# context_DMR_list_file <- "analysis/tiles_CT_GA_filtered_4C_2x_12pctA_DSS_DMRs_samples_6_DMR_lists/DMRs_CG_tiles.tsv"
# max_rows_to_include = 40000
# outFolder = paste0(data_path, "_heatmaps")
# h = 10
# w = 6

###########
  outFolder = paste0(data_path, "_heatmaps")
dir.create(outFolder, showWarnings = F)

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
  # filter(!Sample %in% c("Average_GD", "Average_ST")) %>%
  select(Sample, tile, ratio) %>%
  group_by(tile) %>%
  mutate(ratio_scaled = round(ratio/max(ratio)*100,digits = 1)) %>%
  ungroup()

# Heatmaps
# 1. All data heatmap
# 1. All data heatmap normalised
# 2. Only tiles with data in all/specified number of samples
# 3. Tiles with data normalised

# 2. Only tiles with data in all/specified number of samples
big_data_summary_wide <- big_data_summary %>%
  select(Sample, tile, ratio) %>%
  spread(key = Sample, value = ratio) %>%
   na.omit() 

# check if there are more than 50,000 rows, if so subsample to 50,000
frac_50K <- max_rows_to_include/length(big_data_summary_wide$tile)

big_data_summary_wide_subset <- big_data_summary_wide

if(frac_50K < 1){
set.seed(27)
big_data_summary_wide_subset <- big_data_summary_wide_subset %>% sample_frac(frac_50K)
}

heatmap_data <- data.frame(big_data_summary_wide_subset)
row.names(heatmap_data) <- heatmap_data[,1]
heatmap_data <- data.matrix(heatmap_data[,-c(1)])

# probably put the samples in a better order?

pheat_colour = rev(brewer.pal(11,"RdBu"))

# distance.row = dist(heatmap_data, method = "euclidean")
# cluster.row = hclust(distance.row, method = "ward.D")

p <- pheatmap(heatmap_data, silent = T, cluster_rows=TRUE, cluster_cols=T, color = pheat_colour, show_rownames = F,  show_colnames = TRUE, fontsize = 10, na_col = "purple")
save_pheatmap_pdf_size(p, paste0(outFolder, "/", context_filter, "_", "ratio", "_pheatmap", "_max_", max_rows_to_include, ".pdf"), height = h/2.54, width = w/2.54)
save_pheatmap_png_size_cm(p, paste0(outFolder, "/", context_filter, "_", "ratio", "_pheatmap", "_max_", max_rows_to_include, ".png"), height = h, width = w)

# save col order?
# col_order_heatmap <- p$tree_col$order
}
