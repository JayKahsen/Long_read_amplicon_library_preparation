################################################################################
message("1 # SETUP")
################################################################################

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(normalizePath(rstudioapi::getSourceEditorContext()$path)))
  }
  getwd()
}

source(file.path(get_script_dir(), "_globalStuff.R"))

qLoad("data_tables/meta_ALL.csv")
meta <- meta_ALL %>%
  filter(data_class == "doubleton") %>%
  mutate(
    sample_name = sample_name_old,
    short_name = short_name_old
  )
meta_names <- names(meta)

script_title <- "1 QIIME tables to working files"
script_description <- "
Creates the doubleton-only long-read working data products from the original
QIIME2 exports. The script writes matrix_names.csv, raw matrices, rarefied
matrices, sample_reads.csv, and final ASV reference files trimmed to the ASVs
present in the rarefied ASV matrix. It also rewrites meta.csv and a formatted
meta.xlsx using only the samples retained in the rarefied data.
"

original_path <- "data_tables/original"
raw_matrix_path <- "data_tables/raw_matrix"
rarefied_matrix_path <- "data_tables/rarefied_matrix"

asv_taxa_all_file <- "data_tables/ASV_taxa_all.csv"
asv_taxa_final_file <- "data_tables/ASV_taxa.csv"
asv_sequence_all_file <- "data_tables/ASV_sequence_all.csv"
asv_sequence_final_file <- "data_tables/ASV_sequence.csv"
meta_final_file <- "data_tables/meta.csv"
meta_excel_file <- "data_tables/meta.xlsx"
taxa_collapse_levels <- c("Phylum", "Genus")

selected_groups <- c("Skin", "Feces", "Soil")
experiment_keep <- c("G4", "G6", "G7")
data_set_order <- data_sets <- "Set1"

minimum_counts <- 1000
keep_number <- 1
use_auto_minimum_count <- "yes"

################################################################################
message("40 # END SETUP")
################################################################################


################################################################################
message("42 # MATRIX REGISTRY")
################################################################################

create_directory(raw_matrix_path)
create_directory(rarefied_matrix_path)

names(taxa_plural) <- taxa_levels

matrix_names <- expand.grid(taxa_levs = taxa_levels, data_sets = data_sets) %>%
  mutate(
    taxa_plural = case_when(
      taxa_levs %in% names(taxa_plural) ~ taxa_plural[taxa_levs],
      TRUE ~ taxa_levs
    ),
    raw_path = paste0(raw_matrix_path, "/raw_matrix_", taxa_levs, ".csv"),
    rarefied_path = paste0(rarefied_matrix_path, "/rarefied_matrix_", taxa_levs, "_", data_sets, ".csv"),
    file_path = rarefied_path
  ) %>%
  arrange(taxa_levs)

write.csv(matrix_names, "data_tables/matrix_names.csv", row.names = FALSE)


################################################################################
message("70 # RAW MATRICES")
################################################################################

message("74 # Build ASV raw matrix")

raw_matrix_asv <- read.delim(
  file.path(original_path, "Doubleton", "Table_Features_PostRemoval.tsv"),
  check.names = FALSE,
  row.names = 1,
  skip = 1
) %>%
  t() %>%
  as.data.frame()

raw_matrix_asv1 <- raw_matrix_asv %>%
  rownames_to_column(var = "sample_ID") %>%
  left_join(meta) %>%
  select(any_of(meta_names), everything()) %>%
  filter(primer_set == "Standard") %>%
  filter(sample_type %in% selected_groups) %>%
  filter(experiment %in% experiment_keep) %>%
  filter(!is.na(sample_name)) %>%
  imPlode_sample_name() %>%
  select(which(colSums(.) > 0))

sample_reads_raw <- raw_matrix_asv1 %>%
  mutate(raw_reads = rowSums(.)) %>%
  xPlode_sample_name() %>%
  select(sample_name, raw_reads)

raw_matrix_asv1 %>%
  t() %>%
  as.data.frame() %>%
  write.csv(paste0(raw_matrix_path, "/raw_matrix_ASV.csv"))

################################################################################
message("130 # RAREFIED MATRICES")
################################################################################

auto_min_counts <- function(sample_count_dfx, keep_number) {
  sample_count_dfx %>%
    group_by(short_name) %>%
    arrange(desc(sample_count)) %>%
    slice(1:keep_number) %>%
    ungroup() %>%
    pull(sample_count) %>%
    min()
}

rarefy_function <- function(dfx, min_sample_count) {
  df <- dfx %>%
    filter(sample_count >= min_sample_count) %>%
    select(-sample_count) %>%
    imPlode_sample_name() %>%
    as.matrix()

  set.seed(20250314)

  rrarefy(df, sample = min_sample_count) %>%
    as.data.frame()
}

collapse_asv_matrix <- function(asv_matrix_df, asv_taxa_df, taxa_level) {
  taxa_lookup <- asv_taxa_df[[taxa_level]][match(rownames(asv_matrix_df), asv_taxa_df$ASV)]
  taxa_lookup[is.na(taxa_lookup) | taxa_lookup == ""] <- "Unassigned"

  rowsum(as.matrix(asv_matrix_df), group = taxa_lookup, reorder = TRUE) %>%
    as.data.frame()
}

matrix_df_asv <- read.csv("data_tables/raw_matrix/raw_matrix_ASV.csv", check.names = FALSE, row.names = 1) %>%
  t() %>%
  as.data.frame()

data_set_df_asv <- matrix_df_asv %>%
  mutate(sample_count = rowSums(.)) %>%
  xPlode_sample_name() %>%
  filter(sample_type %in% selected_groups) %>%
  ungroup()

rarefied_df_asv <- NULL

for (sample_type_name in unique(data_set_df_asv$sample_type)) {
  qPrint(sample_type_name)

  subset_df <- data_set_df_asv %>%
    filter(sample_type == sample_type_name)

  min_sample_count <- minimum_counts
  if (use_auto_minimum_count == "yes") {
    min_sample_count <- auto_min_counts(subset_df, keep_number)
  }

  if (sample_type_name == "Skin") min_sample_count <- 35000
  if (sample_type_name == "Feces") min_sample_count <- 25000
  if (sample_type_name == "Soil") min_sample_count <- 10000

  rarefied_subset <- rarefy_function(subset_df, min_sample_count)

  rarefied_df_asv <- bind_rows(rarefied_df_asv, rarefied_subset) %>%
    mutate(across(everything(), ~ replace_na(., 0)))
}

rarefied_df_asv %>%
  t() %>%
  as.data.frame() %>%
  write.csv(paste0(rarefied_matrix_path, "/rarefied_matrix_ASV_Set1.csv"))

sample_reads_rarefied <- rarefied_df_asv %>%
  mutate(rarefied_reads = rowSums(.)) %>%
  xPlode_sample_name() %>%
  select(sample_name, rarefied_reads)

sample_reads_df <- sample_reads_raw %>%
  left_join(sample_reads_rarefied, by = "sample_name") %>%
  arrange(sample_name)

write.csv(sample_reads_df, "data_tables/sample_reads.csv", row.names = FALSE)

retained_meta <- meta %>%
  filter(sample_name %in% rownames(rarefied_df_asv)) %>%
  mutate(sample_name = factor(sample_name, levels = rownames(rarefied_df_asv))) %>%
  arrange(sample_name) %>%
  mutate(sample_name = as.character(sample_name))

write.csv(retained_meta, meta_final_file, row.names = FALSE)
write_formatted_excel(retained_meta, meta_excel_file, sheet_name = "meta")


################################################################################
message("210 # FINAL ASV REFERENCES + COLLAPSED MATRICES")
################################################################################

rarefied_matrix_asv_df <- read.csv("data_tables/rarefied_matrix/rarefied_matrix_ASV_Set1.csv", check.names = FALSE, row.names = 1)

final_asvs <- rownames(rarefied_matrix_asv_df)

asv_taxa <- read.csv(asv_taxa_all_file, check.names = FALSE) %>%
  filter(ASV %in% final_asvs) %>%
  distinct(ASV, .keep_all = TRUE)
write.csv(asv_taxa, asv_taxa_final_file, row.names = FALSE)

if (file.exists(asv_sequence_all_file)) {
  asv_sequence <- read.csv(asv_sequence_all_file, check.names = FALSE) %>%
    filter(ASV %in% final_asvs) %>%
    distinct(ASV, .keep_all = TRUE)
  write.csv(asv_sequence, asv_sequence_final_file, row.names = FALSE)
}

raw_matrix_asv_df <- read.csv("data_tables/raw_matrix/raw_matrix_ASV.csv", check.names = FALSE, row.names = 1)

for (taxa_level in taxa_collapse_levels) {
  qPrint(taxa_level)

  raw_collapsed <- collapse_asv_matrix(raw_matrix_asv_df, asv_taxa, taxa_level)
  rarefied_collapsed <- collapse_asv_matrix(rarefied_matrix_asv_df, asv_taxa, taxa_level)

  write.csv(raw_collapsed, paste0("data_tables/raw_matrix/raw_matrix_", taxa_level, ".csv"))
  write.csv(rarefied_collapsed, paste0("data_tables/rarefied_matrix/rarefied_matrix_", taxa_level, "_Set1.csv"))
}


################################################################################
message(paste("280 # FINISHED", script_title))
################################################################################
