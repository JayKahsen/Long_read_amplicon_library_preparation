###############################################################################
message("1 # SETUP")
###############################################################################

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

script_title <- "global"
script_description <- "
Shared setup for the long-read amplicon balancing workflow.
Loads helper functions, required libraries, data tables, palette objects,
factor orders, and common plotting wrappers used by the run-first and figure scripts.
"

set.seed(20240725)
options(scipen = 999)
options(repos = c(CRAN = "https://cloud.r-project.org"))
.libPaths(unique(c("C:/_Repositories/R_libs", .libPaths())))

ext_list <- c("png", "pdf", "svg")
dpi <- 300

wd <- get_script_dir()
project_folder <- basename(wd)
setwd(wd)

source(file.path(wd, "helperJ.R"))

required_libraries <- c(
  "tidyverse",
  "vegan",
  "compositions",
  "patchwork",
  "cowplot",
  "ggh4x",
  "ggtext",
  "markdown",
  "litedown",
  "svglite",
  "openxlsx"
)

load_libraries(required_libraries)
qPrint(project_folder)

###############################################################################
message("44 # END SETUP")
###############################################################################

###############################################################################
message("48 # LOAD FILES")
###############################################################################

qLoad("data_tables/meta.csv")
qLoad("data_tables/matrix_names.csv")
qLoad("data_tables/ASV_taxa.csv")

if (exists("meta")) {
  meta_names <- names(meta)
  qPrint(meta_names)
}

###############################################################################
message("60 # END LOAD FILES")
###############################################################################

###############################################################################
message("64 # VARIABLES")
###############################################################################

taxa_levels <- c("Phylum", "Genus", "ASV")
taxa_plural <- c("Phyla", "Genera", "ASVs")

if (exists("matrix_names")) {
  data_set_order <- data_sets <- unique(matrix_names$data_sets)
}

set_output()

###############################################################################
message("76 # END VARIABLES")
###############################################################################

###############################################################################
message("80 # ORDERS")
###############################################################################

plate_order <- paste0("plate_", 1:7)

temperature_order <- c("47°C", "52°C", "57°C")

primer_set_order <- primer_order <- c("Standard", "Truncated53", "Truncated50")

sample_type_order <- c("Soil", "Feces", "Skin", "NTC")
sample_type_order <- sample_type_order[1:3]

experiment_order    <- c("G4", "G6", "G7", "G9", "G10")
mastermix_order     <- c("AMPED", "QB")
normalization_order <- c("fixed cycles", "slope", "fold change", "targeted fluorescence")
data_class_order    <- c("singleton", "doubleton")

custom_name    <- ""
plot_set_order <- c("standard", "normalization", "data_class")

cycle            <- 1:24
plot_cycle_order <- paste0("cycle_", str_pad(cycle, width = 2, side = "left", pad = "0"))
plot_cycle_order

###############################################################################
message("103 # END ORDERS")
###############################################################################

###############################################################################
message("107 # PALETTES")
###############################################################################

palette_gray   <- rev(c("#BBBBBB", "#999999", "#777777", "#555555", "#333333"))
palette_green  <- rev(c("#D2EDB2", "#89C189", "#60A740", "#107810", "#004B00"))
palette_gold   <- rev(c("#FFDA65", "#FFCB25", "#DAA520", "#A37B18", "#7D6220"))
palette_purple <- rev(c("#D2AEFA", "#A267E8", "#7A1FD2", "#5B17A1", "#3C1070"))
palette_red    <- rev(c("#FFCCCC", "#E47C7C", "#C94141", "#993333", "#6A2B2B"))
palette_blue   <- rev(c("#C8D6FF", "#92AEFF", "#5C86FF", "#3A62BE", "#20407F"))

palette_color <- c(
  "gray1" = "gray80",
  "gray2" = "gray70",
  
  "general_background" = "#FffaFa",
  
  "panel_background" = "gray90",
  "plot_outline"     = "gray30",
  
  "cumulative" = palette_gray[4],
  "work_flow"  = palette_common[1],
  
  "47°C" = palette_blue[3],
  "52°C" = palette_gold[3],
  "57°C" = palette_red[3],
  
  "Skin"  = palette_purple[1],
  "Feces" = palette_gold[1],
  "Soil"  = palette_gray[1],
  
  "fixed cycles"          = palette_red[1],
  "targeted fluorescence" = palette_green[1],
  
  "doubleton" = "darkblue",
  "singleton" = "darkred",
  "Doubleton" = "darkblue",
  "Singleton" = "darkred",
  
  "ns"  = "gray50",
  "sig" = "red",
  
  "Archaea"   = "#FF4040",
  "Bacteria"  = "#000080",
  "Eukaryote" = "darkgreen",
  "unknown"   = "gray30",
  
  is.na = "gray",
  
  "NTC" = "gray75",
  
  palette_helper_color
)

palette_label <- c(
  "cumulative" = "Remaining Read %",
  "work_flow"  = "Work Flow",
  "targeted fluorescence" = "Targeted Fluorescence",
  
  "chimeric_reads"    = "Chimeric Reads",
  "contaminant_reads" = "Contaminant Reads",
  "denoised_reads"    = "Denoised Reads",
  "fastq_reads"       = "Fastq Reads",
  "filtered_reads"    = "Filtered Reads",
  "pentaton_reads"    = "Pentaton Reads",
  "trimmed_reads"     = "Trimmed Reads",
  "usable_reads"      = "Usable Reads",
  
  "p_chimeric_reads" = "% Chimera",
  "p_usable_reads"   = "% Usable Reads",
  
  "doubleton" = "Doubleton",
  "singleton" = "Singleton",
  "Doubleton" = "Doubleton",
  "Singleton" = "Singleton",
  
  "fixed cycles" = "Fixed Cycles",
  "fixed_cycles" = "Fixed Cycles",
  
  "A" = "A",
  "C" = "C",
  "G" = "G",
  "T" = "T",
  
  "highlight_on"  = "Real data",
  "highlight_off" = "Primer",
  
  "g__Cutibacterium" = "Cutibacterium",
  "no" = "not detected",
  
  "G11"         = "G11",
  "Standard"    = "Standard",
  "Truncated50" = "Truncated50",
  "Truncated53" = "Truncated53",
  
  "42°C" = "42°C",
  "47°C" = "47°C",
  "52°C" = "52°C",
  "57°C" = "57°C",
  
  "Skin"  = "Skin",
  "Feces" = "Feces",
  "Soil"  = "Soil",
  "NTC"   = "NTC",
  
  "AMPED" = "AMPED",
  "QB"    = "QB",
  
  "Slope" = "Slope",
  "Targeted Fluorescence 1500" = "Targeted",
  "Fold Change x3" = "Fold",
  
  "slope"    = "Slope",
  "1500TF"   = "Targeted",
  "x3FC"     = "Fold",
  "24cycles" = "Fixed",
  
  "unique_features" = "Unique ASVs",
  "chimeras"        = "% Chimera",
  
  "percent_Muri"    = "% Muribaculaceae",
  "percent_Lachno"  = "% Lachnospiraceae",
  "percent_Propi"   = "% Cutibacterium",
  "percent_Archaea" = "% Archaea",
  "percent_good"    = "% Good",
  
  "richness" = "Richness",
  "evenness" = "Evenness",
  "shannon"  = "Shannon",
  "n"        = "Reads",
  
  "fastq_reads"            = "fastq_reads",
  "trimmed_reads"          = "trimmed_reads",
  "filtered_reads"         = "filtered_reads",
  "denoised_reads"         = "denoised_reads",
  "chimeric_reads"         = "chimeric_reads",
  "contaminant_table_reads"= "contaminant_table_reads",
  "contaminant_reads"      = "contaminant_reads",
  "usable_reads"           = "usable_reads",
  
  "p_fastq_reads"             = "p_fastq_reads",
  "p_trimmed_reads"           = "p_trimmed_reads",
  "p_filtered_reads"          = "p_filtered_reads",
  "p_denoised_reads"          = "p_denoised_reads",
  "p_chimeric_reads"          = "p_chimeric_reads",
  "p_contaminant_reads"       = "p_contaminant_reads",
  "p_contaminant_table_reads" = "p_contaminant_table_reads",
  "p_usable_reads"            = "p_usable_reads",
  
  "c_fastq_reads"             = "c_fastq_reads",
  "c_trimmed_reads"           = "c_trimmed_reads",
  "c_filtered_reads"          = "c_filtered_reads",
  "c_denoised_reads"          = "c_denoised_reads",
  "c_chimeric_reads"          = "c_chimeric_reads",
  "c_contaminant_reads"       = "c_contaminant_reads",
  "c_contaminant_table_reads" = "c_contaminant_table_reads",
  "c_usable_reads"            = "c_usable_reads",
  
  palette_helper_label
)

palette_shape <- c(
  "doubleton" = 22,
  "singleton" = 21,
  "Doubleton" = 22,
  "Singleton" = 21,
  
  "Standard" = 22,
  
  "targeted fluorescence" = 23,
  "fixed cycles"          = 22,
  "targeted_fluorescence" = 23,
  "fixed_cycles"          = 22,
  
  "42°C" = 21,
  "47°C" = 21,
  "52°C" = 24,
  "57°C" = 23,
  
  "Targeted Fluorescence 1500" = 3,
  "Fold Change x3" = 0,
  
  "default" = 1
)

palette_size <- c(
  "doubleton" = 5,
  "singleton" = 5,
  "Doubleton" = 5,
  "Singleton" = 5,
  
  "targeted fluorescence" = 6,
  "fixed cycles"          = 5,
  "targeted_fluorescence" = 6,
  "fixed_cycles"          = 5,
  
  "Standard"    = 3,
  "Truncated50" = 3,
  "Truncated53" = 3,
  
  "42°C" = 5,
  "47°C" = 5,
  "52°C" = 5,
  "57°C" = 6,
  
  "AMPED" = 3,
  "QB"    = 3,
  
  "slope"    = 5,
  "1500TF"   = 5,
  "x3FC"     = 5,
  "24cycles" = 5,
  
  "Truncated"  = 3,
  "TruncatedR" = 3,
  "TruncatedF" = 3,
  
  "default" = 3
)

palette_linetype <- c(
  "Standard"    = "solid",
  "Truncated50" = "dotted",
  "Truncated53" = "dashed",
  
  "Zymo" = "dashed"
)

###############################################################################
message("334 # END PALETTES")
###############################################################################

###############################################################################
message("338 # THEMES")
###############################################################################

margin_size <- 10

theme_global <- function(base_size = 11, magnify = 1) {
  theme_gray(base_size = base_size) +
    theme(
      plot.margin = margin(margin_size, margin_size, margin_size, margin_size),
      plot.background = element_rect(fill = palette_color["general_background"], color = palette_color["general_background"]),
      panel.background = element_rect(fill = palette_color["panel_background"]),
      legend.background = element_rect(fill = palette_color["general_background"]),
      axis.title = element_text(size = (base_size + 2) * magnify, face = "bold"),
      axis.text = element_text(size = base_size * magnify, face = "bold"),
      legend.title = element_text(size = (base_size + 2) * magnify, face = "bold"),
      legend.text = element_text(size = base_size * magnify),
      strip.text = element_text(size = (base_size + 2) * magnify, face = "bold"),
      plot.title = element_text(size = (base_size + 6) * magnify, face = "bold"),
      plot.subtitle = element_text(size = (base_size + 2) * magnify),
      plot.caption = element_text(size = (base_size + 2) * magnify, face = "bold")
    )
}

theme_common <- theme_global(base_size = 11) + theme()

theme_plot <- theme(
  plot.background = element_rect(fill = palette_color["general_background"], color = palette_color["plot_outline"], linewidth = 2)
)

gPlot <- function(p) {
  p <- p +
    scale_fill_manual(values = palette_color, labels = palette_label) +
    scale_color_manual(values = palette_color, labels = palette_label) +
    scale_linetype_manual(values = palette_linetype, labels = palette_label) +
    scale_shape_manual(values = palette_shape, labels = palette_label)
  print(p)
  p
}

###############################################################################
message("378 # END THEMES")
###############################################################################

###############################################################################
message(paste("382 # FINISHED", script_title))
###############################################################################

cite_R(required_libraries)


