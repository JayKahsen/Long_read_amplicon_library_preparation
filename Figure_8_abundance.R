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

source(file.path(get_script_dir(), "_globalStuff.R"))

script_title <- "abundance"
script_description <- "
Creates Figure 8 abundance outputs from the rarefied Phylum matrix.
Runs the phylum-level comparisons and writes the combined abundance figures.
"

set_output()

plot_set_order <- c("normalization")

###############################################################################
message("24 # END SETUP")
###############################################################################

###############################################################################
message("28 # RUNS")
###############################################################################

for (plot_set in plot_set_order) { qPrint(plot_set)

  ###############################################################################
  message("36 # Plot-set outputs")
  ###############################################################################
  
  df_matrix <- NULL
  
  if (plot_set == "standard") {
    normalization_order <- c("fixed cycles")
    data_class_order    <- c("doubleton")
    normalization_name  <- paste(normalization_order, collapse = "_")
    data_class_name     <- paste(data_class_order, collapse = "_")
  }
  
  if (plot_set == "normalization") {
    normalization_order <- c("fixed cycles", "targeted fluorescence")
    data_class_order    <- c("doubleton")
    normalization_name  <- paste(normalization_order, collapse = "_")
    data_class_name     <- paste(data_class_order, collapse = "_")
    figure_title        <- "Figure 8 The Effects of Auto-Normalization on Taxonomic Profile"
  }
  
  if (plot_set == "data_class") {
    normalization_order <- c("targeted fluorescence")
    data_class_order    <- c("singleton", "doubleton")
    normalization_name  <- paste(normalization_order, collapse = "_")
    data_class_name     <- paste(data_class_order, collapse = "_")
  }
  
  ###############################################################################
  message("58 # Dataset selection")
  ###############################################################################
  
  starting_r <- 1
  ending_r   <- 1
  testing    <- "yes"
  testing_r  <- 1
  
  selected_data_class      <- "doubleton"
  selected_data_class_name <- paste(selected_data_class, collapse = "_")
  zoom <- 2
  
  ###############################################################################
  message("68 # Data processing settings")
  ###############################################################################
  
  use_custom_labels <- "no"
  
  ###############################################################################
  message("74 # Parameter sets")
  ###############################################################################
  
  parameter_sets <- list(
    set1 = list(
      filter1_group  = "data_class",
      y_axis_group   = "normalization",
      x_facet_group  = "temperature",
      y_facet_group  = "sample_type",
      number_size    = 1
    ),
    set3 = list(
      filter1_group  = "normalization",
      y_axis_group   = "data_class",
      x_facet_group  = "temperature",
      y_facet_group  = "sample_type",
      number_size    = 1
    )
  )
  
  ################################################################################
  message("102 # MAIN LOOP: PROCESS DATASETS")
  ################################################################################
  
  Loop_Group_order <- Loop_Group <- "run_once"
  if (exists("Run_Group")) {
    Loop_Group       <- Run_Group
    Loop_Group_order <- Run_Group_order
  }
  
  for (lp in Loop_Group_order) {
    
    if (!exists("ending_r")) {
      ending_r <- nrow(matrix_names)
      message("\tsetting ending r to nrows matrix names")
    }
    if (!exists("starting_r")) {
      starting_r <- 1
      message("\tsetting starting r to 1")
    }
    if (testing == "yes") {
      r <- starting_r <- ending_r <- testing_r
      message("\trunning r= ", testing_r)
    }
    
    for (r in starting_r:ending_r) {
      
      #........................ dataset details .................................#
      
      taxa_levs   <- matrix_names[r, "taxa_levs"]
      data_set    <- matrix_names[r, "data_sets"]
      taxa_plural <- matrix_names[r, "taxa_plural"]
      
      if (plot_set != "none") {
        data_set       <- plot_set
        data_set_order <- plot_set_order
      }
      
      # Assign analysis parameters dynamically
      list2env(parameter_sets[[paste0("set", match(data_set, data_set_order))]], envir = .GlobalEnv)
      
      groups <- c("filter1_group", "y_axis_group", "x_facet_group", "y_facet_group")
      for (var in groups) {
        group_value <- get(var)
        assign(paste0(var, "_order"), get(paste0(group_value, "_order")))
      }
      
      grouping_columns <- c(y_axis_group, x_facet_group, y_facet_group)
      test_across      <- c(x_facet_group, y_facet_group)
      
      run_suffix <- ""
      if (exists("Run_Group")) run_suffix <- paste0("_", lp)
      
      filter1_names <- paste(filter1_group_order, collapse = "_")
      y_axis_names  <- paste(y_axis_group_order,  collapse = "_")
      x_facet_names <- paste(x_facet_group_order, collapse = "_")
      
      p_title <- paste0(
        taxa_levs, "_", script_title, custom_name, "_",
        filter1_names, "_", y_axis_names, "_", x_facet_names, "_",
        data_set, run_suffix
      )
      
      qPrint(p_title)
      
      #==============================================================================#
      message("124 # load and filter data")
      #==============================================================================#
      
      df_matrix1 <- read.csv(matrix_names[r, "file_path"], check.names = FALSE, row.names = 1) %>%
        t() %>% as.data.frame() %>%
        xPlode_sample_name() %>%
        filter(.data[[filter1_group]]  %in% filter1_group_order) %>%
        filter(.data[[y_facet_group]]  %in% y_facet_group_order) %>%
        filter(.data[[x_facet_group]]  %in% x_facet_group_order) %>%
        filter(.data[[y_axis_group]]   %in% y_axis_group_order)
      
      if (exists("Run_Group")) {
        df_matrix1 <- df_matrix1 %>% filter(.data[[Loop_Group]] %in% lp)
      }
      
      df_matrix <- df_matrix1 %>%
        imPlode_sample_name() %>%
        mutate_all(as.numeric) %>%
        select(which(colSums(.) > 0))
      
      feature_names <- names(df_matrix)
      
      #============================== Statistical tests ==========================#
      
      kw_group_results <- Kruskal_Mann_Whitney_Test(
        df = df_matrix,
        testing_group = y_axis_group,
        category_group = test_across,
        collapse_data_used = TRUE
      )
      
      write.csv(kw_group_results, paste0("output_data/kw_", plot_set, ".csv"), row.names = FALSE)
      
      #==============================================================================#
      message("137 # data processing and normalization")
      #==============================================================================#
      
      feature_labels_df  <- feature_labels(y_facet_group)
      y_axis_group_order <- rev(y_axis_group_order)
      
      df_plot <- df_matrix %>%
        xPlode_sample_name() %>%
        pivot_longer(cols = any_of(feature_names), names_to = "feature", values_to = "counts") %>%
        group_by(feature, across(any_of(y_facet_group))) %>%
        mutate(feature_counts_y_facet = sum(counts)) %>%
        filter(feature_counts_y_facet > 0) %>%
        group_by(sample_name) %>%
        mutate(sample_counts = sum(counts)) %>%
        group_by(across(any_of(y_facet_group))) %>%
        mutate(counts_y_facet_mean = mean(sample_counts)) %>%
        ungroup() %>%
        mutate(rel_abun = counts / sample_counts) %>%
        group_by(feature, across(any_of(grouping_columns))) %>%
        summarize(
          counts_y_facet_mean         = mean(counts_y_facet_mean),
          feature_counts_group_mean   = mean(counts),
          counts_sample_mean          = mean(sample_counts),
          rel_abun_group_mean         = mean(rel_abun),
          feature_counts_group_sd     = sd(counts),
          .groups = "drop"
        ) %>%
        mutate(
          feature_counts_group_sd_min = feature_counts_group_mean - feature_counts_group_sd,
          feature_counts_group_sd_max = feature_counts_group_mean + feature_counts_group_sd,
          log10_value_rel = log10(rel_abun_group_mean * counts_y_facet_mean + 1e-6),
          log10_value     = log10(feature_counts_group_mean + 1),
          min_log10       = log10(pmax(feature_counts_group_sd_min, 1e-6)),
          max_log10       = log10(feature_counts_group_sd_max + 1e-6)
        ) %>%
        left_join(
          feature_labels_df %>% dplyr::select(-c(taxa, group_mean_rel_abun)),
          by = c("feature", y_facet_group)
        ) %>%
        mutate(!!sym(y_facet_group) := factor(!!sym(y_facet_group), levels = rev(y_facet_group_order))) %>%
        mutate(!!sym(y_axis_group)  := factor(!!sym(y_axis_group),  levels = y_axis_group_order)) %>%
        mutate(!!sym(x_facet_group) := factor(!!sym(x_facet_group), levels = x_facet_group_order)) %>%
        arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group)))) %>%
        ungroup()
      
      df_plot_mn <- df_plot %>%
        group_by(feature, across(any_of(test_across))) %>%
        mutate(rel_abun_panel_mean = mean(rel_abun_group_mean)) %>%
        ungroup() %>%
        mutate(rel_abun_norm = rel_abun_group_mean / rel_abun_panel_mean) %>%
        left_join(kw_group_results, by = c("feature", test_across)) %>%
        ungroup() %>%
        mutate(log2_value = log2(rel_abun_norm + 1e-6)) %>%
        mutate(!!sym(y_facet_group) := factor(!!sym(y_facet_group), levels = rev(y_facet_group_order))) %>%
        mutate(!!sym(y_axis_group)  := factor(!!sym(y_axis_group),  levels = y_axis_group_order)) %>%
        mutate(!!sym(x_facet_group) := factor(!!sym(x_facet_group), levels = x_facet_group_order)) %>%
        arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group)))) %>%
        ungroup()
      
      #==============================================================================#
      message("164 # generate plots")
      #==============================================================================#
      
      if (taxa_levs == "Phylum") { }
      
      margin_size          <- 10
      caption_text_size    <- 20
      axis_text_size       <- 16
      axis_title_text_size <- 18
      strip_text_size      <- 20
      legend_text_size     <- 20
      rel_x_text_size      <- 1.7
      
      common_theme <- theme(
        strip.text = element_text(face = "bold", size = strip_text_size),
        strip.text.y = element_text(face = "bold", size = strip_text_size),
        plot.margin = margin(margin_size, margin_size, margin_size, margin_size),
        plot.background = element_rect(fill = palette_color["general_background"], color = palette_color["general_background"]),
        panel.background = element_rect(fill = "gray90"),
        legend.background = element_rect(fill = palette_color["general_background"]),
        legend.text = element_text(size = legend_text_size),
        axis.title.y = element_text(face = "bold", size = axis_title_text_size, margin = margin(r = 10)),
        axis.title.x = element_text(face = "bold", size = axis_title_text_size, margin = margin(t = 10)),
        axis.text.x = element_markdown(face = "bold", size = rel(rel_x_text_size)),
        plot.caption = element_text(hjust = 0.5)
      )
      
      gPlot <- function(p) {
        # Example: p <- gPlot(p)
        p <- p +
          common_theme +
          scale_fill_manual(values = palette_color, guide = guide_legend(reverse = TRUE), labels = palette_label) +
          scale_color_manual(values = palette_color, guide = guide_legend(reverse = TRUE), labels = palette_label) +
          labs(title = "", fill = "", color = "", x = "", y = "")
        print(p)
        p
      }
      
      #==============================================================================#
      message("207 # strip backgrounds with outlines")
      #==============================================================================#
      
      unique_y_labels <- unique(df_plot[[y_facet_group]])
      outline_color_y <- "black"
      s <- paste0(
        'backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',
        paste0(unique_y_labels, collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),
        '"]))'
      )
      print(s)
      eval(parse(text = s))
      
      if (use_custom_labels == "yes") {
        custom_y_labels <- c("slope","1500TF","x3FC","slope","1500TF","x3FC","24cycles",rep("Soil",3),rep("Zymo",3),"NTC")
        unique_y_labels <- custom_y_labels
        s <- paste0(
          'backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',
          paste0(unique_y_labels, collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),
          '"]))'
        )
        print(s)
        eval(parse(text = s))
      }
      
      y_strip <- as.data.frame(unique_y_labels) %>%
        mutate(text_color = "white") %>%
        mutate(face = "bold") %>%
        mutate(text_size = strip_text_size) %>%
        ungroup()
      
      unique_x_labels <- unique(df_plot[[x_facet_group]])
      outline_color_x <- "black"
      s <- paste0(
        'backgrounds_x <- list(element_rect(color=outline_color_x,fill = palette_color["',
        paste0(unique_x_labels, collapse='"]),element_rect(color=outline_color_x,fill = palette_color["'),
        '"]))'
      )
      print(s)
      eval(parse(text = s))
      
      x_strip <- as.data.frame(unique_x_labels) %>%
        mutate(text_color = "white") %>%
        mutate(face = "bold") %>%
        mutate(text_size = strip_text_size) %>%
        ungroup()
      
      #==============================================================================#
      message("248 # plot variables")
      #==============================================================================#
      
      x_min <- 2
      x_alt <- 1
      
      plot_st <- ggplot(df_plot, aes(x = log10_value, y = plot_order)) +
        geom_col(aes(fill = !!sym(y_axis_group)), width = 0.6, position = position_dodge(width = 0.8)) +
        geom_errorbarh(
          aes(xmin = log10_value, xmax = max_log10, group = !!sym(y_axis_group)),
          position = position_dodge(width = 0.8),
          height = 0.2
        ) +
        scale_y_discrete(breaks = custom_breaks, labels = custom_labels) +
        facet_grid2(
          formula(paste(y_facet_group, "~", x_facet_group)),
          strip = strip_themed(
            background_y = backgrounds_y,
            text_y = elem_list_text(face = y_strip$face, size = y_strip$text_size, color = y_strip$text_color),
            background_x = backgrounds_x,
            text_x = elem_list_text(face = x_strip$face, size = x_strip$text_size, color = x_strip$text_color)
          ),
          scale = "free",
          space = "free"
        ) +
        scale_x_continuous(labels = log10_labels_bold, breaks = log10_breaks) +
        coord_cartesian(xlim = c(0, NA))
      
      plot_st <- gPlot(plot_st)
      
      plot_mn <- ggplot(df_plot_mn, aes(x = log2_value, y = plot_order, fill = !!sym(y_axis_group))) +
        geom_col(width = 0.8, position = position_dodge(width = 0.8)) +
        geom_text(aes(x = Inf, label = adjusted_significance, color = adjusted_significance), hjust = 1.1, fontface = "bold") +
        geom_text(aes(x = 1.5, label = significance, color = significance), hjust = 0) +
        geom_vline(aes(xintercept = 0), color = "magenta", linetype = "dashed") +
        facet_grid2(
          formula(paste(y_facet_group, "~", x_facet_group)),
          strip = strip_themed(
            background_y = backgrounds_y,
            text_y = elem_list_text(face = y_strip$face, size = y_strip$text_size, color = y_strip$text_color),
            background_x = backgrounds_x,
            text_x = elem_list_text(face = x_strip$face, size = x_strip$text_size, color = x_strip$text_color)
          ),
          scale = "free",
          space = "free"
        ) +
        scale_x_continuous(labels = fold_change_labels_bold, breaks = fold_change_breaks) +
        scale_y_discrete(breaks = custom_breaks, labels = custom_labels, position = "right") +
        coord_cartesian(xlim = c(-zoom, zoom)) +
        theme(legend.position = "none")
      
      if (taxa_levs == "Phylum") {
        
        y_labels <- setNames(df_plot$feature_label2, df_plot$plot_order)
        
        plot_st <- plot_st +
          geom_text(
            aes(
              label = signif(feature_counts_group_mean, 3),
              x = ifelse(group_mean_feature_counts < x_min, x_alt, log10_value),
              color = ifelse(group_mean_feature_counts < x_min, "black", "white"),
              group = !!sym(y_axis_group)
            ),
            size = number_size,
            hjust = 1.1,
            fontface = "bold",
            show.legend = FALSE,
            position = position_dodge(width = 0.8)
          ) +
          scale_y_discrete(labels = y_labels) +
          theme(axis.text.y = element_markdown(face = "bold", size = axis_text_size))
        
        plot_mn <- plot_mn +
          scale_y_discrete(labels = y_labels, position = "right") +
          theme(
            axis.text.y.right = element_markdown(face = "bold", size = axis_text_size),
            legend.position = "none"
          )
      }
      
      plot_st <- plot_st + theme(plot.tag = element_text(size = 36, face = "bold"), legend.position = "bottom")
      plot_mn <- gPlot(plot_mn) + theme(plot.tag = element_text(size = 36, face = "bold"))
      
      #==============================================================================#
      message("321 # combine + save")
      #==============================================================================#
      
      combined_plot <- plot_st + plot_mn +
        plot_layout(nrow = 1, widths = c(1, 1), guides = "collect") &
        plot_annotation(
          tag_levels = "A",
          caption = "ns = not significant; * p < 0.05; ** p < 0.01; *** p < 0.001",
          theme = theme(
            plot.caption = element_text(size = caption_text_size, hjust = 0.5, face = "bold"),
            plot.background = element_rect(fill = palette_color["general_background"], color = "gray30", linewidth = 2),
            legend.position = "bottom",
            legend.justification = "center"
          )
        )
      
      print(combined_plot)
      
      plot_width  <- 36
      plot_height <- 24
      qSave(figure_title, plot = combined_plot, ext = ext_list)
      
    } # end loop over r
  } # end loop over lp
} # end loop over plot_set

###############################################################################
message("446 # END RUNS")
###############################################################################

###############################################################################
message(paste("450 # FINISHED", script_title))
###############################################################################

#============================= meta summary table ==============================#

m <- meta %>%
  filter(data_class == "singleton") %>%
  filter(primer_set == "Standard") %>%
  filter(temperature != "42°C") %>%
  filter(sample_type != "NTC") %>%
  filter(normalization == "targeted fluorescence") %>%
  select(sample_type, temperature, cycle) %>%
  pivot_wider(names_from = temperature, values_from = cycle) %>%
  ungroup()

m
