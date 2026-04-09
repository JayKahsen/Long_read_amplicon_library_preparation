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

script_title <- "taxonomy"
script_description <- "
Creates Figure 2 taxonomy profile panels from the rarefied Phylum matrix.
Writes the consolidated figure outputs for the standard doubleton workflow.
"

plot_set_order <- "standard"

set_output()

###############################################################################
message("24 # END SETUP")
###############################################################################

###############################################################################
message("28 # RUNS")
###############################################################################
for (plot_set in plot_set_order) {

  ###############################################################################
  message("22 # Plot-set outputs")
  ###############################################################################
  
  if (plot_set == "standard") {
    normalization_order <- c("fixed cycles")
    data_class_order    <- c("doubleton")
    normalization_name  <- paste(normalization_order, collapse = "_")
    data_class_name     <- paste(data_class_order, collapse = "_")
    figure_title        <- "Figure 2 Taxonomy Profiles for Skin Feces and Soil"
  }
  
  if (plot_set == "normalization") {
    normalization_order <- c("fixed cycles", "targeted fluorescence")
    data_class_order    <- c("doubleton")
    normalization_name  <- paste(normalization_order, collapse = "_")
    data_class_name     <- paste(data_class_order, collapse = "_")
  }
  
  if (plot_set == "data_class") {
    normalization_order <- c("targeted fluorescence")
    data_class_order    <- c("singleton", "doubleton")
    normalization_name  <- paste(normalization_order, collapse = "_")
    data_class_name     <- paste(data_class_order, collapse = "_")
  }
  
  ###############################################################################
  message("44 # Dataset selection")
  ###############################################################################
  
  starting_r <- 1
  ending_r   <- 1
  testing    <- "yes"
  testing_r  <- 1
  selected_data_class      <- "doubleton"
  selected_data_class_name <- paste(selected_data_class, collapse = "_")
  
  ###############################################################################
  message("56 # Data processing / variable sets")
  ###############################################################################
  
  use_custom_labels    <- "no"
  min_rel_abun_percent <- 0.1
  
  parameter_sets <- list(
    set1 = list(
      filter1_group  = "data_class",
      y_axis_group   = "temperature",
      x_facet_group  = "normalization",
      y_facet_group  = "sample_type",
      number_size    = 1
    ),
    set2 = list(
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
  
  ###############################################################################
  message("82 # MAIN LOOP: PROCESS DATASETS")
  ###############################################################################
  
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
      message(paste0("\trunning r= ", testing_r))
    }
    
    for (r in starting_r:ending_r) {
      
      #..............................................................................#
      message("110 # dataset details")
      #..............................................................................#
      
      taxa_levs   <- matrix_names[r, "taxa_levs"]
      data_set    <- matrix_names[r, "data_sets"]
      taxa_plural <- matrix_names[r, "taxa_plural"]
      
      if (plot_set != "none") {
        # Adjusting for full loop
        data_set       <- plot_set
        data_set_order <- plot_set_order
      }
      
      # Assign analysis parameters dynamically
      list2env(parameter_sets[[paste0("set", match(data_set, data_set_order))]], envir = .GlobalEnv)
      
      groups <- c("filter1_group", "y_axis_group", "x_facet_group", "y_facet_group")
      for (var in groups) {
        group_value <- get(var)  # e.g., "data_class"
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
      message("144 # Load and Filter Data")
      #==============================================================================#
      
      #..............................................................................#
      message("147 # load and filter data")
      #..............................................................................#
      
      df_matrix1 <- read.csv(matrix_names[r, "file_path"], check.names = FALSE, row.names = 1) %>%
        t() %>% as.data.frame() %>%
        xPlode_sample_name() %>%
        filter(.data[[filter1_group]] %in% filter1_group_order) %>%
        filter(.data[[y_facet_group]] %in% y_facet_group_order) %>%
        filter(.data[[x_facet_group]] %in% x_facet_group_order) %>%
        filter(.data[[y_axis_group]]  %in% y_axis_group_order)
      
      if (exists("Run_Group")) {
        df_matrix1 <- df_matrix1 %>% filter(.data[[Loop_Group]] %in% lp)
      }
      
      df_matrix <- df_matrix1 %>%
        imPlode_sample_name() %>%
        mutate_all(as.numeric) %>%
        select(which(colSums(.) > 0))
      
      #==============================================================================#
      message("170 # Combine unassigned taxonomy (optional)")
      #==============================================================================#
      
      combining_unassigned <- "yes"
      if (combining_unassigned == "yes") {
        
        unassigned_cols <- names(df_matrix)[str_detect(names(df_matrix), ";__$|p__$|;__$")]
        
        if (length(unassigned_cols) > 0) {
          df_matrix$Unassigned <- rowSums(df_matrix[, unassigned_cols, drop = FALSE])
          df_matrix <- df_matrix[, !names(df_matrix) %in% unassigned_cols, drop = FALSE]
        }
      }
      
      dfx_matrix <- df_matrix %>% xPlode_sample_name()
      
      feature_names <- names(df_matrix)
      
      #==============================================================================#
      message("187 # Data Processing and Normalization")
      #==============================================================================#
      
      #..............................................................................#
      message("190 # data processing and normalization")
      #..............................................................................#
      
      feature_labels_df <- feature_labels(y_facet_group)
      y_axis_group_order <- rev(y_axis_group_order)
      
      df <- dfx_matrix %>%
        pivot_longer(cols = any_of(feature_names), names_to = "feature", values_to = "counts") %>%
        group_by(feature, across(any_of(y_facet_group))) %>%
        mutate(grouped_feature_counts = sum(counts)) %>%
        filter(grouped_feature_counts > 0) %>%
        group_by(sample_name) %>%
        mutate(sample_counts = sum(counts)) %>%
        ungroup() %>%
        mutate(rel_abun = counts / sample_counts) %>%
        group_by(feature, across(any_of(y_facet_group))) %>%
        summarize(
          group_mean_counts = mean(counts),
          percent_group_mean_rel_abun = 100 * mean(rel_abun),
          .groups = "drop"
        ) %>%
        ungroup()
      
      df_plot <- df %>%
        filter(percent_group_mean_rel_abun > min_rel_abun_percent) %>%
        group_by(across(any_of(y_facet_group))) %>%
        mutate(summed_percent_group_mean_rel_abun = sum(percent_group_mean_rel_abun)) %>%
        ungroup() %>%
        mutate(group_percentage_abun = 100 * percent_group_mean_rel_abun / summed_percent_group_mean_rel_abun) %>%
        left_join(feature_labels_df %>% select(-c(taxa)), by = c("feature", y_facet_group)) %>%
        arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group)))) %>%
        mutate(!!sym(y_facet_group) := factor(!!sym(y_facet_group), levels = y_facet_group_order)) %>%
        ungroup() %>%
        arrange(desc(group_percentage_abun)) %>%
        mutate(feature_label = factor(feature_label, levels = unique(feature_label))) %>%
        ungroup()
      
      phylum_order <- levels(df_plot$feature_label)
      feature_colors <- setNames(palette_common[seq_along(phylum_order)], phylum_order)

      if ("Unassigned" %in% names(feature_colors)) {
        light_gray_index <- which(feature_colors %in% c(palette_gray[1], palette_gray[2]))
        if (length(light_gray_index) > 0) {
          light_gray_name <- names(feature_colors)[light_gray_index[1]]
          light_gray_value <- feature_colors[light_gray_name]
          feature_colors[light_gray_name] <- feature_colors["Unassigned"]
          feature_colors["Unassigned"] <- light_gray_value
        }
      }
      
      #==============================================================================#
      message("226 # Generate Plots")
      #==============================================================================#
      
      #..............................................................................#
      message("229 # generate plots")
      #..............................................................................#
      
      margin_size        <- 10
      annotate_text_size <- 6  
      strip_text_size    <- 20 
      legend_text_size   <- 20
      
      common_theme <- theme(
        strip.text   = element_text(face = "bold", size = strip_text_size),
        strip.text.y = element_text(face = "bold", size = strip_text_size),
        plot.margin  = margin(margin_size, margin_size, margin_size, margin_size),
        plot.background  = element_rect(fill = palette_color["general_background"], color = palette_color["general_background"]),
        panel.background = element_rect(fill = palette_color["general_background"], color = palette_color["general_background"]),
        legend.background = element_rect(fill = palette_color["general_background"]),
        legend.text = element_markdown(size = legend_text_size),
        axis.title  = element_blank(),
        axis.text   = element_blank(),
        axis.ticks  = element_blank(),
        panel.grid  = element_blank(),
        plot.caption = element_text(hjust = 0.5)
      )
      
      gPlot <- function(p) {
        # Example: p <- gPlot(p)
        p <- p +
          common_theme +
          scale_fill_manual(
            name = "",
            values = feature_colors,
            guide = guide_legend(reverse = TRUE)
          )
        print(p)
        p
      }
      
      #==============================================================================#
      message("262 # legend plot")
      #==============================================================================#
      
      #..............................................................................#
      message("265 # legend plot")
      #..............................................................................#
      
      plot_legend <- ggplot(data = df_plot, aes(x = !!sym(y_facet_group), y = group_percentage_abun, fill = feature_label)) +
        geom_col(position = position_stack(reverse = TRUE)) +
        scale_fill_manual(name = "", values = feature_colors) +
        common_theme +
        theme(
          legend.position = "right",
          plot.margin = margin(0, 0, 0, 0)
        ) +
        guides(fill = guide_legend(ncol = 1, reverse = TRUE))
      
      #==============================================================================#
      message("279 # plots by facet (loop)")
      #==============================================================================#
      
      #..............................................................................#
      message("282 # begin plot loop by facets")
      #..............................................................................#
      
      plot_list <- list()
      ratio_vector <- numeric()
      
      for (facet in y_facet_group_order) {
        
        df_plot1 <- df_plot %>%
          filter(.data[[y_facet_group]] %in% facet)
        
        unique_y_labels <- facet
        outline_color_y <- "black"
        
        s <- paste0(
          'backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',
          paste0(unique_y_labels, collapse = '"]),element_rect(color=outline_color_y,fill = palette_color["'),
          '"]))'
        )
        print(s)
        eval(parse(t = s))
        
        if (use_custom_labels == "yes") {
          custom_y_labels <- c("slope","1500TF","x3FC","slope","1500TF","x3FC","24cycles",rep("Soil",3),rep("Zymo",3),"NTC")
          unique_y_labels <- custom_y_labels
          
          s <- paste0(
            'backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',
            paste0(unique_y_labels, collapse = '"]),element_rect(color=outline_color_y,fill = palette_color["'),
            '"]))'
          )
          print(s)
          eval(parse(t = s))
        }
        
        y_strip <- as.data.frame(unique_y_labels) %>%
          mutate(text_color = "white") %>%
          mutate(face = "bold") %>%
          mutate(text_size = strip_text_size) %>%
          ungroup()
        
        df_plot2 <- df_plot %>%
          select(any_of(y_facet_group), group_percentage_abun, feature_label, plot_order) %>%
          distinct() %>%
          mutate(!!sym(y_facet_group) := factor(!!sym(y_facet_group), levels = y_facet_group_order)) %>%
          arrange(across(any_of(c(y_facet_group, x_facet_group, y_axis_group)))) %>%
          filter(.data[[y_facet_group]] %in% facet) %>%
          arrange(desc(group_percentage_abun)) %>%
          mutate(plot_order = factor(plot_order, levels = unique(plot_order))) %>%
          ungroup()
        
        n <- nrow(df_plot2)
        
        plot_title <- paste0(taxa_levs, "_taxonomy_bar_", normalization_name, "_", data_class_name, "_", data_set, run_suffix)
        
        p <- ggplot(data = df_plot2, aes(x = !!sym(y_facet_group), y = group_percentage_abun, fill = feature_label)) +
          geom_col(position = position_stack(reverse = TRUE)) +
          facet_grid2(
            formula(paste(y_facet_group, "~", ".")),
            strip = strip_themed(
              background_y = backgrounds_y,
              text_y = elem_list_text(
                face  = y_strip$face,
                size  = y_strip$text_size,
                color = y_strip$text_color
              )
            ),
            scale = "free",
            space = "free"
          ) +
          scale_fill_manual(values = feature_colors) +
          labs(title = "", fill = "", x = "", y = "") +
          common_theme +
          theme(legend.position = "none")
        
        print(p)
        
        plot_list[[facet]] <- p
        ratio_vector[facet] <- n
        
      } # end loop over facet
      
      type_plots <- wrap_plots(rev(plot_list), ncol = 1, heights = rev(ratio_vector))
      print(type_plots)
      
      #==============================================================================#
      message("338 # assemble multi-panel + legend")
      #==============================================================================#
      
      left_plot    <- plot_list[["Soil"]]
      top_right    <- plot_list[["Skin"]]
      bottom_right <- plot_list[["Feces"]]
      
      right_plots <- top_right / bottom_right +
        plot_layout(heights = ratio_vector[3:2])
      
      #==============================================================================#
      message("349 # legend extraction")
      #==============================================================================#
      
      legend_grob <- get_legend(plot_legend)
      
      legend_plot <- ggplot() +
        theme_void() +
        theme(
          plot.background  = element_rect(fill = palette_color["general_background"], color = NA),
          panel.background = element_rect(fill = palette_color["general_background"], color = NA),
          plot.margin = margin(0, 0, 0, 0)
        ) +
        annotation_custom(legend_grob)
      
      left_plot    <- left_plot + theme(legend.position = "none")
      top_right    <- top_right + theme(legend.position = "none")
      bottom_right <- bottom_right + theme(legend.position = "none")
      
      right_plots <- top_right / bottom_right +
        plot_layout(heights = ratio_vector[3:2])
      
      final_layout <- (right_plots | left_plot | legend_plot) +
        plot_layout(widths = c(1, 1, 1)) +
        plot_annotation(
          caption = paste("Showing", taxa_plural, "with mean relative abundance >", min_rel_abun_percent, "%"),
          theme = theme(
            plot.caption = element_text(hjust = 0.5),
            plot.background = element_rect(
              fill = palette_color["general_background"],
              color = "gray30",
              linewidth = 2
            )
          )
        )
      
      print(final_layout)
      
      #==============================================================================#
      message("383 # combine + save")
      #==============================================================================#
      
      #..............................................................................#
      message("386 # combine + save")
      #..............................................................................#
      
      plot_width  <- plot_height <- 12
      qSave(figure_title, plot = final_layout, ext = ext_list)
      
    } # end loop over r
  } # end loop over lp
} # end loop over plot_set

###############################################################################
message("397 # END RUNS")
###############################################################################

###############################################################################
message(paste("401 # FINISHED", script_title))
###############################################################################
