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

script_title <- "pipeline"
script_description <- "
Creates Figure 5 workflow outputs from the rarefied ASV matrix.
Builds read-retention summaries, compositional metrics, and combined figure outputs using Welch tests.
"
variant_suffix <- ""

tidy_test_result <- function(test_obj) {
  data.frame(
    p.value = test_obj$p.value,
    statistic = unname(test_obj$statistic)[1],
    stringsAsFactors = FALSE
  )
}

set_output()

plot_set_order      <- c("normalization")
make_new_data_files <- "yes"
testing             <- "yes"
testing_r           <- 3
starting_r          <- 3
ending_r            <- 3

use_custom_labels <- "no"

sample_type_order <- c("cumulative","Soil", "Feces", "Skin")
###############################################################################
message("28 # END SETUP")
###############################################################################

###############################################################################
message("32 # RUNS")
###############################################################################

for (plot_set in plot_set_order) {
  qPrint(plot_set)
  
  ###############################################################################
  message("30 # Plot-set outputs")
  ###############################################################################
  
  folder <- plot_set
  
  if (plot_set == "standard") {
    normalization_order <- c("fixed cycles")
    data_class_order    <- c("doubleton")
    normalization_name  <- paste(normalization_order, collapse = "_")
    data_class_name     <- paste(data_class_order, collapse = "_")
    figure_title        <- "Figure X"
  }
  
  if (plot_set == "normalization") {
    normalization_order <- c("fixed cycles", "targeted fluorescence")
    data_class_order    <- c("doubleton")
    normalization_name  <- paste(normalization_order, collapse = "_")
    data_class_name     <- paste(data_class_order, collapse = "_")
    figure_title        <- paste0("Figure 5 Chimera Formation and Data Loss from Over Cycling", variant_suffix)
    summary_file_name   <- paste0("Figure_5_summary", variant_suffix, ".xlsx")
  }
  
  if (plot_set == "data_class") {
    normalization_order <- c("targeted fluorescence")
    data_class_order    <- c("singleton", "doubleton")
    normalization_name  <- paste(normalization_order, collapse = "_")
    data_class_name     <- paste(data_class_order, collapse = "_")
    figure_title        <- "Figure X"
  }
  
  ###############################################################################
  message("56 # Parameter sets")
  ###############################################################################
  
  parameter_sets <- list(
    set1 = list(
      filter1_group  = "data_class",
      x_axis_group   = "normalization",
      x_facet_group  = "temperature",
      plot_group     = "sample_type"
    ),
    set3 = list(
      filter1_group  = "normalization",
      x_axis_group   = "data_class",
      x_facet_group  = "temperature",
      plot_group     = "sample_type"
    )
  )
  
  ###############################################################################
  message("72 # Metrics setup")
  ###############################################################################
  
  message("74 # metrics set up")
  
  metric_order1 <- c(
    "percent_Muri", "percent_Lachno", "percent_Propi", "percent_Archaea",
    "percent_good", "ideal_score", "richness", "evenness", "shannon", "n"
  )
  
  meta_SD_columns <- c(
    "c_fastq_reads","c_trimmed_reads","c_filtered_reads","c_denoised_reads",
    "c_chimeric_reads","c_pentaton_reads","c_contaminant_reads","c_usable_reads",
    "p_contaminant_reads","p_pentaton_reads"
  )
  
  meta_SD <- meta %>%
    mutate(
      c_fastq_reads       = fastq_reads / fastq_reads * 100,
      c_trimmed_reads     = input / fastq_reads * 100,
      c_filtered_reads    = filtered / fastq_reads * 100,
      c_denoised_reads    = denoised / fastq_reads * 100,
      c_chimeric_reads    = non.chimeric / fastq_reads * 100,
      c_pentaton_reads    = reads_with_contaminants / fastq_reads * 100,
      c_contaminant_reads = sample_counts / fastq_reads * 100,
      c_usable_reads      = usable_reads / fastq_reads * 100,
      p_pentaton_reads    = (non.chimeric - reads_with_contaminants) / fastq_reads * 100,
      p_contaminant_reads = (reads_with_contaminants - usable_reads) / fastq_reads * 100
    ) %>%
    select(sample_name, any_of(meta_SD_columns))
  
  meta=meta %>% 
    select(-any_of(meta_SD_columns))
  meta_names=names(meta)
  
  metric_reads <- rev(c(
    "fastq_reads","trimmed_reads","filtered_reads","denoised_reads",
    "chimeric_reads","contaminant_table_reads","contaminant_reads","usable_reads"
  ))
  
  metric_p_reads <- rev(c(
    "p_fastq_reads","p_trimmed_reads","p_filtered_reads","p_denoised_reads",
    "p_chimeric_reads","p_pentaton_reads","p_contaminant_reads","p_usable_reads"
  ))
  
  p_metric_numerator_map <- c(
    "p_fastq_reads"       = "fastq_reads",
    "p_trimmed_reads"     = "trimmed_reads",
    "p_filtered_reads"    = "filtered_reads",
    "p_denoised_reads"    = "denoised_reads",
    "p_chimeric_reads"    = "chimeric_reads",
    "p_pentaton_reads"    = "contaminant_table_reads",
    "p_contaminant_reads" = "contaminant_reads",
    "p_usable_reads"      = "usable_reads"
  )
  
  metric_cumulative <- c(
    "c_fastq_reads","c_trimmed_reads","c_filtered_reads","c_denoised_reads",
    "c_chimeric_reads","c_pentaton_reads","c_contaminant_reads","c_usable_reads"
  )
  
  metric_order <- c(metric_order1, metric_reads, metric_p_reads, rev(metric_cumulative))
  
  selected_metric_groups <- c("metric_cumulative")  # c("metric_reads","metric_p_reads","metric_cumulative")
  
  test <- metric_order
  names(test) <- metric_order
  nvPrint(test)
  
  #==============================================================================#
  message("122 # Main Loop: Process Datasets")
  #==============================================================================#
  
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
      starting_r <- testing_r
      ending_r   <- testing_r
      message(paste0("\trunning r= ", testing_r))
    }
    
    for (r in starting_r:ending_r) {
      
      #..............................................................................#
      message("148 # dataset details")
      #..............................................................................#
      
      taxa_levs   <- matrix_names[r, "taxa_levs"]
      data_set    <- matrix_names[r, "data_sets"]
      taxa_plural <- matrix_names[r, "taxa_plural"]
      
      if (plot_set != "none") {
        data_set       <- plot_set
        data_set_order <- plot_set_order
      } else {
        data_set_order <- unique(matrix_names[, "data_sets"])
      }
      
      # Assign analysis parameters dynamically
      list2env(parameter_sets[[paste0("set", match(data_set, data_set_order))]], envir = .GlobalEnv)
      
      groups <- c("filter1_group", "x_axis_group", "x_facet_group", "plot_group")
      for (var in groups) {
        group_value <- get(var)
        assign(paste0(var, "_order"), get(paste0(group_value, "_order")))
      }
      
      grouping_columns <- c(x_axis_group, x_facet_group, plot_group)
      test_across      <- c(x_facet_group, plot_group)
      group_within     <- c(x_axis_group, plot_group)
      
      run_suffix <- ""
      if (exists("Run_Group")) {
        run_suffix <- paste0("_", lp)
      }
      
      filter1_names <- paste(filter1_group_order, collapse = "_")
      x_axis_names  <- paste(x_axis_group_order,  collapse = "_")
      x_facet_names <- paste(x_facet_group_order, collapse = "_")
      
      p_title <- paste0(
        taxa_levs, "_", script_title, custom_name, "_",
        filter1_names, "_", x_axis_names, "_", x_facet_names, "_",
        data_set, run_suffix
      )
      
      qPrint(p_title)
      
      #==============================================================================#
      message("195 # Load and Filter Data")
      #==============================================================================#
      
      if (make_new_data_files == "yes") {
        
        message("199 # load and filter data")
        
        matrix_df1 <- read.csv(matrix_names[r, "file_path"], check.names = FALSE, row.names = 1) %>%
          t() %>%
          as.data.frame() %>%
          xPlode_sample_name() %>%
          filter(.data[[filter1_group]]  %in% filter1_group_order) %>%
          filter(.data[[plot_group]]    %in% plot_group_order) %>%
          filter(.data[[x_axis_group]]  %in% x_axis_group_order) %>%
          filter(.data[[x_facet_group]] %in% x_facet_group_order)
        
        if (exists("Run_Group")) {
          matrix_df1 <- matrix_df1 %>% filter(.data[[Loop_Group]] %in% lp)
          run_suffix <- paste0("_", lp)
        }
        
        p_title <- paste0(
          taxa_levs, "_", script_title, custom_name, "_",
          normalization_name, "_", data_class_name, "_",
          data_set, run_suffix
        )
        
        qPrint(p_title)
        
        matrix_df <- matrix_df1 %>%
          imPlode_sample_name() %>%
          mutate_all(as.numeric) %>%
          select(which(colSums(.) > 0))
        
        feature_names <- names(matrix_df)
        
        ###############################################################################
        message("228 # Taxa percentages and Zymo score")
        ###############################################################################
        
        message("230 # get percent Archaea and other taxa of interest")
        
        df2 <- matrix_df %>%
          xPlode_sample_name() %>%
          pivot_longer(col = all_of(feature_names), names_to = "feature", values_to = "counts") %>%
          group_by(feature, across(any_of(test_across))) %>%
          mutate(feature_grouped_counts = sum(counts)) %>%
          filter(feature_grouped_counts > 0) %>%
          mutate(taxa = feature) %>%
          group_by(sample_name) %>%
          mutate(sample_counts = sum(counts)) %>%
          ungroup()
        
        if (matrix_names[r, "taxa_levs"] == "ASV") {
          df2 <- df2 %>%
            rename(ASV = taxa) %>%
            left_join(ASV_taxa %>% select(ASV, taxa, Family)) %>%
            ungroup()
        }
        
        Archaea_df <- extract_taxa_percentage(df2, name = "Archaea", target = "archaea")
        Propi_df   <- extract_taxa_percentage(df2, name = "Propi",   target = "g__Cutibacterium")
        Lachno_df  <- extract_taxa_percentage(df2, name = "Lachno",  target = "f__Lachnospiraceae")
        Muri_df    <- extract_taxa_percentage(df2, name = "Muri",    target = "f__Muribaculaceae")
        
        message("254 # for Zymo Testing Ideal Score")
        
        #==============================================================================#
        message("257 # Generate alpha diversity metrics")
        #==============================================================================#
        
        metrics_df <- matrix_df %>%
          as.data.frame() %>%
          xPlode_sample_name() %>%
          rowwise() %>%
          mutate(
            mn           = mean(c_across(any_of(feature_names))),
            sd           = sd(c_across(any_of(feature_names))),
            smp_n        = sum(c_across(any_of(feature_names))),
            smp_shannon  = qShannon(c_across(any_of(feature_names))),
            smp_evenness = qEvenness(c_across(any_of(feature_names))),
            smp_richness = qRichness(c_across(any_of(feature_names)))
          ) %>%
          ungroup() %>%
          select(-any_of(feature_names)) %>%
          full_join(Archaea_df, by = "sample_name") %>%
          full_join(Propi_df,   by = "sample_name") %>%
          full_join(Lachno_df,  by = "sample_name") %>%
          full_join(Muri_df,    by = "sample_name") %>%
          left_join(meta_SD, by = "sample_name") %>%
          rename_with(~ paste0("smp_", .x), .cols = any_of(c(metric_reads, metric_p_reads, metric_cumulative))) %>%
          {
            metrics_df_tmp <- .
            for (p_metric in names(p_metric_numerator_map)) {
              numerator_col <- paste0("smp_", p_metric_numerator_map[[p_metric]])
              percent_col <- paste0("smp_", p_metric)
              metrics_df_tmp[[percent_col]] <- 100 * metrics_df_tmp[[numerator_col]] / metrics_df_tmp[["smp_fastq_reads"]]
            }
            metrics_df_tmp
          } %>%
          group_by(across(any_of(grouping_columns))) %>%
          mutate(across(starts_with("smp_"), ~ mean(., na.rm = TRUE), .names = "{sub('smp_', 'mean_', .col)}")) %>%
          rename_with(~ gsub("smp_", "", .), starts_with("ratio_smp_")) %>%
          ungroup() %>%
          select(where(~ !all(is.na(.))))
        
        write.csv(metrics_df, paste0(output_data, p_title, "_metrics_df.csv"), row.names = FALSE)
      }
      
      if (make_new_data_files != "yes") {
        metrics_df <- read.csv(paste0(output_data, p_title, "_metrics_df.csv"), check.names = FALSE)
      }
      
      #==============================================================================#
      message("309 # Long format and filters")
      #==============================================================================#
      
      result_df <- data_df <- metrics_df %>%
        pivot_longer(cols = starts_with("smp_"), names_to = "metric") %>%
        filter(!(sample_type %in% c("Skin")  & (metric %in% c("smp_percent_Muri",  "smp_percent_Archaea")))) %>%
        filter(!(sample_type %in% c("Soil")  & (metric %in% c("smp_percent_Propi", "smp_percent_Muri")))) %>%
        filter(!(sample_type %in% c("Feces") & (metric %in% c("smp_percent_Propi", "smp_percent_Archaea")))) %>%
        filter(metric != "smp_percent_Lachno") %>%
        group_by(across(any_of(grouping_columns)), metric) %>%
        mutate(mean_value = mean(value, na.rm = TRUE)) %>%
        ungroup()
      
      #==============================================================================#
      message("327 # Wilcoxon test if 2 levels")
      #==============================================================================#
      
      message("329 # T-test if length(x_axis_group)==2")
      
      if (length(unique(data_df[[x_axis_group]])) == 2) {
        
        t_test_results <- data_df %>%
          filter(!is.nan(value)) %>%
          group_by(across(any_of(test_across)), metric) %>%
          summarise(
            t_test_summary = {
              unique_subgroups <- unique(.data[[x_axis_group]])
              
              if (length(unique_subgroups) != 2) {
                list(data.frame(
                  p.value = NA,
                  statistic = NA,
                  fully_separated = NA
                ))
              } else {
                subgroup_counts <- map(unique_subgroups, ~ n_distinct(.data$value[.data[[x_axis_group]] == .x]))
                
                if (any(unlist(subgroup_counts) == 1)) {
                  list(data.frame(
                    p.value = NA,
                    statistic = NA,
                    fully_separated = NA
                  ))
                } else {
                  subgroup_values <- map(unique_subgroups, ~ .data$value[.data[[x_axis_group]] == .x])
                  group1_values <- subgroup_values[[1]]
                  group2_values <- subgroup_values[[2]]
                  fully_separated_flag <- max(group1_values, na.rm = TRUE) < min(group2_values, na.rm = TRUE) |
                    max(group2_values, na.rm = TRUE) < min(group1_values, na.rm = TRUE)
                  
                  list(
                    tidy_test_result(t.test(value ~ .data[[x_axis_group]], var.equal = FALSE)) %>%
                      mutate(fully_separated = fully_separated_flag)
                  )
                }
              }
            },
            .groups = "drop"
          ) %>%
          unnest(t_test_summary)
        
        t_test_df <- t_test_results %>%
          select(any_of(test_across), metric, p.value, fully_separated)
        
        join_by <- setNames(c(test_across, "metric"), c(test_across, "metric"))
        
        result_df <- data_df %>%
          left_join(t_test_df, by = join_by) %>%
          rename(p_value = p.value) %>%
          mutate(metric = factor(metric, levels = rev(paste0("smp_", metric_order))))
      }
      
      #==============================================================================#
      message("386 # Plotting helpers and theme")
      #==============================================================================#
      
      message("388 # looping plots by plot group")
      
      p_value_size    <- 3
      geom_text_size  <- 3
      strip_text_size <- 11
      x_text_size     <- 9
      
      theme_common <- theme_global(base_size = 11) +
        theme(
          plot.tag       = element_text(size = 36, face = "bold"),
          legend.position = "bottom",
          axis.title     = element_blank(),
          axis.text.x    = element_blank(),
          axis.ticks.x   = element_blank(),
          axis.text.y    = element_text(face = "plain")
        )
      
      
      gPlot <- function(p) {
        p <- p +
          scale_color_manual(values = palette_color, labels = palette_label) +
          scale_fill_manual(values  = palette_color, labels = palette_label) +
          theme_common +
          guides(
            color = guide_legend(order = 1),
            fill  = guide_legend(
              order        = 2,
              override.aes = list(shape = 22, size = 6, color = "black", alpha = 1)
            )
          ) +
          labs(color = "", fill = "", title = "", x = "", y = "")
        p
      }
      #==============================================================================#
      message("465 # Plotting loop")
      #==============================================================================#
      
      message("467 # begin plotting loop")
      
      for (smg in selected_metric_groups) {
        
        selected_metric_names <- get(smg)
        
        plot_list <- list()
        plot_export_list <- list()
        list_category <- unique(result_df[[plot_group]])
        
        for (t in list_category) {
          
          qPrint(unique(result_df[[plot_group]]))
          qPrint(t)
          
          df <- result_df %>%
            filter(.data[[plot_group]] == t) %>%
            mutate(metric = factor(metric, levels = rev(paste0("smp_", metric_order))))
          
          #..............................................................................#
          message("487 # strip backgrounds with outlines")
          #..............................................................................#
          
          unique_x_labels <- unique(df[[plot_group]])
          
          outline_color_x <- "black"
          s <- paste0(
            'backgrounds_x <- list(element_rect(color=outline_color_x,fill = palette_color["',
            paste0(unique_x_labels, collapse = '"]),element_rect(color=outline_color_x,fill = palette_color["'),
            '"]))'
          )
          eval(parse(text = s))
          
          unique_y_labels <- c(rep(unique_x_labels, (length(unique(df$metric)) / length(unique_x_labels))))
          
          outline_color_y <- "black"
          s <- paste0(
            'backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',
            paste0(unique_y_labels, collapse = '"]),element_rect(color=outline_color_y,fill = palette_color["'),
            '"]))'
          )
          eval(parse(text = s))
          
          y_strip <- as.data.frame(unique_y_labels) %>%
            mutate(text_color = "white", face = "bold", text_size = strip_text_size) %>%
            ungroup()
          
          x_strip <- as.data.frame(unique_x_labels) %>%
            mutate(text_color = "white", face = "bold", text_size = strip_text_size) %>%
            ungroup()
          
          #..............................................................................#
          message("525 # Plotting")
          #..............................................................................#
          
          df_plot <- df %>%
            mutate(metric = gsub("smp_", "", metric)) %>%
            mutate(metric = factor(metric, levels = rev(metric_order))) %>%
            mutate(mean_value_label = signif(mean_value, 3)) %>%
            mutate(mean_value_label = if_else(mean_value_label > 1000, round(mean_value, 0), mean_value_label)) %>%
            mutate(mean_value_label = paste(mean_value_label, "%")) %>%
            mutate(!!x_axis_group   := factor(.data[[x_axis_group]],   levels = x_axis_group_order)) %>%
            mutate(!!x_facet_group  := factor(.data[[x_facet_group]],  levels = x_facet_group_order)) %>%
            mutate(!!plot_group     := factor(.data[[plot_group]],     levels = plot_group_order)) %>%
            ungroup() %>%
            arrange(metric)
          
          df_cumulative_summary <- df_plot %>%
            filter(metric %in% c(metric_cumulative)) %>%
            select(metric, value, everything()) %>%
            mutate(metric = gsub("c_", "", metric)) %>%
            group_by(metric, across(any_of(grouping_columns))) %>%
            summarize(mean_value = mean(value), .groups = "drop") %>%
            mutate(metric = ifelse(metric == "chimerireads", "chimeric_reads", metric)) %>%
            mutate(cumulative = "cumulative") %>%
            mutate(metric = factor(metric, levels = unique(metric))) %>% 
            mutate(cumulative = factor(cumulative, levels = sample_type_order))
          
          df_fastq <- df_plot %>%
            filter(metric %in% c("fastq_reads")) %>%
            group_by(across(any_of(x_facet_group))) %>%
            mutate(max_reads = max(value)) %>%
            ungroup() %>%
            mutate(relative_fastq_reads = 100 * value / max_reads) %>%
            select(metric, value, max_reads, relative_fastq_reads, any_of(grouping_columns), everything()) %>%
            ungroup()
          
          df_fastq_summary <- df_fastq %>%
            group_by(metric, across(any_of(grouping_columns))) %>%
            summarize(mean_value = mean(value), .groups = "drop") %>%
            mutate(mean_value_label = paste(round(mean_value), "reads"))
          
          df_percentage <- df_plot %>%
            filter(metric %in% metric_p_reads) %>%
            select(metric, value, everything()) %>%
            mutate(metric = gsub("p_", "", metric)) %>%
            ungroup() %>%
            mutate(metric = factor(metric, levels = unique(metric)))
          
          percentage_sample_detail <- metrics_df %>%
            select(
              sample_name,
              any_of(c(plot_group, x_facet_group, x_axis_group)),
              smp_fastq_reads,
              any_of(paste0("smp_", metric_p_reads)),
              any_of(paste0("smp_", metric_reads))
            ) %>%
            pivot_longer(
              cols = any_of(paste0("smp_", metric_p_reads)),
              names_to = "percent_metric",
              values_to = "point_value_percent"
            ) %>%
            mutate(
              metric = gsub("^smp_p_", "", percent_metric),
              reads_metric_col = case_when(
                metric == "fastq_reads" ~ "smp_fastq_reads",
                metric == "trimmed_reads" ~ "smp_trimmed_reads",
                metric == "filtered_reads" ~ "smp_filtered_reads",
                metric == "denoised_reads" ~ "smp_denoised_reads",
                metric == "chimeric_reads" ~ "smp_chimeric_reads",
                metric == "pentaton_reads" ~ "smp_contaminant_table_reads",
                metric == "contaminant_reads" ~ "smp_contaminant_reads",
                metric == "usable_reads" ~ "smp_usable_reads",
                TRUE ~ NA_character_
              )
            ) %>%
            rowwise() %>%
            mutate(point_value_reads = get(reads_metric_col)) %>%
            ungroup() %>%
            transmute(
              sample_name = sample_name,
              !!plot_group := .data[[plot_group]],
              !!x_facet_group := .data[[x_facet_group]],
              !!x_axis_group := .data[[x_axis_group]],
              metric = metric,
              point_value_percent = point_value_percent,
              point_value_reads = point_value_reads,
              sample_fastq_reads = smp_fastq_reads
            )
          
          fastq_sample_detail <- df_fastq %>%
            transmute(
              sample_name = sample_name,
              !!plot_group := .data[[plot_group]],
              !!x_facet_group := .data[[x_facet_group]],
              !!x_axis_group := .data[[x_axis_group]],
              metric = "fastq_reads",
              point_value_percent = relative_fastq_reads,
              point_value_reads = value,
              sample_fastq_reads = value
            )
          
          sample_point_base <- bind_rows(
            percentage_sample_detail,
            fastq_sample_detail
          )
          
          vector_summary_by_normalization <- sample_point_base %>%
            mutate(
              p_metric = factor(paste0("p_", metric), levels = metric_p_reads),
              metric_label = palette_label[as.character(p_metric)]
            ) %>%
            group_by(normalization = .data[[x_axis_group]], p_metric, metric_label) %>%
            summarize(
              n = dplyr::n(),
              mean_percent = mean(point_value_percent, na.rm = TRUE),
              sd_percent = sd(point_value_percent, na.rm = TRUE),
              min_percent = min(point_value_percent, na.rm = TRUE),
              max_percent = max(point_value_percent, na.rm = TRUE),
              range_percent = paste0(signif(min_percent, 4), " to ", signif(max_percent, 4)),
              .groups = "drop"
            ) %>%
            arrange(normalization, p_metric)

          vector_summary_by_normalization_sample_type <- sample_point_base %>%
            mutate(
              p_metric = factor(paste0("p_", metric), levels = metric_p_reads),
              metric_label = palette_label[as.character(p_metric)]
            ) %>%
            group_by(
              normalization = .data[[x_axis_group]],
              sample_type = .data[[plot_group]],
              p_metric,
              metric_label
            ) %>%
            summarize(
              n = dplyr::n(),
              mean_percent = mean(point_value_percent, na.rm = TRUE),
              sd_percent = sd(point_value_percent, na.rm = TRUE),
              min_percent = min(point_value_percent, na.rm = TRUE),
              max_percent = max(point_value_percent, na.rm = TRUE),
              range_percent = paste0(signif(min_percent, 4), " to ", signif(max_percent, 4)),
              .groups = "drop"
            ) %>%
            arrange(normalization, sample_type, p_metric)

          p_vector_values <- sample_point_base %>%
            mutate(
              p_metric = factor(paste0("p_", metric), levels = metric_p_reads),
              metric_label = palette_label[as.character(p_metric)]
            ) %>%
            select(
              normalization = all_of(x_axis_group),
              sample_type = all_of(plot_group),
              temperature = all_of(x_facet_group),
              sample_name,
              p_metric,
              metric_label,
              point_value_percent,
              point_value_reads,
              sample_fastq_reads
            ) %>%
            arrange(normalization, sample_type, temperature, p_metric, sample_name)

          p_vector_values_wide <- sample_point_base %>%
            mutate(
              p_metric = paste0("p_", metric)
            ) %>%
            select(
              normalization = all_of(x_axis_group),
              sample_type = all_of(plot_group),
              temperature = all_of(x_facet_group),
              sample_name,
              p_metric,
              point_value_percent
            ) %>%
            pivot_wider(
              names_from = p_metric,
              values_from = point_value_percent,
              values_fn = ~ .x[[1]]
            ) %>%
            arrange(normalization, sample_type, temperature, sample_name)

          p_vector_test_base <- sample_point_base %>%
            mutate(
              sample_type = factor(.data[[plot_group]], levels = plot_group_order),
              normalization = factor(.data[[x_axis_group]], levels = x_axis_group_order),
              p_metric = factor(paste0("p_", metric), levels = metric_p_reads),
              metric_label = palette_label[as.character(p_metric)]
            ) %>%
            select(sample_type, normalization, p_metric, metric_label, point_value_percent)

          run_percent_test <- function(df) {
            normalize_decimal_numeric <- function(x) {
              if (is.na(x)) return(NA_real_)
              as.numeric(format(unname(x)[1], scientific = FALSE, trim = TRUE, digits = 15))
            }
            norm_levels <- as.character(unique(df$normalization))
            if (length(norm_levels) != 2) {
              return(tibble(
                group_1 = NA_character_,
                group_2 = NA_character_,
                n_group_1 = NA_integer_,
                n_group_2 = NA_integer_,
                mean_group_1 = NA_real_,
                mean_group_2 = NA_real_,
                statistic = NA_real_,
                p_value = NA_real_
              ))
            }

            group_1 <- norm_levels[1]
            group_2 <- norm_levels[2]
            values_1 <- df %>% filter(normalization == group_1) %>% pull(point_value_percent)
            values_2 <- df %>% filter(normalization == group_2) %>% pull(point_value_percent)

            test_tbl <- tryCatch(
              {
                tidy_test_result(t.test(point_value_percent ~ normalization, data = df, var.equal = FALSE))
              },
              error = function(e) tibble(p.value = NA_real_, statistic = NA_real_)
            )

            tibble(
              group_1 = group_1,
              group_2 = group_2,
              n_group_1 = length(values_1),
              n_group_2 = length(values_2),
              mean_group_1 = mean(values_1, na.rm = TRUE),
              mean_group_2 = mean(values_2, na.rm = TRUE),
              statistic = normalize_decimal_numeric(test_tbl$statistic[[1]]),
              p_value = normalize_decimal_numeric(test_tbl$p.value[[1]])
            )
          }

          p_vector_tests_by_type <- p_vector_test_base %>%
            group_by(sample_type, p_metric, metric_label) %>%
            group_modify(~ run_percent_test(.x)) %>%
            ungroup() %>%
            arrange(sample_type, p_metric) %>%
            rename(
              welch_statistic = statistic,
              welch_p_value = p_value
            )

          round_p_value <- function(x) {
            ifelse(
              is.na(x),
              NA_real_,
              round(x, 8)
            )
          }

          round_p_value_columns <- function(df) {
            p_cols <- grepl("(^p$)|(^p_value$)|(_p_value$)|(\\.p\\.value$)|(^p\\.value$)|(^adj\\.p\\.value$)|(^anova\\.p\\.value$)",
              names(df), ignore.case = TRUE)
            if (any(p_cols)) {
              df[p_cols] <- lapply(df[p_cols], function(col) {
                if (is.numeric(col)) round_p_value(col) else col
              })
            }
            df
          }

          p_vector_tests_by_norm <- p_vector_test_base %>%
            group_by(p_metric, metric_label) %>%
            group_modify(~ run_percent_test(.x)) %>%
            ungroup() %>%
            arrange(p_metric) %>%
            rename(
              welch_statistic = statistic,
              welch_p_value = p_value
            ) %>%
            mutate(
              welch_p_value = round_p_value(welch_p_value)
            ) %>%
            as_tibble() %>%
            mutate(across(everything(), ~ if (is.factor(.x)) as.character(.x) else .x)) %>%
            as.data.frame(stringsAsFactors = FALSE)
          
          sample_point_export <- build_boxplot_export(
            df = sample_point_base,
            x_axis_group = x_axis_group,
            x_facet_group = x_facet_group,
            y_facet_group = "metric",
            extra_groups = plot_group,
            value_col = "point_value_percent",
            value_name = "point_value_percent"
          ) %>%
            left_join(
              build_boxplot_export(
                df = sample_point_base,
                x_axis_group = x_axis_group,
                x_facet_group = x_facet_group,
                y_facet_group = "metric",
                extra_groups = plot_group,
                value_col = "point_value_reads",
                value_name = "point_value_reads"
              ) %>%
                select(-panel, -boxplot),
              by = c("x_faceting", "y_faceting", "x_axis", plot_group, "sample_name")
            ) %>%
            left_join(
              build_boxplot_export(
                df = sample_point_base,
                x_axis_group = x_axis_group,
                x_facet_group = x_facet_group,
                y_facet_group = "metric",
                extra_groups = plot_group,
                value_col = "sample_fastq_reads",
                value_name = "sample_fastq_reads"
              ) %>%
                select(-panel, -boxplot),
              by = c("x_faceting", "y_faceting", "x_axis", plot_group, "sample_name")
            )
          
          boxplot_export <- bind_rows(
            build_boxplot_export(
              df = df_percentage %>%
                select(sample_name, any_of(c(plot_group, x_facet_group, "metric", x_axis_group)), value),
              x_axis_group = x_axis_group,
              x_facet_group = x_facet_group,
              y_facet_group = "metric",
              extra_groups = plot_group,
              value_col = "value",
              value_name = "boxplot_value"
            ),
            build_boxplot_export(
              df = df_fastq %>%
                select(sample_name, any_of(c(plot_group, x_facet_group, "metric", x_axis_group)), relative_fastq_reads) %>%
                rename(value = relative_fastq_reads),
              x_axis_group = x_axis_group,
              x_facet_group = x_facet_group,
              y_facet_group = "metric",
              extra_groups = plot_group,
              value_col = "value",
              value_name = "boxplot_value"
            )
          )
          
          cumulative_export_sample <- build_boxplot_export(
            df = df_plot %>%
              filter(metric %in% c(metric_cumulative)) %>%
              mutate(metric = gsub("c_", "", metric)) %>%
              mutate(metric = ifelse(metric == "chimerireads", "chimeric_reads", metric)) %>%
              select(sample_name, any_of(c(plot_group, x_facet_group, "metric", x_axis_group)), value) %>%
              distinct(),
            x_axis_group = x_axis_group,
            x_facet_group = x_facet_group,
            y_facet_group = "metric",
            extra_groups = plot_group,
            value_col = "value",
            value_name = "cumulative_value"
          ) %>%
            select(-panel, -boxplot)
          
          cumulative_export_group <- build_group_value_export(
            df = df_cumulative_summary %>%
              select(any_of(c(plot_group, x_facet_group, "metric", x_axis_group)), mean_value),
            x_axis_group = x_axis_group,
            x_facet_group = x_facet_group,
            y_facet_group = "metric",
            extra_groups = plot_group,
            value_col = "mean_value",
            value_name = "col_value"
          ) %>%
            select(-panel, -boxplot)
          
          plot_export_list[[t]] <- boxplot_export %>%
            left_join(
              sample_point_export,
              by = c("panel", "boxplot", "x_faceting", "y_faceting", "x_axis", plot_group, "sample_name")
            ) %>%
            left_join(
              cumulative_export_sample,
              by = c("x_faceting", "y_faceting", "x_axis", plot_group, "sample_name")
            ) %>%
            left_join(
              cumulative_export_group,
              by = c("x_faceting", "y_faceting", "x_axis", plot_group)
            )
          
          fill_breaks <- c("cumulative", t)
          
          linetype_legend_df <- tibble(
            sig_group = factor(
              c("sig", "fully_separated", "ns"),
              levels = c("sig", "fully_separated", "ns")
            ),
            x_dummy = factor(x_axis_group_order[1], levels = x_axis_group_order),
            xend_dummy = factor(x_axis_group_order[1], levels = x_axis_group_order),
            y_dummy = 0,
            yend_dummy = 0
          )
          
          
          p <- ggplot(df_percentage, aes(x = get(x_axis_group), y = value, color = get(x_axis_group))) +
            
            #................ adding cumulative/remaining reads to legend .................# 
            
            geom_point(data = df_cumulative_summary,
                       aes(x = get(x_axis_group), y = mean_value,fill = cumulative),
                       alpha = 0,
                       show.legend = TRUE
            ) +
            
            #........................ adding sample_type to legend ........................# 
            
            geom_point(
              aes(fill = get(plot_group)),
              alpha = 0,
              show.legend = TRUE
            ) +
            
            geom_col(
              data = df_cumulative_summary,
              aes(x = get(x_axis_group), y = mean_value,),
              fill = palette_color['cumulative'],color = palette_color['cumulative']
              #,size=7
            ) +
            geom_boxplot(
              data        = df_percentage %>% filter(metric != "fastq_reads"),
              width       = .6,
              outlier.size = .5
            ) +
            geom_boxplot(
              data        = df_fastq,
              aes(y       = relative_fastq_reads),
              width       = .6,
              outlier.size = .5
            ) +
            geom_text(
              data  = df_percentage %>% filter(metric != "fastq_reads"),
              aes(label = paste0("  ", mean_value_label)),
              y     = 0,
              vjust = .5,
              hjust = 0,
              color = "black",
              show.legend = FALSE,
              size  = geom_text_size,
              angle = 90
            ) +
            geom_text(
              data  = df_fastq_summary,
              aes(label = paste0("  ", mean_value_label)),
              y     = 0,
              vjust = .5,
              hjust = 0,
              color = "black",
              show.legend = FALSE,
              size  = geom_text_size,
              angle = 90
            ) +
            geom_text(
              data  = df_cumulative_summary,
              aes(y = 90, label = paste(round(mean_value, 0), "%")),
              vjust = .5,
              hjust = .5,
              size  = geom_text_size,
              fontface = "bold",
              color    = "gray30"
            ) +
            facet_grid2(
              get(x_facet_group) ~ metric,
              strip = strip_themed(
                background_x = backgrounds_x,
                background_y = backgrounds_y,
                text_y       = elem_list_text(face = y_strip$face, color = y_strip$text_color),
                text_x       = elem_list_text(face = x_strip$face, color = x_strip$text_color)
              ),
              labeller = labeller(.cols = palette_label, .rows = palette_label),
              scale    = "free"
            ) +
            scale_y_continuous(limits = c(0, NA)) +
            scale_x_discrete(labels = function(x) {
              color <- palette_color[x]
              label <- palette_label[x]
              paste0("<span style='color:", color, "'>", label, "</span>")
            }) +
            theme(
              axis.text.x = element_markdown(
                angle = 90,
                hjust = 1,
                vjust = 0.5,
                face  = "bold",
                size  = x_text_size
              )
            ) +
            labs(title = paste(taxa_levs, t, smg), x = "", y = "")
          p
          gPlot(p)
          q=gPlot(p)
          #..............................................................................#
          message("605 # Add p-values if 2 levels")
          #..............................................................................#
          
          if (length(unique(data_df[[x_axis_group]])) == 2) {
            
            df_plot2 <- df_percentage %>%
              mutate(p_label      = paste0("p = ", format(p_value, scientific = TRUE, digits = 3))) %>%
              mutate(significance = ifelse(p_value < .05, "sig", "ns")) %>%
              
              mutate(
                comparison_characteristic = case_when(
                  p_value < .05 ~ "sig",
                  coalesce(fully_separated, FALSE) ~ "fully_separated",
                  TRUE ~ "ns"
                )) %>%
              
              mutate(p_label      = ifelse(metric == "fastq_reads", "", p_label)) %>%
              group_by(metric, across(any_of(plot_group))) %>%
              mutate(y_p_label = 70) %>%
              ungroup() %>%
              mutate(!!x_axis_group  := factor(.data[[x_axis_group]],  levels = x_axis_group_order)) %>%
              mutate(!!x_facet_group := factor(.data[[x_facet_group]], levels = x_facet_group_order)) %>%
              mutate(!!plot_group    := factor(.data[[plot_group]],    levels = plot_group_order)) %>%
              filter(metric!='fastq_reads') %>% 
              ungroup()
            
            p <- q +
              geom_label(
                data        = df_plot2 %>% filter(significance == "sig"),
                aes(x = 1.5, y = y_p_label, label = p_label),
                color       = "red",
                inherit.aes = FALSE,
                fill        = scales::alpha("gray95", 0.1),
                show.legend = FALSE,
                hjust       = 0.5,
                vjust       = 1,
                size        = p_value_size,
                fontface    = "bold"
              ) +
              geom_label(
                data        = df_plot2 %>% filter(significance == "ns" | is.na(significance)),
                aes(x = 1.5, y = y_p_label, label = p_label),
                color       = "gray50",
                inherit.aes = FALSE,
                fill        = scales::alpha("gray90", 0.05),
                show.legend = FALSE,
                hjust       = 0.5,
                vjust       = 1,
                size        = p_value_size,
                fontface    = "bold"
              )
            
          }
          
          plot_list[[t]] <- gPlot(p)
        } # end loop over t (list_category)
        
        #==============================================================================#
        message("653 # Combine and save plots")
        #==============================================================================#
        
        message("655 # save plots")
        qPrint(smg)
        
        combined_plot <- wrap_plots(plot_list, ncol = 1) +
          plot_annotation(tag_levels = "A", theme = theme_plot)
        
        plot_export_df <- bind_rows(plot_export_list) %>%
          arrange(panel, x_axis, sample_name)
        
        summary_sheets <- list(
          plot_data = plot_export_df,
          vector_by_norm = vector_summary_by_normalization,
          vector_by_norm_type = vector_summary_by_normalization_sample_type,
          p_vector_values = p_vector_values,
          p_vector_values_wide = p_vector_values_wide,
          p_vector_tests_by_norm = p_vector_tests_by_norm,
          p_vector_tests_by_type = p_vector_tests_by_type
        )
        normalize_list_column <- function(col) {
          if (!is.list(col)) return(col)
          all_numeric_like <- all(vapply(col, function(x) {
            length(x) == 0 || all(is.na(x)) || is.numeric(x)
          }, logical(1)))
          if (all_numeric_like) {
            return(vapply(col, function(x) {
              if (length(x) == 0 || all(is.na(x))) {
                NA_real_
              } else {
                as.numeric(x[[1]])
              }
            }, numeric(1)))
          }
          vapply(col, function(x) {
            if (length(x) == 0 || all(is.na(x))) {
              NA_character_
            } else {
              paste(as.character(x), collapse = ", ")
            }
          }, character(1))
        }
        summary_sheets <- lapply(summary_sheets, function(x) {
          x %>%
            as_tibble() %>%
            mutate(across(where(is.list), normalize_list_column)) %>%
            mutate(across(everything(), ~ if (is.factor(.x)) as.character(.x) else .x)) %>%
            round_p_value_columns()
        })
        write_formatted_workbook(
          summary_sheets,
          file.path(output_data, summary_file_name)
        )
        
        plot_width  <- 18
        plot_height <- 18
        
        qSave(figure_title, plot = combined_plot, ext = ext_list)
      } # end loop over smg (selected_metric_groups)
    } # end loop over r (starting_r:ending_r)
  } # end loop over lp (Loop_Group_order)
} # end loop over plot_set (plot_set_order)

###############################################################################
message("673 # END RUNS")
###############################################################################

###############################################################################
message(paste("677 # FINISHED", script_title))
###############################################################################
