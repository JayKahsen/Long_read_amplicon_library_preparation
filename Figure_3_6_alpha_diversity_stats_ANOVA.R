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

script_title <- "alpha_diversity"
script_description <- "
Creates Figure 3 and Figure 6 alpha-diversity outputs from the rarefied ASV matrix.
Builds summary metrics, runs the active statistical tests, and writes figure outputs.
"

plot_set_order <- c("standard", "normalization")

selected_metric <- metric_order <- c(
  "percent_Muri", "percent_Propi", "percent_Archaea",
  "richness", "evenness", "shannon", "n"
)

###############################################################################
message("28 # END SETUP")
###############################################################################

###############################################################################
message("32 # RUNS")
###############################################################################

for(plot_set in plot_set_order){

  ###############################################################################
  message('Plot-Set Outputs')
  ###############################################################################
  
  if(plot_set == 'standard'){
    normalization_order = c('fixed cycles')
    data_class_order    = c('doubleton')
    normalization_name  = paste(normalization_order, collapse = '_')
    data_class_name     = paste(data_class_order, collapse = '_')
    
    metric_filter  = setdiff(metric_order, c())
    selected_metric = c(
      'percent_Muri','percent_Lachno',
      'percent_Propi','percent_Archaea',
      'richness','evenness','shannon','n'
    )
    
    figure_title = 'Figure 3 Alpha Diversity Metrics by Sample Type and Annealing Temperature'
  }
  
  if(plot_set == 'normalization'){
    normalization_order = c('fixed cycles','targeted fluorescence')
    data_class_order    = c('doubleton')
    normalization_name  = paste(normalization_order, collapse = '_')
    data_class_name     = paste(data_class_order, collapse = '_')
    
    metric_filter  = setdiff(metric_order, c())
    selected_metric = c(
      'percent_Muri','percent_Lachno',
      'percent_Propi','percent_Archaea',
      'richness','evenness','shannon','n'
    )
    
    figure_title = 'Figure 6 Alpha Diversity Effects from Auto-Normalization'
  }
  
  if(plot_set == 'data_class'){
    normalization_order = c('targeted fluorescence')
    data_class_order    = c('singleton','doubleton')
    normalization_name  = paste(normalization_order, collapse = '_')
    data_class_name     = paste(data_class_order, collapse = '_')
    
    metric_filter  = setdiff(metric_order, c())
    selected_metric = c(
      'unique_features','percent_Muri','percent_Lachno',
      'percent_Propi','percent_Archaea',
      'richness','evenness','shannon','n'
    )
  }
  
  ###############################################################################
  message('Dataset Selection')
  ###############################################################################
  
  make_new_data_files <- "yes"

  starting_r <- 3
  ending_r   <- 3
  testing <- 3

  ###############################################################################
  message('Data Processing Settings')
  ###############################################################################
  
  use_custom_labels <- 'no'
  
  ###############################################################################
  message('Parameter Sets')
  ###############################################################################
  
  parameter_sets <- list(
    set1 = list(
      filter1_group = 'data_class',
      filter2_group = 'normalization',
      x_axis_group  = 'temperature',
      x_facet_group = 'sample_type',
      plot_group    = 'sample_type'
    ),
    set2 = list(
      filter1_group = 'data_class',
      filter2_group = 'data_class',
      x_axis_group  = 'normalization',
      x_facet_group = 'temperature',
      plot_group = 'sample_type'  
    ),
    set3 = list(
      filter1_group = 'normalization',
      filter2_group = 'normalization',
      x_axis_group  = 'data_class',
      x_facet_group = 'temperature',
      plot_group = 'sample_type'  
    )
  )
  
  ################################################################################
  message('113 # MAIN LOOP: DATASETS')
  ################################################################################
  
  Loop_Group_order = Loop_Group = 'run_once'
  if(exists('Run_Group')){
    Loop_Group = Run_Group
    Loop_Group_order = Run_Group_order
  }
  
  for(lp in Loop_Group_order){
    
    if(!exists('ending_r')){
      ending_r = nrow(matrix_names)
      message('\tsetting ending r to nrows matrix names')
    }
    if(!exists('starting_r')){
      starting_r = 1
      message('\tsetting starting r to 1')
    }
    if(exists('testing')){
      r = starting_r = ending_r = testing
      message('\trunning r= ', testing)
    }
    
    for (r in starting_r:ending_r) {
      
      #........................ dataset details .................................#
      
      taxa_levs   <- matrix_names[r, 'taxa_levs']
      data_set    <- matrix_names[r, 'data_sets']
      taxa_plural <- matrix_names[r, 'taxa_plural']
      
      if(plot_set != 'none'){
        data_set = plot_set
        data_set_order = plot_set_order
      }
      
      # Assign analysis parameters dynamically
      list2env(parameter_sets[[paste0('set', match(data_set, data_set_order))]], envir = .GlobalEnv)
      
      groups <- c('filter1_group','filter2_group','x_axis_group','x_facet_group','plot_group')
      for (var in groups) {
        group_value <- get(var)
        assign(paste0(var, "_order"), get(paste0(group_value, "_order")))
      }
      
      grouping_columns = unique(c(x_axis_group, x_facet_group, plot_group))
      test_across      = unique(c(x_facet_group, plot_group))
      group_within     = unique(grouping_columns)
      
      run_suffix = ''
      if(exists('Run_Group')){
        run_suffix = paste0('_', lp)
      }
      
      filter1_names = (paste(filter1_group_order, collapse = '_'))
      filter2_names = (paste(filter2_group_order, collapse = '_'))
      x_axis_names  = (paste(x_axis_group_order,  collapse = '_'))
      x_facet_names = (paste(x_facet_group_order, collapse = '_'))
      
      p_title <- paste0(
        taxa_levs, '_', script_title, custom_name, '_',
        filter1_names, '_', x_axis_names, '_', x_facet_names, '_',
        data_set, run_suffix
      )
      
      qPrint(p_title)
      
      #==============================================================================#
      message('Build / Load Metrics')
      #==============================================================================#
      
      if (make_new_data_files == 'yes'){
        
        #..............................................................................#
        message('163 # load and filter data')
        #..............................................................................#
        
        matrix_df1 = read.csv(matrix_names[r,'file_path'], check.names = FALSE, row.names = 1) %>%
          t() %>% as.data.frame() %>%
          xPlode_sample_name() %>%
          filter(.data[[filter1_group]] %in% filter1_group_order) %>%
          filter(.data[[filter2_group]] %in% filter2_group_order) %>%
          filter(.data[[plot_group]]   %in% plot_group_order) %>%
          filter(.data[[x_axis_group]] %in% x_axis_group_order) %>%
          filter(.data[[x_facet_group]] %in% x_facet_group_order)
        
        if(exists('Run_Group')){
          matrix_df1 = matrix_df1 %>% filter(.data[[Loop_Group]] %in% lp)
        }
        
        matrix_df = matrix_df1 %>%
          imPlode_sample_name() %>%
          mutate_all(as.numeric) %>%
          select(which(colSums(.) > 0))
        
        feature_names = names(matrix_df)
        
        #==============================================================================#
        message('Taxa Of Interest + Alpha Diversity')
        #==============================================================================#
        
        #..............................................................................#
        message('187 # get percent archaea and other taxa of interest')
        #..............................................................................#
        
        df2 = matrix_df %>%
          xPlode_sample_name() %>%
          pivot_longer(col = all_of(feature_names), names_to = 'feature', values_to = 'counts') %>%
          group_by(feature, across(any_of(test_across))) %>%
          mutate(feature_grouped_counts = sum(counts)) %>%
          filter(feature_grouped_counts > 0) %>%
          mutate(taxa = feature) %>%
          group_by(sample_name) %>%
          mutate(sample_counts = sum(counts)) %>%
          ungroup()
        
        num_unique_features_df = df2 %>%
          select(sample_name, feature, counts, any_of(group_within)) %>%
          filter(counts > 0) %>%
          group_by(across(all_of(group_within)), feature) %>%
          filter(n_distinct(sample_name) == 1) %>%
          ungroup() %>%
          group_by(across(all_of(group_within)), sample_name) %>%
          summarise(smp_unique_features = n(), .groups = "drop") %>%
          select(sample_name, smp_unique_features) %>%
          distinct()
        
        if(matrix_names[r,"taxa_levs"] == 'ASV'){
          df2 = df2 %>%
            rename(ASV = taxa) %>%
            left_join(ASV_taxa %>% select(ASV, taxa, Family)) %>%
            ungroup()
        }
        
        Archaea_df <- extract_taxa_percentage(df2, name = "Archaea", target = "archaea")
        Propi_df   <- extract_taxa_percentage(df2, name = "Propi",   target = "g__Cutibacterium")
        Muri_df    <- extract_taxa_percentage(df2, "Muri", "f__Muribaculaceae")
        
        #..............................................................................#
        message('230 # for Zymo testing ideal score')
        #..............................................................................#
        
        #..............................................................................#
        message('264 # generate alpha diversity metrics')
        #..............................................................................#
        
        metrics_df <- matrix_df %>%
          as.data.frame() %>%
          xPlode_sample_name() %>%
          rowwise() %>%
          mutate(
            mn  = mean(c_across(any_of(feature_names))),
            sd  = sd(c_across(any_of(feature_names))),
            smp_n = sum(c_across(any_of(feature_names))),
            smp_shannon   = qShannon(c_across(any_of(feature_names))),
            smp_evenness  = qEvenness(c_across(any_of(feature_names))),
            smp_richness  = qRichness(c_across(any_of(feature_names)))
          ) %>%
          select(-any_of(feature_names)) %>%
          left_join(Archaea_df) %>%
          left_join(Propi_df) %>%
          left_join(Muri_df) %>%
          group_by(across(any_of(grouping_columns))) %>%
          mutate(across(starts_with("smp_"), ~mean(., na.rm = TRUE), .names = "{sub('smp_', 'mean_', .col)}")) %>%
          rename_with(~gsub("smp_", "", .), starts_with("ratio_smp_")) %>%
          ungroup() %>%
          select(where(~!all(is.na(.))))
        
        write.csv(metrics_df, paste0(output_data, p_title, '_metrics_df.csv'), row.names = FALSE)
        
      } # end if(make_new_data_files == 'yes')
      
      if (make_new_data_files != 'yes'){
        metrics_df = read.csv(paste0(output_data, p_title, '_metrics_df.csv'), check.names = FALSE)
      }
      
      #==============================================================================#
      message('Reshape Metrics')
      #==============================================================================#
      
      result_df = data_df = metrics_df %>%
        pivot_longer(cols = starts_with("smp_"), names_to = "metric") %>%
        filter(!(sample_type %in% c('Skin')  & (metric %in% c('smp_percent_Muri','smp_percent_Archaea')))) %>%
        filter(!(sample_type %in% c('Soil')  & (metric %in% c('smp_percent_Propi','smp_percent_Muri')))) %>%
        filter(!(sample_type %in% c('Feces') & (metric %in% c('smp_percent_Propi','smp_percent_Archaea')))) %>%
        filter(metric != 'smp_percent_Lachno') %>%
        group_by(across(any_of(grouping_columns)), metric) %>%
        mutate(mean_value = mean(value, na.rm = TRUE)) %>%
        ungroup()
      
      #==============================================================================#
      message('319 # T-test if length(x_axis_group)==2s')
      #==============================================================================#
      
      unique_values <- unique(data_df[[x_axis_group]])
      
      if (length(unique(data_df[[x_axis_group]])) == 2) {
        
        t_test_results <- data_df %>%
          filter(!is.nan(value)) %>%
          group_by(across(any_of(test_across)), metric) %>%
          summarise(
            t_test_summary = {
              unique_subgroups <- unique(.data[[x_axis_group]])
              if (length(unique_subgroups) != 2) {
                list(data.frame(p.value = NA, statistic = NA))
              } else {
                subgroup_counts <- map(unique_subgroups, ~ n_distinct(.data$value[.data[[x_axis_group]] == .x]))
                if (any(unlist(subgroup_counts) == 1)) {
                  list(data.frame(p.value = NA, statistic = NA))
                } else {
                  list(broom::tidy(wilcox.test(value ~ .data[[x_axis_group]])))
                }
              }
            },
            .groups = 'drop'
          ) %>%
          unnest(t_test_summary)
        
        t_test_df <- t_test_results %>%
          select(any_of(test_across), metric, p.value)
        
        join_by <- setNames(c(test_across, "metric"), c(test_across, "metric"))
        
        result_df <- data_df %>%
          left_join(t_test_df, by = join_by) %>%
          rename(p_value = p.value) %>%
          mutate(metric = factor(metric, levels = rev(paste0("smp_", metric_order))))
      }
      
      #==============================================================================#
      message('Stats: Three Groups')
      #==============================================================================#
      
      unique_values <- unique(data_df[[x_axis_group]])
      
      if (length(unique(data_df[[x_axis_group]])) == 3) {
        
        message("368 # ANOVA Tukey's if length(x_axis_group)==3")
        
        data_df <- data_df %>%
          mutate(!!x_axis_group := factor(.data[[x_axis_group]], levels = x_axis_group_order))
        
        anova_results <- data_df %>%
          filter(!is.nan(value)) %>%
          group_by(across(any_of(c(test_across, "metric")))) %>%
          reframe(
            anova_result = {
              n_groups <- n_distinct(.data[[x_axis_group]])
              all_n <- table(.data[[x_axis_group]])
              valid <- all(all_n > 1)
              
              if (n_groups >= 2 && valid) {
                fit <- tryCatch(aov(value ~ .data[[x_axis_group]], data = cur_data()), error = function(e) NULL)
                if (!is.null(fit)) broom::tidy(fit) else data.frame(term = NA, df = NA, statistic = NA, p.value = NA)
              } else {
                data.frame(term = NA, df = NA, statistic = NA, p.value = NA)
              }
            }
          ) %>%
          ungroup() %>%
          unnest(anova_result) %>%
          mutate(term = ifelse(grepl("^\\.data\\[\\[", term), x_axis_group, term))
        
        anova_results
        
        tukey_results <- data_df %>%
          filter(!is.nan(value)) %>%
          group_by(across(any_of(c(test_across, "metric")))) %>%
          reframe(
            tukey_result = {
              fit <- tryCatch(aov(value ~ .data[[x_axis_group]], data = cur_data()), error = function(e) NULL)
              if (!is.null(fit)) {
                tukey_df <- as.data.frame(TukeyHSD(fit)[[1]])
                colnames(tukey_df) <- c("estimate","conf.low","conf.high","adj.p.value")
                tukey_df$term <- x_axis_group
                tukey_df$comparison <- rownames(tukey_df)
                rownames(tukey_df) <- NULL
                tukey_df
              } else {
                data.frame(
                  estimate = NA, conf.low = NA, conf.high = NA, adj.p.value = NA,
                  term = deparse(substitute(x_axis_group)), comparison = NA
                )
              }
            },
            .groups = "drop"
          ) %>%
          unnest(tukey_result) %>%
          mutate(
            adj.significance = case_when(
              is.na(adj.p.value)  ~ NA_character_,
              adj.p.value < 0.001 ~ "***",
              adj.p.value < 0.01  ~ "**",
              adj.p.value < 0.05  ~ "*",
              TRUE                ~ "ns"
            )
          )
        
        df_residual = anova_results %>%
          filter(term == 'Residuals') %>%
          mutate(residual_df = df) %>%
          select(any_of(test_across), residual_df) %>%
          distinct()
        
        anova_join = anova_results %>%
          filter(!term == 'Residuals') %>%
          left_join(df_residual) %>%
          mutate(anova.p.value = p.value) %>%
          mutate(anova.sig = ifelse(anova.p.value < .05, TRUE, FALSE)) %>%
          mutate(anova_df = df) %>%
          select(any_of(test_across), metric, anova.p.value, anova.sig, anova_df, residual_df)
        
        final_results <- tukey_results %>%
          left_join(anova_join, by = c(test_across, "metric"))
        
        head(final_results)
        
        results <- summary_brackets(
          df = data_df,
          axis_column = x_axis_group,
          dodge_column = NULL,
          value_column = "value",
          facet_columns = c(test_across, 'metric'),
          adj_y_max_columns = 'metric',
          testing_col = x_axis_group,
          test_type = NULL,
          adjust_method = "BH",
          dodge_width = 0.8,
          width = 0.8,
          log = FALSE,
          space_adj = 1.5,
          bracket_space_adj = 3,
          tip_adj = 2,
          run_tests = TRUE,
          sample_id = 'sample_name',
          clear_ns_ANOVA = FALSE,
          clear_ns_message = FALSE
        )
        
        df_segments <- results$segments
        df_segments
        
        df_tests <- results$tests %>%
          mutate(test_adj.p.value = adj.p.value) %>%
          select(any_of(test_across), metric, group1, group2, label_position, label_value, test_adj.p.value)
        
        df_check1 = results$tests %>%
          mutate(comparison = paste0(group2, '-', group1)) %>%
          select(any_of(names(final_results))) %>%
          mutate(df = 'tests')
        
        df_check = df_check1 %>%
          mutate(df = 'tests') %>%
          bind_rows(final_results %>% mutate(df = 'final')) %>%
          arrange(comparison, sample_type, metric)
        
        df_check
        
        df_tests_anova = final_results %>%
          separate(comparison, into = c("g1", "g2"), sep = "-", remove = FALSE) %>%
          mutate(
            g1 = factor(g1, levels = x_axis_group_order),
            g2 = factor(g2, levels = x_axis_group_order),
            group1 = if_else(as.integer(g1) < as.integer(g2), as.character(g1), as.character(g2)),
            group2 = if_else(as.integer(g1) < as.integer(g2), as.character(g2), as.character(g1))
          ) %>%
          select(-g1, -g2) %>%
          filter(anova.sig) %>%
          left_join(df_tests)
        
        df_segments <- results$segments %>%
          left_join(df_tests_anova %>% select(any_of(test_across), metric, anova.sig, group1, group2)) %>%
          filter(anova.sig) %>%
          ungroup()
      }
      
      ################################################################################
      message('500 # PLOTTING SETUP')
      ################################################################################
      
      magnify = 1
      first_height  = 14
      first_width   = 14
      second_height = 14
      second_width  = 10
      
      p_value_size   = 3 * magnify
      geom_text_size = 3 * magnify
      ind_text_size  = 3 * magnify
      strip_text_size = 11 * magnify
      x_text_size    = 9 * magnify
      
      theme_common <- theme_global(base_size = 11, magnify = 1) +
        theme(
          plot.tag = element_text(size = 36, face = "bold"),
          legend.position = "bottom",
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(face = 'plain')
        )
      
      gPlot <- function(p) {
        p = p +
          scale_color_manual(values = palette_color, labels = palette_label) +
          scale_fill_manual(values = palette_color, labels = palette_label) +
          theme_common +
          guides(
            title = '',
            color = guide_legend(order = 1),
            fill  = guide_legend(
              order = 2,
              override.aes = list(
                shape = 22,
                size  = 5,
                color = "black",
                alpha = 1
              )
            )
          ) +
          labs(title = '', color = '', fill = '', x = '', y = '')
        
        print(p)
      }
      
      ################################################################################
      message('592 # PLOTTING LOOP')
      ################################################################################
      
      plot_list <- list()
      list_category = unique(result_df[[plot_group]])
      
      for (t in list_category) {
        
        qPrint(unique(result_df[[plot_group]]))
        qPrint(t)
        
        df = result_df %>%
          filter(.data[[plot_group]] == t) %>%
          mutate(metric = factor(metric, levels = rev(paste0("smp_", metric_order))))
        
        #..............................................................................#
        message('609 # strip backgrounds with outlines')
        #..............................................................................#
        
        unique_x_labels <- unique(df[[plot_group]])
        unique_x_labels
        
        outline_color_x = 'black'
        
        s = paste0(
          'backgrounds_x <- list(element_rect(color=outline_color_x,fill = palette_color["',
          paste0(unique_x_labels, collapse='"]),element_rect(color=outline_color_x,fill = palette_color["'),
          '"]))'
        )
        print(s)
        eval(parse(t = s))
        backgrounds_x
        
        unique_y_labels = c(rep(unique_x_labels, (length(unique(df$metric)) / length(unique_x_labels))))
        unique_y_labels
        
        outline_color_y = 'black'
        
        s = paste0(
          'backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',
          paste0(unique_y_labels, collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),
          '"]))'
        )
        eval(parse(t = s))
        
        y_strip = as.data.frame(unique_y_labels) %>%
          mutate(text_color = 'white') %>%
          mutate(face = 'bold') %>%
          mutate(text_size = strip_text_size) %>%
          ungroup()
        
        x_strip = as.data.frame(unique_x_labels) %>%
          mutate(text_color = 'white') %>%
          mutate(face = 'bold') %>%
          mutate(text_size = strip_text_size) %>%
          ungroup()
        
        #..............................................................................#
           message('646 # plotting')
        #..............................................................................#
        
        df_segments_plot = df_segments %>%
          mutate(metric = gsub("smp_", "", metric)) %>%
          filter(.data[[plot_group]] == t) %>%
          filter(metric %in% selected_metric) %>%
          mutate(metric = factor(metric, levels = rev(metric_order))) %>%
          ungroup()
        
        df_tests_anova_plot = df_tests_anova %>%
          mutate(metric = gsub("smp_", "", metric)) %>%
          filter(metric %in% selected_metric) %>%
          filter(.data[[plot_group]] == t) %>%
          mutate(metric = factor(metric, levels = rev(metric_order))) %>%
          ungroup()
        
        summary_df <- df %>%
          select(sample_name, any_of(grouping_columns), metric, value, cycle) %>%
          filter(metric == 'smp_n') %>%
          distinct() %>%
          group_by(across(any_of(grouping_columns))) %>%
          mutate(max_reads = max(value)) %>%
          mutate(y_label = max_reads/2) %>%
          mutate(count = n()) %>%
          mutate(metric = 'n') %>%
          mutate(mean_cycles = mean(cycle)) %>%
          select(-c(sample_name, value)) %>%
          distinct() %>%
          filter(metric %in% selected_metric) %>%
          mutate(metric = factor(metric, levels = rev(metric_order))) %>%
          mutate(!!x_axis_group := factor(.data[[x_axis_group]], levels = x_axis_group_order)) %>%
          mutate(!!x_facet_group := factor(.data[[x_facet_group]], levels = x_facet_group_order)) %>%
          mutate(!!plot_group := factor(.data[[plot_group]], levels = plot_group_order)) %>%
          filter(metric %in% metric_filter)
        
        plot_df = df %>%
          mutate(metric = gsub("smp_", "", metric)) %>%
          filter(metric %in% selected_metric) %>%
          mutate(metric = factor(metric, levels = rev(metric_order))) %>%
          mutate(mean_value_label = signif(mean_value, 3)) %>%
          mutate(mean_value_label = if_else(mean_value_label > 1000, round(mean_value, 0), mean_value_label)) %>%
          mutate(!!x_axis_group := factor(.data[[x_axis_group]], levels = x_axis_group_order)) %>%
          mutate(!!x_facet_group := factor(.data[[x_facet_group]], levels = x_facet_group_order)) %>%
          mutate(!!plot_group := factor(.data[[plot_group]], levels = plot_group_order)) %>%
          ungroup() %>%
          filter(metric %in% metric_filter)
        
        x_label = (1 + length(x_axis_group_order)) / 2
        plot_title = paste(t, '\n', p_title)
        
        p = ggplot(plot_df, aes(x = get(x_axis_group), y = value, color = get(x_axis_group))) +
          geom_point(aes(fill = get(plot_group)), alpha = 0) +
          geom_boxplot(width = .6, outlier.size = .5) +
          geom_text(aes(label = paste0('  ', mean_value_label)), y = 0, vjust = .5, hjust = 0,
                    color = 'black', show.legend = FALSE, size = geom_text_size, angle = 90) +
          geom_text(data = summary_df, aes(y = y_label, label = paste('n =', count)),
                    vjust = -4, hjust = .5, size = geom_text_size, fontface = 'bold', color = 'black') +
          geom_text(data = summary_df, aes(x = x_label, y = y_label, label = 'avg cycles'),
                    vjust = 0, hjust = .5, size = geom_text_size, fontface = 'bold', color = 'gray30') +
          geom_text(data = summary_df, aes(y = y_label, label = round(mean_cycles, 1)),
                    vjust = 2, hjust = .5, size = geom_text_size, fontface = 'bold', color = 'gray30') +
          facet_grid2(
            metric ~ get(x_facet_group),
            strip = strip_themed(
              background_x = backgrounds_x,
              background_y = backgrounds_y,
              text_y = elem_list_text(face = y_strip$face, color = y_strip$text_color),
              text_x = elem_list_text(face = x_strip$face, color = x_strip$text_color)
            ),
            labeller = labeller(.cols = palette_label, .rows = palette_label),
            scale = "free"
          ) +
          scale_y_continuous(limits = c(0, NA)) +
          scale_x_discrete(labels = function(x) {
            color <- palette_color[x]
            label <- palette_label[x]
            paste0("<span style='color:", color, "'>", label, "</span>")
          }) +
          theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5, face = 'bold', size = x_text_size)) +
          labs(title = plot_title, x = '', y = '')
        
        p
        
        #==============================================================================#
        message('704 # Annotate P-Values: Two Groups')
        #==============================================================================#
        
        if (length(unique(data_df[[x_axis_group]])) == 2) {
          
          plot_df2 = df %>%
            mutate(metric = gsub("smp_", "", metric)) %>%
            filter(metric %in% selected_metric) %>%
            mutate(metric = factor(metric, levels = rev(metric_order))) %>%
            mutate(p_label = paste0('p = ', format(p_value, scientific = TRUE, digits = 3))) %>%
            mutate(significance = ifelse(p_value < .05, 'sig', 'ns')) %>%
            mutate(p_label = ifelse(metric == 'n', "", p_label)) %>%
            mutate(p_label_colored = paste0("<span style='color:", palette_color[significance], "'>", p_label, "</span>")) %>%
            mutate(mean_value_label = signif(mean_value, 3)) %>%
            mutate(mean_value_label = if_else(mean_value_label > 1000, round(mean_value, 0), mean_value_label)) %>%
            group_by(across(any_of(plot_group)), metric) %>%
            mutate(y_p_label = max(value)/2) %>%
            mutate(!!x_axis_group := factor(.data[[x_axis_group]], levels = x_axis_group_order)) %>%
            mutate(!!x_facet_group := factor(.data[[x_facet_group]], levels = x_facet_group_order)) %>%
            mutate(!!plot_group := factor(.data[[plot_group]], levels = plot_group_order)) %>%
            ungroup() %>%
            filter(metric %in% metric_filter) %>%
            filter(metric != 'n')
          
          q = p +
            geom_richtext(
              data = plot_df2 %>% filter(significance == 'sig'),
              aes(x = 1.5, y = y_p_label, label = p_label_colored),
              color = 'red',
              inherit.aes = FALSE,
              fill = scales::alpha("gray95", 0.1),
              show.legend = FALSE,
              hjust = 0.5,
              vjust = 1,
              size = p_value_size,
              fontface = "bold"
            ) +
            geom_richtext(
              data = plot_df2 %>% filter(significance == 'ns' | is.na(significance)),
              aes(x = 1.5, y = y_p_label, label = p_label),
              color = 'gray50',
              inherit.aes = FALSE,
              fill = scales::alpha("gray90", 0.05),
              show.legend = FALSE,
              hjust = 0.5,
              vjust = 1,
              size = p_value_size,
              fontface = "bold"
            )
          
          q
          p = q
        }
        
        #==============================================================================#
        message('759 # Annotate Brackets: Three Groups')
        #==============================================================================#
        
        if (length(unique(data_df[[x_axis_group]])) == 3) {
          
          print(p)
          
          q = p +
            geom_segment(
              data = df_segments_plot,
              aes(x = x, xend = xend, y = y, yend = yend),
              inherit.aes = FALSE
            ) +
            geom_text(
              data = df_tests_anova_plot,
              aes(x = label_position, y = label_value, label = adj.significance),
              inherit.aes = FALSE,
              size = 3,
              vjust = -.5
            ) +
            geom_text(
              data = df_tests_anova_plot,
              aes(x = label_position, y = label_value, label = scales::scientific(adj.p.value)),
              inherit.aes = FALSE,
              size = 3,
              vjust = 1.5
            ) +
            scale_y_continuous(
              limits = c(0, NA),
              expand = expansion(mult = c(0, 0.2))
            )
          
          print(q)
          p = q
        }
        
        plot_list[[t]] <- gPlot(p)
        
      } # end loop over t (list_category)
      
      #..............................................................................#
   message('800 # save plots')
      #..............................................................................#
      
      combined_plot <- wrap_plots(plot_list, nrow = 1) +
        plot_annotation(
          tag_levels = 'A',
          caption = "ns = not significant; * p < 0.05; ** p < 0.01; *** p < 0.001",
          theme = theme_plot
        )
      
      combined_plot
      
      plot_width = 14
      plot_height = 14
      qSave(figure_title, plot = combined_plot, ext = ext_list)
      
    } # end loop over r (starting_r:ending_r)
  } # end loop over lp (Loop_Group_order)
} # end loop over plot_set (plot_set_order)

###############################################################################
message("817 # END RUNS")
###############################################################################

###############################################################################
message(paste("821 # FINISHED", script_title))
###############################################################################
