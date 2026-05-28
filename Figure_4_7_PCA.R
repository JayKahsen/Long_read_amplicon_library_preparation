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

script_title <- "PCA"
script_description <- "
Creates Figure 4 and Figure 7 PCA outputs from the rarefied ASV matrix.
Calculates ordinations, permDISP, and PERMANOVA summaries for the active workflows.
"

plot_set_order <- c("standard", "normalization")

###############################################################################
message("24 # END SETUP")
###############################################################################

###############################################################################
message("28 # RUNS")
###############################################################################

for(plot_set in plot_set_order){ qPrint(plot_set)

  ###############################################################################
  message('Plot-Set Parameters')
  ###############################################################################
  
  if(plot_set == 'standard'){
    normalization_order = c('fixed cycles')
    data_class_order    = c('doubleton')
    normalization_name  = paste(normalization_order, collapse = '_')
    data_class_name     = paste(data_class_order, collapse = '_')
    figure_title        = 'Figure 4 PCA by Sample Type and Annealing Temperature'
  }
  
  if(plot_set == 'normalization'){
    normalization_order = c('fixed cycles','targeted fluorescence')
    data_class_order    = c('doubleton')
    normalization_name  = paste(normalization_order, collapse = '_')
    data_class_name     = paste(data_class_order, collapse = '_')
    figure_title        = 'Figure 7 The Effects of Auto-Normalization on Community'
  }
  
  if(plot_set == 'data_class'){
    normalization_order = c('targeted fluorescence')
    data_class_order    = c('singleton','doubleton')
    normalization_name  = paste(normalization_order, collapse = '_')
    data_class_name     = paste(data_class_order, collapse = '_')
  }
  
  ###############################################################################
  message('Dataset Selection')
  ###############################################################################
  
  starting_r <- 3
  ending_r   <- 3
  testing    = 'yes'
  testing_r  = 3
  
  number_of_permutations = 999
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
      shape_group   = 'temperature',
      color_group   = 'normalization',
      plot_group    = 'sample_type'
    ),
    set2 = list(
      filter1_group = 'data_class',
      shape_group   = 'normalization',
      color_group   = 'temperature',
      plot_group    = 'sample_type'
    ),
    set3 = list(
      filter1_group = 'normalization',
      shape_group   = 'data_class',
      color_group   = 'temperature',
      plot_group    = 'sample_type'
    )
  )
  
  ################################################################################
  message('107 # MAIN LOOP: PROCESS DATASETS')
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
    if(testing == 'yes'){
      r = starting_r = ending_r = testing_r
      message('\trunning r= ', testing_r)
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
      
      groups <- c('filter1_group','shape_group','color_group','plot_group')
      for (var in groups) {
        group_value <- get(var)
        assign(paste0(var, "_order"), get(paste0(group_value, "_order")))
      }
      
      grouping_columns = c(shape_group, color_group, plot_group)
      test_across      = c(color_group, plot_group)
      testing_group    = shape_group
      
      run_suffix = ''
      if(exists('Run_Group')){
        run_suffix = paste0('_', lp)
      }
      
      filter1_names = (paste(filter1_group_order, collapse = '_'))
      shape_names   = (paste(shape_group_order,  collapse = '_'))
      color_names   = (paste(color_group_order,  collapse = '_'))
      
      p_title <- paste0(
        taxa_levs, '_', script_title, custom_name, '_',
        filter1_names, '_', shape_names, '_', color_names, '_',
        data_set, run_suffix
      )
      
      qPrint(p_title)
      
      ################################################################################
      message('152 # LOAD AND FILTER DATA')
      ################################################################################
      
      matrix_df1 = read.csv(matrix_names[r,'file_path'], check.names = FALSE, row.names = 1) %>%
        t() %>% as.data.frame() %>%
        xPlode_sample_name() %>%
        filter(.data[[filter1_group]] %in% filter1_group_order) %>%
        filter(.data[[plot_group]]   %in% plot_group_order) %>%
        filter(.data[[shape_group]]  %in% shape_group_order) %>%
        filter(.data[[color_group]]  %in% color_group_order)
      
      if(exists('Run_Group')){
        matrix_df1 = matrix_df1 %>% filter(.data[[Loop_Group]] %in% lp)
      }
      
      matrix_df = matrix_df1 %>%
        imPlode_sample_name() %>%
        mutate_all(as.numeric) %>%
        select(which(colSums(.) > 0))
      
      feature_names = names(matrix_df)
      
      #==============================================================================#
      message('Build matrix_dfx')
      #==============================================================================#
      
      matrix_dfx = matrix_df %>%
        as.data.frame() %>%
        xPlode_sample_name() %>%
        ungroup()
      
      ################################################################################
      message('188 # PER-PLOT_GROUP PCA LOOP')
      ################################################################################
      
      plot_groups = unique(matrix_dfx[[plot_group]])
      qPrint(plot_groups)
      
      pcs_df = loadings_df = percent_explained_df = data.frame()
      
      df_permdisp = df_permanova = NULL
      
      p_vector = numeric()
      test_group_vector = f_vector = disp_p_vector = disp_f_vector = p_vector
      
      #==============================================================================#
      message('200 # Begin Plotting Group Loop For Data')
      #==============================================================================#
      
      for (pg in plot_groups){
        
        qPrint(pg)
        
        plot_group_df = matrix_dfx %>%
          filter(.data[[plot_group]] == pg) %>%
          imPlode_sample_name() %>%
          select(where(~ sum(.) > 0))
        
        df_meta = plot_group_df %>%
          xPlode_sample_name()
        
        #..............................................................................#
        message('215 # distance matrix, mds, points, percentage explained, assemble types')
        #..............................................................................#
        
        pseudo_count    = 0.5 * (min(plot_group_df[plot_group_df > 0]))
        pseudo_count_df = plot_group_df + pseudo_count
        normalized_df   = pseudo_count_df / rowSums(pseudo_count_df)
        data            = as.data.frame(clr(normalized_df))
        
        results = pca_result <- prcomp(data, scale. = FALSE)
        summary(results)
        
        pca_result_df = as.data.frame(pca_result$x) %>% ungroup()
        pca_result_df = as.data.frame(pca_result$rotation) %>% ungroup()
        
        pcs <- as.data.frame(pca_result$x) %>%
          xPlode_sample_name() %>%
          rename(x = 'PC1', y = 'PC2') %>%
          select(any_of(meta_names), x, y)
        
        percent_explained <- format(round(results$sdev^2 / sum(results$sdev^2) * 100, 1), nsmall = 1, trim = TRUE) %>%
          as.data.frame() %>%
          slice(1:2) %>%
          t() %>%
          as.data.frame() %>%
          `colnames<-`(c('x_lab', 'y_lab')) %>%
          mutate(!!plot_group := pg)
        
        percent_explained_df = rbind(percent_explained_df, percent_explained)
        pcs_df               = rbind(pcs_df, pcs)
        
        ################################################################################
        message('253 # PERMANOVA + PERMDISP')
        ################################################################################
        
        pca_scores <- pca_result$x[, 1:2]
        
        dist_matrix <- dist(data)
        
        #==============================================================================#
        message('266 # permANOVA as NMDS begin')
        #==============================================================================#
        
        if(!exists('df_meta')){ df_meta = meta_df }
        
        if(!exists('df_meta')){df_meta=meta_df}
        run_adonis <- function(formula, method = NULL) {
          set.seed(123)
          formula_string <- paste('dist_matrix ~', formula)
          formula_str <- as.formula(formula_string)
          print(formula_str)
          
          # Run adonis2 with renamed argument
          res <- adonis2(formula_str, data = df_meta, permutations = number_of_permutations, by = method)
          
          if (is.null(method)) method <- NA
          
          # Tidy the results
          as.data.frame(res) %>%
            rownames_to_column(var = "Term") %>%
            rename(
              Df = Df,
              SumOfSqs = SumOfSqs,
              R2 = R2,
              F_value = F,
              p_value = `Pr(>F)`
            ) %>%
            mutate(
              formula = formula,
              method = method,
              test = "permANOVA"
            ) %>% 
            select(method,formula,Term,everything())
        }
        
        # Function to run PERMDISP and return tidy result
        run_permdisp <- function(dist_matrix, grouping_var) {
          # Example usage: run_permdisp(dist_matrix, "Treatment")
          
          set.seed(123)
          
          # Compute betadisper object using median
          bd <- betadisper(dist_matrix, group = df_meta[[grouping_var]], type = "median")
          
          # Run permutation test
          perm_disp <- permutest(bd, permutations = number_of_permutations)
          
          # Extract distances to median
          dist_to_median <- bd$distances
          
          # Compute group-wise average distances
          group_means <- tapply(dist_to_median, bd$group, mean, na.rm = TRUE)
          group_means <- round(group_means, 2)
          
          # Create formatted string
          average_dist_median <- paste0(names(group_means), "=", group_means, collapse = "___")
          
          # Format output table
          result_df <- as.data.frame(perm_disp$tab) %>%
            rownames_to_column(var = "Term") %>%
            mutate(
              test = "permDISP",
              grouping_var = grouping_var,
              average_dist_median = average_dist_median
            )
          
          return(result_df)
        }
        
        is_valid_group <- function(varname, df_meta) {
          var <- df_meta[[varname]]
          var <- var[!is.na(var)]
          length(unique(var)) > 1
        }
        
        perm_group <- unique(c(
          if (is_valid_group(color_group, df_meta))  color_group,
          if (is_valid_group(shape_group, df_meta))  shape_group,
          if (is_valid_group(filter1_group, df_meta)) filter1_group
        ))
        
        panov_margin <- if (length(perm_group) > 0)
          run_adonis(formula = paste(perm_group, collapse = " + "), method = 'margin')
        else
          NULL
        
        panov_terms <- NULL
        
        if (is_valid_group(color_group, df_meta) && is_valid_group(shape_group, df_meta)) {
          panov_terms <- rbind(panov_terms,
                               run_adonis(formula = paste(color_group, "*", shape_group), method = 'terms'))
        }
        
        if (is_valid_group(color_group, df_meta) && is_valid_group(filter1_group, df_meta)) {
          panov_terms <- rbind(panov_terms,
                               run_adonis(formula = paste(color_group, "*", filter1_group), method = 'terms'))
        }
        
        if (is_valid_group(shape_group, df_meta) && is_valid_group(filter1_group, df_meta)) {
          panov_terms <- rbind(panov_terms,
                               run_adonis(formula = paste(shape_group, "*", filter1_group), method = 'terms'))
        }
        
        panov_onedf <- if (length(perm_group) > 0)
          run_adonis(formula = paste(perm_group, collapse = " + "), method = 'onedf')
        else
          NULL
        
        df_permanova1 <- dplyr::bind_rows(panov_margin, panov_terms, panov_onedf) %>%
          mutate(!!plot_group := pg)
        
        disp_color   <- if (is_valid_group(color_group, df_meta))  run_permdisp(dist_matrix, color_group) else NULL
        disp_shape   <- if (is_valid_group(shape_group, df_meta))  run_permdisp(dist_matrix, shape_group) else NULL
        disp_filter1 <- if (is_valid_group(filter1_group, df_meta)) run_permdisp(dist_matrix, filter1_group) else NULL
        
        df_permdisp1 <- list(disp_color, disp_shape, disp_filter1) %>%
          compact() %>%
          bind_rows() %>%
          mutate(!!plot_group := pg)
        
        print(df_permdisp1)
        print(df_permanova1)
        
        df_permanova = bind_rows(df_permanova, df_permanova1)
        df_permdisp  = bind_rows(df_permdisp,  df_permdisp1)
        
        #==============================================================================#
        message('396 # permANOVA as NMDS end')
        #==============================================================================#
        
      } # end loop over pg (plot_groups)
      
      ################################################################################
      message('COMMON THEME + PLOTTING HELPERS')
      ################################################################################
      
      annotate_text_size = 6
      vjust_margin = 1.5
      text_size = 18
      title_size = text_size + 2
      strip_text_size = title_size
      caption_size = text_size
      
      margin_size = 20
      loading_text_size = 2.1
      
      common_theme = theme(
        plot.background = element_rect(fill = palette_color['general_background']),
        panel.background = element_rect(fill = palette_color["panel_background"]),
        panel.grid.major.x = element_line(color = palette_color["gridline_major"], linewidth = 0.4),
        panel.grid.major.y = element_line(color = palette_color["gridline_major"], linewidth = 0.4),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.background = element_rect(fill = palette_color['general_background']),
        legend.text = element_text(size = title_size),
      axis.title.x = element_text(size = title_size, margin = margin(t = 10)),
      axis.title.y = element_text(size = title_size),
        axis.text = element_text(size = text_size),
        plot.caption = element_text(hjust = 0, size = caption_size)
      )
      
      gPlot <- function(p) {
        p = p +
          scale_color_manual(values = palette_color, labels = palette_label) +
          scale_size_manual(values = palette_size, labels = palette_label) +
          scale_shape_manual(values = palette_shape, labels = palette_label) +
          scale_fill_manual(values = palette_color, labels = palette_label) +
          common_theme +
          labs(title = '', color = '', shape = '')
        p
      }
      
      tPlot <- function(p) {
        p = p +
          theme_common +
          theme(
            strip.background = element_rect(fill = "gray40", color = "black"),
            strip.text = element_text(color = "white", face = 'bold'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0, "pt"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x  = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_blank(),
            legend.position = 'none'
          ) +
          scale_fill_manual(values = palette_color, labels = palette_label, guide = "none") +
          scale_x_discrete(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0))
        print(p)
        p
      }
      
      #..............................................................................#
      message('438 # pre loop initializing')
      #..............................................................................#
      
      plot_list_bottom = plot_list_richness = plot_list_log10 = plot_list = plot_list_anova = plot_list_disp <- list()
      plot_list = plot_list_test <- list()
      
      df_permdisp1 = df_permdisp %>%
        rename('SumOfSqs' = 'Sum Sq', 'F_value' = 'F', 'p_value' = 'Pr(>F)', 'formula' = 'grouping_var') %>%
        mutate(Term = case_when(!Term %in% c('Residual','Total','Residuals') ~ formula, TRUE ~ Term))
      
      df_tests <- bind_rows(df_permdisp1, df_permanova) %>%
        select(any_of(plot_group), test, method, formula, Term, Df, R2, F_value, p_value, everything())
      
      print(df_tests)
      
      ################################################################################
      message('443 # PLOTTING LOOP (PG)')
      ################################################################################
      
      for (pg in plot_groups){ qPrint(pg)
        
        #..............................................................................#
        message('529 # annotations as nmds begin')
        #..............................................................................#
        
        col_order=c('test','method','formula','Term','average_dist_median','R2','Df','F_value','p_value')
        
        df_tests1=df_tests %>%   
          filter(.data[[plot_group]] == pg)%>%
          mutate(rownames = row_number())%>%
          filter(!Term %in% c('Residual','Total','Residuals'))%>% 
          mutate(F_value=round(F_value,2))%>%
          mutate(R2=round(R2,2))%>%
          mutate(p_value=round(p_value,3))%>%
          mutate(short_name=paste(test,method,formula,sep='_')) %>% 
          mutate(test_name=paste(short_name,rownames,sep='_')) %>% 
          mutate(across(everything(), as.character)) %>% 
          select(any_of(col_order),everything()) %>% 
          ungroup()
        
        tests_text = df_tests1 %>%
          select(any_of(col_order)) %>%
          mutate(
            label = ifelse(
              test == 'permDISP',
              paste0(test, ' ', testing, ' (', formula, '):   ADM=(', average_dist_median, ') Df=', Df, ' F=', F_value, ' p=', p_value),
              paste0(test, ' ', testing, ' (', formula, ') ', Term, ':   R2=', R2, ' Df=', Df, ' F=', F_value, ' p=', p_value)
            )
          ) %>%
          pull(label) %>%
          paste(collapse = "\n")
        
        #..............................................................................#
        message('450 # annotations as nmds end')
        #..............................................................................#
        
        col_order=c('method','formula','Term','average_dist_median','R2','Df','F_value','p_value')
        col_disp=c('formula','Term','average_dist_median','R2','Df','F_value','p_value')
        col_anova=c('method','formula','Term','R2','Df','F_value','p_value')
        
        df_tests_plot1 <- df_tests1 %>% 
          mutate(average_dist_median = str_replace_all(average_dist_median, "___", "\n")) %>% 
          mutate(formula = paste0('(~',formula,')')) %>% 
          mutate(filter_method=method) %>% 
          filter(!is.na(p_value)) %>% 
          pivot_longer(any_of(col_order), names_to = "col", values_to = "val")%>%
          filter(!is.na(val)) %>% 
          mutate(widths = case_when(
            col == "average_dist_median" ~ 10,
            col == "formula" ~ 5,
            col == "Term" ~ 5,
            col == "method" ~ 2,
            col == "test" ~ 2,
            TRUE~2
          ))%>%
          mutate(col = factor(col, levels = col_order))%>%
          mutate(test_name = factor(test_name, levels = unique(test_name))) %>% 
          group_by(short_name) %>%
          mutate(group_color = as.integer(as.factor(short_name)) %% 2 + 1)
        
        alt_colors <- rep(c("gray60","gray80"), 10)
        
        #..............................................................................#
        message('519 # first plot permdisp')
        #..............................................................................#
        
        df_disp_plot=df_tests_plot1 %>% 
          filter(test=='permDISP') %>% 
          filter(col %in% col_disp) %>% 
          ungroup() %>% 
          mutate(short_name=factor(short_name,levels=unique(short_name))) %>% 
          group_by(short_name) %>% 
          mutate(fill = ifelse(as.integer(short_name) %% 2 == 1,'gray1','gray2'))%>% 
          #filter(!col%in% c('average_dist_median','formula')) %>% 
          filter(!col%in% c('formula')) %>% 
          ungroup() %>% 
          arrange(desc(short_name)) %>%
          mutate(short_name=factor(short_name,levels=unique(short_name)))
        
        p = ggplot(df_disp_plot, aes(x = col, y = short_name, fill = fill)) +
          geom_tile(aes(width = widths), color = "black") +
          geom_text(aes(label = paste(val))) +
          facet_grid2(test ~ col, scales = "free", space = "free")
        
        plot_disp = tPlot(p)
        
        plot_list_disp[[pg]] <- plot_disp
        print(plot_disp)
        
        #..............................................................................#
        message('544 # repeat for permanova')
        #..............................................................................#
        
        df_anova_plot = df_tests_plot1 %>%
          filter(test == 'permANOVA') %>%
          filter(col %in% col_anova) %>%
          ungroup() %>%
          mutate(short_name = factor(short_name, levels = unique(short_name))) %>%
          group_by(short_name) %>%
          mutate(fill = ifelse(as.integer(short_name) %% 2 == 1, 'gray1', 'gray2')) %>%
          filter(!col %in% c('method','formula')) %>%
          ungroup() %>%
          mutate(test_main = 'permANOVA') %>%
          mutate(test_block = case_when(
            filter_method == 'margin' ~ 'by Margin',
            filter_method == 'terms' ~ 'by Term',
            TRUE ~ 'test type'
          ))
        
        if(plot_set == 'standard'){      df_anova_plot = df_anova_plot %>% filter(filter_method != 'onedf') }
        if(plot_set == 'normalization'){ df_anova_plot = df_anova_plot %>% filter(filter_method != 'onedf') }
        
        p = ggplot(df_anova_plot, aes(x = col, y = rev(test_name), fill = fill)) +
          geom_tile(aes(width = widths), color = "black") +
          geom_text(aes(label = paste(val))) +
          facet_nested(rows = vars(test_main, test_block), cols = vars(col), scales = "free", space = "free")
        
        plot_anova = tPlot(p)
        
        plot_list_anova[[pg]] <- plot_anova
        print(plot_anova)
        
        #..............................................................................#
        message('569 # annotations as nmds end')
        #..............................................................................#
        
        x_lab = percent_explained_df %>% filter(.data[[plot_group]] == pg) %>% pull(x_lab)
        y_lab = percent_explained_df %>% filter(.data[[plot_group]] == pg) %>% pull(y_lab)
        
        df_plot = pcs_df %>%
          mutate(!!sym(filter1_group) := factor(!!sym(filter1_group), levels = filter1_group_order)) %>%
          mutate(!!sym(shape_group)   := factor(!!sym(shape_group),   levels = shape_group_order)) %>%
          mutate(!!sym(color_group)   := factor(!!sym(color_group),   levels = color_group_order)) %>%
          filter(.data[[plot_group]] == pg)
        
        ################################################################################
        message('568 # STRIP BACKGROUNDS WITH OUTLINES')
        ################################################################################
        
        unique_y_labels <- rev(unique(df_plot[[plot_group]]))
        
        outline_color_y = 'black'
        
        s = paste0(
          'backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',
          paste0(unique_y_labels, collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),
          '"]))'
        )
        print(s)
        eval(parse(t = s))
        
        if (use_custom_labels == 'yes'){
          custom_y_labels = c("slope","1500TF","x3FC","slope","1500TF","x3FC",'24cycles',rep('Soil',3),rep('Zymo',3),'NTC')
          unique_y_labels = custom_y_labels
          s = paste0(
            'backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',
            paste0(unique_y_labels, collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),
            '"]))'
          )
          print(s)
          eval(parse(t = s))
        }
        
        y_strip = as.data.frame(unique_y_labels) %>%
          mutate(text_color = ifelse(unique_y_labels == 'Soil', 'white', 'black')) %>%
          mutate(text_color = 'white') %>%  # FLAG: overwrites previous line
          mutate(face = 'bold') %>%
          mutate(text_size = strip_text_size) %>%
          ungroup()
        
        unique_x_labels = unique(df_plot[[plot_group]])
        
        outline_color_x = 'black'
        
        s = paste0(
          'backgrounds_x <- list(element_rect(color=outline_color_x,fill = palette_color["',
          paste0(unique_x_labels, collapse='"]),element_rect(color=outline_color_x,fill = palette_color["'),
          '"]))'
        )
        print(s)
        eval(parse(t = s))
        
        if (use_custom_labels == 'yes'){
          custom_x_labels = c("slope","1500TF","x3FC","slope","1500TF","x3FC",'24cycles',rep('Soil',3),rep('Zymo',3),'NTC')
          custom_y_labels
          unique_y_labels = custom_y_labels
          s = paste0(
            'backgrounds_y <- list(element_rect(color=outline_color_y,fill = palette_color["',
            paste0(unique_y_labels, collapse='"]),element_rect(color=outline_color_y,fill = palette_color["'),
            '"]))'
          )
          print(s)
          eval(parse(t = s))
        }
        
        x_strip = as.data.frame(unique_x_labels) %>%
          mutate(text_color = ifelse(unique_x_labels == 'Soil', 'white', 'black')) %>%
          mutate(text_color = 'white') %>%  # FLAG: overwrites previous line
          mutate(face = 'bold') %>%
          mutate(text_size = strip_text_size) %>%
          ungroup()
        
        ################################################################################
        message('653 # PCA PLOTTING')
        ################################################################################
        
        p = ggplot(df_plot, aes(x = x, y = y)) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          geom_vline(xintercept = 0, linetype = "dashed") +
          geom_point(aes(fill = !!sym(color_group), shape = !!sym(shape_group), size = !!sym(shape_group))) +
          xlab(paste('PC1:\t', x_lab, '% Explained')) +
          ylab(paste('PC2:\t', y_lab, '% Explained')) +
          facet_wrap2(
            formula(paste("~", plot_group)),
            strip.position = "right",
            strip = strip_themed(
              background_x = backgrounds_x,
              text_x = elem_list_text(face = x_strip$face, size = x_strip$text_size, color = x_strip$text_color),
              background_y = backgrounds_y,
              text_y = elem_list_text(face = y_strip$face, size = y_strip$text_size, color = y_strip$text_color)
            ),
            scale = 'free'
          ) +
          guides(
            size = "none",
            color = guide_legend(title = 'Ellipse Color:'),
            fill  = guide_legend(title = 'Point Color:', override.aes = list(shape = 22, size = 5)),
            shape = guide_legend(title = 'Point Shape:', override.aes = list(size = 5))
          ) +
          ggtitle(paste(p_title, pg))
        
        print(p)
        
        if(plot_set == 'standard'){
          p = ggplot(df_plot, aes(x = x, y = y)) +
            geom_hline(yintercept = 0, linetype = "dashed") +
            geom_vline(xintercept = 0, linetype = "dashed") +
            geom_point(aes(fill = !!sym(shape_group), shape = !!sym(shape_group), size = !!sym(shape_group))) +
            xlab(paste('PC1:\t', x_lab, '% Explained')) +
            ylab(paste('PC2:\t', y_lab, '% Explained')) +
            facet_wrap2(
              formula(paste("~", plot_group)),
              strip.position = "right",
              strip = strip_themed(
                background_x = backgrounds_x,
                text_x = elem_list_text(face = x_strip$face, size = x_strip$text_size, color = x_strip$text_color),
                background_y = backgrounds_y,
                text_y = elem_list_text(face = y_strip$face, size = y_strip$text_size, color = y_strip$text_color)
              ),
              scale = 'free'
            ) +
            guides(
              size = "none",
              color = guide_legend(title = 'Ellipse Color:'),
              fill  = guide_legend(title = 'Point Color:', override.aes = list(shape = 22, size = 5)),
              shape = guide_legend(title = 'Point Shape:', override.aes = list(size = 5))
            ) +
            ggtitle(paste(p_title, pg))
        }
        
        if(plot_set == 'standard'){
          p <- p + stat_ellipse(aes(group = !!sym(shape_group), color = !!sym(shape_group)), linewidth = 1.2, level = 0.95)
        } else {
          
          if (shape_group %in% colnames(df_plot) && length(unique(df_plot[[shape_group]])) > 1) {
            p <- p + stat_ellipse(aes(group = !!sym(shape_group), color = !!sym(shape_group)), linewidth = 5, level = 0.95, alpha = .3)
          }
          
          if (color_group %in% colnames(df_plot) && length(unique(df_plot[[color_group]])) > 1) {
            p <- p + stat_ellipse(aes(group = !!sym(color_group), color = !!sym(color_group)), linewidth = 1.2, level = 0.95)
          }
        }
        
        print(p)
        
        gPlot(p)
        
        plot_list[[pg]] <- gPlot(p)
        plot_list_bottom[[pg]] <- gPlot(p) + theme(legend.position = 'bottom')
        
      } # end loop over pg (plot_groups)
      
      ################################################################################
      message('691 # COMBINE + SAVE')
      ################################################################################
      
      height_ratios = c(15, 2, 3)
      height_adj = 2
      if(plot_set == 'standard'){      height_ratios = c(13, 2, 2); height_adj = 4/3 }
      if(plot_set == 'normalization'){ height_ratios = c(13, 3, 4); height_adj = 3/2 }
      
      #==============================================================================#
      message('881 # Combined Plots')
      #==============================================================================#
      
      top_row <- reduce(plot_list_bottom, `+`) +
        plot_layout(guides = "collect") &
        theme(
          legend.position = "bottom",
          legend.title = element_text(size = rel(1.8), face = 'bold')
        )
      
      middle_row <- reduce(plot_list_disp, `+`)
      bottom_row <- reduce(plot_list_anova, `+`)
      
      patch_combined_plots <- wrap_elements(top_row) /
        wrap_elements(middle_row) /
        wrap_elements(bottom_row) +
        plot_layout(heights = height_ratios)
      
      print(patch_combined_plots)
      
      plot_width  = 12 * length(plot_list) / max(1, floor(length(plot_list)/2))
      plot_height = 12 * max(1, floor(length(plot_list)/2)) * height_adj
      
      qSave(figure_title, plot = patch_combined_plots, ext = ext_list)
      
      #==============================================================================#
      message('881 # Combined Test')
      #==============================================================================#
      
      combined_test_list <- Map(function(top, bottom) {
        wrap_elements(top) / wrap_elements(bottom) +
          plot_layout(heights = c(1, 2))
      }, plot_list_disp, plot_list_anova)
      
      patch_test_plots <- wrap_plots(combined_test_list)
      print(patch_test_plots)
      
      #==============================================================================#
      message('881 # Plot List')
      #==============================================================================#
      
      patch_plots = wrap_plots(plot_list)
      
      #==============================================================================#
      message('881 # CSV')
      #==============================================================================#
      
    } # end loop over r (starting_r:ending_r)
  } # end loop over lp (Loop_Group_order)
  
} # end loop over plot_set (plot_set_order)

###############################################################################
message("813 # END RUNS")
###############################################################################

###############################################################################
message(paste("817 # FINISHED", script_title))
###############################################################################
