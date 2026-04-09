################################################################################
#source(file.path(dirname(wd), "helperJ.R"))
library(stringr)
library(tidyverse)
library(readxl)

################################################################################
# removes row names, adds meta data, adjusts data frame
################################################################################
xPlode_sample_name <- function(df) {
  message('xPloding')
  rows_df0=nrow(df)
  bad_cols <- intersect(colnames(df), meta_names)
  if (length(bad_cols) > 0) { stop(paste('duplicate columns found\n',bad_cols))}
  
  df_out=df %>%
    as.data.frame() %>%
    rownames_to_column(var='sample_name')%>%
    left_join(meta, by = "sample_name") %>%
    select(any_of(meta_names),everything())%>%
    ungroup()
  
  if(nrow(df_out)!=rows_df0){stop('meta has duplicate rows')}
  return(df_out)
}

################################################################################
# adds row names,removes meta data
################################################################################
imPlode_sample_name <- function(df) {
  message('xPloding')
  df %>%
    as.data.frame() %>%
    column_to_rownames(var='sample_name')%>%
    select(-any_of(meta_names))
}

################################################################################
# Load package libraries
################################################################################
load_libraries <- function(required_libraries, quietly = TRUE) {
  # Example: load_libraries(required_libraries)
  stopifnot(is.character(required_libraries))
  
  for (lib in required_libraries) {
    if (!requireNamespace(lib, quietly = TRUE)) {
      install.packages(lib, dependencies = TRUE)
    }
    suppressPackageStartupMessages(
      library(lib, character.only = TRUE, quietly = quietly, warn.conflicts = FALSE)
    )
  }
  
  print(paste0("Loaded libraries: ", paste(required_libraries, collapse = ", ")))
  invisible(required_libraries)
}

################################################################################
# Citation text for libraries and base R
################################################################################
cite_R <- function(libs = NULL, file_name = "R_package_citations.txt") {
  # Example:
  # cite_R()                      # uses default package list
  # cite_R(c("vegan", "ggplot2")) # custom list
  
  if (is.null(libs)) {
    libs <- c(
      "tidyverse", "ggplot2", "seqinr", "vegan", "gridExtra",
      "compositions", "ggh4x", "scales", "grid", "ggtext", "patchwork"
    )
  }
  
  # Always include R citation automatically
  libs <- c("base", libs) |> unique()
  
  # Clear existing file
  cat("", file = file_name)
  
  for (pkg in libs) {
    
    # -------------------------
    # Handle BASE R citation
    # -------------------------
    if (tolower(pkg) == "base") {
      raw <- capture.output(citation())
      
      # Find boundaries of first citation block (the real R citation)
      start <- grep("^  R Core Team", raw)
      next_block <- grep("^  \\w", raw)[-1]  # next top-level block
      
      if (length(start) > 0 && length(next_block) > 0) {
        end <- min(next_block[next_block > start[1]]) - 1
        cit <- raw[start[1]:end]
      } else {
        cit <- raw
      }
      
    } else {
      cit <- capture.output(citation(pkg))
    }
    
    # Remove BibTeX entries
    cutoff <- grep("A BibTeX entry", cit)
    if (length(cutoff) > 0) cit <- cit[1:(cutoff[1] - 1)]
    
    # Write header
    cat(
      "\n#############################################\n",
      "Citation for", pkg, "\n",
      "#############################################\n",
      file = file_name, append = TRUE
    )
    
    # Write citation
    cat(paste(cit, collapse = "\n"), "\n\n", 
        file = file_name, append = TRUE)
  }
  
  # Screen output
  cat("\n✅ All citations saved to:", normalizePath(file_name), "\n\n")
  cat(readLines(file_name), sep = "\n")
}

##################################################################################
# dataframe,vector of column names to be factorized levels will be in same order
#
qMatrix=function(df,x,y,d,f='mean'){
  df=df
  ux=unique(df[,x])
  uy=unique(df[,y])
  m=matrix(0,length(uy),length(ux))
  row.names(m)=uy
  colnames(m)=ux
  s=paste0("ag=aggregate(",d,"~",x,"+",y,",data=df,",f,")")
  print(s)
  eval(parse(t=s))
  for (c in ux){
    for(r in uy){
      m[r,c]=ag[ag[,y]==r & ag[,x]==c,d]
      
    }}
  m
}
################################################################################
######### simple aggregate will add more cName, factorize ####################################
################################################################################
# dataframe,named vector or list?
#qAg=function(df='df',iv='iv',dv='dv',fun='mean'){

sAg=function(df,values,categories,fn="mean"){
  #df=df;av=values;avs=categories;fn=fn;
  s=paste0("aggregate(cbind(",paste(values,collapse=","),")~",paste(categories,collapse="+"),",data=df,",fn,")")
  print(s)
  eval(parse(t=s))
  
  
}

################################################################################
######### quick aggregate will add more cName, factorize ####################################
################################################################################
# dataframe,named vector or list?
#qAg=function(df='df',iv='iv',dv='dv',fun='mean'){

qAg=function(df,av,avs,fn="mean"){
  df='df';av='av';avs='avs';fn='fn'
  df=df;av=av;avs=avs;fn=fn
  s=paste0("aGdf=aggregate(cbind(",paste(av,collapse=","),")~",paste(avs,collapse="+"),",data=df,",fn,")")
  print(s)
  eval(parse(t=s))
  for(d in dv){
    s=paste0("aGdf$",d,"=as.character(aGdf$",d,")")
    eval(parse(t=s))
  }
  
  for (r in 1:nrow(aGdf)){
    aGdf$cNames[r]=paste(aGdf[r,1:length(dv)],collapse=".")
    aGdf$Jitx[r]=1-(jitter(1,30))
    aGdf$Jity[r]=1-(jitter(1,30))
  }
  dv=c(dv,"cNames")
  aGdf=qFactordf(aGdf,dv)
  aGdf
}
################################################################################
#log 2
# natural_log((count+1)/(geometric_mean (count+1)))
################################################################################
qLog2=function(x,s='off'){
  if(s=='on'){print('qLog2 vector');print(x)}
  
  #Test for negatives
  if (any(x < 0)) {
    stop("Data frame contains negative values. CLR transformation requires non-negative values.")
  }
  if(sum(x)==0){return(NA)}
  log_base_2=log2(x/mean(x[x > 0]))
  
  # print(paste('mean(x)',mean(x)))
  print(paste('mean(x[x > 0]))',mean(x[x > 0])))
  print(paste('log_base_2=log2(x/mean(x[x > 0]))',log_base_2=log2(x/mean(x[x > 0]))))
  
  return(log_base_2)
  if(sum(x)==0){return(NAN)}
}
################################################################################
################################################################################
#centered log transformation
# natural_log((count+1)/(geometric_mean (count+1)))
################################################################################
# ignores zeros
qRCLR=function(x,s='off'){
  if(s=='on'){print('qCLR vector');print(x)}
  
  #Test for negatives
  if (any(x < 0)) {
    stop("Data frame contains negative values. CLR transformation requires non-negative values.")
  }
  centered_log=x
  centered_log[centered_log>0] <- log(x[x>0]) - mean(log(x[x>0]))
  return(centered_log)
  
}
################################################################################
# ignores zeros
qCLR=function(x,s='off'){
  if(s=='on'){print('qCLR vector');print(x)}
  
  #Test for negatives
  if (any(x < 0)) {
    stop("Data frame contains negative values. CLR transformation requires non-negative values.")
  }
  centered_log=x
  centered_log[centered_log>0] <- log(x[x>0]) - mean(log(x[x>0]))
  return(centered_log)
  
}

################################################################################
################################################################################
# PermANOVA Test
# perm_ANOVA_p=qPermANOVA(.,sample_name,FP,ASV,counts,s='on',r='p')
################################################################################
qPermANOVA <- function(df, row_names, grouping, features, values, r = 'p', s = 'off') {
  
  if (s == 'on') {
    print('PermANOVA row_names')
    print(row_names)
    print('PermANOVA grouping')
    print(grouping)
    print('PermANOVA features')
    print(features)
    print('PermANOVA values')
    print(values)
  }
  
  raw_df <- df
  
  df <- as.data.frame(cbind(row_names, grouping, features, values))
  colnames(df) <- c('row_name', 'grp', 'feature', 'value')
  
  pre_matrix <- df %>%
    pivot_wider(names_from = feature, values_from = value)
  
  distance_matrix <- pre_matrix %>%
    select(-grp) %>%
    mutate_at(-1, as.numeric) %>%
    column_to_rownames(var = 'row_name')
  
  tryCatch({
    result <- adonis2(distance_matrix ~ pre_matrix$grp, permutations = 999)
    if (r == 'p') {
      return(result$Pr[1])
    }
    if (r == 'R2') {
      return(result$R2[1])
    }
    if (r == 'Df') {
      return(result$Df[1])
    }
  }, error = function(e) {
    if(s=='on'){print(paste("Error:", e$message))}
    return(NA)
  })
}

################################################################################
ev=function(x){eval(parse(text=x))}

pasteP=function(x){paste(x,collapse="+")}
pasteC=function(x){paste(x,collapse=",")}
pasteQ=function(x){noquote(paste0(x))} 
nQ=function(x){noquote(x)} 
pasteQC=function(x){noquote(paste(x,collapse=","))}
evP=function(x){eval(parse(text=paste0(x)))}

###############################################################################

qPrint <- function(var = NULL) {
  # Get the name of the variable from the environment if not provided
  if (is.null(var)) {
    env <- parent.frame()
    objects <- ls(envir = env)
    vars <- objects[sapply(objects, function(x) !is.function(get(x, envir = env)))]
    
    # Check if there are any variables
    if (length(vars) > 0) {
      last_var_name <- tail(vars, 1)
      var_name=last_var_name
      var=get(var_name, envir = env)
    } else {
      # If no variables found, check the global environment
      env <- .GlobalEnv
      objects <- ls(envir = env)
      vars <- objects[sapply(objects, function(x) !is.function(get(x, envir = env)))]
      
      if (length(vars) > 0) {
        last_var_name <- tail(vars, 1)
        var_name=last_var_name
        var=get(var_name, envir = env)
      } else {
        cat("No variables found in the global environment.\n")
        return(NULL)  # Return NULL if no variables are found
      }
    }
  }else{
    var_name <- deparse(substitute(var))
  }
  
  # Convert the variable to a string, handling vectors appropriately
  if (is.vector(var)) {
    value_str <- paste(var, collapse = ", ")
  } else {
    value_str <- toString(var)
  }
  
  # Remove slashes from the string
  value_str <- gsub("/", "", value_str)
  
  # Create the output string
  output_string <- sprintf('%s = %s', var_name, value_str)
  
  # Print the output
  library(crayon)
  cat(bgGreen$black$bold(paste0('\n\t', output_string, '\t\n')))
  
  
}

#qPrint()  # Prints the last variable in the calling environment
################################################################################
# Define the tPrint function
tPrint <- function(message='') {
  library(crayon)
  # Get the current system time
  current_time <- Sys.time()
  
  # Print the message followed by the current time
  cat(bgBlue$black$bold(paste0(message, " at ", format(current_time, "%Y-%m-%d %H:%M:%S")), "\n"))
}
#black,red,green,yellow,blue,magenta,cyan,white
# Example usage
tPrint("start")


################################################################################
# Kruskal-Wallace Test
# qKW(df, df$clr_counts, df$FP)
# df%>%
#   group_by(ASV)%>%
#   summarize(KW=qKW(.,clr_counts,FP))
################################################################################
qKW=function(df,x,grouping,s='off'){
  if(s=='on'){print('Kruskal-Wallace vector');print(x)}
  #Test for negatives
  KW=kruskal.test(x ~ grouping, data = df)
  return(KW$p.value)
}
################################################################################
#return index from population vector 
#input a vector representing populations
################################################################################
# population vector
qShannon=function(x,s='off'){
  if(s=='on'){print('shannon vector');print(x)}
  x=x[x!=0];if(length(x)==0){return(NA)}
  s=sum(x);-sum(x/s*log(x/s))}

qEvenness=function(x,s='off'){
  if(s=='on'){print('evenness vector');print(x)}
  x=x[x!=0];if(length(x)==0){return(NA)};if(length(x)==1){return(1)}
  s=sum(x);-sum(x/s*log(x/s))/log(length(x))}

# y = reference population ratios. default is even distribution
# if distribution uneven need to make sure reference aligns
qIdealScore <- function(x, y = rep(1 / length(x), length(x)), s = 'off') {
  if (length(y) > length(x)) {
    x <- c(x, rep(0, length(y) - length(x)))
  }
  if (s == 'on') {
    print('Ideal Score vector')
    print(x)
  }
  sum(abs((x / sum(x)) - (y / sum(y)))) * 100
}

qBray <- function(x, y, s='off') {
  if (s == 'on') {
    print('Input vectors:')
    print(x)
    print(y)
  }
  sum(abs(x - y)) / sum(x + y)
}

qEU <- function(x, y, s='off') {
  if (s == 'on') {
    print('Input vectors:')
    print(x)
    print(y)
  }
  sum((x - y)^2)^0.5
}


################################################################################
inverse_simpson <- function(x){
  n <- sum(x)
  1/(sum(x * (x-1) / (n * (n-1)))) 
}

qRichness<- function(x){
  sum(x>0)
} # compare to vegan specnumber(x) 1 - sum((x/n)^2) 

shannon<-function(x){
  n <- sum(x)
  relative_abundance = ra=x[x>0]/n
  -sum(ra*log(ra))
} # compare to vegan diversity(x,index='shannon') 


simpson<- function(x){
  n <- sum(x)
  sum(x * (x-1) / (n * (n-1))) 
} # compare to vegan diversity(x,index='simpson') 1 - sum((x/n)^2
##################################################################
hash_line <- function(line_symbol = "#", text_string = NULL, width = 80, use_spaces = TRUE) {
  # Example:
  # hash_line("=")                       # full line of "="
  # hash_line("=", "Load Libraries")     # banner with centered text
  
  # If no text: print full line
  if (is.null(text_string)) {
    line <- paste0("#", paste(rep(line_symbol, width - 2), collapse = ""), "#")
    cat(line, "\n")
    return(invisible(line))
  }
  
  # If text supplied: centered banner
  text <- if (use_spaces) paste0(" ", text_string, " ") else text_string
  
  text_len <- nchar(text)
  remaining <- width - text_len - 2  # two border "#"
  
  left_pad  <- floor(remaining / 2)
  right_pad <- ceiling(remaining / 2)
  
  banner <- paste0(
    "#",
    paste(rep(line_symbol, left_pad), collapse = ""),
    text,
    paste(rep(line_symbol, right_pad), collapse = ""),
    "#"
  )
  
  cat(banner, "\n")
  invisible(banner)
}




################################################################################
################################################################################
log_runtime <- function(start_time, end_time, function_name, log_file = "runtime_log.txt") {
  runtime <- end_time - start_time
  log_entry <- paste(Sys.time(), "Function:", function_name, "Runtime:", runtime, "seconds", "\n")
  
  # Append the log entry to the specified log file
  write(log_entry, file = log_file, append = TRUE)
}
################################################################################
# write file to csv (names vector, filepath, prefix)
################################################################################
qWrite <- function(df, name = NULL, folder = NULL) {
  # Example: qWrite(loadings_long_df)
  # Example: qWrite(loadings_long_df, "my_name")
  
  if (is.null(name)) {
    name <- deparse(substitute(df))
  }
  
  if (is.null(folder)) {
    folder <- get("output_data", envir = parent.frame())
  }
  
  filename <- paste0(folder, name, ".csv")
  qPrint(filename)
  write.csv(df, filename, row.names = FALSE)
}


################################################################################

################################################################################
qSave <- function(x,
                  plot = NULL,
                  folder = NULL,
                  ext = "png",
                  dpi = NULL,
                  plot_width = NULL,
                  plot_height = NULL,
                  show_plot = TRUE) {
  # Example: qSave("hambone", ext = c("png","pdf","svg"))
  
  ext <- tolower(gsub("^\\.", "", ext))
  ext <- as.character(ext)
  
  # Case 1: x is a plot object
  if (inherits(x, "ggplot")) {
    plot <- x
    name <- deparse(substitute(x))
  }
  
  # Case 2: x is a name (character)
  if (is.character(x)) {
    name <- x
    plot <- if (is.null(plot)) ggplot2::last_plot() else plot
  }
  
  if (is.null(plot) || is.null(name)) {
    stop("qSave() requires either a ggplot object or a character name.")
  }
  
  # Resolve output folder
  if (is.null(folder)) {
    folder <- get("output_plot", envir = parent.frame())
  }
  
  # Pull defaults from environment if present
  if (is.null(plot_width)  && exists("plot_width",  envir = parent.frame()))
    plot_width  <- get("plot_width",  envir = parent.frame())
  if (is.null(plot_height) && exists("plot_height", envir = parent.frame()))
    plot_height <- get("plot_height", envir = parent.frame())
  if (is.null(dpi) && exists("dpi", envir = parent.frame()))
    dpi <- get("dpi", envir = parent.frame())
  
  out_files <- character(0)
  
  for (ext_i in ext) {
    
    filename <- paste0(folder, name, ".", ext_i)
    qPrint(filename)
    out_files <- c(out_files, filename)
    
    args <- list(
      filename = filename,
      plot = plot
    )
    
    if (!is.null(plot_width))  args$width  <- plot_width
    if (!is.null(plot_height)) args$height <- plot_height
    
    # Raster formats (dpi only applies there)
    if ((ext_i %in% c("png", "jpg", "jpeg", "tiff", "bmp")) && !is.null(dpi)) {
      args$dpi <- dpi
    }
    
    # Vector formats
    if (ext_i == "svg") args$device <- "svg"
    if (ext_i == "pdf") args$device <- "pdf"
    
    do.call(ggplot2::ggsave, args)
  }
  
  if (isTRUE(show_plot)) print(plot)
  
  invisible(out_files)
}


################################################################################
qHead <- function(data, n_cols = 6, n_rows = 3, max_colname_length = 20) {
  # Check if input is a data frame or matrix
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input must be a data frame or matrix.")
  }
  
  # Truncate column names and ensure uniqueness
  colnames(data) <- str_trunc(colnames(data), width = max_colname_length, side = "right", ellipsis = "...")
  
  # Ensure unique column names by adding a prefix if needed
  colnames(data) <- make.unique(colnames(data), sep = "_")
  
  # Print a subset of the data (up to n_cols and n_rows)
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  
  # Print the data
  print(head(data %>% select(1:n_cols), n_rows))
}
################################################################################
# mean centering a data frame by category sub groups
################################################################################
mean_normalizing_by_category <- function(df, category_col) {
  # df with columns = features rows = observations row.names=identifier
  data_cols <- names(df)
  
  # Apply xPlode_sample_name to remove rownames and add meta data tied to them
  exploded_df <- df %>% xPlode_sample_name()
  
  # Identify the meta_data columns and data columns
  meta_data_cols <- setdiff(names(exploded_df), data_cols)
  
  mean_centered_df <- exploded_df %>%
    split(.[[category_col]]) %>%
    lapply(function(sub_df) {
      # Separate meta_data and data columns
      sub_meta_data <- sub_df[, meta_data_cols, drop = FALSE]
      sub_data <- sub_df[, data_cols, drop = FALSE]
      
      # Calculate the column-wise (feature) means
      sub_mean_standard <- colMeans(sub_data)
      
      sub_mean_centered_data <- sub_data %>%
        t() %>%
        `/`(sub_mean_standard) %>%
        t() %>%
        as.data.frame() %>%
        mutate(across(everything(), ~ replace_na(.x, 0))) # handle 0 count columns
      
      
      # Combine meta_data with normalized data columns
      sub_mean_centered_df <- cbind(sub_meta_data, sub_mean_centered_data)
      
      return(sub_mean_centered_df)
    }) %>%
    do.call(rbind, .) %>%
    {rownames(.) <- NULL; .} %>%  # Remove row names
    imPlode_sample_name() # Apply imPlode_sample_name to reverse xPlode
  
  return(mean_centered_df)
}
################################################################################
# Define a function to perform NMDS with increasing k
nmds_with_increasing_k <- function(data, max_tries = 100, max_iterations = 999) {
  
  #example: results <- nmds_with_increasing_k(data)
  
  # Compute distance matrix
  dist_matrix <- vegdist(data, method = "bray")
  
  # Initialize variables to store results
  stress_values <- c()
  k_values <- c()
  
  # Start with k = 1 and increase until convergence or max_tries reached
  k <- 1
  converged <- FALSE
  
  while (!converged && k <= max_tries) {
    # Perform PCA for initialization
    qPrint(k)
    set.seed(123)
    pca_result <- prcomp(dist_matrix, scale. = FALSE)
    init <- pca_result$x[, 1:k]
    
    # Perform NMDS with current k
    nmds_result <- metaMDS(dist_matrix, k = k, 
                           maxit = 999, 
                           trymax = 250,
                           autotransform = FALSE, init = init,wascores = TRUE)
    
    # Check convergence
    if (nmds_result$converged) {
      converged <- TRUE
      # Store stress value and k
      stress_values <- c(stress_values, nmds_result$stress)
      k_values <- c(k_values, k)
    } else {
      # Store non-converged stress (NA) and k
      stress_values <- c(stress_values, NA)
      k_values <- c(k_values, k)
    }
    
    # Increase k
    k <- k + 1
    qPrint(k)
  }
  
  # Return a data frame with k values and corresponding stress
  results_df <- data.frame(K = k_values, Stress = stress_values)
  qPrint(results_df)
  write.csv(results_df,'nmds_k_results.csv')
  return(nmds_result)
}
################################################################################
# xAgg aggregate an exploded matrix by sub dividing 
################################################################################

xAgg <- function(dfx, aggregate_columns) {
  # dfx=df1; aggregate_columns='primer_set'
  # Example usage: result_df <- xAgg(dfx, c('Lat_Location','Type','Species'))
  qPrint(aggregate_columns)
  # Check if all specified columns exist in the data frame
  missing_cols <- aggregate_columns[!aggregate_columns %in% colnames(dfx)]
  
  # If there are missing columns, return an error
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are not found in the data frame:", paste(missing_cols, collapse = ", ")))
  }
  
  # Remove excess metadata columns that are not part of the aggregation
  dfx <- dfx %>% select(-any_of(setdiff(meta_names, aggregate_columns)))
  
  # Get unique combinations of the aggregation columns
  group_combinations <- dfx %>%
    select(any_of(aggregate_columns)) %>%
    distinct()
  
  # Initialize an empty data frame to hold the summarized results
  dfx_summarized <- data.frame(matrix(ncol = ncol(dfx) - length(aggregate_columns), 
                                      nrow = nrow(group_combinations)))
  colnames(dfx_summarized) <- names(dfx)[!names(dfx) %in% aggregate_columns]  # Set column names to ASV names
  
  # Loop through each unique group combination to calculate sums
  for (combination_row in 1:nrow(group_combinations)) {
    group_name <- paste(group_combinations[combination_row, ], collapse = '_')
    tPrint(group_name)
    
    if(ncol(group_combinations)==1){
      # Subset the data frame to include only the current group combination
      ag_group <- as.data.frame(group_combinations[combination_row, ]) %>%
        setNames(aggregate_columns)%>%
        left_join(dfx) %>%
        select(-all_of(aggregate_columns))
    }else{
      
      # Subset the data frame to include only the current group combination
      ag_group <- group_combinations[combination_row, ]%>%
        left_join(dfx) %>%
        select(-all_of(aggregate_columns))
    }  
    
    # Calculate column sums for the current group and store them in the summarized data frame
    dfx_summarized[combination_row, ] <- colSums(ag_group, na.rm = TRUE)
  }
  
  # Combine the group combinations with their corresponding summarized data
  dfx <- cbind(group_combinations, dfx_summarized)
  
  return(dfx)  # Return the resulting data frame
}
################################################################################
# xAgg_prevalence aggregate an exploded matrix by subdividing (only columns that meet prevalence) 
################################################################################

xAgg_prevalence <- function(dfx, aggregate_columns,threshold) {
  # dfx=df1; aggregate_columns='primer_set'
  # Example usage: result_df <- xAgg(dfx, c('Lat_Location','Type','Species'))
  qPrint(aggregate_columns)
  # Check if all specified columns exist in the data frame
  missing_cols <- aggregate_columns[!aggregate_columns %in% colnames(dfx)]
  
  # If there are missing columns, return an error
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are not found in the data frame:", paste(missing_cols, collapse = ", ")))
  }
  
  # Remove excess metadata columns that are not part of the aggregation
  dfx <- dfx %>% select(-any_of(setdiff(meta_names, aggregate_columns)))
  
  # Get unique combinations of the aggregation columns
  group_combinations <- dfx %>%
    select(any_of(aggregate_columns)) %>%
    distinct()
  
  # Initialize an empty data frame to hold the summarized results
  dfx_summarized <- data.frame(matrix(ncol = ncol(dfx) - length(aggregate_columns), 
                                      nrow = nrow(group_combinations)))
  colnames(dfx_summarized) <- names(dfx)[!names(dfx) %in% aggregate_columns]  # Set column names to ASV names
  
  # Loop through each unique group combination to calculate sums
  for (combination_row in 1:nrow(group_combinations)) {
    group_name <- paste(group_combinations[combination_row, ], collapse = '_')
    tPrint(group_name)
    
    if(ncol(group_combinations)==1){
      # Subset the data frame to include only the current group combination
      ag_group <- as.data.frame(group_combinations[combination_row, ]) %>%
        setNames(aggregate_columns)%>%
        left_join(dfx) %>%
        select(-all_of(aggregate_columns))
      
      prevalence <- colSums(ag_group > 0) / nrow(ag_group)
      ag_group[, prevalence < threshold] <- 0
    }else{
      
      # Subset the data frame to include only the current group combination
      ag_group <- group_combinations[combination_row, ]%>%
        left_join(dfx) %>%
        select(-all_of(aggregate_columns))
      
      prevalence <- colSums(ag_group > 0) / nrow(ag_group)
      ag_group[, prevalence < threshold] <- 0
    }  
    
    # Calculate column sums for the current group and store them in the summarized data frame
    dfx_summarized[combination_row, ] <- colSums(ag_group, na.rm = TRUE)
  }
  
  # Combine the group combinations with their corresponding summarized data
  dfx <- cbind(group_combinations, dfx_summarized)
  
  return(dfx)  # Return the resulting data frame
}
###############################################################################
# xAgg_mean aggregate an exploded matrix by subsampling by mean
################################################################################

xAgg_mean <- function(dfx, aggregate_columns) {
  
  # Example usage: result_df <- xAgg(dfx, c('Lat_Location','Type','Species'))
  
  # Check if all specified columns exist in the data frame
  missing_cols <- aggregate_columns[!aggregate_columns %in% colnames(dfx)]
  
  # If there are missing columns, return an error
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are not found in the data frame:", paste(missing_cols, collapse = ", ")))
  }
  
  # Remove excess metadata columns that are not part of the aggregation
  dfx1 <- dfx %>% select(-any_of(setdiff(meta_names, aggregate_columns)))
  
  # Get unique combinations of the aggregation columns
  group_combinations <- dfx1 %>%
    select(aggregate_columns) %>%
    distinct()
  
  # Initialize an empty data frame to hold the summarized results
  dfx_summarized <- data.frame(matrix(ncol = ncol(dfx1) - length(aggregate_columns), 
                                      nrow = nrow(group_combinations)))
  colnames(dfx_summarized) <- names(dfx1)[!names(dfx1) %in% aggregate_columns]  # Set column names to ASV names
  
  # Loop through each unique group combination to calculate sums
  for (combination_row in 1:nrow(group_combinations)) {
    group_name <- paste(group_combinations[combination_row, ], collapse = '_')
    tPrint(group_name)
    
    if(ncol(group_combinations)==1){
      # Subset the data frame to include only the current group combination
      ag_group <- as.data.frame(group_combinations[combination_row, ]) %>%
        setNames(aggregate_columns)%>%
        left_join(dfx1) %>%
        select(-all_of(aggregate_columns))
    }else{
      
      # Subset the data frame to include only the current group combination
      ag_group <- group_combinations[combination_row, ]%>%
        left_join(dfx1) %>%
        select(-all_of(aggregate_columns))
    }  
    
    # Calculate column sums for the current group and store them in the summarized data frame
    dfx_summarized[combination_row, ] <- colMeans(ag_group, na.rm = TRUE)
  }
  
  # Combine the group combinations with their corresponding summarized data
  dfx2 <- cbind(group_combinations, dfx_summarized)
  
  return(dfx2)  # Return the resulting data frame
}
##############################################################################
# generate simplified taxon labels
#df = df_matrix;group = 'sample_type'
##############################################################################
feature_labels <- function(group = 'feature', df = df_matrix,averaging_group='sample_name') {
  # Ensure 'feature' is included in the group list
  feature_names <- names(df)
  if (!('feature' %in% group)) {
    group <- c(group, 'feature')
  }
  
  # Function to generate the new feature label based on taxonomy
  generate_new_feature_label <- function(feature) {
    # Example: generate_new_feature_label("k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__;g__Lactobacillus;s__uncultured")
    
    # Split the feature taxonomic chain
    feature_split <- strsplit(feature, ";")[[1]]
    num_segments <- length(feature_split)
    
    # Start from the tail and walk backward to find the first valid piece
    idx <- num_segments
    while (
      idx > 1 &&
      (grepl("__$", feature_split[idx]) ||
       grepl("uncultured", feature_split[idx], ignore.case = TRUE) ||
       grepl("unclassified", feature_split[idx], ignore.case = TRUE))
    ) {
      idx <- idx - 1
    }
    
    # Always include the valid segment
    start_idx <- idx
    
    # If the next one is uncultured/unclassified, include one step deeper too
    if (idx < num_segments) {
      next_seg <- feature_split[idx + 1]
      if (grepl("uncultured", next_seg, ignore.case = TRUE) ||
          grepl("unclassified", next_seg, ignore.case = TRUE)) {
        start_idx <- idx
        num_segments <- idx + 1
      } else {
        num_segments <- idx
      }
    } else {
      num_segments <- idx
    }
    
    # Return the formatted label from the selected segment onward
    return(paste(feature_split[start_idx:num_segments], collapse = ";"))
  }
  
  
  # Process the input data frame
  df_processed <- df %>%
    as.data.frame() %>%
    xPlode_sample_name() %>% 
    pivot_longer(cols = all_of(feature_names), names_to = 'feature', values_to = 'counts')
  
  if (taxa_levs == 'ASV') {
    df_processed <- df_processed %>%
      mutate(ASV=feature)%>%
      left_join(ASV_taxa %>% select(ASV, taxa), by = "ASV")
  } else {
    df_processed <- df_processed %>%
      mutate(taxa = feature)
  }
  
  df_processed1<-df_processed%>%
    mutate(
      Domain = sapply(taxa, extract_domain),  # Assuming extract_domain() is defined elsewhere
      feature_label_original = sapply(taxa, generate_new_feature_label),
      feature_label = if_else(feature_label_original == '__', 'unclassified', feature_label_original),
      feature_label = sub("^p__", "", feature_label)  # Clean up 'p__' prefix if needed
    ) %>%
    
    group_by(feature)%>%
    mutate(feature_counts=sum(counts))%>%
    ungroup() %>%
    mutate(total_counts=sum(counts))%>%
    mutate(total_rel_abun=feature_counts/total_counts)%>%
    
    group_by(sample_name) %>%
    mutate(sample_counts = sum(counts)) %>%
    ungroup() %>%
    mutate(rel_abun=counts/sample_counts)%>%
    
    group_by(across(all_of(c(group,averaging_group))), feature_label,feature_label_original, Domain,taxa,total_rel_abun) %>%
    summarize(average_rel_abun = mean(rel_abun),
              mean_feature_counts=mean(feature_counts),
              mean_total_counts=mean(total_counts),
              mean_total_rel_abun=mean(total_rel_abun),
              summarized_rows = n(), 
              .groups = "drop") %>%
    ####
    #group_by(feature, feature_label, Domain) %>%
    #mutate(group_rel_abun = rel_abun) %>%
    group_by(across(all_of(group)), feature_label,feature_label_original, Domain,taxa,total_rel_abun) %>%
    summarize(group_mean_rel_abun = mean(average_rel_abun),
              group_mean_feature_counts = mean(mean_feature_counts),
              group_mean_total_counts = mean( mean_total_counts),
              group_mean_total_rel_abun = mean(mean_total_rel_abun),
              mean_prevous_summarized_rows=mean(summarized_rows),
              summarized_rows = n(), 
              .groups = "drop") %>%
    ungroup() %>%
    mutate(
      rel_abun_label = paste0(signif(100 * group_mean_rel_abun, 2), '%'),
      total_rel_abun_label = paste0(signif(100 * total_rel_abun, 2), '%'),
      feature_label2 = paste0("*", feature_label, "* ", rel_abun_label),
      feature_label3 = paste0(feature_label, ' ', signif(100 * group_mean_rel_abun, 2), '%'),
      total_feature_label = paste0("*", feature_label, "* ", total_rel_abun_label),
    ) %>%
    arrange(desc(group_mean_rel_abun)) %>%
    mutate(plot_order_test = factor(feature_label3, levels = rev(unique(feature_label3))))%>%
    mutate(plot_order = factor(feature_label2, levels = rev(unique(feature_label2))))
  
  return(df_processed1)
}
##############################################################################
# Kruskal Wallis or Mann_Whitney_Test NEW
#df = df_matrix; testing_group = y_axis_group; category_group = test_across
##############################################################################
Kruskal_Mann_Whitney_Test <- function(df = df_matrix, testing_group = "y_axis_group", category_group = "test_across", collapse_data_used = FALSE) {
  # Example usage: kw_group_results = Kruskal_Mann_Whitney_Test(df_matrix, plotting_set)
  
  dfx <- df %>%
    xPlode_sample_name() %>%
    ungroup()
  
  kw_group_results <- data.frame()
  
  unique_categories <- dfx %>%
    distinct(across(all_of(category_group))) %>%
    split(seq(nrow(.)))
  
  for (category in unique_categories) {
    category_df <- as.data.frame(category)
    category_str <- paste(category_df, collapse = "_")
    category_vals <- as.list(category_df[1, , drop = TRUE])

    keep_rows <- rep(TRUE, nrow(dfx))
    for (i in seq_along(category_group)) {
      keep_rows <- keep_rows & dfx[[category_group[i]]] == category_vals[[i]]
    }

    category_data1 <- dfx[keep_rows, , drop = FALSE]
    rownames(category_data1) <- NULL
    
    # Prepare category columns for skipped rows
    category_columns <- as.list(category_df)
    names(category_columns) <- category_group
    
    if (nrow(category_data1) < 2) {
      warning(paste("Skipping category", category_str, "due to insufficient rows"))
      kw_group_results <- bind_rows(kw_group_results,
                                    tibble(
                                      feature = NA,
                                      estimate = NA,
                                      statistic = NA,
                                      p.value = 1,
                                      testing = testing_group,
                                      test_type = NA,
                                      adjusted_pvalue = 1,
                                      significance = "99",
                                      adjusted_significance = "99",
                                      adj_check = NA,
                                      !!!category_columns  # unquote-splice category columns
                                    )
      )
      next
    }
    
    category_data <- category_data1 %>%
      imPlode_sample_name() %>%
      select(where(~ sum(.) != 0)) %>%
      # {
      #   pseudo_count <- min(.[. > 0], na.rm = TRUE) / 2
      #   . + pseudo_count
      # } %>%
      t() %>% as.data.frame() %>%
      mutate(across(everything(), ~ . / sum(.))) %>%
      t() %>%
      #clr() %>%
      as.data.frame() %>%
      xPlode_sample_name()
    
    check <- category_data %>%
      select(all_of(testing_group)) %>%
      distinct() %>%
      pull()
    
    cat("Number of unique categories in", testing_group, "for", category_str, ":", length(check), "\n")
    
    if (length(check) < 2) {
      warning(paste("Skipping category", category_str, "due to insufficient unique categories"))
      kw_group_results <- bind_rows(kw_group_results,
                                    tibble(
                                      feature = NA,
                                      estimate = NA,
                                      statistic = NA,
                                      p.value = 1,
                                      testing = testing_group,
                                      test_type = NA,
                                      adjusted_pvalue = 1,
                                      significance = "ns",
                                      adjusted_significance = "ns",
                                      adj_check = NA,
                                      !!!category_columns
                                    )
      )
      next
    }
    
    test_type <- if (length(check) == 2) {
      "Welch’s t-test"
    } else {
      "Kruskal-Wallis"
    }
    test_func <- if (length(check) == 2) wilcox.test else kruskal.test
    
    for (group in testing_group) {
      group_var <- category_data[[group]]
      
      if (length(unique(group_var)) > 1) {
        # new
        feature_data <- category_data %>% select(-any_of(meta_names))
        
        # changing
        # test_results <- category_data %>%
        #   select(-any_of(meta_names)) %>%
        test_results <- feature_data %>%
          map_df(function(.x) {
            tryCatch(
              broom::tidy(test_func(.x ~ group_var)),
              error = function(e) {
                tibble(
                  estimate = NA,
                  statistic = NA,
                  p.value = 1
                )
              }
            )
          }, .id = "feature") %>%
          
          mutate(data_used = if (collapse_data_used) {
            map_chr(.$feature, function(f) {
              vals <- feature_data[[f]]
              grouped_vals <- split(vals, group_var)
              paste(
                names(grouped_vals),
                map(grouped_vals, ~ paste0(signif(.x, 4), collapse = ", ")),
                sep = ": ",
                collapse = " | "
              )
            })
          } else {
            NA_character_
          }) %>%
          mutate(
            testing = group,
            test_type = test_type,
            adjusted_pvalue = p.adjust(p.value, method = "BH"),
            significance = case_when(
              p.value < 0.001 ~ "***",
              p.value < 0.01 ~ "**",
              p.value < 0.05 ~ "*",
              TRUE ~ "ns"
            ),
            adjusted_significance = case_when(
              adjusted_pvalue < 0.001 ~ "***",
              adjusted_pvalue < 0.01 ~ "**",
              adjusted_pvalue < 0.05 ~ "*",
              TRUE ~ "ns"
            ),
            adj_check = adjusted_pvalue / p.value
          ) %>%
          bind_cols(tibble(!!!category_columns))
        
        
        kw_group_results <- bind_rows(kw_group_results, test_results)
      } else {
        warning(paste("Skipping group", group, "in category", category_str, "due to only one unique category."))
        kw_group_results <- bind_rows(kw_group_results,
                                      tibble(
                                        feature = NA,
                                        estimate = NA,
                                        statistic = NA,
                                        p.value = 1,
                                        testing = group,
                                        test_type = NA,
                                        adjusted_pvalue = 1,
                                        significance = "ns",
                                        adjusted_significance = "ns",
                                        adj_check = NA,
                                        !!!category_columns
                                      )
        )
      }
    }
  }
  
  return(kw_group_results)
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# end? ANCOMBC
#global_formula <- paste(c(testing_group, category_group), collapse = " + ")
################################################################################
ANCOM_BC_Test <- function(df, testing_group, category_group=NA, adjust_method = "BH",testing_formula=NA,contrast=NULL,tax_level = NULL) {
  # Example:
  # ANCOM_group_results <- ANCOM_BC_Test(df = df_matrix, testing_group = y_axis_group, category_group = test_across)
  
  #testing
  # category_group=NA; adjust_method = "BH";testing_formula=NA;contrast=NULL
  # df = df_matrix; testing_group = y_axis_group; category_group = test_across
  # tax_level='Genus'
  
  #table(sample_data(physeq_obj)$extraction_method)
  
  
  library(dplyr)
  library(tidyr)
  library(phyloseq)
  library(ANCOMBC)
  
  # if (inherits(testing_formula, "formula")) {
  #   testing_formula1 <- testing_formula
  # } else if (is.character(testing_formula)) {
  #   if (!is.na(testing_formula)) {
  #     testing_formula1 <- as.formula(paste("~",testing_formula))
  #   } else {
  #     testing_formula1 <- as.formula(paste("~", testing_group))
  #   }
  # } else {
  #   testing_formula1 <- as.formula(paste("~", testing_group))
  # }
  
  
  if (!is.na(testing_formula)) {
    testing_formula1 <- testing_formula
  } else {
    testing_formula1 <- testing_group
  }
  
  dfx <- df %>%
    xPlode_sample_name() %>%
    select(-setdiff(meta_names,c('sample_name',category_group,testing_group))) %>% 
    ungroup()
  
  
  if (length(category_group) == 1 && is.na(category_group)) {
    dfx$default_grouping_category='none'
    category_group='default_grouping_category'
  }
  
  ANCOM_group_results <- data.frame()
  
  # Extract unique categories for grouping
  unique_categories <- dfx %>%
    distinct(across(all_of(category_group))) %>%
    split(seq(nrow(.)))
  
  # Iterate over each category
  for (category in unique_categories) {
    qPrint(category)
    category_str <- paste(category, collapse = "_")
    
    # Create a list of category values
    category_vals <- as.list(category)
    
    # Filter data for the current category
    category_data1 <- dfx %>%
      filter(if_any(all_of(category_group), ~ . %in% category)) %>%
      select(-any_of(c("default_grouping_category"))) %>% 
      imPlode_sample_name() %>%
      select(where(~ sum(.) != 0)) %>%
      xPlode_sample_name()
    
    OTU <- category_data1 %>% imPlode_sample_name()
    
    # Prepare taxa data
    TAX1 <- OTU %>%
      t() %>%
      as.data.frame() %>%
      mutate(Taxa = rownames(.)) %>%
      select(Taxa)
    
    # Handle ASV-specific processing
    if (taxa_levs == 'ASV') {
      if (!exists("ASV_taxa")) {
        ASV_taxa <- read.csv('data_tables/raw_ASV_taxa.csv')  # Or update with the actual path
      }
      TAX1 <- TAX1 %>%
        rename(ASV = Taxa) %>%
        left_join(ASV_taxa, by = "ASV")
    }
    
    # Taxonomy column splitting
    TAX <- TAX1 %>%
      mutate(Taxa = str_replace_all(Taxa, " ", "")) %>%  # Remove all spaces
      separate(Taxa, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               sep = ";", fill = "right") %>%
      mutate(across(everything(), ~ sub("^.*__", "", .)))
    
    # Extract and prepare metadata
    META <- category_data1 %>%
      select(any_of(meta_names)) %>%
      mutate(row_name=sample_name) %>%
      column_to_rownames(var='row_name')
    
    # Create phyloseq object
    physeq_obj <- phyloseq(
      otu_table(as.matrix(OTU), taxa_are_rows = FALSE),
      tax_table(as.matrix(TAX)),
      sample_data(META)
    )
    # strsplit(testing_formula1, split = "\\s*\\+\\s*")
    # strsplit(testing_formula, split = "\\s*\\+\\s*")
    # Run ANCOMBC2 on the phyloseq object
    ancom_out <- ancombc2(
      data = physeq_obj,
      tax_level = tax_level,
      fix_formula = testing_formula1, # can use formula
      p_adj_method = adjust_method,
      prv_cut = 0.10, # minimum prevalence
      lib_cut = 1000, # minimum sample count
      s0_perc = 0.05, # sensitivity adjustment?
      group = testing_group,
      alpha = 0.05, # significance limit 
      n_cl = 1, # nodes to be forked?
      verbose = TRUE,
      global = FALSE, 
      #lfc_control = list(contrast = contrast),
      iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE),
      em_control = list(tol = 1e-05, max_iter = 100),
      lme_control = lme4::lmerControl(),
      mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
      trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)
    )
    names(META)
    # Extract the results
    # Determine available prefixes in the column names
    available_parts <- unique(gsub("_.+", "", colnames(ancom_out$res)[-1]))  # remove 'taxon'
    
    # Build the pivot pattern dynamically
    pattern <- paste0("(", paste(available_parts, collapse = "|"), ")_(.*)")
    
    # Pivot the result
    res_df <- ancom_out$res %>%
      as.data.frame() %>%
      pivot_longer(
        cols = -taxon,
        names_to = c(".value", "location"),
        names_pattern = pattern
      ) %>%
      {
        # If 'q' is present, rename to adjust_method
        if ("q" %in% available_parts) {
          rename(., !!paste0(adjust_method, "_value") := q)
        } else {
          .
        }
      } %>%
      mutate(
        significance = if ("p" %in% available_parts) ifelse(p < 0.05, "yes", "no") else NA,
        adjusted_significance = if ("q" %in% available_parts) ifelse(!!sym(paste0(adjust_method, "_value")) < 0.05, "yes", "no") else NA,
        test_type = "ANCOM-BC",
        testing = testing_group,
        !!paste(category_group, collapse = "_") := category_str,
        !!!category_vals
      ) %>%
      select(
        taxon,
        any_of(c(paste0(adjust_method, "_value"), "adjusted_significance", "lfc", "location",
                 "p", "significance", "se", "W", "diff", "passed_ss")),
        test_type,
        testing,
        all_of(category_group),
        !!paste(category_group, collapse = "_")
      )
    
    # Combine the results from each category
    ANCOM_group_results <- bind_rows(ANCOM_group_results, res_df)
  }
  ANCOM_group_results=ANCOM_group_results %>% 
    mutate(feature=taxon,
           p.value=p,
           adjusted_pvalue=BH_value) %>% 
    mutate(significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )) %>%
    mutate(adjusted_significance = case_when(
      adjusted_pvalue < 0.001 ~ "***",
      adjusted_pvalue < 0.01 ~ "**",
      adjusted_pvalue < 0.05 ~ "*",
      TRUE ~ "ns"
    )) %>%
    mutate(adj_check = adjusted_pvalue / p.value) %>%
    ungroup()
  
  return(ANCOM_group_results)
}



################################################################################
# end?
################################################################################
sc_log_labels <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else if (val == 1*sc) {
      return(expression(10^0))
    } else {val=val/sc
    return(as.expression(bquote(10^.(round(log10(val))))))
    }
  })
}
################################################################################
################################################################################
# end?
################################################################################
library(scales)

scale_x_log10_shift <- function(shift = 0, ...) {
  # Create a shifted log10 transformation
  scale_x_continuous(
    trans = trans_new(
      "log10_shifted",
      transform = function(x) log10(x) + shift,   # add shift to move right
      inverse   = function(x) 10^(x - shift)      # subtract shift to map back
    ),
    ...
  )
}



################################################################################
################################################################################
# end?
################################################################################
log10_labels_percent <- function(x) {
  # Converts:
  # 1, 10, 100 → "1", "10", "100"
  # 0.01, 0.1, 1 → ".01", ".1", "1"
  
  out <- vapply(x, function(val) {
    if (is.na(val)) return(NA_character_)
    if(exists('log_shift')){val=val/log_shift}
    
    if (val < 1) {
      s <- format(val, digits = 6, nsmall = 3, scientific = FALSE, trim = TRUE)
      s <- sub("^0", "", s)
      s <- sub("0+$", "", s)
      s <- sub("\\.$", "", s)
    } else {
      s <- formatC(val, format = "f", digits = 0)
    }
    s
  }, character(1))
  
  out
}


################################################################################
# end?
################################################################################
library(scales)

# transformation with adjustable factor
log10_decimal <- function(factor = 1) {
  offset <- log10(factor)
  trans_new(
    name      = paste0("log10_factor_", factor),
    transform = function(y) log10(y) + offset,
    inverse   = function(y) 10^(y - offset),
    breaks    = log_breaks(base = 10),
    format    = label_number(accuracy = 1)   # 1 sig digit
  )
}

#   scale_x_continuous(trans = log10_decimal(factor = 100)) +

################################################################################
# end?
################################################################################
sc_labels <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else {val=val/sc
    return(val)
    }
  })
}
################################################################################
# end?
################################################################################
log_labels2 <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else if (val == 1) {
      return(expression(2^0))
    } else {val=val
    return(as.expression(bquote(2^.(round(log2(val))))))
    }
  })
}
################################################################################
# end?
################################################################################
custom_strip_text <- function(strip_colors) {
  theme(
    strip.text = element_text(color = function(label) strip_colors[label])
  )
}
################################################################################
# end?
################################################################################
adjust_labels <- function(x) {
  scales::trans_format("log", math_format(exp(.x)))(x)
}
custom_breaks <- function(x) {
  max_number=max(as.numeric(x))
  breaks=unique(x)
  return(breaks)
}
################################################################################
# end?
################################################################################
custom_labels <- function(x) {
  labels=c((rep('',(length(x)-1))),paste(length(x),taxa_plural))
  
  return(labels)
}
################################################################################
# end?
################################################################################
custom_breaks_labels <- function(x) {
  unique_x <- unique(x)  # Ensure breaks are unique values
  breaks <- tail(unique_x, 1)  # Get the last unique value if needed
  
  qPrint(breaks)
  qPrint(labels)
  
  labels <- length(unique_x)
  assign("labels_from_breaks", labels, envir = .GlobalEnv)
  qPrint(breaks)
  qPrint(labels)
  qPrint(labels_from_breaks)
  return(breaks)
}
################################################################################

################################################################################
extract_domain <- function(taxa) {
  if (grepl("archaea", taxa, ignore.case = TRUE)) {
    return("Archaea")
  } else if (grepl("bacteria", taxa, ignore.case = TRUE)) {
    return("Bacteria")
  } else if (grepl("eukaryot", taxa, ignore.case = TRUE)) {
    return("Eukaryote")
  } else {
    return("unknown")
  }
}
################################################################################

################################################################################
# Taxonomy extraction function (same as before)
extract_taxa <- function(taxa, target) {
  #extract_taxa(taxa='Muri,target='f__Muribaculaceae')
  if (grepl(target, taxa, ignore.case = TRUE)) {
    return(target)
  } else {
    return("no")
  }
}

# Generalized function to extract and calculate percentage for any given taxonomic level
extract_taxa_percentage <- function(long_dfx, name, target) {
  #extract_taxa_percentage(df=df,taxa='Muri',target='f__Muribaculaceae')
  # long_dfx=df_matrix%>%
  #   dfx=df_matrix%>%
  #     xPlode_sample_name()%>%
  #     pivot_longer(col=all_of(feature_names),names_to='feature',values_to='counts')%>%
  #     group_by(feature,across(any_of(test_across)))%>%
  #     mutate(feature_grouped_counts=sum(counts))%>%
  #     filter(feature_grouped_counts>0)%>%
  #     mutate(taxa=feature)%>%
  #     group_by(sample_name) %>%
  #     mutate(sample_counts=sum(counts))%>%
  #     ungroup()
  #   
  #   if(matrix_names[r,"taxa_levs"]=='ASV'){
  #     
  #     long_dfx=long_dfx%>%
  #       rename(ASV=taxa)%>%
  #       left_join(ASV_taxa%>%select(ASV,taxa,Family))%>%
  #       ungroup()
  #   }
  long_dfx%>%
    mutate(found = sapply(taxa, extract_taxa, target = target)) %>%
    # select(sample_name,feature,counts,found,sample_counts)%>%
    # distinct()%>%
    group_by(sample_name)%>%
    summarise(
      taxa_counts = sum(counts[found == target]),
      sample_counts=sum(counts),
      # sample_counts = mean(sample_counts),
      .groups = 'drop'
    ) %>%
    mutate(taxa_percentage = 100 * taxa_counts / sample_counts) %>%
    ungroup() %>%
    rename(!!paste0("smp_percent_", name) := taxa_percentage) %>%
    select(sample_name, !!paste0("smp_percent_", name))
}

################################################################################

################################################################################
fold_change_log2 <- function(x, y) {
  mean_x <- mean(x) + 1
  mean_y <- mean(y) + 1
  fold_change <- mean_x / mean_y
  fold_change_log2 <- log2(fold_change)
  return(fold_change_log2)
}
################################################################################
# checks for files existence before trying to read
################################################################################
qLoad <- function(file_path) {
  # Check if the file exists
  if (!file.exists(file_path)) {
    message(paste('\t',file_path, "does not exist."))
  }else{
    
    # Extract the filename without extension
    file_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Determine file type and read accordingly
    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
      data <- read.csv(file_path)
    } else if (grepl("\\.xlsx$", file_path, ignore.case = TRUE)) {
      data <- read_excel(file_path)
    } else {
      stop("\tUnsupported file type. Please provide a CSV or XLSX file.")
    }
    # Assign the data frame to a variable with the filename
    assign(file_name, data, envir = .GlobalEnv)
    
    message("\tFile successfully loaded:\t", basename(file_path))
  }
}

write_formatted_excel <- function(df, file_path, sheet_name = "Sheet1") {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet = sheet_name, x = df, withFilter = TRUE)

  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center"
  )

  body_style <- openxlsx::createStyle(
    halign = "center",
    valign = "center"
  )

  openxlsx::addStyle(
    wb, sheet = sheet_name, style = header_style,
    rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE
  )

  if (nrow(df) > 0) {
    openxlsx::addStyle(
      wb, sheet = sheet_name, style = body_style,
      rows = 2:(nrow(df) + 1), cols = seq_len(ncol(df)),
      gridExpand = TRUE, stack = TRUE
    )
  }

  openxlsx::freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet = sheet_name, cols = seq_len(ncol(df)), widths = "auto")
  openxlsx::saveWorkbook(wb, file = file_path, overwrite = TRUE)
}
################################################################################
# creates directory from filepath... checks existence first
################################################################################
create_directory <- function(dir_path) {
  # Check if the directory exists
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Directory created: ", dir_path)
  } else {
    message("Directory already exists: ", dir_path)
  }
}
################################################################################

################################################################################
vPrint <- function(x) {
  cat("c(", paste0("'", x, "'", collapse = ","), ")")
}

################################################################################
nvPrint <- function(x) {
  #stopifnot(is.vector(x) && !is.null(names(x)))  
  cat(paste0("'", names(x), "'='", x, "'", collapse = ",\n"), "\n")
}
################################################################################
qnvPrint <- function(x) {
  names(x)=x
  #stopifnot(is.vector(x) && !is.null(names(x)))  
  cat(paste0("'", names(x), "'='", x, "'", collapse = ",\n"), "\n")
}
################################################################################
# Function to get the reverse complement of a DNA sequence
reverse_complement <- function(seq) {
  as.character(reverseComplement(DNAString(seq)))
}
################################################################################
log_breaks <- 10^seq(0, 12, by = 2)
log_breaks <- c(0,10^seq(0, 12, by = 1))
log10_breaks <-seq(0, 12, by = 1)

log2_breaks <-seq(0, 12, by = 1)
fold_change_breaks <-seq(-12, 12, by = 1)


log10_labels_bold <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else {
      return(paste0("10^", val))
    }
  })
}
fold_change_labels_bold <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else {
      return(paste0("2^", val))
    }
  })
}
sc_log10_labels_bold <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else {
      return(paste0("10^", val - sc))
    }
  })
}

sc_log10_labels <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else {
      val=val-sc
      return(as.expression(bquote(10^.(round(val)))))
    }
  })
}



log_labels <- function(x) {
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else if (val == 1) {
      return(expression(10^0))
    } else {
      return(as.expression(bquote(10^.(round(log10(val))))))
    }
  })
}

log_labels_bold <- function(x) {
  # Requires theme(axis.text.y = element_markdown())
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else {
      exponent <- round(log10(val))
      return(paste0("10^", exponent))
    }
  })
}
log2_labels_bold <- function(x) {
  # Requires theme(axis.text.y = element_markdown())
  sapply(x, function(val) {
    if (is.na(val)) {
      return(NA)  # Return NA if the value is missing
    } else {
      exponent <- round(log2(val))
      return(paste0("2^", exponent))
    }
  })
}
################################################################################
# Function for significance asterisks
sig_stars <- function(p) {
  case_when(
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE       ~ "ns"
  )
}


################################################################################
################################################################################
# function to run summarize data, run significance tests, and supply bracket and plotting info
################################################################################
summary_brackets <- function(
    df=df_plot1,
    axis_column = 'axis_name',
    dodge_column = NULL,
    value_column = "value",
    facet_columns = 'facet_name',
    testing_col='axis_name',
    test_type = NULL,
    adjust_method = "BH",# does this apply?
    dodge_width = 0.8,
    width = 0.8,
    log = FALSE,
    space_adj = 1,
    bracket_space_adj=1,
    tip_adj = 1,
    run_tests = TRUE,
    clear_ns_ANOVA=FALSE,
    clear_ns_message=FALSE,
    sample_id='sample_id',
    adj_y_max_columns=NULL,
    return = c("summary", "segments", "tests")
) {
  
  
  # test
  
  # df=data_df; axis_column=x_axis_group; dodge_column=NULL; value_column="value";
  # facet_columns=c(test_across,'metric'); testing_col=x_axis_group; test_type=NULL;
  # adjust_method="BH"; dodge_width=0.8; width=0.8; log=FALSE; space_adj=1;
  # bracket_space_adj=1; tip_adj=1; run_tests=TRUE; clear_ns_ANOVA=FALSE; clear_ns_message=FALSE;
  # sample_id="sample_id"; adj_y_max_columns='metric';
  
  # example
  
  #   
  # results <- summary_brackets(
  #   df = df_plot3,
  #   axis_column = "sample_type",
  #   dodge_column = "extraction_method",
  #   value_column = "value",
  #   facet_columns = NA,
  #   testing_col= "extraction_method",
  #   test_type = NULL,
  #   adjust_method = "BH",
  #   dodge_width = 0.8,
  #   width = 0.8,
  #   log = TRUE,
  #   space_adj = 1,
  # bracket_space_adj=1,
  #   tip_adj = 1,
  #   run_tests = TRUE,
  # clear_ns_ANOVA='FALSE',
  # clear_ns_message='FALSE',
  # sample_id='sample_id',
  #  )
  # df_summary  <- results$summary
  # df_segments <- results$segments  
  # df_tests    <- results$tests
  # geom_segment(
  #   data = df_segments,
  #   aes(x = x, xend = xend, y = y, yend = yend),
  #   # linewidth = 0.8,
  #   inherit.aes = FALSE
  # ) +
  #   geom_text(
  #     data = df_tests,
  #     aes(
  #       x = label_position,
  #       y = label_value ,  
  #       label = sig_label
  #     ),
  #     inherit.aes = FALSE,
  #     size = 3,vjust=-.5
  #   )+
  #   geom_text(
  #     data = df_tests,
  #     aes(
  #       x = label_position,
  #       y = label_value ,  
  #       label = p_value_label
  #     ),
  #     inherit.aes = FALSE,
  #     size = 3,vjust=1.5
  #   )+
  # geom_errorbar(
  #   data = df_summary,
  #   aes(
  #     x = axis_name,
  #     ymin = mean_minus_sd,
  #     ymax = mean_plus_sd,
  #     #group = dodge column
  #   ),
  #   position = position_dodge(width = 0.8),
  #   width = 0.1, size = 1, inherit.aes = FALSE, alpha = 0.5
  # )+
  #   geom_errorbar(  # Mean segment
  #     data = df_summary,
  #     aes(
  #       x = axis_name,
  #       ymin = mean_value,
  #       ymax = mean_value,
  #       #group = dodge column
  #     ),
  #     position = position_dodge(width = 0.8),
  #     width = 0.4, size = 1.5, inherit.aes = FALSE, alpha = 0.3
  #   )+
  
  
  library(purrr)
  
  if (is.na(facet_columns[1])) {
    df$facet_dummy <- "none"
    facet_columns <- "facet_dummy"
  }
  
  axis_sym   <- sym(axis_column)
  value_sym  <- sym(value_column)
  facet_syms <- syms(facet_columns)
  
  
  
  axis_order_name <- paste0(axis_column, "_order")
  testing_order_name <- paste0(testing_col, "_order")
  
  # Order for axis_column based on first occurrence
  
  axis_order <- df %>%
    arrange(.data[[axis_column]]) %>%
    pull({{ axis_column }}) %>%
    unique()
  
  
  # Order for testing_col based on first occurrence
  
  testing_order <- df %>%
    arrange(.data[[testing_col]]) %>%
    pull({{ testing_col }}) %>%
    unique()
  
  # testing_levels <- setNames(seq_along(testing_order), testing_order)
  
  dodge_defined <- !is.null(dodge_column) && !is.na(dodge_column)
  
  if (dodge_defined) {
    dodge_sym <- sym(dodge_column)
    summary_group_cols <- c(facet_columns, axis_column, testing_col)
    test_across=c(facet_columns,axis_column)
  } else {
    summary_group_cols <- c(facet_columns, testing_col)
    dodge_width <- 0
    test_across=unique(c(facet_columns))
  }
  
  df_group_testing_order <- df %>%
    arrange(.data[[testing_col]]) %>%
    group_by(across(all_of(test_across))) %>%
    summarize(group_testing_order = list(unique(.data[[testing_col]]))) %>%
    ungroup() 
  
  df=df %>% 
    select(any_of(c(sample_id,test_across,testing_col,value_column))) %>% 
    distinct()
  
  
  df_summary <- df %>%
    group_by(across(all_of(summary_group_cols))) %>%
    summarise(
      mean_value = mean(.data[[value_column]], na.rm = TRUE),
      sd_value = sd(.data[[value_column]], na.rm = TRUE),
      n = sum(!is.na(.data[[value_column]])),
      mean_minus_sd = mean_value - sd_value,
      mean_plus_sd = mean_value + sd_value,
      max_value = max(.data[[value_column]], na.rm = TRUE),  # <-- added here
      .groups = "drop" )%>% 
    ungroup()
  
  df_wide <- df_summary %>%
    select(all_of(c(test_across, testing_col, 'max_value'))) %>%
    pivot_wider(
      names_from = all_of(testing_col),
      values_from = 'max_value'
    )
  
  
  # Get group-specific testing order
  df_group_testing_order <- df %>%
    arrange(.data[[testing_col]]) %>%
    group_by(across(all_of(test_across))) %>%
    summarize(group_testing_order = list(unique(.data[[testing_col]])), .groups = "drop")
  
  df_combs <- df %>%
    group_by(across(any_of(test_across))) %>%
    summarise(
      groups = list(na.omit(unique(.data[[testing_col]]))),
      .groups = "drop"
    ) %>%
    mutate(
      pairs = map(groups, ~ {
        if (length(.x) < 2) return(tibble(
          group1 = character(),
          group2 = character(),
          group1_level = integer(),
          group2_level = integer()
        ))
        
        combn(sort(.x), 2, simplify = FALSE) %>%
          map(~ {
            g1 <- .x[1]
            g2 <- .x[2]
            tibble(
              group1 = g1,
              group2 = g2,
              group1_level = which(testing_order == g1),
              group2_level = which(testing_order == g2)
            )
          }) %>%
          bind_rows()
      })
    ) %>%
    select(-groups) %>%
    unnest(pairs) %>%
    left_join(
      df %>%
        group_by(across(any_of(test_across)), !!sym(testing_col)) %>%
        summarise(group_values = list(.data[[value_column]]), .groups = "drop"),
      by = c(test_across, "group1" = testing_col)
    ) %>%
    rename(group1_values = group_values) %>%
    left_join(
      df %>%
        group_by(across(any_of(test_across)), !!sym(testing_col)) %>%
        summarise(group_values = list(.data[[value_column]]), .groups = "drop"),
      by = c(test_across, "group2" = testing_col)
    ) %>%
    rename(group2_values = group_values) %>% 
    left_join(df_group_testing_order) %>% 
    rowwise() %>%
    mutate(
      group1_level = which(group1 == group_testing_order),
      group2_level = which(group2 == group_testing_order)
    ) %>%
    ungroup()
  
  df_combs
  
  
  # Add axis level if dodge is defined (only one value per row so no pairwise comparison)
  if (dodge_defined) {
    df_combs <- df_combs %>%
      mutate(axis_level = match(.data[[axis_column]], axis_order))
  }
  
  # get max values and span size
  df_segments0 <- df_combs %>%
    #  group_by(across(any_of(test_across)))
    mutate(
      seg_min = group1_level,
      seg_max = group2_level,
      group_level_min = pmin(group1_level, group2_level),
      group_level_max = pmax(group1_level, group2_level),
      span_size = group_level_max - group_level_min
    ) %>%
    left_join(df_wide)%>% 
    rowwise() %>%
    mutate(
      max_value = max(group1_values,group2_values)) %>% 
    # 
    #   
    #         unlist(cur_data()[testing_order[group_level_min:group_level_max]]),
    #         na.rm = TRUE
    #       )
    #     ) %>%
    group_by(across(any_of(adj_y_max_columns))) %>% 
    mutate(maximum_value=max(max_value)) %>% 
    ungroup() %>% 
    mutate(
      
      value_adj = if (log) {
        2 * space_adj
      } else {
        value_adj = 0.05 * maximum_value * space_adj
      },
      seg_value = if (log) {
        max_value * value_adj * span_size * bracket_space_adj
      } else {
        max_value + value_adj * span_size * bracket_space_adj
      })
  
  # if (!dodge_defined) {
  
  # Arrange by span and group positions to set a hierarchy
  df_segments <- df_segments0 %>%
    arrange(span_size,seg_value) %>% 
    mutate(final_seg_value = seg_value)
  
  # Loop to resolve overlapping brackets
  for (i in seq_len(nrow(df_segments))) { #error here
    for (j in seq_len(i - 1)) {
      # Test for horizontal overlap
      overlap <- !(df_segments$group_level_max[i] < df_segments$group_level_min[j] ||
                     df_segments$group_level_min[i] > df_segments$group_level_max[j])
      
      # Test for vertical proximity (too close)
      if(log){
        too_close <- abs(df_segments$final_seg_value[i]/df_segments$final_seg_value[j]) < df_segments$value_adj[i]
      }else{
        too_close <- abs(df_segments$final_seg_value[i] - df_segments$final_seg_value[j]) < df_segments$value_adj[i]
      }
      # If overlapping and too close, shift the current seg_value up
      if (overlap && too_close) {
        df_segments$final_seg_value[i] <- df_segments$final_seg_value[j] + df_segments$value_adj[i]
      }
    }
  }
  
  # Final segment height and tip_end position
  df_segments <- df_segments %>%
    mutate(
      seg_value = final_seg_value,
      tip_end = if (log) {
        seg_value/(value_adj * tip_adj* .6)
      } else {
        seg_value - (value_adj * tip_adj / 4)
      } ) %>%
    select(-final_seg_value)
  
  if (dodge_defined) {
    df_segments <- df_segments %>%
      mutate(
        seg_min=axis_level-dodge_width/4,
        seg_max=axis_level+dodge_width/4)
  }
  
  df_brackets_long <- df_segments %>%
    select(all_of(test_across), group1, group2, seg_min, seg_max, seg_value, tip_end) %>%
    pivot_longer(cols = c("seg_min", "seg_max"), names_to = "pos", values_to = "x_bracket") %>%
    mutate(
      part = ifelse(pos == "seg_min", "left", "right")
    ) %>%
    transmute(
      !!!syms(test_across),
      group1,
      group2,
      part,
      x = x_bracket,
      xend = x_bracket,
      y = tip_end,
      yend = seg_value
    ) %>%
    bind_rows(
      df_segments %>%
        transmute(
          !!!syms(test_across),
          group1,
          group2,
          part = "top",
          x = seg_min,
          xend = seg_max,
          y = seg_value,
          yend = seg_value
        )
    )
  
  #------------------------------------------------------------------------------#
  # start significance testing (brackets function)
  #------------------------------------------------------------------------------#
  
  # functions
  
  safe_test <- function(x, y) {
    tryCatch({
      test <- t.test(x, y)
      list(
        estimate = unname(test$estimate[[1]] - test$estimate[[2]]),
        statistic = unname(test$statistic),
        p.value = unname(test$p.value),
        conf.low = unname(test$conf.int[1]),
        conf.high = unname(test$conf.int[2])
      )
    }, error = function(e) {
      warning("Test failed: ", conditionMessage(e))
      list(
        estimate = NA_real_,
        statistic = NA_real_,
        p.value = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_
      )
    })
  }
  
  ## ============================================================
  ## Helper functions for Shapiro + tests + significance
  ## ============================================================
  
  run_single_shapiro <- function(v) {
    # v may be a vector or list-column element
    if (is.list(v)) v <- unlist(v, use.names = FALSE)
    v <- v[is.finite(v)]
    if (length(v) < 3) return(NA_real_)
    tryCatch(
      shapiro.test(v)$p.value,
      error = function(e) NA_real_
    )
  }
  
  run_welch_full <- function(x, y) {
    # x, y are vectors or list elements
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    if (is.list(y)) y <- unlist(y, use.names = FALSE)
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    out <- tryCatch(
      t.test(x, y),
      error = function(e) NULL
    )
    if (is.null(out)) {
      return(list(
        p         = NA_real_,
        estimate  = NA_real_,
        stat      = NA_real_,
        conf.low  = NA_real_,
        conf.high = NA_real_
      ))
    }
    list(
      p         = out$p.value,
      estimate  = unname(out$estimate[[1]] - out$estimate[[2]]),
      stat      = unname(out$statistic),
      conf.low  = out$conf.int[1],
      conf.high = out$conf.int[2]
    )
  }
  
  run_mann <- function(x, y) {
    # x, y are vectors or list elements
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    if (is.list(y)) y <- unlist(y, use.names = FALSE)
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    suppressWarnings(
      tryCatch(
        wilcox.test(x, y)$p.value,
        error = function(e) NA_real_
      )
    )
  }
  
  sig_labeler <- function(p) {
    dplyr::case_when(
      is.na(p)      ~ NA_character_,
      p < 0.001     ~ "***",
      p < 0.01      ~ "**",
      p < 0.05      ~ "*",
      TRUE          ~ "ns"
    )
  }
  
  # end functions
  # initialize
  
  library(broom)
  
  test_results=ttest_results=final_results=NULL
  
  if(run_tests){
    
    df_tests <- df %>%
      group_by(across(all_of(test_across))) %>%
      mutate(n_groups = n_distinct(.data[[testing_col]])) %>%
      ungroup()
    
    df_for_ANOVA <- df_tests %>% filter(n_groups > 2)
    
    final_results <- tibble() # default empty
    
    #----------------------------------- ANOVA ------------------------------------#
    
    if (nrow(df_for_ANOVA) > 0) {
      qPrint('running ANova')
      anova_results <- df_for_ANOVA %>%
        #  mutate(axis_name = as.factor(axis_column)) %>% # mistake?
        mutate(!!axis_column := as.factor(!!sym(axis_column))) %>% 
        group_by(across(all_of(test_across))) %>%
        summarise(
          anova_result = list({
            n_groups <- n_distinct(.data[[testing_col]])
            all_n <- table(.data[[testing_col]])
            valid <- sum(all_n > 1) >= 2   # at least 2 groups with ≥2 replicates
            
            if (n_groups >= 2 && valid) {
              fit <- tryCatch(
                aov(stats::as.formula(paste(value_column, "~", testing_col)), data = cur_data()),
                error = function(e) NULL
              )
              if (!is.null(fit)) {
                res <- broom::tidy(fit)
                if (all(is.na(res$statistic))) {
                  tibble(term = testing_col, df = NA_real_, statistic = 88, p.value = 88)
                } else {
                  res
                }
              } else {
                tibble(term = testing_col, df = NA_real_, statistic = 77, p.value = 77)
              }
            } else {
              tibble(term = testing_col, df = NA_real_, statistic = 77, p.value = 77)
            }
          }),
          .groups = "drop"
        ) %>%
        unnest(anova_result)
      anova_results
      
      qPrint('running Tukeys')
      tukey_results <- df_for_ANOVA %>%
        filter(!is.nan(.data[[value_column]])) %>%
        group_by(across(all_of(test_across))) %>%
        reframe(
          tukey_result = {
            fit <- tryCatch(
              aov(stats::as.formula(paste(value_column, "~", testing_col)), data = cur_data()),
              error = function(e) NULL
            )
            if (!is.null(fit)) {
              tukey_df <- as.data.frame(TukeyHSD(fit)[[1]])
              colnames(tukey_df) <- c("estimate", "conf.low", "conf.high", "adj.p.value")
              tukey_df$term <- testing_col
              tukey_df$comparison <- rownames(tukey_df)
              rownames(tukey_df) <- NULL
              tukey_df
            } else {
              tibble(
                estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, adj.p.value = NA_real_,
                term = testing_col, comparison = NA_character_
              )
            }
          }
        ) %>%
        unnest(tukey_result)
      
      
      df_residual=anova_results %>% 
        filter(term=='Residuals') %>%
        mutate(residual_df=df) %>% 
        select(any_of(test_across),residual_df) %>% 
        distinct()
      
      anova_join=anova_results%>%
        filter(!term=='Residuals') %>% 
        left_join(df_residual) %>% 
        mutate(anova.p.value = p.value) %>%
        mutate(anova_df=df)%>%
        select(any_of(test_across), anova.p.value,anova_df,residual_df)
      
      
      final_results <- tukey_results %>%
        left_join(anova_join, by = test_across) %>%
        mutate(comparison2=comparison) %>% 
        separate(comparison2, into = c("group2", "group1"), sep = "-") %>% 
        left_join(df_combs) %>% 
        mutate(adjust_method='Tukey',
               test_type = 'Tukey') %>% # adjusting for new ttest
        mutate(
          adj.significance = case_when(
            is.na(adj.p.value)      ~ NA_character_,
            adj.p.value < 0.001     ~ "***",
            adj.p.value < 0.01      ~ "**",
            adj.p.value < 0.05      ~ "*",
            TRUE                    ~ "ns"
          )
        )  %>% 
        mutate(sig_label = adj.significance) %>% 
        mutate(
          sig_label = ifelse(is.na(anova.p.value), adj.significance,
                             ifelse(anova.p.value > 1, "error", adj.significance)),
          p_value_label = ifelse(is.na(anova.p.value), as.character(round(adj.p.value,4)),
                                 ifelse(anova.p.value > 1, "error", as.character(round(adj.p.value,4))))
        ) %>% 
        mutate(p_value_label = formatC(adj.p.value, format = "e", digits = 2)) %>% 
        ungroup()
    }
    
    df_for_ttest <- df_tests %>% filter(n_groups < 3)
    
    #----------------------------------- ttest ------------------------------------#
    
    
    ## ============================================================
    ## T-test subsystem for summary_brackets()
    ## Input: df_for_ttest, df_combs, test_across, facet_columns,
    ##        value_column, log (logical)
    ## Output: ttest_results
    ## ============================================================
    
    ttest_results <- NULL
    
    if (nrow(df_for_ttest) > 0) {
      
      ## ------------------------------------------
      ## 1) Build shapiro_table from df_combs
      ##    (one row per comparison)
      ## ------------------------------------------
      comb_keys <- df_for_ttest %>%
        dplyr::select(dplyr::any_of(test_across)) %>%
        dplyr::distinct()
      
      shapiro_table <- df_combs %>%
        dplyr::semi_join(comb_keys, by = test_across) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          # raw + log vectors
          g1_raw = list(as.numeric(unlist(group1_values))),
          g2_raw = list(as.numeric(unlist(group2_values))),
          g1_log = list(log10(as.numeric(unlist(group1_values)) + 1)),
          g2_log = list(log10(as.numeric(unlist(group2_values)) + 1)),
          
          # Shapiro per group
          shapiro_raw_1 = run_single_shapiro(g1_raw),
          shapiro_raw_2 = run_single_shapiro(g2_raw),
          shapiro_log_1 = run_single_shapiro(g1_log),
          shapiro_log_2 = run_single_shapiro(g2_log),
          
          # Comparison-level Shapiro (max p, but NA if both NA)
          shapiro_raw = {
            vals <- c(shapiro_raw_1, shapiro_raw_2)
            if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
          },
          shapiro_log = {
            vals <- c(shapiro_log_1, shapiro_log_2)
            if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
          }
        ) %>%
        dplyr::ungroup()
      
      ## ------------------------------------------
      ## 2) Attach tests (Welch & Mann, raw + log)
      ## ------------------------------------------
      ttest_results <- shapiro_table %>%
        dplyr::select(
          dplyr::any_of(c(test_across,'axis_level')),
          group1, group2,
          group1_level, group2_level,
          group1_values, group2_values,
          group_testing_order, 
          #axis_level,
          g1_raw, g2_raw, g1_log, g2_log,
          shapiro_raw, shapiro_log
        ) %>%
        dplyr::distinct() %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          # set up list-columns for tests
          raw1 = list(g1_raw),
          raw2 = list(g2_raw),
          log1 = list(g1_log),
          log2 = list(g2_log),
          
          # Welch RAW
          welch_raw_full = list(run_welch_full(raw1, raw2)),
          welch_raw_val  = welch_raw_full$p,
          welch_raw_estimate = welch_raw_full$estimate,
          welch_raw_stat     = welch_raw_full$stat,
          welch_raw_conf_low  = welch_raw_full$conf.low,
          welch_raw_conf_high = welch_raw_full$conf.high,
          
          # Mann RAW
          mann_raw_val = run_mann(raw1, raw2),
          
          # Welch LOG
          welch_log_full = list(run_welch_full(log1, log2)),
          welch_log_val  = welch_log_full$p,
          welch_log_estimate = welch_log_full$estimate,
          welch_log_stat     = welch_log_full$stat,
          welch_log_conf_low  = welch_log_full$conf.low,
          welch_log_conf_high = welch_log_full$conf.high,
          
          # Mann LOG
          mann_log_val = run_mann(log1, log2)
        ) %>%
        dplyr::ungroup()
      
      ## ------------------------------------------
      ## 3) BH on ALL test types (per facet family)
      ## ------------------------------------------
      ttest_results <- ttest_results %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(facet_columns))) %>%
        dplyr::mutate(
          n_tests = dplyr::n(),
          
          welch_raw_BH_val = p.adjust(welch_raw_val, method = "BH"),
          mann_raw_BH_val  = p.adjust(mann_raw_val,  method = "BH"),
          welch_log_BH_val = p.adjust(welch_log_val, method = "BH"),
          mann_log_BH_val  = p.adjust(mann_log_val,  method = "BH")
        ) %>%
        dplyr::ungroup()
      
      ## ------------------------------------------
      ## 4) Significance on all value + BH columns
      ## ------------------------------------------
      ttest_results <- ttest_results %>%
        dplyr::mutate(
          welch_raw_sig = sig_labeler(welch_raw_val),
          mann_raw_sig  = sig_labeler(mann_raw_val),
          welch_log_sig = sig_labeler(welch_log_val),
          mann_log_sig  = sig_labeler(mann_log_val),
          
          welch_raw_BH_sig = sig_labeler(welch_raw_BH_val),
          mann_raw_BH_sig  = sig_labeler(mann_raw_BH_val),
          welch_log_BH_sig = sig_labeler(welch_log_BH_val),
          mann_log_BH_sig  = sig_labeler(mann_log_BH_val)
        )
      
      ## ------------------------------------------
      ## 5) Facet-level default test selection
      ##    (normality: NA = non-normal)
      ## ------------------------------------------
      facet_default <- ttest_results %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(facet_columns))) %>%
        dplyr::summarise(
          facet_shapiro_raw_ok = all(!is.na(shapiro_raw) & shapiro_raw > 0.05),
          facet_shapiro_log_ok = all(!is.na(shapiro_log) & shapiro_log > 0.05),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          default_test_label = dplyr::case_when(
            isTRUE(log)  & facet_shapiro_log_ok ~ "LW_BH",
            isTRUE(log)                         ~ "LM_BH",
            !isTRUE(log) & facet_shapiro_raw_ok ~ "W_BH",
            TRUE                                ~ "M_BH"
          )
        )
      
      ttest_results <- ttest_results %>%
        dplyr::left_join(facet_default, by = facet_columns)
      
      ## ------------------------------------------
      ## 6) Build default_val + default_sig
      ## ------------------------------------------
      ttest_results <- ttest_results %>%
        dplyr::mutate(
          selected_BH =
            dplyr::case_when(
              default_test_label == "W_BH"  ~ welch_raw_BH_val,
              default_test_label == "M_BH"  ~ mann_raw_BH_val,
              default_test_label == "LW_BH" ~ welch_log_BH_val,
              default_test_label == "LM_BH" ~ mann_log_BH_val
            ),
          selected_sig =
            dplyr::case_when(
              default_test_label == "W_BH"  ~ welch_raw_BH_sig,
              default_test_label == "M_BH"  ~ mann_raw_BH_sig,
              default_test_label == "LW_BH" ~ welch_log_BH_sig,
              default_test_label == "LM_BH" ~ mann_log_BH_sig
            ),
          
          default_val = paste0(
            default_test_label, "(", n_tests, ")=",
            formatC(selected_BH, format = "e", digits = 2)
          ),
          default_sig = paste0(
            default_test_label, "(", n_tests, ")=",
            selected_sig
          )
        )
      
      ## ------------------------------------------
      ## 7) Final column ordering
      ## ------------------------------------------
      ttest_results <- ttest_results %>%
        dplyr::select(
          dplyr::any_of(c(test_across,'axis_level')),
          group1, group2,
          group1_level, group2_level,
          group_testing_order, 
          #axis_level,
          
          # Defaults
          default_test_label, default_val, default_sig,
          
          # Shapiro
          shapiro_raw, shapiro_log,
          
          # Raw tests
          welch_raw_val, welch_raw_sig, welch_raw_BH_val, welch_raw_BH_sig,
          mann_raw_val,  mann_raw_sig,  mann_raw_BH_val,  mann_raw_BH_sig,
          
          # Log tests
          welch_log_val, welch_log_sig, welch_log_BH_val, welch_log_BH_sig,
          mann_log_val,  mann_log_sig,  mann_log_BH_val,  mann_log_BH_sig,
          
          n_tests,
          
          # Values
          group1_values, group2_values,
          
          # Welch stats (raw)
          welch_raw_estimate, welch_raw_stat,
          welch_raw_conf_low, welch_raw_conf_high,
          
          # Welch stats (log)
          welch_log_estimate, welch_log_stat,
          welch_log_conf_low, welch_log_conf_high,
          
          dplyr::everything()
        ) %>% 
        mutate(
          p_value_label = default_val,
          sig_label = default_sig
        )
    }
    
    
    qPrint('joining tests')
    
    test_results1=bind_rows(ttest_results,final_results)%>%
      left_join(df_segments%>%select(any_of(c(summary_group_cols,'seg_value','seg_min','seg_max','group1','group2'))))%>%
      mutate(
        label_position = (seg_min + seg_max) / 2,
        label_value = seg_value
      ) %>%
      
      select(-c(seg_min,seg_max,seg_value))%>%
      
      select(all_of(test_across), any_of(c('group1', 'group2','adj.p.value','adj.significance','adjust_method','p.value','significance',
                                           'estimate', 'statistic')),everything())%>%
      
      
      ungroup()
    names(test_results1)
    
    if( clear_ns_ANOVA){ # remove ns Tukuys?
      qPrint('clear_ns_ANOVA')
      df_ns_names=test_results1 %>% 
        filter(anova.p.value>=.05) %>% 
        select(any_of(test_across)) %>% 
        mutate(remove='yes') %>% 
        distinct()
      
      df_brackets_long=df_brackets_long %>% 
        left_join(df_ns_names) %>%
        filter(is.na(remove))
      
      names(test_results1)
      
      test_results <- test_results1 %>% 
        left_join(df_ns_names) %>%
        group_by(across(any_of(test_across))) %>% 
        mutate(
          label_position = ifelse(!is.na(remove), mean(label_position), label_position),
          label_value = ifelse(!is.na(remove), max(label_value), label_value),
          # change 20251219
          p_value_label = ifelse(!is.na(remove),as.character(round(anova.p.value, 2)),p_value_label),
          sig_label = ifelse(!is.na(remove), 'ANOVA ns', sig_label),
          # change 20251219
          #  p_value_label = ifelse(!is.na(remove) & clear_ns_message, '', p_value_label),
          sig_label = ifelse(!is.na(remove) & clear_ns_message, '', sig_label)
        ) %>%
        ungroup()
      
      df_summary
      names(df_brackets_long)
      names( test_results)
      
    }else{test_results=test_results1}
  }
  
  
  output <- list()
  # if ("summary" %in% return) output$summary <- df_summary
  # if ("segments" %in% return) output$segments <- df_brackets_long
  # if ("tests" %in% return) output$tests <- test_results
  output$summary <- df_summary
  output$segments <- df_brackets_long
  output$tests <- test_results
  return(output)
}
################################################################################
# end summary brackets
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
run_adonis <- function(formula, by = NULL, method = NULL) {
  news='2665 adonis';qPrint(news);qPrint(news);qPrint(news)
  # resolve 'by'
  if (is.null(by) && is.null(method)) {
    stop("run_adonis: either 'by' or 'method' must be provided")
  }
  if (is.null(by)) {
    by <- method
  }
  
  set.seed(123)
  
  formula_string <- paste("dist_matrix ~", formula)
  formula_str <- stats::as.formula(formula_string)
  print(formula_str)
  
  res <- adonis2(
    formula_str,
    data = df_meta,
    permutations = 999,
    by = by
  )
  
  as.data.frame(res) %>%
    rownames_to_column(var = "Term") %>%
    rename(
      Df = Df,
      SumOfSqs = SumOfSqs,
      R2 = R2,
      F_value = F,
      p_value = `Pr(>F)`
    ) %>%
    mutate(formula = formula, by = by)
}

################################################################################
run_permdisp <- function(dist_matrix, grouping_var) {
  news='2701 permdisp';qPrint(news);qPrint(news);qPrint(news)
  set.seed(123)
  # grouping_var: string, e.g. "color_group"
  
  # Compute betadisper object
  bd <- betadisper(dist_matrix, group = df_meta[[grouping_var]])
  
  # Run permutation test
  perm_disp <- permutest(bd, permutations = 999)
  
  # Extract table from result
  df <- as.data.frame(perm_disp$tab) %>%
    rownames_to_column(var = "Term") %>%
    mutate(grouping_var = grouping_var)
  
  return(df)
}


################################################################################
assign_palette <- function(vec) {
  # Example: print_palette_assignments(c("Clipper", "Field Control"))
  cat(sapply(seq_along(vec), function(i) {
    paste0("'", vec[i], "' = palette_common[", i, "]")
  }), sep = ",\n")
}
################################################################################
assign_default <- function(vec, default) {
  # Example: assign_default(c("A", "B", "C"), "#888888")
  unique_vals <- unique(vec)
  
  cat(sapply(unique_vals, function(v) {
    key <- if (is.numeric(v)) v else paste0("'", v, "'")
    val <- if (is.numeric(default)) default else paste0("'", default, "'")
    paste0(key, " = ", val)
  }), sep = ",\n")
}
################################################################################
build_palette <- function(color_vector) {
  # Example of use:
  # full_palette <- build_palette(c("blue", "red", "gold", "purple", "green", "gray"))
  
  # Get palette names and assign them to a list
  palettes <- list(
    blue = palette_blue,
    red = palette_red,
    gold = palette_gold,
    purple = palette_purple,
    green = palette_green,
    gray = palette_gray
  )
  
  # Number of shades (assuming all palettes have same length)
  n_shades <- length(palette_blue)
  
  # Create final palette by looping through shade levels
  out <- unlist(lapply(1:n_shades, function(i) {
    sapply(color_vector, function(col) palettes[[col]][i])
  }))
  
  return(out)
}
################################################################################
################################################################################
# Example of use:
# diagonal_palette(c("blue", "red", "gold", "purple", "green", "gray"))
diagonal_palette <- function(color_vector) {
  palettes <- list(
    blue = palette_blue,
    red = palette_red,
    gold = palette_gold,
    purple = palette_purple,
    green = palette_green,
    gray = palette_gray
  )
  
  n_colors <- length(color_vector)
  n_shades <- length(palette_blue)  # assuming consistent length
  
  out <- character(0)
  
  for (i in seq_len(n_shades)) {
    for (j in seq_along(color_vector)) {
      shade_level <- (7-(i + j)) %% n_shades + 1  # wrap around
      color_name <- color_vector[j]
      out <- c(out, palettes[[color_name]][shade_level])
    }
  }
  
  return(out)
}
################################################################################
# filter data set by group prevalence and relative abundan
# dfx_matrix includes df_matrix + 'sample_name'+grouping columns
################################################################################
filter_by_prevalence_and_abundance <- function(dfx_matrix, grouping_columns, minimum_prevalence, minimum_abundance) {
  #testing
  #dfx_matrix=dfx_matrix;grouping_columns=grouping_columns;minimum_prevalence=.9; minimum_abundance=.0001
  
  # Example usage:
  # filtered_df <- filter_by_prevalence_and_abundance(dfx_matrix,
  #                  grouping_columns = c("sample_type", "plant_type"),
  #                  minimum_prevalence = 0.9, minimum_abundance = 0.00001)
  
  # Split data by unique group combination
  
  grouped_list <- split(dfx_matrix, dfx_matrix[, grouping_columns], drop = TRUE)
  
  # Filter each group
  filtered_list <- lapply(grouped_list, function(group_df) {
    count_data <- group_df[, !(colnames(group_df) %in% c('sample_name',grouping_columns)), drop = FALSE]
    
    # 1. Prevalence filter: zero taxa not meeting threshold
    taxa_prevalence <- colSums(count_data > 0) / nrow(count_data)
    low_prev_taxa <- names(taxa_prevalence)[taxa_prevalence < minimum_prevalence]
    count_data[, low_prev_taxa] <- 0
    
    # 2. Relative abundance filter: zero values below threshold
    row_totals <- rowSums(count_data)
    rel_abundance <- sweep(count_data, 1, row_totals, FUN = "/")
    count_data[rel_abundance < minimum_abundance] <- 0
    
    # Return count data with group metadata
    cbind(group_df[, c('sample_name',grouping_columns), drop = FALSE], count_data)
  })
  
  # Combine all groups back together
  combined_filtered_df <- do.call(rbind, filtered_list) 
  rownames(combined_filtered_df) <- NULL
  # combined_filtered_df<-combined_filtered_df%>% 
  #   column_to_rownames(var='sample_name')
  # combined_filtered_df$sample_name <- rownames(combined_filtered_df)
  #  rownames(combined_filtered_df) <- NULL
  
  return(combined_filtered_df)
}
################################################################################
################################################################################
################################################################################
combine_keep_top10 <- function(df = df_matrix, averaging_columns = c("location"), name_prefix = "Consolidated") {
  # Example usage:
  # combine_keep_top10(df_matrix, averaging_columns = c("location", "sample_type"))
  
  df_top10 <- df %>%
    # explode metadata from sample names
    xPlode_sample_name() %>%
    pivot_longer(
      cols = all_of(colnames(df)),
      names_to = "feature",
      values_to = "counts"
    ) %>%
    # per-sample totals
    group_by(sample_name) %>%
    mutate(sample_counts = sum(counts, na.rm = TRUE)) %>%
    ungroup() %>%
    # relative abundance
    mutate(rel_abun = counts / sample_counts) %>%
    # group-level mean (e.g. mean per location)
    group_by(across(all_of(averaging_columns)), feature) %>%
    mutate(group_average = mean(rel_abun, na.rm = TRUE)) %>%
    ungroup() %>%
    # overall mean (each group contributes equally)
    group_by(feature) %>%
    mutate(feature_average = mean(group_average, na.rm = TRUE)) %>%
    ungroup() %>% 
    distinct(feature,feature_average) %>% 
    arrange(desc(feature_average)) %>% 
    slice_head(n = 10)
  
  less_than <- min(round(df_top10$feature_average / 0.01, 1), na.rm = TRUE)
  consolidated_name <- paste0(name_prefix, " < ", less_than, "%")
  
  # Get top feature names
  top_features <- df_top10$feature
  
  # Split and recombine counts
  feature_names <- colnames(df)
  top_df <- df[, feature_names %in% top_features, drop = FALSE]
  other_df <- df[, !(feature_names %in% top_features), drop = FALSE]
  
  # Sum all low-abundance features into one column
  consolidated_counts <- rowSums(other_df, na.rm = TRUE)
  
  # Combine top features + consolidated column
  df_out <- cbind(top_df, consolidated_counts)
  colnames(df_out)[ncol(df_out)] <- consolidated_name
  
  return(df_out)
  
}


################################################################################
################################################################################
combine_keep_topN <- function(
    df = df_matrix,
    averaging_columns = c("location"),
    name_prefix = "Consolidated",
    top_n = 10
) {
  
  df_topN <- df %>%
    xPlode_sample_name() %>%
    pivot_longer(
      cols = all_of(colnames(df)),
      names_to = "feature",
      values_to = "counts"
    ) %>%
    group_by(sample_name) %>%
    mutate(sample_counts = sum(counts, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(rel_abun = counts / sample_counts) %>%
    group_by(across(all_of(averaging_columns)), feature) %>%
    mutate(group_average = mean(rel_abun, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(feature) %>%
    mutate(feature_average = mean(group_average, na.rm = TRUE)) %>%
    ungroup() %>%
    distinct(feature, feature_average) %>%
    arrange(desc(feature_average)) %>%
    slice_head(n = top_n)
  
  less_than <- min(round(df_topN$feature_average / 0.01, 1), na.rm = TRUE)
  consolidated_name <- paste0(name_prefix, " < ", less_than, "%")
  
  top_features <- df_topN$feature
  
  feature_names <- colnames(df)
  top_df <- df[, feature_names %in% top_features, drop = FALSE]
  other_df <- df[, !(feature_names %in% top_features), drop = FALSE]
  
  consolidated_counts <- rowSums(other_df, na.rm = TRUE)
  
  df_out <- cbind(top_df, consolidated_counts)
  colnames(df_out)[ncol(df_out)] <- consolidated_name
  
  df_out
}

################################################################################
combine_low_abundance <- function(df = df_matrix, threshold = 0.01, name = 'Low_Abundance') {
  # Example usage:
  # combined_df <- combine_low_abundance(df_matrix, threshold = 0.01, name = "Other")
  
  df_rel <- df / rowSums(df, na.rm = TRUE)
  feature_means <- colMeans(df_rel, na.rm = TRUE)
  
  low_abundance_features <- names(feature_means[feature_means < threshold])
  high_abundance_features <- names(feature_means[feature_means >= threshold])
  
  df_return <- df_rel[, high_abundance_features, drop = FALSE]
  
  # Taxonomic name format for combined low-abundance group
  #combined_name <- paste0(c("k_", "p_", "c_", "o_", "f_", "g_"), name, collapse = ";")
  combined_name=name
  df_return[[combined_name]] <- rowSums(df_rel[, low_abundance_features, drop = FALSE], na.rm = TRUE)
  
  return(df_return)
}
################################################################################
combine_other <- function(df, patterns, starts_with = FALSE, ignore_case = TRUE, fun = sum) {
  # Example usage:
  # combine_other(df, patterns = c("Gene","score"), starts_with = FALSE)
  
  # build regex depending on starts_with
  regex_patterns <- if (starts_with) {
    paste0("^", patterns)   # anchor at beginning
  } else {
    patterns                 # anywhere in name
  }
  
  keep_idx <- grep(paste(regex_patterns, collapse = "|"), 
                   names(df), 
                   ignore.case = ignore_case)
  keep_cols <- names(df)[keep_idx]
  
  other_cols <- setdiff(names(df), keep_cols)
  
  df %>%
    mutate(
      Other = apply(across(all_of(other_cols)), 1, function(x) fun(x, na.rm = TRUE))
    ) %>%
    select(all_of(keep_cols), Other)
}


################################################################################
library(dplyr)

qSummary <- function(df, cols, groups = NULL) {
  
  cols <- rlang::syms(cols) # turn strings into symbols
  
  if (!is.null(groups)) {
    group_syms <- rlang::syms(groups)
    df <- df %>% group_by(!!!group_syms)
  }
  
  if (length(cols) == 1) {
    col <- cols[[1]]
    out <- df %>%
      summarise(
        mean = mean(!!col, na.rm = TRUE),
        sd   = sd(!!col, na.rm = TRUE),
        min  = min(!!col, na.rm = TRUE),
        max  = max(!!col, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    out <- df %>%
      summarise(
        across(
          .cols = !!!cols,
          .fns = list(
            mean = ~mean(.x, na.rm = TRUE),
            sd   = ~sd(.x, na.rm = TRUE),
            min  = ~min(.x, na.rm = TRUE),
            max  = ~max(.x, na.rm = TRUE)
          ),
          .names = "{.col}_{.fn}"
        ),
        .groups = "drop"
      )
  }
  
  return(out)
}


################################################################################
QC_range <- function(df, value_vec, grouping_vec, output_prefix = "QC_summary_") {
  # Example usage:
  # value_vec=c('ideal_score','total_counts')
  # grouping_vec=c('axis_name','facet_name')
  # QC_range(df=df2,value_vec,grouping_vec)
  
  # QC_range(df, value_vec = c("pH","moisture"), grouping_vec = c("Location","Treatment"))
  
  library(dplyr)
  library(gridExtra)
  library(patchwork)
  library(cowplot)
  
  value_vec    <- intersect(value_vec, colnames(df))
  grouping_vec <- intersect(grouping_vec, colnames(df))
  
  plots <- list()
  total_rows = 0
  
  for (val in value_vec) {
    summary_tbl <- df %>%
      group_by(across(all_of(grouping_vec))) %>%
      summarise(
        n    = sum(!is.na(.data[[val]])),
        min  = min(.data[[val]], na.rm = TRUE),
        max  = max(.data[[val]], na.rm = TRUE),
        mean = mean(.data[[val]], na.rm = TRUE),
        .groups = "drop"
      )
    
    # title grob rotated vertically
    title <- textGrob(paste("QC summary for", val),
                      gp = gpar(fontsize = 12, fontface = "bold"),
                      rot = 90,
                      just = "centre")
    title_plot <- ggdraw() + draw_grob(title)
    
    # table grob
    tbl_grob <- tableGrob(summary_tbl, rows = NULL, theme = ttheme_default(base_size = 10))
    tbl_plot <- ggdraw() + draw_grob(tbl_grob)
    
    # combine title (left) + table (right) side by side
    combined_side <- title_plot + tbl_plot + plot_layout(widths = c(1, 6))
    
    plots[[val]] <- combined_side
    total_rows = total_rows + nrow(summary_tbl)
  }
  
  # stack all tables vertically
  combined <- wrap_plots(plots, ncol = 1)
  
  grouping_name <- paste(grouping_vec, collapse = "_by_")
  filename <- paste0(output_prefix, grouping_name, ".pdf")
  
  ggsave(
    filename = filename,
    plot     = combined,
    width    = 8.5,
    height   = min(total_rows * .3, 49)
  )
}
################################################################################
run_adonis <- function(formula, method = NULL, strata_var = NULL) {
  news='3094 Adonis';qPrint(news);qPrint(news);qPrint(news)
  # Example:
  # run_adonis("sample_time", method = "margin", strata_var = "location")
  
  set.seed(123)
  
  formula_string <- paste("dist_matrix ~", formula)
  formula_str <- as.formula(formula_string)
  print(formula_str)
  
  strata_vec <- NULL
  if (!is.null(strata_var)) {
    if (!strata_var %in% colnames(df_meta)) stop(paste("strata_var", strata_var, "not found in df_meta"))
    strata_vec <- df_meta[[strata_var]]
  }
  
  res <- adonis2(
    formula_str,
    data = df_meta,
    permutations = 1000,
    by = method,
    strata = strata_vec
  )
  
  if (is.null(method)) method <- NA
  
  as.data.frame(res) %>%
    rownames_to_column(var = "Term") %>%
    rename(
      SumOfSqs = SumOfSqs,
      R2 = R2,
      F_value = F,
      p_value = `Pr(>F)`
    ) %>%
    mutate(
      formula = formula,
      method = method,
      strata_var = ifelse(is.null(strata_var), NA, strata_var),
      test = "permANOVA",
      Df_residual = Df[Term == "Residual"],
      Df_total = Df[Term == "Total"]
    ) %>%
    mutate(Variable=Term) %>% 
    select(method, strata_var, formula, Variable, everything())
}


################################################################################

# Function to run PERMDISP and return tidy result
run_permdisp <- function(dist_matrix, grouping_var) {
  news='3145 permdisp';qPrint(news);qPrint(news);qPrint(news)
  # Example usage: run_permdisp(dist_matrix, "Treatment")
  
  set.seed(123)
  
  # Compute betadisper object using median
  bd <- betadisper(dist_matrix, group = df_meta[[grouping_var]], type = "median")
  
  # Run permutation test
  perm_disp <- permutest(bd, permutations = 999)
  
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
    )%>% 
    mutate(Df_residual = Df[Term == "Residuals"]) %>% 
    mutate(Variable=Term) %>% 
    select(everything())
  
  return(result_df)
}
################################################################################
################################################################################
add_tags <- function(plot,
                     tags,
                     x_values = 0.5,
                     y_values = 0.5,
                     color_values = "black",
                     size_values = 8,
                     fontface = "bold",
                     shape_values = 0,
                     shape_size_values = 2,
                     shape_color_values = "gray70",
                     alpha_values = 1,
                     grid = TRUE,
                     grid_major_by = 0.1,
                     grid_minor_by = 0.02,
                     grid_label_positions = c(0.3, 0.6)) {
  # Example:
  # add_tags(final_layout, tags = c("A","B"),
  #          x_values = c(0.2,0.8), y_values = 0.9)
  
  library(ggplot2)
  library(patchwork)
  
  n <- length(tags)
  
  # Invisible placeholder if shape_values == 0
  add_points='yes'
  if (length(shape_values) == 1 && shape_values == 0) {
    shape_values=1;alpha_values=0
  }
  
  # Helper for recycling vectors
  recycle_to <- function(x, n) rep(x, length.out = n)
  x_values           <- recycle_to(x_values, n)
  y_values           <- recycle_to(y_values, n)
  color_values       <- recycle_to(color_values, n)
  size_values        <- recycle_to(size_values, n)
  shape_values       <- recycle_to(shape_values, n)
  shape_size_values  <- recycle_to(shape_size_values, n)
  shape_color_values <- recycle_to(shape_color_values, n)
  alpha_values       <- recycle_to(alpha_values, n)
  
  # Tag and optional points
  
  tag_layer <- ggplot()+
    geom_point(aes(x = x_values, y = y_values),
               shape = shape_values,
               alpha = alpha_values,
               size = shape_size_values,
               color = ifelse(shape_values %in% 21:25,color_values, shape_color_values),
               fill = ifelse(shape_values %in% 21:25, shape_color_values, NA))+
    geom_text(aes(x = x_values, y = y_values, label = tags),
              color = color_values,
              size = size_values,
              fontface = fontface,
              vjust = 0.5, hjust = 0.5) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
  
  # Optional grid overlay
  if (grid) {
    overlay=tag_layer +
      geom_vline(xintercept = seq(0, 1, by = grid_minor_by/10),
                 color = "gray90", linewidth = 0.3, alpha = 0.7) +
      geom_vline(xintercept = seq(0, 1, by = grid_minor_by),
                 color = "gray70", linewidth = 0.3, alpha = 0.7) +
      geom_vline(xintercept = seq(0, 1, by = grid_major_by),
                 color = "gray30", linewidth = 0.6, alpha = 0.7) +
      
      geom_vline(xintercept = .5,color = "navy", linewidth = 0.6, alpha = 1,linetype='dashed') +
      
      geom_hline(yintercept = seq(0, 1, by = grid_minor_by/10),
                 color = "gray90", linewidth = 0.3, alpha = 0.7) +
      
      geom_hline(yintercept = seq(0, 1, by = grid_minor_by),
                 color = "gray70", linewidth = 0.3, alpha = 0.7) +
      geom_hline(yintercept = seq(0, 1, by = grid_major_by),
                 color = "gray30", linewidth = 0.6, alpha = 0.7) +
      geom_hline(yintercept = .5,color = "navy", linewidth = 0.6, alpha = 1,linetype='dashed') +
      
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_void()
    
  } else {
    overlay <- tag_layer
  }
  
  final_plot <- wrap_elements(full = plot) +
    inset_element(overlay,
                  left = 0, bottom = 0,
                  right = 1, top = 1,
                  align_to = "full")
  
  return(final_plot)
}







################################################################################

################################################################################

################################################################################
################################################################################

################################################################################
################################################################################

################################################################################
################################################################################

################################################################################
################################################################################
qPrint('# adding stuff for scripts')
################################################################################
qPrint('# palettes')
################################################################################

color_vector=c('blue','red','gold','purple','green','gray')

# Define color groups
palette_gray <- (c('#BBBBBB', '#999999', '#777777', '#555555','#333333'))
palette_gray <- (c('#E0E0E0', '#C0C0C0', '#A0A0A0', '#606060','#333333'))
palette_green <- (c('#D2EDB2', '#89C189', '#60A740', '#107810', '#004B00'))
palette_gold  <- (c('#FFDA65', '#FFCB25', '#DAA520', '#A37B18', '#7D6220'))
palette_purple <- (c('#D2AEFA', '#A267E8', '#7A1FD2', '#5B17A1', '#3C1070'))
palette_red <- (c('#FFCCCC', '#E47C7C', '#C94141', '#993333', '#6A2B2B'))
palette_blue <- (c('#C8D6FF', '#92AEFF', '#5C86FF', '#3A62BE', '#20407F'))

df_palette <- tibble(
  palette = c(
    rep("gray",5),rep("green",5),rep("gold",5),rep("purple",5),rep("red",5),rep("blue",5)
  ),
  shade = rep(1:5, times = 6),
  fill = c(palette_gray,palette_green,palette_gold,palette_purple,palette_red,palette_blue
  )
) %>%
  mutate(
    palette = factor(palette,levels = c("gray","green","gold","purple","red","blue")),
    label = paste0(palette, "[", shade, "]")
  )

palette_plot <- df_palette %>%
  ggplot(aes(x = shade, y = palette)) +
  geom_tile(aes(fill = fill), color = "black") +
  geom_text(aes(label = label), size = 3) +
  scale_fill_identity() +
  coord_equal() +
  theme_void()
palette_plot

df_palette$label <- paste0(df_palette$palette, "[", df_palette$shade, "]")
# palette_common <- build_palette(color_vector)
# palette_common <- diagonal_palette(color_vector)
palette_common_names=c('blue[5]','red[4]','gold[3]','purple[2]','green[3]','gray[4]',
                       'blue[3]','red[2]','gold[5]','purple[4]','green[5]','gray[2]',
                       'blue[1]','red[3]','gold[4]','purple[5]','green[4]','gray[3]',
                       'blue[2]','red[1]','gold[2]','purple[1]','green[2]','gray[1]',
                       'blue[4]','red[5]','gold[1]','purple[3]','green[1]','gray[5]')

palette_common=c(palette_blue[5],palette_red[4],palette_gold[3],palette_purple[2],palette_green[3],palette_gray[4],
                 palette_blue[3],palette_red[2],palette_gold[5],palette_purple[4],palette_green[5],palette_gray[2],
                 palette_blue[1],palette_red[3],palette_gold[4],palette_purple[5],palette_green[4],palette_gray[3],
                 palette_blue[2],palette_red[1],palette_gold[2],palette_purple[1],palette_green[2],palette_gray[1],
                 palette_blue[4],palette_red[5],palette_gold[1],palette_purple[3],palette_green[1],palette_gray[5])
df_palette <- data.frame(
  label = palette_common_names,
  fill = palette_common,
  x = rep(1:6, each = 5),
  y = rep(5:1, times= 6)
)

common_palettte_plot=ggplot(df_palette, aes(x, y)) +
  geom_tile(aes(fill = fill), color = "black") +
  geom_text(aes(label = label), size = 3) +
  scale_fill_identity() +
  coord_equal() +
  theme_void()

common_palettte_plot

palette(palette_common)

palette_common




palette_gradient=c( '#20407F','#5B17A1','#C94141','#DAA520','#E2E21C')
#F0D88F
palette_gradient <- c(
  "#1B2ACC",  # Vivid dark blue
  "#9137CB",  # Medium bright purple
  "#D62728",  # Bright red (colorblind-friendly)
  "#FF7F0E",  # Bright orange (colorblind-friendly)
  "#FFD700"   # Gold (bright, more colorblind-safe than pale yellow)
)
palette_gradient2=palette_gradient[c(1,3)]

palette_helper_color=c(          
  'gray1'=palette_gray[1],
  'gray1'=palette_gray[2],
  'gray3'=palette_gray[3],
  'ns'=palette_gray[2],
  'sig'=palette_red[3],
  '*'=palette_red[3],
  '**'=palette_red[3],
  '***'=palette_red[3],
  '****'=palette_red[3],
  "Archaea" = palette_red[3],
  "Bacteria" = palette_blue[5],
  "Eukaryote"=palette_green[5],
  "unknown"=palette_gray[3],
  is.na=palette_gray[1],
  'white'='white',
  'black'='black',
  'default'='white')

palette_helper_label=c(          
  '*'='*',
  '**'='**',
  '***'='***',
  '****'='****',
  'ns'='ns',
  'sig'='sig',
  
  'white'='white',
  'gray1'='gray1',
  'gray2'='gray2',
  'gray3'='gray3',
  'black'='black',
  'default'='white')

# Example: shape_pal_default <- shape_pal_default
palette_shape_fill <- c(
  21, 22, 23, 24, 25,   # fillable first
  16, 17, 15, 3, 4, 8,  # common non-fillable
  7, 9, 10, 12, 13, 14, # less common
  0, 1, 2, 5, 6, 11     # outline shapes
)


# +scale_shape_manual(values = palette_shape_fill)

################################################################################
qPrint('# variables')
################################################################################
current_date <- format(Sys.Date(), "%Y%m%d")

################################################################################
# 2026-01-27
set_palette_common <- function() {
  options(
    ggplot2.discrete.colour = palette_common,
    ggplot2.discrete.fill   = palette_common
  )
}

# Example usage:
set_palette_common()

################################################################################
set_output <- function(subfolder = NULL) {
  # Example: set_outputs(script_title)
  
  base_data <- "output_data"
  base_plot <- "output_plot"
  
  output_data <- if (is.null(subfolder)) {
    paste0(base_data, "/")
  } else {
    paste0(base_data, "/", subfolder, "/")
  }
  
  output_plot <- if (is.null(subfolder)) {
    paste0(base_plot, "/")
  } else {
    paste0(base_plot, "/", subfolder, "/")
  }
  
  dir.create(output_data, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_plot, recursive = TRUE, showWarnings = FALSE)
  
  assign("output_data", output_data, envir = parent.frame())
  assign("output_plot", output_plot, envir = parent.frame())
}


################################################################################
################################################################################
# stat_count(
#   aes(x=bar_color,
#       
#       label = after_stat(paste0("n = ", count))
#   ),
#   inherit.aes = FALSE,
#   geom = "text",
#   vjust = 0,
#   fontface = "bold",
#   size = 3.5
# )
