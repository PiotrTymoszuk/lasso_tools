# This toolbox script contains functions for LASSO (or other model type provided by glmnet package)
# modeling taking a data frame as input data and functions used for selection of the optimal variable set
# based on bootstrap subsets of the given data

# tools -----

  require(tidyverse)
  require(glmnet)
  require(furrr)

# globals -----

  lasso_globs <- list()

  lasso_globs$common_text <- element_text(size = 8, 
                                          face = 'plain', 
                                          color = 'black')
  
  lasso_globs$common_margin <- ggplot2::margin(t = 5, 
                                               l = 5, 
                                               r = 5, 
                                               unit = 'mm')
  
  lasso_globs$common_theme <- theme_classic() + theme(axis.text = lasso_globs$common_text, 
                                                      axis.title = lasso_globs$common_text, 
                                                      plot.title = element_text(size = 10, 
                                                                                face = 'bold'), 
                                                      plot.subtitle = lasso_globs$common_text, 
                                                      plot.tag = element_text(size = 8, 
                                                                              face = 'plain', 
                                                                              color = 'black', 
                                                                              hjust = 0, 
                                                                              vjust = 1), 
                                                      plot.tag.position = 'bottom', 
                                                      legend.text = lasso_globs$common_text, 
                                                      legend.title = lasso_globs$common_text, 
                                                      strip.text = lasso_globs$common_text,
                                                      plot.margin = lasso_globs$common_margin)
  
# helper functions -----

  make_splits <- function(inp_tbl, nrow_train, n_splits = 100) {
    
    ## creates a series of train/test data sets out of the given table
    ## the size of the training data set can be specified by the user
    ## for optimal performance and compatibility with the knn prediction
    ## and accuracy testing functions, a data frame with specified row names
    ## is required
    
    if(is.null(rownames(inp_tbl))) {
      
      stop('The function requires a data frame with specified rownames')
      
    }
    
    tbl_len <- nrow(inp_tbl)
    
    split_ids <- paste('split', 1:n_splits, sep = '_') %>% 
      map(function(x) sample(1:tbl_len, nrow_train)) %>% 
      set_names(paste('split', 1:n_splits, sep = '_'))
    
    split_lst <- split_ids %>% 
      map(function(x) list(train = inp_tbl[x, ], 
                           test = inp_tbl[!(1:tbl_len) %in% x, ]))
    
    return(split_lst)
    
}

  count_occurrence <- function(branched_lst, element) {
    
    ## for the given list of vectors with non-repeating elements,
    ## the occurrence of the given unique vector element in the list
    ## is returned
    
    counting_res <- branched_lst %>% 
      map(function(x) element %in% x) %>% 
      reduce(sum)
    
    counting_res <- tibble(element = element, 
                           n = counting_res)
    
    return(counting_res)
    
}

  make_occurr_plot <- function(count_tbl, cutoff = 0.95, 
                               plot_title = NULL, plot_subtitle = NULL, plot_tag = NULL) {
    
    ## draws a bar plot with the occurrence of the given element
    
    count_plot <- count_tbl %>% 
      mutate(selected = ifelse(frac_models >= cutoff, 
                               'yes', 
                               'no')) %>% 
      ggplot(aes(x = frac_models, 
                 y = reorder(co_variate, frac_models), 
                 fill = selected)) + 
      geom_bar(stat = 'identity', 
               color = 'black') + 
      geom_vline(xintercept = cutoff, 
                 linetype = 'dashed') + 
      geom_text(aes(label = signif(frac_models, 2)), 
                size = 2.75, 
                hjust = 0, 
                nudge_x = 0.025) + 
      scale_fill_manual(values = c('no' = 'gray80', 
                                   'yes' = 'coral3'), 
                        name = 'Included\nin the best set') + 
      lasso_globs$common_theme + 
      theme(axis.title.y = element_blank()) +
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           x = 'Fraction of the bootstraps')
        
    
    return(count_plot)
    
  } 
  
  point_plot_ <- function(data, x_var, y_var, x_lab = x_var, y_lab = y_var, 
                          plot_title = NULL, smooth = T, silent = T, ...) {
    
    ## draws a simple point plot for diagnostic purposes, takes the output of get_qc_tbl() as data argument
    ## color-codes model missfits
    
    ## table for plotting 
    
    data <- data %>% 
      mutate(misslab = ifelse(.candidate_missfit == 'yes',
                              .rownames, 
                              NA))
    
    ## fill colors
    
    fill_colors <- c(no = 'cornflowerblue', 
                     yes = 'firebrick4')
    
    ## point plot
    
    point_plot <- data %>% 
      ggplot(aes(x = .data[[x_var]], 
                 y = .data[[y_var]], 
                 fill = .candidate_missfit)) + 
      geom_point(size = 2, 
                 shape = 21) + 
      geom_text_repel(aes(label = misslab), 
                      show.legend = F) + 
      scale_fill_manual(values = fill_colors, 
                        name = 'Candidate missfit') + 
      labs(x = x_lab, 
           y = y_lab, 
           title = plot_title)
    
    if(smooth) {
      
      if(silent) {
        
        suppressWarnings(point_plot <- point_plot + 
                           geom_smooth(show.legend = F, 
                                       color = 'black', 
                                       fill = 'dodgerblue2', ...))
        
      } else {
        
        point_plot <- point_plot + 
          geom_smooth(show.legend = F, 
                      color = 'black', 
                      fill = 'dodgerblue2', ...)
        
      }
      
    }
    
    return(point_plot)
    
  }
  
  calc_expected_ <- function(inp_data, observed) {
    
    ## calculates expected normal distribution of a variable observed
    ## credits to: https://stackoverflow.com/questions/43217104/coloring-points-in-a-geom-qq-plot
    
    
    inp_data <- inp_data[order(inp_data[[observed]]), ] %>% 
      mutate(.expect.norm = qnorm(ppoints(nrow(.))))
    
    return(inp_data)
    
  }
  
# glmnet modeling wrapper for data frame input data -----
  
  model_cv_lasso <- function(inp_tbl, response, variables = names(inp_tbl)[names(inp_tbl) != response], 
                             family = 'poisson', alpha = 1, ...) {
    
    ## a customized wrapped around cv.glmnet() for handling data frames
    ## categorical covariates are turned into 'dummy' variables by model.matrix()
    
    
    ## model components
    
    y <- inp_tbl[[response]]
    
    x <- model.matrix(~ ., 
                      inp_tbl[, variables])
    
    ## cross-validated model selection
    
    cv_model <- cv.glmnet(x = x, 
                          y = y, 
                          family = family, 
                          alpha = alpha, ...)
    
    return(cv_model)    
    
  }
  
# Prediction and QC functions handling the data frame input ------  
  
  get_qc_tbl_lasso <- function(cv_object, new_data_tbl, 
                               response, 
                               variables = names(inp_tbl)[names(inp_tbl) != response], 
                               lambda_crit = 'lambda.1se', 
                               type = 'response', ...) {
    
    ## a customized wrapper around predict.glmnet that handles data frames
    ## returns predictions, true values and residuals
    
    x <- model.matrix(~ ., 
                      new_data_tbl[, variables])
    
    true_vals <- new_data_tbl[[response]]
    
    predictions <- predict(object = cv_object, 
                           newx = x, 
                           s = lambda_crit, 
                           type = type)
    
    predictions <- tibble(.rownames = rownames(predictions), 
                          .fitted = predictions[, 1], 
                          y = true_vals) %>% 
      mutate(.resid = y - .fitted, 
             .sq.resid = .resid ^ 2, 
             .std.resid = scale(.resid)[, 1], 
             .sq.std.resid = .std.resid ^ 2, 
             .candidate_missfit = ifelse(abs(.std.resid) > qnorm(0.975), 
                                         'yes', 
                                         'no'))
    
    predictions <- calc_expected_(predictions, '.std.resid')
    
    return(predictions)    
    
  } 
  
  get_qc_plots_lasso <- function(cv_object, new_data_tbl, 
                                 response, 
                                 variables = names(inp_tbl)[names(inp_tbl) != response], 
                                 lambda_crit = 'lambda.1se', ...) {
    
    ## draws standard model qc plots with missfits labeled by obs. numbers
    ## ... are additional arguments passed to get_qc_tbl_lasso()
    
    ## QC table
    
    qc_tbl <- get_qc_tbl_lasso(cv_object = cv_object, 
                               new_data_tbl = new_data_tbl, 
                               response = response, 
                               variables = variables, 
                               lambda_crit = lambda_crit, 
                               type = 'response')
    
    ## QC plots
    
    qc_plotting_lst <- list(x_var = c('.fitted', '.fitted', '.fitted', '.expect.norm'), 
                            y_var = c('.resid', '.std.resid', '.sq.std.resid', '.std.resid'), 
                            plot_title = c('Residuals vs. fitted', 
                                           'Standardized residuals vs. fitted', 
                                           'Sqared residuals vs. fitted', 
                                           'QQ standardized residuals vs expected normal'),
                            method = c('loess', 'loess', 'loess', 'lm'), 
                            smooth = c(T, T, T, T))
    
    qc_plots <- qc_plotting_lst %>% 
      pmap(point_plot_, 
           data = qc_tbl) %>% 
      set_names(c('resid_fitted', 
                  'std.resid_fitted', 
                  'sq.resid_fitted', 
                  'qq.std.resid'))
    
    return(qc_plots)
    
  }
  
# variable selection functions ------
  
  select_var_set_ <- function(boot_tbls, response, 
                              variables = names(boot_tbls[[1]])[names(boot_tbls[[1]]) != response], 
                              family = 'poisson', alpha = 1, lambda_crit = 'lambda.1se', 
                              frac_cutoff = 0.95, .parallel = F, ...) {
    
    ## selects the optimal co-variate pool from the initial set ('variables')
    ## and the given bootstrap subsets of a table unsing the glmnet() lasso technique. 
    ## Lambda_crit defines a criterion for the variable set selection in each bootstrap
    ## frac_cutoff defines the minimum fraction of the bootstrap model where the given
    ## co-variate was selected as a criterion for the inclusion in the final co-variate dataset
    ## Values: the vectors of non-zero co-efficient in the bootstrap models (co-eff selection by lambda_crit), 
    ## occurrence table with the counts of each unique co-variate within the non-zero co-effcicient set
    ## occurrence plot and the best variable set selected with the frac_cutoff paramater
    
    start_time <- Sys.time()
    
    message(paste('Looking for the optimal variable set with '), 
            length(boot_tbls), 
            ' bootstrap tables')
    
    ## serial testing
    
    if(.parallel) {
      
      require(furrr)
      
      plan('multisession')
      
      res_models <- boot_tbls %>% 
        future_map(model_cv_lasso, 
                   response = response, 
                   variables = variables, 
                   family = family, 
                   alpha = alpha, 
                   .options = furrr_options(seed = T), 
                   ...)
      
      plan('sequential')
      
    } else {
      
      res_models <- boot_tbls %>% 
        map(model_cv_lasso, 
            response = response, 
            variables = variables, 
            family = family, 
            alpha = alpha, ...)
      
    }
    
    ## selection of non-zero coefficients
    ## coefficent names are extracted from the results 
    ## via a regex made of the input variable vector
    
    coef_lst <- res_models %>% 
      map(coef, 
          s = lambda_crit) %>% 
      map(as.matrix) %>% 
      map(as.data.frame) %>% 
      map(function(x) filter(x, x[[1]] != 0)) %>% 
      map(rownames)
    
    extr_regex <- paste(variables, 
                        collapse = '|')
    
    coef_lst <- coef_lst %>% 
      map(stri_extract, 
          regex = extr_regex) %>%
      map(function(x) x[!is.na(x)]) %>% 
      map(unique)
    
    ## counting the occurrence of each unique co-variate in the bootstrap models
    ## calculating the fraction of all models, where ach of the co-variates was selected
    
    all_vars <- coef_lst %>% 
      reduce(c) %>% 
      unique
    
    occurrence_counts <- all_vars %>% 
      map_dfr(count_occurrence, 
              branched_lst = coef_lst) %>% 
      mutate(n_total = length(coef_lst), 
             frac_models = n/n_total) %>% 
      set_names(c('co_variate', 
                  'n', 
                  'n_total', 
                  'frac_models')) %>% 
      arrange( - frac_models)
    
    ## selection of the best set of the co-variates
    
    best_cov_set <- occurrence_counts %>% 
      filter(frac_models >= frac_cutoff) %>% 
      .$co_variate
    
    ## occurrence plot
    
    plot_tag = paste('\nBootstraps: n = ', 
                     length(boot_tbls), 
                     '\ncutoff = ', 
                     frac_cutoff, 
                     '\nbest co-variates: n = ', 
                     length(best_cov_set), 
                     '\ninitial co-variate set: n = ', 
                     length(variables), 
                     sep = '')
    
    occurrence_plot <- occurrence_counts %>% 
      make_occurr_plot(cutoff = frac_cutoff, 
                       plot_title = 'Co-variate selection', 
                       plot_subtitle = 'Co-variate occurrence in the bootstrap models',
                       plot_tag = plot_tag)
    
    
    message(paste('Elapsed:', Sys.time() - start_time))
    
    return(list(coef_lst = coef_lst, 
                occurrence_counts = occurrence_counts, 
                occurrence_plot = occurrence_plot, 
                best_cov_set = best_cov_set))
    
  }
  
  select_var_set <- function(inp_tbl, response, 
                             variables = names(inp_tbl)[names(inp_tbl) != response], 
                             nrow_boot = floor(nrow(inp_tbl) * 0.8),  n_boots = 100, 
                             family = 'poisson', alpha = 1, lambda_crit = 'lambda.1se', 
                             frac_cutoff = 0.95, .parallel = F, seed = NULL, ...) {
    
    ## The key function of the project: selects the optimal variable set
    ## based on the non-zero estimates of LASSO regression done with the random 
    ## subsets (bootstrap) of the input data. The optimal co-variate set is finally
    ## fed into the LASSO model applied to the input data.
    ## Arguments:
    ## inp_tbl: input data, ideally with named rows,
    ## response: modeling response, 
    ## variables: vector with the names of the variables included in the initial co-variate pool
    ## nrow_boot: size of the bootstrap subsets
    ## n_boots: number of the bootstraps
    ## family: glm family
    ## alpha: glmnet alpha parameter
    ## lambda_crit: defines a criterion for the non-zero co-efficient co-variate set selection in each bootstrap
    ## frac_cutoff: occurrence cutoff for the variable selection from the bootstrap model list
    ## .parallel: should the analysis be run in parallel
    
    ## bootstrap sets
    
    if(!is.null(seed)) {
      
      set.seed(seed = seed)
      
    }
    
    boot_list <- make_splits(inp_tbl = inp_tbl, 
                             nrow_train = nrow_boot, 
                             n_splits = n_boots)
    
    ## variable selection
    
    var_sel_results <- boot_list %>% 
      map(function(x) x$train) %>% 
      select_var_set_(boot_tbls = ., 
                      response = response, 
                      variables = variables, 
                      family = family, 
                      alpha = alpha, 
                      lambda_crit = lambda_crit, 
                      frac_cutoff = frac_cutoff, 
                      .parallel = .parallel, ...)
    
    ## fitting to the initial data
    
    final_fit <- model_cv_lasso(inp_tbl = inp_tbl, 
                                response = response, 
                                variables = var_sel_results$best_cov_set, 
                                family = family, 
                                alpha = alpha)
    
    ## obtaining the model measures: fit, error and pseudo-Rsq
    
    final_fit_measures <- assess.glmnet(final_fit,
                                        newx = model.matrix(~., inp_tbl[, var_sel_results$best_cov_set]), 
                                        newy = psych_lasso$analysis_tbl[[response]], 
                                        s = lambda_crit)
    
    final_fit_rsq_tbl <- tibble(lambda = final_fit$lambda, 
                                rsq = 1 - final_fit$cvm/var(psych_lasso$analysis_tbl[[response]])) ## credits to: https://myweb.uiowa.edu/pbreheny/7600/s16/notes/2-22.pdf
    
    
    final_fit_measures$rsq <- final_fit_rsq_tbl %>% 
      filter(lambda == final_fit[[lambda_crit]]) %>% 
      .$rsq %>% 
      set_names('rsq')
    
    attr(final_fit_measures$rsq, 'measure') = 'pseudo R squared'
    
    ## coefficents of the final fit, obtaining the co_variate names
    ## and levels with a regex search
    
    extr_regex <- paste(variables, 
                        collapse = '|')
    
    coef_tbl <- final_fit %>% 
      coef(s = lambda_crit) %>% 
      as.matrix %>% 
      as.data.frame %>% 
      set_names('estimate')
    
    coef_tbl <- coef_tbl %>% 
      rownames_to_column('coefficient_name') %>% 
      mutate(covariate = stri_extract(coefficient_name, 
                                      regex = extr_regex), 
             level = stri_replace(coefficient_name, 
                                  regex = extr_regex, 
                                  replacement = '')) %>% 
      as_tibble %>% 
      select(coefficient_name, 
             covariate, 
             level, 
             estimate) %>%
      filter(!(stri_detect(coefficient_name, fixed = 'Intercept') & estimate == 0)) %>% 
      mutate(coefficient_name = ifelse(stri_detect(coefficient_name, fixed = 'Intercept'), 
                                       'Intercept', 
                                       coefficient_name))
    
    ## final co-variate set
    
    final_var_set <- coef_tbl %>% 
      filter(!is.na(covariate), 
             estimate != 0) %>% 
      .$covariate %>% 
      unique
 
    ## output list
    
    result_list <- c(var_sel_results, 
                     list(final_fit = final_fit, 
                          final_fit_measures = final_fit_measures, 
                          final_fit_coef = coef_tbl, 
                          final_var_set = final_var_set))
    
    return(result_list)
    
  }
  
# END ------