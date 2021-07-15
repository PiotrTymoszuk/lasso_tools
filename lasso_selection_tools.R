# This toolbox script contains functions for LASSO (or other model type provided by glmnet package)
# modeling taking a data frame as input data and functions used for selection of the optimal variable set
# and inference based on bootstrap subsets of the given data.
# At the moment, there may be a bit of problem with handling NAs in the input data sets.
# Please make sure, your input data contain complete observations only

# tools -----

  require(tidyverse)
  require(glmnet)
  require(furrr)
  require(ggrepel)
  require(coxed)
  require(stringi)

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
    
    ## if the nrow_train equals to nrow of the input table
    ## a simple bootstrap is performed
    
    if(is.null(rownames(inp_tbl))) {
      
      stop('The function requires a data frame with specified rownames')
      
    }
    
    tbl_len <- nrow(inp_tbl)
    
    split_names <- paste('split', 1:n_splits, sep = '_')
    
    if(nrow_train < nrow(inp_tbl)) {
      
      split_ids <- split_names %>% 
        map(function(x) sample(1:tbl_len, nrow_train, replace = F))
      
    } else {
      
      split_ids <- split_names %>% 
        map(function(x) sample(1:tbl_len, nrow_train, replace = T))
      
    }
    
    split_ids <- split_ids %>% 
      set_names(split_names)
    
    split_lst <- split_ids %>% 
      map(function(x) list(train = inp_tbl[x, ], 
                           test = inp_tbl[!(1:tbl_len) %in% x, ]))
    
    return(split_lst)
    
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
    
    inp_tbl <- inp_tbl %>% 
      select(all_of(c(response, variables))) %>% 
      filter(complete.cases(.))
    
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
  
  model_lasso <- function(inp_tbl, response, variables = names(inp_tbl)[names(inp_tbl) != response], 
                          family = 'poisson', alpha = 1, ...) {
    
    ## a customized wrapped around glmnet() for handling data frames
    ## categorical covariates are turned into 'dummy' variables by model.matrix()
    
    
    ## model components
    
    inp_tbl <- inp_tbl %>% 
      select(all_of(c(response, variables))) %>% 
      filter(complete.cases(.))
    
    y <- inp_tbl[[response]]
    
    x <- model.matrix(~ ., 
                      inp_tbl[, variables])
    
    ## cross-validated model selection
    
    net_model <- glmnet(x = x, 
                        y = y, 
                        family = family, 
                        alpha = alpha, ...)
    
    return(net_model)    
    
  }
  
# Prediction and QC functions handling the data frame input ------  
  
  get_qc_tbl_lasso <- function(cv_object, new_data_tbl, 
                               response, 
                               transf_fun = NULL, 
                               variables = names(inp_tbl)[names(inp_tbl) != response], 
                               lambda_crit = 'lambda.1se', 
                               type = 'response', ...) {
    
    ## a customized wrapper around predict.glmnet that handles data frames
    ## returns predictions, true values and residuals
    
    if(is.null(transf_fun)) {
      
      transf_fun <- function(x) x
      
    }
    
    
    x <- model.matrix(~ ., 
                      new_data_tbl[, variables])
    
    true_vals <- transf_fun(new_data_tbl[[response]])
    
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
                                 transf_fun = NULL, 
                                 variables = names(inp_tbl)[names(inp_tbl) != response], 
                                 lambda_crit = 'lambda.1se', 
                                 type = 'response', ...) {
    
    ## draws standard model qc plots with missfits labeled by obs. numbers
    ## ... are additional arguments passed to get_qc_tbl_lasso()
    
    ## QC table
    
    qc_tbl <- get_qc_tbl_lasso(cv_object = cv_object, 
                               new_data_tbl = new_data_tbl, 
                               response = response, 
                               variables = variables, 
                               transf_fun, 
                               lambda_crit = lambda_crit, 
                               type = type, ...)
    
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
  
# bootstrap modeling function ------
  
  boot_lasso_ <- function(boot_tbls, response, 
                          variables = names(boot_tbls[[1]])[names(boot_tbls[[1]]) != response], 
                          mod_fun = model_lasso, 
                          family = 'poisson', alpha = 1, lambda_crit = 'lambda.1se',
                          .parallel = F,  ...) {
    
    ## selects the optimal co-variate pool from the initial set ('variables')
    ## and the given bootstrap subsets of a table using the glmnet() lasso technique. 
    ## Lambda_crit defines a criterion for the variable set selection in each bootstrap
    ## co-variate was selected as a criterion for the inclusion in the final co-variate dataset
    ## Values: a table of co-efficients in the bootstrap models (co-eff selection by lambda_crit)
    
    start_time <- Sys.time()
    
    message(paste('Looking for the optimal variable set with '), 
            length(boot_tbls), 
            ' bootstrap tables')
    
    ## serial testing
    
    if(.parallel) {
      
      require(furrr)
      
      plan('multisession')
      
      res_models <- boot_tbls %>% 
        future_map(mod_fun, 
                   response = response, 
                   variables = variables, 
                   family = family, 
                   alpha = alpha, 
                   .options = furrr_options(seed = T), 
                   ...)
      
      plan('sequential')
      
    } else {
      
      res_models <- boot_tbls %>% 
        map(mod_fun, 
            response = response, 
            variables = variables, 
            family = family, 
            alpha = alpha, ...)
      
    }
    
    ## selection of non-zero coefficients
    ## co-efficent names are extracted from the results 
    ## via a regex made of the input variable vector
    
    coef_lst <- res_models %>% 
      map(coef, 
          s = lambda_crit) %>% 
      map(as.matrix) %>% 
      map(as.data.frame)
    
    message(paste('Elapsed:', Sys.time() - start_time))
    
    return(coef_lst %>% 
             map(rownames_to_column, 
                 'coeff') %>% 
             map(as_tibble))

  }
  
# inference GLM net modeling -----
  
  model_boot_lasso <- function(inp_tbl, response, 
                               variables = names(inp_tbl)[names(inp_tbl) != response], 
                               nrow_boot = nrow(inp_tbl),  n_boots = 100, 
                               family = 'poisson', alpha = 1, lambda_crit = 'lambda.1se', 
                               ci_method = 'BCA', .parallel = F, seed = NULL, ...) {
    
    ## The key function of the project: fits a LASSO model to the data and obtains estimate
    ## inference statistics (95% CI) by bootstrap. P value for the coeffcients are calculated
    ## by Wilcoxon test
    ## Arguments:
    ## inp_tbl: input data, ideally with named rows,
    ## response: modeling response, 
    ## variables: vector with the names of the variables included in the initial co-variate pool
    ## nrow_boot: size of the bootstrap subsets
    ## n_boots: number of the bootstraps
    ## family: glm family
    ## alpha: glmnet alpha parameter
    ## lambda_crit: defines a criterion for the non-zero co-efficient co-variate set selection in each bootstrap
    ## .parallel: should the analysis be run in parallel
    ## The bootstrap CI may be calculated by BCA (coxed) or percentile method
    
    ## fitting to the initial data, obtaining the coefficients and merging with the bootstrap stats
    
    final_fit <- model_cv_lasso(inp_tbl = inp_tbl, 
                                response = response, 
                                variables = variables, 
                                family = family, 
                                alpha = alpha, ...)
    
    final_coeffs <- coef(final_fit, lambda_crit) %>% 
      as.matrix %>% 
      as.data.frame %>% 
      rownames_to_column('coeff') %>% 
      as_tibble %>% 
      set_names(c('coeff', 'estimate'))
    
    ## fixed lambda used for bootstrap models
    
    lambda_fix <- final_fit[[lambda_crit]]
    
    ## bootstrap sets
    
    if(!is.null(seed)) {
      
      set.seed(seed = seed)
      
    }
    
    boot_list <- make_splits(inp_tbl = inp_tbl, 
                             nrow_train = nrow_boot, 
                             n_splits = n_boots)
    
    ## coefficients in the bootstrap models
    
    boot_coefs <- boot_list %>% 
      map(function(x) x$train) %>% 
      boot_lasso_(boot_tbls = ., 
                  response = response, 
                  variables = variables, 
                  family = family, 
                  alpha = alpha, 
                  lambda_crit = lambda_fix, 
                  mod_fun = model_lasso, 
                  .parallel = .parallel, ...)
    
    coeff_tbl <- boot_coefs %>% 
      reduce(left_join, 
             by = 'coeff') %>% 
      column_to_rownames('coeff') %>% 
      t %>% 
      as_tibble
    
    ## obtaining the coefficient bootstrap percentiles and 95% CI by the BCA method provided by coxed
    
    if(ci_method == 'BCA') {
      
      coeff_stats <- coeff_tbl %>% 
        map2_dfr(.,
                 names(.), 
                 function(x, y) tibble(coeff = y, 
                                       lower_ci = bca(x)[1], 
                                       upper_ci = bca(x)[2]))
      
    } else {
      
      coeff_stats <- coeff_tbl %>% 
        map2_dfr(.,
                 names(.), 
                 function(x, y) tibble(coeff = y, 
                                       lower_ci = quantile(x, 0.025), 
                                       upper_ci = quantile(x, 0.975)))
      
    }
    
    coeff_tests <- coeff_tbl %>% 
      map(function(x) suppressWarnings(wilcox.test(x))) %>% 
      map2_dfr(., names(.),
               function(x, y) tibble(coeff = y, 
                                     p_value = x$p.value))
    
    coef_tbl <- left_join(final_coeffs, 
                           coeff_stats, 
                           by = 'coeff') %>% 
      left_join(coeff_tests, 
                by = 'coeff') %>% 
      mutate(response = response)
    
    ## obtaining the model measures: n_cases, fit, error and pseudo-Rsq
    
    final_fit_measures <- assess.glmnet(final_fit,
                                        newx = model.matrix(~., inp_tbl[, variables]), 
                                        newy = inp_tbl[[response]], 
                                        s = lambda_crit)
    
    final_fit_rsq_tbl <- tibble(lambda = final_fit$lambda, 
                                rsq = 1 - final_fit$cvm/var(inp_tbl[[response]])) ## credits to: https://myweb.uiowa.edu/pbreheny/7600/s16/notes/2-22.pdf
    
    final_fit_measures$n <- nrow(inp_tbl)
    
    attr(final_fit_measures$n, 'measure') = 'n observations'
    
    final_fit_measures$rsq <- final_fit_rsq_tbl %>% 
      filter(lambda == final_fit[[lambda_crit]]) %>% 
      .$rsq %>% 
      set_names('rsq')
    
    attr(final_fit_measures$rsq, 'measure') = 'pseudo R squared'
    
    ## coefficents of the final fit, obtaining the co_variate names
    ## and levels with a regex search
    
    extr_regex <- paste(variables, 
                        collapse = '|')

    coef_tbl <- coef_tbl %>% 
      mutate(covariate = stri_extract(coeff, 
                                      regex = extr_regex), 
             level = stri_replace(coeff, 
                                  regex = extr_regex, 
                                  replacement = ''), 
             coefficient_name = coeff)  %>%
      filter(!(stri_detect(coefficient_name, fixed = 'Intercept') & estimate == 0)) %>% 
      mutate(coefficient_name = ifelse(stri_detect(coefficient_name, fixed = 'Intercept'), 
                                       'Intercept', 
                                       coefficient_name), 
             response = response) %>% 
      select(response, 
             coefficient_name, 
             covariate, 
             level, 
             estimate,  
             lower_ci, 
             upper_ci,  
             p_value)
    
    ## final co-variate set
    
    final_var_set <- coef_tbl %>% 
      filter(!is.na(covariate), 
             estimate != 0) %>% 
      .$covariate %>% 
      unique
    
    ## output list
    
    result_list <- list(response = response, 
                        fit = final_fit, 
                        fit_measures = final_fit_measures, 
                        fit_coef = coef_tbl, 
                        final_var_set = final_var_set) 
                     
    
    return(result_list)
    
  }
  
# END ------