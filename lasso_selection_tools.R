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
  require(hdi)
  require(OptimalCutpoints)

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
      mutate(misslab = try(ifelse(.candidate_missfit == 'yes',
                                  .rownames, 
                                  NA), silent = T),
             misslab = ifelse(any(class(misslab) == 'try-error'), NA, misslab))
    
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
  
  get_estimates_ <- function(linear_model, transf_fun = NULL, 
                             fallback = F, silent_messages = T, estimate_only = T, ...) {
    
    ## extract coefficients from a linear or generalized linear model (betas)
    ## together with confidence intervals obtained by OLS. The argument trans_fun allows for 
    ## transformation of the coefficients and CI (e.g. to convert them to OR in logistic regression)
    ## ... specifies other arguments to the CI-calculating function confint()
    ## the fallback option: calculates the CI based on the normality assumption: i. e. SEM * critical norm distr
    
    ## transforming function
    
    if(is.null(transf_fun)) {
      
      transf_fun <- function(x) x
      
    }
    
    ## model summary: to get se and p values
    
    mod_summary <- summary(linear_model)
    
    ## model estimates, error and CI, transforming
    
    model_coefs <- coefficients(linear_model)
    
    model_se <- mod_summary$coefficients[, 2]
    
    if(estimate_only) {
      
      ## returns solely the beta estimates
      
     return(tibble(parameter = names(model_coefs), 
                   estimate = unname(model_coefs), 
                   p_value = unname(mod_summary$coefficients[, 4])))
      
    }
    
    if(fallback) {
      
      model_ci <- tibble(lower_ci = model_coefs + qnorm(0.025)*model_se, 
                         upper_ci = model_coefs + qnorm(0.975)*model_se)
      
      
    } else {
      
      if(silent_messages) {
        
        model_ci <- suppressMessages(confint(linear_model, ...))
        
        
      } else {
        
        model_ci <- confint(linear_model, ...)
        
      }
      
      model_ci <- model_ci %>% 
        as_tibble %>% 
        set_names(c('lower_ci', 
                    'upper_ci'))
      
    }
    
    
    
    est_tibble <- model_ci %>% 
      mutate(estimate = unname(model_coefs)) %>% 
      map_dfc(transf_fun) %>% 
      mutate(parameter = names(model_coefs), 
             se = unname(model_se)) %>% 
      select(parameter, 
             estimate, 
             se, 
             lower_ci, 
             upper_ci)
    
    ## p values extracted from model summary
    ## n number of complete observations extracted from the model frame
    
    model_p <- mod_summary$coefficients[, 4]
    
    est_tibble <- est_tibble %>% 
      mutate(p_value = unname(model_p), 
             n_complete = nrow(model.frame(linear_model)))
    
    return(est_tibble)
    
  }
  
  make_ols_ <- function(inp_data, response, variables, weights, family, error_resistant = T, 
                        estimate_only = T) {
    
    ## fits a (gneralized) linear model to the data using OLS
    ## returns a table with model estimates
    
    
    if(error_resistant) {
      
      mod_formula <- try(paste(response, 
                               '~', 
                               paste(variables, collapse = '+')) %>% 
                           as.formula, 
                         silent = T)
      
      model <- try(glm(formula = mod_formula, 
                       data = inp_data, 
                       weights = weights, 
                       family = 'poisson'), silent = T)
      
      est_tbl <- try(get_estimates_(model, 
                                    estimate_only = estimate_only), 
                     silent = T)
      
      if(any(c(class(model), class(est_tbl), class(mod_formula)) == 'try-error')) {
        
        return(NULL)
        
      }
      
    } else {
      
      mod_formula <- paste(response, 
                           '~', 
                           paste(variables, collapse = '+')) %>% 
        as.formula
      
      model <- glm(formula = mod_formula, 
                   data = inp_data, 
                   weights = weights, 
                   family = 'poisson')
      
      est_tbl <- model %>% 
        get_estimates_(estimate_only = estimate_only)
      
    }
    
    return(est_tbl)
    
  }
  
  lotto_p_ <- function(p_val_vector) {
    
    ## aggregates the p value by the quantile method described in http://arxiv.org/abs/1408.4026
    
    quant_vec <- seq(0, 1, by = 0.01)
    
    corr_p_values <- quant_vec %>% 
      map_dbl(function(x) quantile(p_val_vector, x)/x)
    
    agg_p <- min(corr_p_values)
    
    return(agg_p)
    
  }
  
  test_perm_ <- function(est_sign, est_boots) {
    
    ## permutation bootstrap test
    
    est_boots <- est_boots[!is.na(est_boots)]
    
    greater_zero <- (est_sign*est_boots > 0) %>% 
      sum
    
    p <- (length(est_boots) - greater_zero)/length(est_boots)
    
    if(p == 0) {
      
      p <- 1/length(est_boots)
      
    }
    
    return(p)      
    
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
  
  get_measures_cv_lasso <- function(cv_object, 
                                    new_data_tbl, 
                                    response, 
                                    variables, 
                                    family,
                                    alpha = 1, 
                                    weights = NULL, 
                                    lambda_crit = 'lambda.1se') {
    
    ## gets basic model goodness measures out of a cv_model
    ## pseudo RSq should be interpreted as an approximation!
    
    fit_measures <- assess.glmnet(cv_object,
                                  newx = model.matrix(~., new_data_tbl[, variables]), 
                                  newy = new_data_tbl[[response]], 
                                  weights = weights, 
                                  family = family, 
                                  s = lambda_crit)

    fit_measures$n <- nrow(new_data_tbl)
    
    attr(fit_measures$n, 'measure') = 'n observations'

    fit_measures$rsq  <- tibble(lambda = cv_object$glmnet.fit$lambda, 
                                rsq = 1 - cv_object$glmnet.fit$dev.ratio) %>% 
      filter(lambda == cv_object[[lambda_crit]]) %>% 
      .$rsq
    
    attr(fit_measures$rsq, 'measure') = 'pseudo R squared'
    
    return(fit_measures)
    
  }
  
  get_qc_tbl_lasso <- function(cv_object, 
                               new_data_tbl, 
                               response, 
                               transf_fun = NULL, 
                               variables = names(new_data_tbl)[names(new_data_tbl) != response], 
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
                           type = type, ...)
    
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
  
  get_qc_plots_lasso <- function(cv_object, 
                                 new_data_tbl, 
                                 response, 
                                 transf_fun = NULL, 
                                 variables = names(new_data_tbl)[names(new_data_tbl) != response], 
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
  
  get_lasso_cv_stats <- function(lasso_modeling_list, 
                                 lambda_crit = 'lambda.1se', 
                                 cv = TRUE) {
    
    ## an upper level function.
    ## retrieves cross-validation and re-distribution stats from a lasso modeling result list
    
    redist_tbl <- as_tibble(lasso_modeling_list$fit_measures)
    
    if(cv) {
      
      lambda_type <- stri_replace(lambda_crit, fixed = 'lambda.', replacement = '')
      
      md_fit <- lasso_modeling_list$fit
      
      lambda_index <- md_fit$index[lambda_type, 1]
      
      cv_tbl <- tibble(error_type = md_fit$name, 
                       error_cv = md_fit$cvm[lambda_index], 
                       error_cv_lower_ci = md_fit$cvlo[lambda_index], 
                       error_cv_upper_ci = md_fit$cvup[lambda_index])
      
      redist_tbl <- as_tibble(cbind(redist_tbl, cv_tbl))
      
    }
    
    return(redist_tbl)
    
  } 
  
  get_lasso_pred_stats <- function(lasso_qc_tbl, family = 'binomial', tag.healthy = 0, ...) {
    
    ## calculates error measures based on the provided prediction table
    
    pred_tbl <- tibble(n = nrow(lasso_qc_tbl), 
                       mse = 2 * mean((lasso_qc_tbl$y - lasso_qc_tbl$.fitted)^2, na.rm = TRUE), 
                       mae = 2 * mean(abs(lasso_qc_tbl$y - lasso_qc_tbl$.fitted), na.rm = TRUE))
    
    if(family == 'binomial') {
      
      pred_tbl <- pred_tbl %>% 
        mutate(n_cases = nrow(filter(lasso_qc_tbl, y == 1)))
      
      opt_cutpoint_obj <- optimal.cutpoints(.fitted ~ y, 
                                            data = lasso_qc_tbl %>% 
                                              mutate(.fitted = unname(.fitted)) %>% 
                                              as.data.frame, 
                                            methods = 'Youden', 
                                            tag.healthy = tag.healthy, ...)
      
      opt_stats <- summary(opt_cutpoint_obj)
      
      cutpoint_stats <- opt_stats$p.table$Global$Youden %>% 
        as.data.frame %>% 
        t %>% 
        as_tibble
      
      auc_stats <- opt_stats$p.table$Global$AUC_CI %>% 
        stri_replace(., fixed = ",", replacement = '') %>% 
        stri_replace(., fixed = '(', replacement = '') %>% 
        stri_replace(., fixed = ')', replacement = '') %>%
        stri_split(., fixed = ' ', simplify = T) %>% 
        c
      
      cutpoint_stats <- cutpoint_stats %>% 
        mutate(auc = as.numeric(auc_stats[1]), 
               auc_lower_ci = as.numeric(auc_stats[2]), 
               auc_upper_ci = as.numeric(auc_stats[3]))
      
      pred_tbl <- cbind(pred_tbl, cutpoint_stats) %>% 
        as_tibble
      
    }
    
    return(pred_tbl)
    
  }
  
  
# bootstrap modeling function ------
  
  boot_lasso_ <- function(boot_tbls, response, 
                          variables = names(boot_tbls[[1]])[names(boot_tbls[[1]]) != response], 
                          weights = NULL, 
                          mod_fun = model_lasso, 
                          family = 'poisson', alpha = 1, lambda_crit = 'lambda.1se', seed = NULL, 
                          .parallel = F,  ...) {
    
    ## selects the optimal co-variate pool from the initial set ('variables')
    ## and the given bootstrap subsets of a table using the glmnet() lasso technique. 
    ## Lambda_crit defines a criterion for the variable set selection in each bootstrap
    ## co-variate was selected as a criterion for the inclusion in the final co-variate dataset
    ## Values: a table of co-efficients in the bootstrap models (co-eff selection by lambda_crit)
    
    start_time <- Sys.time()
    
    if(!is.null(seed)) {
      
      set.seed(seed = seed)
      
    }
    
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
                   weights = weights, 
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
            weights = weights,
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
  
# Inference GLM net modeling based on bootstrap -----
  
  model_boot_lasso <- function(inp_tbl, response, 
                               variables = names(inp_tbl)[names(inp_tbl) != response], 
                               weights = NULL, inference = T, 
                               nrow_boot = nrow(inp_tbl),  n_boots = 100, 
                               family = 'poisson', alpha = 1, lambda_crit = 'lambda.1se', 
                               ci_method = 'BCA', .parallel = F, seed = NULL, ...) {
    
    ## The key function of the project: fits a LASSO model to the data and, optionally, obtains estimate
    ## inference statistics (95% CI) by bootstrap. 
    ## P value for the coeffcients are calculated by a standard bootstrap test
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
    
    message(paste('Modeling:', response))
    
    ## fitting to the initial data, obtaining the coefficients and merging with the bootstrap stats
    
    final_fit <- model_cv_lasso(inp_tbl = inp_tbl, 
                                response = response, 
                                variables = variables, 
                                weights = weights, 
                                family = family, 
                                alpha = alpha, ...)
    
    ## getting the coeffcient values, variable names and their levels
    
    final_coeffs <- coef(final_fit, lambda_crit) %>% 
      as.matrix %>% 
      as.data.frame %>% 
      rownames_to_column('coeff') %>% 
      as_tibble %>% 
      set_names(c('coeff', 'estimate'))
    
    extr_regex <- paste(variables, 
                        collapse = '|')
    
    final_coeffs <- final_coeffs %>% 
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
             response = response, 
             n_complete = nrow(inp_tbl))
    
    ## obtaining model goodness paramaters
    
    final_fit_measures <- get_measures_cv_lasso(cv_object = final_fit, 
                                                new_data_tbl = inp_tbl, 
                                                response = response, 
                                                variables = variables, 
                                                family = family, 
                                                alpha = alpha, 
                                                weights = weights, 
                                                lambda_crit = lambda_crit)
    
    ## a vector with non-zero coeffcients
    
    final_var_set <- final_coeffs %>% 
      filter(!is.na(covariate), 
             estimate != 0) %>% 
      .$covariate %>% 
      unique
    
    if(!inference) {
      
      result_list <- list(response = response, 
                          fit = final_fit, 
                          fit_measures = final_fit_measures, 
                          fit_coef = final_coeffs, 
                          final_var_set = final_var_set) 
      
      return(result_list)
      
    }
    
    ## model inference, be sure, you want it!
    
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
                  weights = weights, 
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
    
    ## calculation of the p value with the classical permutation test
    
    sign_tbl <- final_coeffs %>% 
      mutate(est_sign = sign(estimate), 
             est_boots = coeff_tbl %>% 
               select(-X.Intercept..1) %>%
               as.list) %>% 
      select(coeff, 
             est_sign, 
             est_boots)
    
    p_boot <- sign_tbl %>% 
      select(est_sign, 
             est_boots) %>% 
      pmap_dbl(test_perm_)
    
    sign_tbl <- sign_tbl %>% 
      mutate(p_boot = p_boot)
    
    ## merging

    coef_tbl <- left_join(final_coeffs, 
                           coeff_stats, 
                           by = 'coeff') %>% 
      left_join(sign_tbl, 
                by = 'coeff') %>% 
      mutate(response = response) %>% 
      select(response, 
             coefficient_name, 
             covariate, 
             level, 
             n_complete, 
             estimate,  
             lower_ci, 
             upper_ci,  
             p_boot)

    ## output list
    
    result_list <- list(response = response, 
                        fit = final_fit, 
                        fit_measures = final_fit_measures, 
                        fit_coef = coef_tbl, 
                        final_var_set = final_var_set) 
    
    return(result_list)
    
  }
  
# Inference GLM net modeling based on the sample splitting -------

  model_split_lasso <- function(inp_tbl, response, 
                                variables = names(inp_tbl)[names(inp_tbl) != response], 
                                weights = NULL, 
                                n_boots = 100, 
                                family = 'poisson', alpha = 1, lambda_crit = 'lambda.1se', 
                                ci_method = 'BCA', 
                                .parallel_lasso = F, .parallel_ols = F, seed = NULL, ...) {
    
    ## The key function of the project: fits a LASSO model to the data and obtains estimate
    ## inference statistics (95% CI) is obtained by sample splitting as described at https://arxiv.org/pdf/1408.4026.pdf.
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
    
    message(paste('Modeling:', response))
    
    message(paste('Looking for the optimal variable set with '), 
            n_boots, 
            ' bootstrap tables')
    
    start_time <- Sys.time()
    
    ## generating the LASSO (training) and OLS splits (test)
    
    if(!is.null(seed)) {
      
      set.seed(seed = seed)
      
    }
    
    nrow_boot <-  floor(nrow(inp_tbl)/2)
    
    boot_list <- make_splits(inp_tbl = inp_tbl %>% 
                               mutate(weights = weights), 
                             nrow_train = nrow_boot, 
                             n_splits = n_boots)
    
    weight_list <- boot_list %>% 
      map(function(x) list(train = x$train$weights, 
                           test = x$test$weights))
    
    ## fitting the lasso models to the training splits, extracting the non-zero coefficients
    
    if(.parallel_lasso) {
      
      plan('multisession')
      
      lasso_models <- list(inp_tbl = map(boot_list, function(x) x$train), 
                           weights = map(weight_list, function(x) x$train)) %>% 
        future_pmap(model_cv_lasso, 
                    response = response, 
                    variables = variables, 
                    family = family, 
                    alpha = alpha, 
                    .options = furrr_options(seed = T), ...)
      
      plan('sequential')
      
    } else {
      
      lasso_models <- list(inp_tbl = map(boot_list, function(x) x$train), 
                           weights = map(weight_list, function(x) x$train)) %>% 
        pmap(model_cv_lasso, 
             response = response, 
             variables = variables, 
             family = family, 
             alpha = alpha, ...)
      
    }
    
    lasso_coefs <- lasso_models %>% 
      map(coef, 
          s = lambda_crit) %>% 
      map(as.matrix) %>% 
      map(as.data.frame) %>% 
      map(set_names, 
          'coeff') %>% 
      map(rownames_to_column, 
          'coeff_name')
    
    lasso_non_zero <- lasso_coefs %>% 
      map(filter, 
          coeff != 0, 
          !stri_detect(coeff_name, fixed = 'Intercept')) %>% 
      map(function(x) x$coeff_name)
    
    ## extracting the initial variable names from the list of non-zero estimates
    
    extr_regex <- paste(variables, 
                        collapse = '|')
    
    sel_variables <- lasso_non_zero %>% 
      map(stri_extract, 
          regex = extr_regex)
    
    ## fitting the OLS models to the test halves of the data splits
    ## with the variable sets selected by LASSO
    
    if(.parallel_ols) {
      
      plan('multisession')
      
      ols_est <- list(inp_data = map(boot_list, function(x) x$test), 
                      variables = sel_variables, 
                      weights = map(weight_list, function(x) x$test)) %>% 
        future_pmap(make_ols_, 
                    response = response, 
                    family = family, 
                    estimate_only = T, 
                    .options = furrr_options(seed = T))
      
      plan('sequential')
      
    } else {
      
      ols_est <- list(inp_data = map(boot_list, function(x) x$test), 
                      variables = sel_variables, 
                      weights = map(weight_list, function(x) x$test)) %>% 
        pmap(make_ols_, 
             response = response, 
             family = family, 
             estimate_only = T)
      
    }
    
    ## obtaining the expected value from the bootstrap estimate samples
    ## and 95%CI: either with the percentile or BCA method
    
    est_tbl <- ols_est %>% 
      compact %>% 
      map(select, 
          parameter, 
          estimate) %>% 
      reduce(full_join, 
             by = 'parameter') %>% 
      column_to_rownames('parameter') %>% 
      t %>% 
      as_tibble

    if(ci_method == 'BCA') {
      
      mod_summ_tbl <- est_tbl %>% 
        map(function(x) tibble(estimate = mean(x, na.rm = T), 
                               lower_ci = bca(x[!is.na(x)])[1], 
                               upper_ci = bca(x[!is.na(x)])[2]))
      
      
    } else {
      
      mod_summ_tbl <- est_tbl %>% 
        map(function(x) tibble(estimate = mean(x, na.rm = T), 
                               lower_ci = quantile(x, 0.025, na.rm = T), 
                               upper_ci = quantile(x, 0.975, na.rm = T)))
      
      
    }
    
    mod_summ_tbl <- mod_summ_tbl %>% 
      map2_dfr(., names(.), 
               function(x, y) mutate(x, coefficient_name = y)) %>% 
      mutate(covariate = stri_extract(coefficient_name, 
                                      regex = extr_regex), 
             covariate = ifelse(coefficient_name == '(Intercept)', 
                                'Intercept', 
                                covariate),
             level = stri_replace_all(coefficient_name, 
                                      regex = extr_regex, 
                                      replacement = ''), 
             level = ifelse(coefficient_name == '(Intercept)', 
                            'Intercept', 
                            level), 
             response = response, 
             n_complete = nrow(inp_tbl))
    
    ## calculation of the aggregated p values by the 'p value lottery' as described by http://arxiv.org/abs/1408.4026
    
    p_val_tbl <- ols_est %>% ## extracting t test p values for OLS estimates
      compact %>% 
      map(select, 
          parameter, 
          p_value) %>% 
      map(mutate, 
          p_value = p.adjust(p_value, 'BH')) %>% 
      reduce(full_join, 
             by = 'parameter') %>% 
      column_to_rownames('parameter') %>% 
      t %>% 
      as_tibble
    
    p_val_tbl <- p_val_tbl %>% ## substituting the p values for the variables not selected in the particular split with 1
      map_dfc(function(x) ifelse(is.na(x), 1, x))
    
    agg_p <- p_val_tbl %>% 
      map(lotto_p_) %>% 
      map2_dfr(., 
               names(.), 
               function(x, y) tibble(coefficient_name = y, 
                                     p_aggreg = unlist(x)))
    
    mod_summ_tbl <- left_join(mod_summ_tbl, 
                              agg_p, 
                              by = 'coefficient_name')
    
    ## calculation of the p values via permutation test, i.e the fractions of OLS splits
    ## where the estimate was greater or lower than 0, correction for multiple testing by 'BH' method
    
    sign_tbl <- mod_summ_tbl %>% 
      mutate(est_sign = sign(estimate), 
             est_boots = est_tbl %>% 
               as.list) %>% 
      select(coefficient_name, 
             est_sign, 
             est_boots)
    
    p_boot <- sign_tbl %>% 
      select(est_sign, 
             est_boots) %>% 
      pmap_dbl(test_perm_)
    
    sign_tbl <- sign_tbl %>% 
      mutate(p_boot = p_boot)
    
    mod_summ_tbl <- left_join(mod_summ_tbl, 
                              sign_tbl[, c('coefficient_name', 'p_boot')], 
                              by = 'coefficient_name')
    
    ## selecting the significant coefficient set
    
    final_var_set <- sign_tbl %>% 
      filter(p_boot < 0.05) %>% 
      .$coefficient_name
  
    ## the output list
    
    result_list <- list(response = response, 
                        fit = NULL, 
                        fit_measures = NULL, 
                        fit_coef = mod_summ_tbl, 
                        final_var_set = final_var_set) 
    
    message(paste('Elapsed:', 
                  Sys.time() - start_time))

    return(result_list)
      
  }
  
# END ------