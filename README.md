# lasso_tools
A medley of tools for LASSO modeling and prediction (data frame handling wrapper), LASSO-driven co-variate selection and generation of QC plots (residuals vs fitted, std residuals, qq plot).

Event though quite controversial, tools for inference of LASSO regression coefficients are provided: (1) bootstrap as implemented in the GitHub hdrm package (with a hook for weighted models, the root algorithm: https://rdrr.io/github/pbreheny/hdrm/man/boot.glmnet.html) and (2) multiple-split as described by Dezeure et al (https://arxiv.org/pdf/1408.4026).
