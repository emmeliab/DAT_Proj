# This code creates a function to create a semi-parametric bootstrap courtesy of Dusty Gannon (dustin.gannon@oregonstate.edu)
# Code is associated with the article at DOI: 10.1093/treephys/tpae153
# Licence information:
# Questions can be directed to Loren Albert (corresponding author) at loren.albert@oregonstate.edu, Emmelia Braun (first-author) at emmelia.braun@oregonstate.edu, or Charles Southwick (first author) at charles.southwick@oregonstate.edu


boot <- function(mfit, max_tries = 10, resp_var, calc_null = FALSE){
    sd_re <- VarCorr(mfit)["(Intercept)", "StdDev"] |> as.numeric()
    # simulate new random effects
    current_re <- mfit$coefficients$random$treeid
    current_fe <- ifelse(calc_null, 0, mfit$coefficients$fixed)
    
    new_re <- rnorm(nrow(current_re), sd = sd_re)
    
    # build new observations
    newdat <- data.frame(
        treeid = row.names(current_re) |> factor(),
        re = new_re
    )
    # subset original data
    dat_orig <- mfit$data[, which(colnames(mfit$data) %in% c("treeid", resp_var))]
    newdat <- merge(newdat, dat_orig, all.y = T)
    
    # now replace resp_var with a bootstrapped version
    newdat[, resp_var] <- current_fe + newdat$re +
        sample(residuals(mfit, type = "response"), size = nrow(newdat), replace = T)
    
    # new refit and return the intercept estimate
    mfit_boot <- tryCatch({
        update(mfit, data = newdat)
    }, error = function(e) {
        # Handle the error
        message("singular fit")
        return(NULL)  # or some other fallback value
    })
    counter <- 1
    
    while(is.null(mfit_boot) & counter <= max_tries) {
        # try again in case the resampling yielded zero variance
        # within one of the trees
        newdat[, resp_var] <- current_fe + newdat$re +
            sample(residuals(mfit, type = "response"), size = nrow(newdat), replace = T)
        
        # new refit and return the intercept estimate
        mfit_boot <- tryCatch({
            update(mfit, data = newdat)
        }, error = function(e) {
            # Handle the error
            message("singular fit")
            return(NULL)  # or some other fallback value
        })
    }
    
    
    if(is.null(mfit_boot)){
        return(NA)
    } else {
        return(mfit_boot$coefficients$fixed)
    }
    
    
}

# # now apply the boot function over 500 iterations
# null_boot <- sapply(
#     1:500,
#     boot,
#     mfit = mfit,
#     resp_var = "vc_diff",
#     calc_null = TRUE
# )
# 
# estims_boot <- sapply(
#     1:500,
#     boot,
#     mfit = mfit,
#     resp_var = "vc_diff",
#     calc_null = FALSE
# )


# mean( abs(mean(estims_boot)) <= null_boot ) * 2
# 
# mean(estims_boot)
# 
# summary(mfit)
