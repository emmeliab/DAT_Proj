
library(tidyverse)
library(nlme)

dat <- read_csv("~/Downloads/lf_diffs_summ_TPU.csv")
dat <- dat %>% mutate(
  treeid = factor(treeid)
)

# fit the current model
mfit <- lme(
  j_diff ~ 1,
  random = ~ 1 | treeid,
  weights = varIdent(form = ~ 1 | treeid),
  data = dat
)

# plot a histogram of the residuals
hist(residuals(mfit, type = "pearson"))

# not as nicely symmetric as we would hope

# bootstrapping simulation to see how much this matters

# function to create a semi-parametric bootstrap
boot <- function(mfit, max_tries = 10, resp_var){
  sd_re <- VarCorr(mfit)["(Intercept)", "StdDev"] |> as.numeric()
  # simulate new random effects
  current_re <- mfit$coefficients$random$treeid
  
  new_re <- rnorm(nrow(current_re), sd = sd_re)
  
  # build new observations
  newdat <- data.frame(
    treeid = row.names(current_re) |> factor(),
    re = new_re
  )
  # subset original data
  dat_orig <- mfit$data[, which(colnames(mfit$data) %in% c("treeid", resp_var))]
  newdat <- merge(newdat, dat_orig, all.y = T)
  
  # now replace j_diff with a bootstrapped version
  newdat[, resp_var] <- mfit$coefficients$fixed + newdat$re +
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
    newdat[, resp_var] <- mfit$coefficients$fixed + newdat$re +
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
  } else{
    return(mfit_boot$coefficients$fixed)
  }
}

# now apply the boot function over 250 iterations
estims_boot <- sapply(
  1:250,
  boot,
  mfit = mfit,
  resp_var = "j_diff"
)


# plot density of bootstrapped estimates against the 
# theoretical sampling distribution
hist(estims_boot, freq = F)

# theoretical density
m <- mfit$coefficients$fixed
se <- sqrt(mfit$varFix[1,1])

lines(
  x = seq(min(estims_boot) - 1, max(estims_boot) + 1, length.out = 100),
  y = dnorm(seq(min(estims_boot) - 1, max(estims_boot) +1, length.out = 100), mean = m, sd = se),
  col = "blue"
)

# bootstrapped estimate
mean(estims_boot)

# bootstrapped CI
quantile(estims_boot, probs = c(0.025, 0.975))

