# demo for data analysis
source('/main_estS.R')
source('/appl_estbeta.R')


# evaluation function
evaluation_func <- function(boot_CI,est,method = 'scad'){
  # post-processing
  significant_columns <- which(boot_CI[1, ] > 0 | boot_CI[2, ] < 0)
  # regions and variable (exposures) names
  colnames(boot_CI) = as.vector(t(outer(colnames(est[[1]]),1:nrow(est[[1]]),paste,sep='_')))
  # results
  # and the corresponding effect estimate
  if(method=='scad'){
    estimate = as.vector(est$scad)[significant_columns]
  }else if (method == 'wang'){
    estimate = as.vector(est[[1]])[significant_columns]
  }else{
    stop('Method out of range')
  }
  return(list(columns = significant_columns,regions = colnames(boot_CI)[significant_columns],
              CI = boot_CI[,significant_columns], est = estimate))
}


# application
start_time = Sys.time()
appl_result = estbeta(yourX,yourY,yourXdemo)
Sys.time()-start_time


seresult = evaluation_func(appl_result$bootstrap$scad[[1]],appl_result$beta)


