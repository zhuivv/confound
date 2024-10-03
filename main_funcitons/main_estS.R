# estimation of S
## model function
model_list <-function(Y,r){
  # pre-value calculation
  require(nFactors);require(POET)
  hatr = nScree(cor(Y))$Components$nkaiser-r
  hatlambda = eigen(cor(Y))$values; hatev = eigen(cor(Y))$vectors
  hatsigma = 1/(ncol(Y))*sum(hatlambda[(hatr+r+1):ncol(Y)])
  delta = sqrt(2/nrow(Y)*log(ncol(Y))*hatsigma)
  gammahat = do.call(cbind,lapply(1:(hatr+r),function(j) sqrt(hatlambda[j])*hatev[,j]))
  # optimization
  library(gurobi)
  ## model
  M = ceiling(max(Y)/10)*10
  model_l0 <- list() # the input parameters working as (w,Y,eta)
  model_l0$obj <- matrix(c(rep(0,hatr+r),rep(0,ncol(Y)),rep(1,ncol(Y))),nrow = 1)
  model_l0$modelsense <- "min"
  model_l0$A <- rbind(cbind(gammahat,-diag(ncol(Y)),diag(ncol(Y))-diag(ncol(Y))),
                      cbind(matrix(0,nrow = 2*ncol(Y),ncol = hatr+r),
                            rbind(diag(-1,ncol(Y)),diag(1,ncol(Y))),
                            rbind(diag(-M,ncol(Y)),diag(-M,ncol(Y)))))
  model_l0$sense <- matrix(c(rep('=',ncol(Y)),rep('<',2*ncol(Y))),ncol = 1)
  model_l0$rhs <- matrix(c(rep(0,ncol(Y)),rep(delta,2*ncol(Y))),ncol = 1)
  model_l0$vtype <- matrix(c(rep('C',hatr+r),rep('C',ncol(Y)),rep('B',ncol(Y))),ncol = 1)
  model_l0$lb <- c(rep(-Inf,hatr+r+ncol(Y)),rep(0,ncol(Y)))
  model_l0$quadcon[[1]]$Qc <- diag(c(rep(1,hatr+r),rep(0,2*ncol(Y))))
  model_l0$quadcon[[1]]$rhs <- 1; model_l0$quadcon[[1]]$sense <- '='
  return(list(model_l0,hatr))
}

# greedy step
greedy_l0 <- function(Y,r,X,rate = 0.8){
  library(gurobi)
  ind_list <- vector('list',length = ncol(X))
  
  model = model_list(Y,r)
  opt_result = gurobi(model[[1]])
  ind_list[[1]] = which(opt_result$x[((model[[2]]+r)+ncol(Y)+1):((model[[2]]+r)+2*ncol(Y))]!=0)
  S0 <- ind_list[[1]]
  # S0_ = 
  for(l in 2:ncol(X)){
    r = r-1
    model = model_list(Y[,-S0],r)
    opt_result = gurobi(model[[1]])
    S0_ = which(opt_result$x[((model[[2]]+r)+ncol(Y[,-S0])+1):((model[[2]]+r)+2*ncol(Y[,-S0]))]!=0)
    ind_list[[l]] = setdiff(1:ncol(Y),S0)[S0_]
    S0 <- unique(unlist(ind_list))
    if((rate*ncol(Y))<length(S0)){
      S0 <- unique(unlist(ind_list[1:(l-1)]))
      break
    }
  }
  return(list(ind_list = ind_list, truest = S0))
}


