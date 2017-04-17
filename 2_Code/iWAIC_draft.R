#####################################################################
### Function to calculate iWAIC as derived by Yunyang Wang (2016) ###
### Code by D.K. Okamoto for Burt et al. urchin paper             ###
### Updated 4/14/2017                                             ###
#####################################################################


### Step 1: for each iteration, for each location, for each point in time, for each stage, sample probability distribution of latent abundance by generating 100 random observations

### step 2: for each observation within each site for each random draw, generate the log-probability

### step 3: integrate: find the log of the mean across random samples

### step 4: calculate WAIC of the integrated likelihood (iWAIC)


# The calculation of waic (integrated or not)  Returns lppd, p_waic_1, p_waic_2, and waic, which we define
# as 2*(lppd - p_waic_2), as recommmended in BDA
WAIC <- function (lp_mat){
  lppd <- sum (log (colMeans(exp(lp_mat))))
  p_waic_1 <- 2*sum (log(colMeans(exp(lp_mat))) - colMeans(lp_mat))
  p_waic_2 <- sum (colVars(lp_mat))
  waic_2 <- -2*lppd + 2*p_waic_2
  return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
}

iWAIC_BCurch <- function (model,n=100){
  
  ### extract the fitted model 
  fit <- extract(fitted_model)
  
  ### loop to calculate random latent estimates 
  ### and calculate log-probability given parameter estimates
  ### and prior latent estimates (forward projection only... 
  ### ... for now, could do a VPA type calulation too)
  
  ### generate blank log-probability array 
  lp_array <- array(0,dim= c(length(data$Obs),n_iter,n_rand))
  
  for (i in 1:n_iter){
    for (j in 1:n_stages){
      for (k in 1:n_sites) {
        for (l in 2:(n_times-1)){
          
          ### assumes deviations from expected survival are uncorrelated... 
          LV_rand[i,j,k,l,] = N[i,j-1,k,l-1]*exp(fit$mort[i,j,k,l]+rnorm(n_rand,0,1/fit$tau[i,j]))+
                        N[i,j-1,k,l-1]*exp(fit$mort[i,j,k,l]+rnorm(n_rand,0,1/fit$tau[i,j]))
          
          ### caluclate the log-probability for each random draw and observation
          for (m in 1:n_rand){
            
            ### generate obervation index for subsetting
            obs_index = with(data,stage== j&site==k&time==l)
            
            ### generate log-probability using appropriate likelihood
            lp_array[obs_index,i,] = dnegbin(subset(data,stage== j&site==k&time==l)$Obs,
                                           mean= LV_rand[i,j,k,l,], 
                                           dispers=fit$dispers[i],log= T) 
          }
        }
      }
    }
  }
  
  ### calculate the matrix log of the mean probability 
  ### for each observation and iteration
  lp_mat <- apply(lp_array,c(1,2),log(mean(exp(x))))
  
  ### calculate iWAIC
  iWAIC <- WAIC(lp_mat)
  
  return(iWAIC)
}
