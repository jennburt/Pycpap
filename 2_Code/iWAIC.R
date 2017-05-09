#####################################################################
### Function to calculate iWAIC as derived by Yunyang Wang (2016) ###
### Code by D.K. Okamoto for Burt et al. urchin paper             ###
### Updated 5/8/2017                                             ###
#####################################################################


### Step 1: for each iteration, for each location, for each point in time, for each stage, sample probability distribution of latent abundance by generating 100 random observations

### step 2: for each observation within each site for each random draw, generate the log-probability

### step 3: integrate: find the log of the mean across random samples

### step 4: calculate WAIC of the integrated likelihood (iWAIC)

library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)

# The calculation of waic (integrated or not)  Returns lppd, p_waic_1, p_waic_2, and waic, which we define
# as 2*(lppd - p_waic_2), as recommmended in BDA
WAIC_fun <- function (lp_mat){
  lppd <- sum (log (colMeans(exp(lp_mat))))
  p_waic_1 <- 2*sum (log(colMeans(exp(lp_mat))) - colMeans(lp_mat))
  p_waic_2 <- sum (colVars(lp_mat))
  waic_2 <- -2*lppd + 2*p_waic_2
  return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1))
}

data <- DatQ%>%
  select(Y,S,Q,sm,med,lg)%>%
  rename(s= sm,m=med,l= lg,time=Y,site=S)%>%
  melt(id.vars= c("time","site","Q"),value.name= "obs",variable.name="class")%>%
  mutate(time= time-min(time)+1)

iWIAC <- function(jags_posterior,data,n_cores,n_reps){
  char_list_N<- c("NsP","NmP","NlP")
  char_list_g <- c("gamma_s","gamma_m","gamma_l")
  char_list_s <- c("sigS","sigY","Dispers")
  index_N <- grep(paste(char_list_N,collapse= "|"),dimnames(out$mcmc[[1]])[[2]])
  index_g <- grep(paste(char_list_g,collapse= "|"),dimnames(out$mcmc[[1]])[[2]])
  index_s <- grep(paste(char_list_s,collapse= "|"),dimnames(out$mcmc[[1]])[[2]])
  
  post <- ldply(out$mcmc)%>%
    add_rownames(var= "iter")%>%
    mutate(iter= as.numeric(as.character(iter)))
  
  ### dispersion should probably vary by group
  ### I think the model shouldn't be double exponential 
  ### additive instantanous mortality is a single exponential 
  ### double means otters and pycnos multiplicatively impact mortality. 
  
  post_N <- post[,c(1,index_N+1)]%>%
    melt(id.vars= "iter",value.name= "N")%>%
    separate(variable,into= c("class","time","site"),sep= c("\\[|\\,|\\]"))%>%
    mutate(class = gsub("N|P","",class))
  
  post_gamma <- post[,c(1,index_g+1)]%>%
    melt(id.vars= "iter",value.name= "gamma")%>%
    separate(variable,into= c("class","time","site"),sep= c("\\[|\\,|\\]"))%>%
    mutate(class = gsub("gamma_","",class))
  
  post_agg <- join(post_N,post_gamma)%>%
    right_join(post[,c(1,index_s+1)])
  ### loop to calculate random latent estimates 
  ### and calculate log-probability given parameter estimates
  ### and prior latent estimates (forward projection only... 
  ### ... for now, could do a VPA type calulation too)
  ### list of repetitions over which to integrate the likelihood
  iter_list <- as.list(1:n_reps)
  ### function to paralellize the likelihood integration
  lp_integrator <- function(x){
    lp_vec <- array(0,dim= c(nrow(data),1))
      for (i in 1:n_times){
        for(j in 1:n_sites){
          for(k in 1:n_classes){
            index <- with(data,time==i&site==j&class==unique(class)[k])
            with(subset(post_agg,time==i&site==j&class==unique(class)[k]&iter==x),
            lp <- sapply(rnorm(n_reps,0,1),function(z) dnbinom(x=data[index,"obs"],
              mu=N*exp(-exp(gamma+z*sigS)),
              size=Dispers)))
            lp_mat[index] <- apply(lp,1,function(x)log(mean(exp(x))))        
        }
      }
    }
  }
  ### generate the integrated likelihood over parallel threads
  lp_list <- mclapply(iter_list,lp_integrator,mc.cores= n_cores,mc.preschedule= FALSE)
  ### calculate iWAIC
  iWAIC <- WAIC(laply(lp_list))
  return(iWAIC)
}
