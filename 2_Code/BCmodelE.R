# Set current directory
dirname = "O:/Documents/MATLAB/SanNic" # Set this to desired directory
setwd(dirname)

## -------- Load packages --------
library(lattice)
library(coda)
library(boot)
## -------- Bayesian model description --------
jagsfile = 'BCUrcE.jags'
## -------- Read and prepare data for Bayesian modeling --------
library(readxl)
DatS <- read_excel("Site_Sum_stats.xlsx")
DatQ <- read_excel("QuadCounts.xlsx")
attach(DatS)
attach(DatQ)

Growth_Sm_Md = 0.517
Growth_Md_Lg = 0.057
R_est = 10
SiteN = SiteNum; 
YearN = Year - min(Year) + 1
Yr = Y -  min(Year) + 1
St = S 
Nsites =  max(SiteN)
Nyears =  max(YearN)
Nobs = length(Yr) 
Ott = matrix(0,nrow = Nsites, ncol = Nyears)
Pyc = matrix(0,nrow = Nsites, ncol = Nyears)
counts = cbind(sm, med, lg)
N0 = cbind(AvgOfsm[1:Nsites], AvgOfmed[1:Nsites], AvgOflg[1:Nsites])

# Deternine Otter and Pycno status (relative abund, 0-1 scale) at each site/year
for (s in 1:Nsites){
   for (y in 1:Nyears){
     Ott[s,y] = mean(Otter[which(SiteN==s & YearN==y)])
     Pyc[s,y] = mean(Pycno[which(SiteN==s & YearN==y)])
   }
} 

# Create data and initials objects for JAGS
jags.data <- list(NS = Nsites, obs = counts, NY = Nyears, NEQ = 10, 
                  Y = Yr, S = St, NObs = Nobs, g1 = Growth_Sm_Md,  
                  g2 = Growth_Md_Lg, O = Ott, P = Pyc, N0 = N0, R=R_est) # 
inits <- function() list(Pobs = runif(1, .1, .5), dbar = runif(1, .5, 2),
                         beta1 = runif(1, .1, .5), beta2 = runif(1, .05, .3),
                         sigS = runif(1, .5, 2), sigY = runif(1, .5, 2),
                         Dispers = runif(1, .2, 1), # R = runif(1, 2, 20), 
                         alphaP = runif(3, .1, 1),
                         alphaO = runif(2, .1, 1)) 
# R = runif(1, 5, 60), 
params <- c("dbar","Pobs","sigS","sigY","Dispers","beta1","beta2",
            "alphaP","alphaO","loglikS","Ns","Nm","Nl","d") # NOTE: removed "R", made constant
nsamples <- 1000
nt <- 5
nb <- 3000

# For parallel (comment out for serial)
library(parallel)
cores = detectCores()
ncore = min(25,cores-1)
cl <- makeCluster(ncore)
nc <- ncore


## -------- Call JAGS from R --------
library(rjags)
library(runjags)
out <- run.jags(data = jags.data, 
                inits = inits, 
                monitor = params, 
                model = jagsfile, 
                n.chains = nc, 
                thin = nt, 
                sample = nsamples, 
                burnin = nb,
                method="rjparallel", cl=cl)

sumstats = summary(out)
vn = row.names(sumstats)
# create matrix of all posterior chains (coda)
post = rbind(out$mcmc[[1]], out$mcmc[[2]])
for (i in 3:nc){
  post = rbind(post, out$mcmc[[i]])
}

# Diagnostic plots---------------------------
plot(out, c('trace','histogram'),
     vars=c("Pobs","Dispers","sigS","sigY","alphaP","alphaO",
            "dbar","beta1","beta2")) # "R",

plot(out, c('trace','histogram'),
     vars=c("d")) # "R",

# "LgHazratPsml",
## --------- Calculate WAIC ----------
# library(loo)
# mc_ll <- post[,paste0("loglikS[",1:Nobs,"]")]
# WAIC = waic(mc_ll)
# LOO = loo(mc_ll)
# print(paste0("WAIC = ", as.character(format(WAIC$waic, digits=4 ))))

save(list = ls(all.names = TRUE),file='BC_UrcE_stoch_Results.RData')