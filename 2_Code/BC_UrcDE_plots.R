# Load MCMC results from JAGS run
setwd("~/Dropbox/Urch-Pyc-Ott Paper")

# load("BC_Urc5_det_Results.Rdata")
load("BC_UrcE_stoch_Results.Rdata")

## -------- Plot some results --------
library(ggplot2)
library(matrixStats)
library(boot)
library(MASS)

rr = c(which(row.names(sumstats)=="alphaO[1]"),which(row.names(sumstats)=="alphaP[1]"),
       which(row.names(sumstats)=="alphaO[2]"),which(row.names(sumstats)=="alphaP[2]"),
       which(row.names(sumstats)=="alphaP[3]"))
d1=data.frame(names=c("Otter-->Medium","Pycno-->Medium","Otter-->Large","Pycno-->Large",
                      "Pycno-->Small"),
              mean=sumstats[rr,4],
              lower=sumstats[rr,1], upper=sumstats[rr,3])
d1$names <- factor(d1$names, levels = d1$names)
plt1 = ggplot() + 
  geom_errorbar(data=d1, mapping=aes(x=names, ymin=lower, ymax=upper), stat="identity", width=0.2, size=1, color="blue") + 
  geom_point(data=d1, mapping=aes(x=names, y=mean),stat="identity", size=4, shape=21, fill="white") +
  xlab('Predator-->Urchin Size Class') +
  ylab('Log(Hazard Ratio)') +
  theme_bw() +
  ggtitle("Predation Rates, Sea Otter vs Pycnopodia")
# 
print(plt1)


# ---Time Series by Site--------------------------------------

PobsMn = sumstats[which(vn=="Pobs"),4]

Treatment_labs = c("Otters all years, High->Lo stars by 2015",
                   "Otters all years, Med->Lo stars by 2015",
                   "Otters all years, Med->Lo stars by 2015",
                   "Otters all years, Med->Lo stars by 2015",
                   "Otters all years, High->Lo stars by 2015",
                   "Otters after 2013, High->Lo stars by 2015",
                   "Otters after 2013, High->Lo stars by 2015",
                   "Otters after 2013, High->Lo stars by 2015",
                   "No Otters, Med->Lo stars by 2015",
                   "No Otters, High->Lo stars by 2015",
                   "No Otters, High->Lo stars by 2015")

for (i in 1:Nsites){
  # Estimate bootstrap mean from sample quadrats
  # (for validation, plot bootstrap mean vals)  
  
  samplemean <- function(x, d) {
    return(mean(x[d]))
  }
  mnqdS = numeric(length=4)
  mnqdM = numeric(length=4)
  mnqdL = numeric(length=4)
  YearsF = as.factor(c(2013,2014,2015,2016))
  for (yy in 2013:2016){
    Yearscol = rep(yy,5000)
    ii = which(DatQ$S==i & DatQ$Y==yy)
    b = boot(DatQ$sm[ii],samplemean,5000)
    stS = as.numeric(b$t)
    mnqdS[yy-2012] = mean(stS)
    ii = which(DatQ$S==i & DatQ$Y==yy)
    b = boot(DatQ$med[ii],samplemean,5000)
    stM = as.numeric(b$t)
    mnqdM[yy-2012] = mean(stM)
    ii = which(DatQ$S==i & DatQ$Y==yy)
    b = boot(DatQ$lg[ii],samplemean,5000)
    stL = as.numeric(b$t)
    mnqdL[yy-2012] = mean(stL)
  }
  Quadmeans = data.frame(names=YearsF,Small=stS,Medium=stM,Large=stL)

  rrS = c(which(row.names(sumstats)==paste0("Ns[",i,",1]")),which(row.names(sumstats)==paste0("Ns[",i,",2]")),
          which(row.names(sumstats)==paste0("Ns[",i,",3]")),which(row.names(sumstats)==paste0("Ns[",i,",4]")))
  rrM = c(which(row.names(sumstats)==paste0("Nm[",i,",1]")),which(row.names(sumstats)==paste0("Nm[",i,",2]")),
          which(row.names(sumstats)==paste0("Nm[",i,",3]")),which(row.names(sumstats)==paste0("Nm[",i,",4]")))
  rrL = c(which(row.names(sumstats)==paste0("Nl[",i,",1]")),which(row.names(sumstats)==paste0("Nl[",i,",2]")),
          which(row.names(sumstats)==paste0("Nl[",i,",3]")),which(row.names(sumstats)==paste0("Nl[",i,",4]")))

  #        which(row.names(sumstats)=="DensPyc"),which(row.names(sumstats)=="DensBoth"))
  d2s=data.frame(names=YearsF,groupN=c("Small"),
                 Mean=sumstats[rrS,4]*PobsMn,
                 Lower=sumstats[rrS,1]*PobsMn,
                 Upper=sumstats[rrS,3]*PobsMn,
                 Quadmean = mnqdS)              
  d2m=data.frame(names=YearsF,groupN=c("Medium"),
                 Mean=sumstats[rrM,4],Lower=sumstats[rrM,1],
                 Upper=sumstats[rrM,3],
                 Quadmean = mnqdM)              
  d2l=data.frame(names=YearsF,groupN=c("Large"),
                 Mean=sumstats[rrL,4],Lower=sumstats[rrL,1],
                 Upper=sumstats[rrL,3],
                 Quadmean = mnqdL)              
  d2 = rbind(d2s,d2m,d2l)
  
  pd <- position_dodge(0.1) # move them .05 to the left and right
  
  plt3 = ggplot(d2, aes(x=names, y=Mean, colour=groupN, group=groupN)) + 
    geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=4, shape=19) +
#               aes(name = "Urchin Size Class", fill = factor(groupN)),
#               show.legend=FALSE) + 
    # scale_fill_manual(values=c("red", "green", "blue")) +
    xlab("Year") +
    ylab("Estimated mean abundance") +
    scale_colour_hue(name="Urchin Size Class",    # Legend label, use darker colors
                     breaks=c("Large", "Medium", "Small"),
                     labels=c("Large", "Medium", "Small"),
                     l=40) +                    # Use darker colors, lightness=40
    ggtitle(paste0("Estimated True Abundance, Site ",i," \n (", Treatment_labs[i]  )) +
    expand_limits(y=0) +                        # Expand y range
    # scale_y_continuous(trans = "log10") +
    theme_bw() +
    theme(legend.justification=c(1,0),
          legend.position=c(.5,.7))               # Position legend in bottom right
  print(plt3)
   
  val = data.frame(year = d2$names, site = as.factor(i), size = d2$groupN,
                   meanpred = d2$Mean,meanquad=d2$Quadmean)
  if (i==1){
    Valdf = val
  }else{
    Valdf = rbind(Valdf,val)
  }
}

# Model validation plot ------------------------------
# Linear model, observed vs predicted mean density
fit = lm(meanquad~meanpred,Valdf)
resid = as.numeric(fit$residuals)
jj = which(abs(resid)<2.5)
# Linear model, excluding two outliers
fit = lm(meanquad~meanpred,Valdf[jj,])
R2 = summary(fit)$r.squared
coef = as.numeric(fit$coefficients)

plt4 = ggplot(Valdf, aes(x=meanpred, y=meanquad, colour=size, group=size)) +
  geom_point(size=4, shape=19) +
  geom_abline(intercept = coef[1], slope = coef[2]) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  xlab("Model Predicted Abundance") +
  ylab("Mean Quadrat Counts") +
  annotate("text", x = 6.5, y = 4, label = paste0("R2 = ", format(R2,digits=3) )) +
  ggtitle("Observed vs. Predicted Density by Size Class, All Sites/Years") 

print(plt4)  

# Violin plots of Hazard rates -------------------
simreps = 5000
reps = dim(post); reps = reps[1]
rs = sample(reps,simreps) ## 
postsamp = post[rs,]


d31 = data.frame(effect = "alphaPsm",
                 value = (postsamp[,which(vn=="alphaP[3]")]) ) 
d32 = data.frame(effect = "alphaPmd",
                 value = (postsamp[,which(vn=="alphaP[1]")]) )
d33 = data.frame(effect = "alphaPlg",
                 value = (postsamp[,which(vn=="alphaP[2]")]) )
d34 = data.frame(effect = "alphaOsm",
                 value = numeric(length=simreps)) 
d35 = data.frame(effect = "alphaOmd",
                 value = (postsamp[,which(vn=="alphaO[1]")]) )
d36 = data.frame(effect = "alphaOlg",
                 value = (postsamp[,which(vn=="alphaO[2]")]) )

d3 = rbind(d31,d32,d33,d34,d35,d36)

plt5 <- ggplot(d3, aes(x=effect, y=value)) + 
  geom_violin() +
  geom_boxplot(width=0.1) +
  xlab("Predator Effect") +
  ylab("Log(Hazard Ratio)") +
  scale_x_discrete(breaks=c("alphaPsm", "alphaPmd", "alphaPlg",
                            "alphaOsm", "alphaOmd", "alphaOlg"),
                   labels=c("Pyc-Sm", "Pyc-Med", "Pyc-Lg",
                            "Ott-Sm", "Ott-Med", "Ott-Lg")) +
  ggtitle("Effect of Predators on Urchin Size Classes")

print(plt5)

# ---Estimate equilibrial states-------------------------

NYEQ = 25
# 4 scenerios: No preds, Both preds, Otters, Pycno
P = c(0, 1, 0, 1); O = c(0, 1, 1, 0);

g1 = Growth_Sm_Md
g2 = Growth_Md_Lg
R = R_est
alphaP = matrix(nrow=simreps,ncol=3); alphaO = matrix(nrow=simreps,ncol=2); 

dbar = postsamp[,which(vn=="dbar")]
beta1 = postsamp[,which(vn=="beta1")]
beta2 = postsamp[,which(vn=="beta2")]
Pobs = postsamp[,which(vn=="Pobs")]
b1 = dbar*beta1
b2 = dbar*beta2 
alphaP[,1] = postsamp[,which(vn=="alphaP[1]")]
alphaP[,2] = postsamp[,which(vn=="alphaP[2]")]
alphaP[,3] = postsamp[,which(vn=="alphaP[3]")]
alphaO[,1] = postsamp[,which(vn=="alphaO[1]")]
alphaO[,2] = postsamp[,which(vn=="alphaO[2]")]

NeqS = matrix(nrow=simreps,ncol=4)
NeqM = matrix(nrow=simreps,ncol=4)
NeqL = matrix(nrow=simreps,ncol=4)
Ns = numeric(length=NYEQ); Ns[1] = 25 
Nm = numeric(length=NYEQ); Nm[1] = 5
Nl = numeric(length=NYEQ); Nl[1] = 5

for (r in 1:simreps){
  for (s in 1:4){
    for (t in 2:NYEQ){
      NsP = Ns[t-1] + R - g1*Ns[t-1] 
      NmP = Nm[t-1] + g1*Ns[t-1] - g2*Nm[t-1] 
      NlP = Nl[t-1] + g2*Nm[t-1]
      Ns[t] = NsP*(exp(-1*exp(-5 + dbar[r] + alphaP[r,3]*P[s])))
      #Ns[t] = NsP*(exp(-1*(dbar[r])))
      Nm[t] = NmP*(exp(-1*exp(-5 + alphaP[r,1]*P[s] + alphaO[r,1]*O[s] + b1[r]))) #*(NmP+NlP)
      Nl[t] = NlP*(exp(-1*exp(-5 + alphaP[r,2]*P[s] + alphaO[r,2]*O[s] + b2[r])))  #*(NmP+NlP)
    }
    NeqS[r,s] = Ns[NYEQ]*Pobs[r]
    NeqM[r,s] = Nm[NYEQ]
    NeqL[r,s] = Nl[NYEQ]
  }
}
# 
d4s=data.frame(names=c("No Pred","Both Pred","Otters only","Pycno only"),
              groupN=c("Small"),
              mean=colMeans(NeqS),
              lower=colQuantiles(NeqS, probs = c(0.05)), 
              upper=colQuantiles(NeqS, probs = c(0.95)))
d4m=data.frame(names=c("No Pred","Both Pred","Otters only","Pycno only"),
               groupN=c("Medium"),
               mean=colMeans(NeqM),
               lower=colQuantiles(NeqM, probs = c(0.05)), 
               upper=colQuantiles(NeqM, probs = c(0.95)))
d4l=data.frame(names=c("No Pred","Both Pred","Otters only","Pycno only"),
               groupN=c("Large"),
               mean=colMeans(NeqL),
               lower=colQuantiles(NeqL, probs = c(0.05)), 
               upper=colQuantiles(NeqL, probs = c(0.95)))
d4 = rbind(d4s,d4m,d4l)

pd <- position_dodge(0.1) # move them .05 to the left and right

plt6 = ggplot(d4, aes(x=names, y=mean, colour=groupN, group=groupN)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1, position=pd) +
#  geom_line(position=pd) +
  geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
  xlab("Predation Scenerio") +
  ylab("Equilibrium Abundance") +
  scale_colour_hue(name="Urchin Size Class",    # Legend label, use darker colors
                   breaks=c("Small", "Medium", "Large"),
                   labels=c("Small", "Medium", "Large"),
                   l=40) +                    # Use darker colors, lightness=40
  ggtitle("Equilibrium Urchin Abundance by Size Class") +
  expand_limits(y=0) +                        # Expand y range
  # scale_y_continuous(trans = "log10") +
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(.9,.55))               # Position legend in bottom right
print(plt6)
