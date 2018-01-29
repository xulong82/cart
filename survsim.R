### implement GUIDE simulations
rm(list=ls())
#local
#dir.create("~/Documents/method/subpop_analysis/GUIDE_sim")
#setwd("~/Documents/method/subpop_analysis/GUIDE_sim")

#starmaster
setwd("/users/wzheng/GUIDE_testing/GUIDE_sim/")
#setwd("~/Documents/method/subpop_analysis/GUIDE_sim/")
dir.create("null_models")
### null models for assessing selection bias:
smpsz = 100
iter = 2500
for (i in 1:iter){
  Y = rbinom(smpsz, 1, 0.5)
  Z = rbinom(smpsz, 1, 0.5)
  X1 = rnorm(smpsz)
  X2 = sample(1:4, smpsz, replace = TRUE, prob = rep(1/4, 4)) #ordinal 4
  X3 = sample(c('a','b','c'), smpsz, replace = TRUE, prob = rep(1/3, 3)) # cat3
  X4 = sample(c('a', 'b', 'c', 'd', 'e', 'f', 'g'), smpsz, replace = TRUE, prob = rep(1/7, 7)) #cat7
write.table(cbind(Y, Z, X1, X2, X3, X4), paste0("null_models/null_", i, ".txt"))  
}

### models to test accuracy 
dir.create("predictive_models")
Z = rbinom(smpsz, 1, 0.5)
for (i in 1:iter){
  Xmat = matrix(0, smpsz, 100)
  for (j in 1:100){
    if (j > 2){
      pi = rbeta(1, 2, 3)
      Xmat[,j] = sample(c(0, 1, 2), smpsz, replace = TRUE, prob = c((1-pi)^2, 2*pi*(1-pi), pi^2))
    }else{
      Xmat[,j] = sample(c(0, 1, 2), smpsz, replace = TRUE, prob = c(0.4, 0.465, 0.135))
    }
  }
  #M1: X1 and X2 are predictive
  Y1 = 0.4 + 0.05 * (Z==1) * (4*(Xmat[,1]!=0) + 3*(Xmat[,2]!=0) + (Xmat[,1]!=0 & Xmat[,2]!=0))
  
  #M2: X1 and X2 are predictive, X3 and X4 are prognostic
  Y2 = 0.3 + 0.2 * ( (2*(Z==1) - 1) * (Xmat[,1]!=0 & Xmat[,2]!=0) + (Xmat[,3]!=0) + (Xmat[,4]!=0))
  
  #M3: X1 and X2 are prognostic
  Y3 = 0.5 + 0.1 * ( 2*( (Z==1) + (Xmat[,3]!=0) + (Xmat[,4]!=0)) - 3 )
write.table(cbind(Z, Xmat, Y1, Y2, Y3), paste0("predictive_models/pred_", i, ".txt"))  
}

### smulate survival data that has similar distn with OAK
# Estimates of the number of events required to demonstrate efficacy with regard to OS were based on the following assumptions:
#   • Event times exponentially distributed
# • A 7.5% 24-month dropout rate assumed for both treatment arms
# • Greater than 95% power for the primary analysis of OS in the ITT and TC1/2/3 or IC1/2/3 of PP and > 80% power for the secondary analysis of OS in the ITT, TC1/2/3 or IC1/2/3, TC2/3 or IC2/3, TC3 or IC3 of SP (detailed HR and minimum detectable difference [MDD] are shown in Table 9 and Table 10)
# • Median survival of 10 months in the docetaxel arm for the ITT and PD-L1 subgroups
# • 65% prevalence rate for TC1/2/3 or IC1/2/3

# p118 table 9 shows that we need 850 patients and 595 events to reach HR = 0.73 (median OS from 10 months in control arm to 13.7 months in treatment arm)

rm(list=ls())
library(igraph)
#setwd("~/Documents/method/subpop_analysis/GUIDE_sim/")
setwd("/users/wzheng/GUIDE_testing/GUIDE_sim/")
dir.create("causal_models")

## simulate censor time so that exp(-l * 24) = 0.075,  7.5% 24-month dropout rate assumed for both treatment arms
smpsz = 500
l = - log(0.075) /24
maxfu = 40

dcd = "/users/dcunha/subpopulationInferenceEngine/Test_Sets/"

for (testset in 1:25){
  for (sim in c("Discrete", "Continuous")){
    tc = pmin(rexp(smpsz, rate = l), maxfu)
    ## simulate survival time from causal network Dan generated
    md <- read.csv(paste0(dcd, "Test_Sets_", sim, "/Test_Sets", testset, "/causal_network_training_data.csv"))
    mdt <- read.csv(paste0(dcd, "Test_Sets_", sim, "/Test_Sets", testset, "/causal_network_testing_data.csv"))
    ms <- read.csv( paste0(dcd, "Test_Sets_", sim, "/Test_Sets", testset, "/true_network.csv"))
    ms = graph_from_data_frame(ms, directed = TRUE)
    D0l = shortest_paths(ms, V(ms)["D0"], V(ms)["CC"])$vpath    ## find the variable on the causal pathway of D0->CC
    lastlayer <- grep("_L4", names(V(ms)), value = TRUE)
    sm.pred <- setdiff(lastlayer, names(unlist(D0l))) # variable on the last layer but not on the causal pathway should have a small coefficient
    lg.pred <- setdiff(lastlayer, sm.pred) # variable on the last layer AND on the causal pathway should have a large coefficient
    
    y <- sapply(list(md, mdt), function(X){-log(log(2)/10) + 1.2 * X[,lg.pred] + 0.1 * apply(X[,sm.pred], 1, sum)}) 
    tevt = apply(y, 2, function(X)rexp(smpsz, rate = exp(-1 * X)))
    md$surv.time = pmin(tevt[,1], tc)
    md$censor = as.integer(tevt[,1] < tc)
    mdt$surv.time = pmin(tevt[,2], tc)
    mdt$censor = as.integer(tevt[,2] < tc) # censored 0, not censored 1
    
    write.csv(md, paste0("causal_models/test", testset, "_", sim, "_training_data.csv"))
    write.csv(mdt, paste0("causal_models/test", testset, "_", sim, "_testing_data.csv"))
    
  }#dscr cntn
}# testset 1-25

## results saved in /users/wzheng/GUIDE_testing/GUIDE_sim/causal_models
