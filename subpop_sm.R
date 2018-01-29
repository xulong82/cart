# ssh starmaster
# cd /home/wzheng/GUIDE_testing
#ln -s /home/dcunha/subpopulationInferenceEngine/Test_Sets_Continuous sim_data_cntn
#ln -s /home/dcunha/subpopulationInferenceEngine/Test_Sets_Discrete sim_data_dscr

# /gnshealthcare/software/R/3.2.5/bin/R

rm(list=ls())
setwd("/home/wzheng/GUIDE_testing")
library(gnsutils)

nsim = 25
simtype = c("dscr", "cntn")
for (j in simtype){
  for (i in 1:nsim){
  ## plot true network graph
  # true_network = read.csv(paste0("/home/wzheng/GUIDE_testing/sim_data_", j, "/Test_Sets", i, "/true_network.csv"), header=TRUE)
  # ig = graph.data.frame(true_network, T)
  # pdf(paste0("true_network_plots/true_network_", j, "_", i, ".pdf"))
  # ppPlotIgraph(ig)
  # dev.off()
  
  ## generate dsc file to describe data input file
  inf = paste0("/home/wzheng/GUIDE_testing/sim_data_", j, "/Test_Sets", i, "/causal_network_training_data.csv")
  inf2 = paste0("/home/wzheng/GUIDE_testing/test_", j, "_", i, "_training_data.csv")
  
  system(paste("cp", inf, inf2)) #GUIDE truncate long file path so has to copy it here
  cat(paste0('\"', inf2, '\"', "\n"), file=paste0("test_", j, "_", i, ".dsc"))  
  cat("NA\n", file=paste0("test_", j, "_", i, ".dsc"), append = TRUE)  
  cat("2\n", file=paste0("test_", j, "_", i, ".dsc"), append = TRUE)  
  
  infc = read.csv(inf, nrows=1)
  nsmp = dim(infc)[2]
  tab = data.frame(colid = 1:nsmp, varnm = colnames(infc), vartyp =  rep("n", nsmp))
  tab = as.matrix(tab)
  tab[grep("^D", tab[,2]), 3] = "R"
  tab[grep("^CC$", tab[,2]), 3] = "d"
  write.table(tab, file=paste0("test_", j, "_", i, ".dsc"), append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)

  ## create .in file
  system(paste0("sed 's/template/test_", j, "_", i, "/g' template/template.in >test_", j, "_", i, ".in"))
  ## run GUIDE
  system(paste0("~/bin/guide <test_", j, "_", i, ".in"))
  
  ## convert tree .tex file to .pdf
  system(paste0("/gnshealthcare/software/texlive/bin/x86_64-linux/latex test_", j, "_", i, "_tree"))
  system(paste0("/gnshealthcare/software/texlive/bin/x86_64-linux/dvips test_", j, "_", i, "_tree"))
  system(paste0("/usr/bin/ps2pdf test_", j, "_", i, "_tree.ps"))
  
  }#for i

}#for j

#######################################
## !!! edit test*_predict.R files in bash
#######################################
## edit input file name
## write prediction results in output files
### too many escapes to wrap in R
### for TYPE in dscr cntn; do for VAR in {1..25}; do sed -i.bak '1i setwd("/home/wzheng/GUIDE_testing/")' test_${TYPE}_${VAR}_predict.R; sed -i.bak "s/new.rdata/test_${TYPE}_${VAR}_testing_data.txt/" test_${TYPE}_${VAR}_predict.R; sed -i.bak '2s/$/\n cat(c("ID", \"\\t\", "node", \"\\t\", "rpred", \"\\n\"), file=\"test_'"${TYPE}_${VAR}"'_results.txt\") \n/' test_${TYPE}_${VAR}_predict.R; tac test_${TYPE}_${VAR}_predict.R | sed '2s/}/cat(c(i, \"\\t\", node, \"\\t\", rpred, \"\\n\"), file=\"test_'"${TYPE}_${VAR}"'_results.txt\", append = TRUE)}/' | tac >test_${TYPE}_${VAR}_predict_new.R; done; done


###########################################
## check out of sample prediction accuracy
###########################################

rsquare = rep(0, length(simtype) * nsim)
for (j in simtype){
  for (i in 1:nsim){
    ## predict new data
    newdf = paste0("/home/wzheng/GUIDE_testing/sim_data_", j, "/Test_Sets", i, "/causal_network_testing_data.csv")
    testdat = read.csv(newdf, colClasses = "character")
    write.table(testdat, file=paste0("test_", j, "_", i, "_testing_data.txt"), row.names=FALSE, col.names = TRUE, sep="\t")
    system(paste0("R --vanilla <test_", j, "_", i, "_predict_new.R")) # run prediction R code
    pred.res = read.table(paste0("test_", j, "_", i, "_results.txt"), header=TRUE, sep="\t")
    jj = as.numeric(j=="cntn")
    if(is.numeric(pred.res$rpred)){
      rsquare[jj*nsim + i] = cor(pred.res$rpred, as.numeric(testdat$CC))^2
    }else{ #problem in predicting new test data
      rsquare[jj*nsim + i] = "NA"
    }

  }#for i
}#for j
round(rsquare, 2)
# [1] 0.00 0.00 0.01 0.00 0.01 0.17 0.17 0.15 0.15 0.18 0.38 0.38 0.39 0.42 0.34
# [16] 0.59 0.59 0.61 0.60 0.62 0.75 0.81 0.79 0.80 0.80 0.01 0.00 0.00 0.00 0.00
# [31] 0.19 0.23 0.21 0.23 0.20 0.41 0.38 0.36 0.36 0.38 0.57 0.60 0.53 0.55 0.56
# [46] 0.76 0.78 0.82 0.77 0.81


###########################################
## check split var selection
###########################################
G1sel = GPsel = rep(0, length(simtype) * nsim)
testid = splitvar = rep("", length(simtype) * nsim)
for (j in simtype){
  for (i in 1:nsim){
    ## read in split var variable
    newdf = paste0("/home/wzheng/GUIDE_testing/test_", j, "_", i, "_split_var.txt")
    testdat = read.table(newdf, colClasses = "character", nrows = 1)
    jj = as.numeric(j=="cntn")
    if(testdat[1,3] == "G1_L1"){
      G1sel[jj*nsim + i] = GPsel[jj*nsim + i] = 1
    }else if(testdat[1,3] == "P0"){
      GPsel[jj*nsim + i] = 1
    }
    testid[jj*nsim + i] = paste0("test_", j, "_", i)
    splitvar[jj*nsim + i] = paste(testdat[1,2:3], collapse=";")
  }#for i
}#for j
G1sel
# [1] 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# [39] 0 1 0 0 1 0 0 0 0 0 1 1
GPsel
# [1] 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0 0 0 1 1 0 1 1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 1
# [39] 0 1 1 0 1 0 1 1 0 0 1 1
sum(G1sel)
# [1] 10
sum(GPsel)
# [1] 18
write.table(data.frame(testid, splitvar, G1sel, GPsel, rsquare), file="GUIDE_preformance.txt", sep="\t")
###########################################
## check classification accuracy, MCC for discrete interaction (local)
###########################################
## starmaster:
 # mv *results.txt out_of_sample_prediction/
 # mv *fitted.txt training_data_recovery/
## local:
 # scp -r starmaster:~/GUIDE_testing/out_of_sample_prediction /Users/wzheng/Documents/method/subpop_analysis/GUIDE_testing/
 # scp -r starmaster:~/GUIDE_testing/training_data_recovery /Users/wzheng/Documents/method/subpop_analysis/GUIDE_testing/

### local environment setup:
rm(list=ls())
setwd("/Users/wzheng/Documents/method/subpop_analysis/GUIDE_testing/")
library(ggplot2)
library(reshape2)
nsim = 25

mcc <- function (actual, predicted)
{
  # handles zero denominator and overflow error on large-ish products in denominator.
  #
  # actual = vector of true outcomes, 1 = Positive, 0 = Negative
  # predicted = vector of predicted outcomes, 1 = Positive, 0 = Negative
  # function returns MCC
  
  TP <- sum(actual == 1 & predicted == 1)
  TN <- sum(actual == 0 & predicted == 0)
  FP <- sum(actual == 0 & predicted == 1)
  FN <- sum(actual == 1 & predicted == 0)
  
  sum1 <- TP+FP; sum2 <-TP+FN ; sum3 <-TN+FP ; sum4 <- TN+FN;
  denom <- as.double(sum1)*sum2*sum3*sum4 # as.double to avoid overflow error on large products
  
  if (any(sum1==0, sum2==0, sum3==0, sum4==0)) {
    denom <- 1
  }
  
  mcc <- ((TP*TN)-(FP*FN)) / sqrt(denom)
  return(mcc)
}

omcc = imcc = rep(0, nsim)

for (i in 1:nsim){
    ## if there are more than 1 terminal node, use the node overlapping most with actual.group as group 1, and the rest as one group
  
    ## read in oos predicted grouping
    grp = read.csv(paste0("sim_data_dscr/Test_Sets", i, "/subgroup_identifier_testing_data.csv"))
    actual.grp = as.integer(grp==0.8) 
    prd.grp = read.table(paste0("out_of_sample_prediction/test_dscr_", i, "_results.txt"), header = TRUE)
    pg = unique(prd.grp$node)
    if(length(pg) > 1) {
      pgid = as.matrix(data.frame(lapply(pg, function(X) as.integer(prd.grp$node == X))))
      concord = apply(pgid, 2, function(X)sum(X==actual.grp))
      g1 = which(concord == max(concord))   
      pred.grp = as.integer(prd.grp$node == pg[g1])
     # table(actual.grp, pred.grp)
      omcc[i] = mcc(actual.ogrp, pred.grp)  
      }#if
    
    ## read in in-sample predicted grouping
    grp = read.csv(paste0("sim_data_dscr/Test_Sets", i, "/subgroup_identifier_training_data.csv"))
    actual.grp = as.integer(grp==0.8) 
    prd.grp = read.table(paste0("training_data_recovery/test_dscr_", i, "_fitted.txt"), header = TRUE)
    pg = unique(prd.grp$node)
    if(length(pg) > 1) {
      pgid = as.matrix(data.frame(lapply(pg, function(X) as.integer(prd.grp$node == X))))
      concord = apply(pgid, 2, function(X)sum(X==actual.grp))
      g1 = which(concord == max(concord))   
      pred.grp = as.integer(prd.grp$node == pg[g1])
      # table(actual.grp, pred.grp)
      imcc[i] = mcc(actual.igrp, pred.grp)  
    }#if
    
  }#for i
round(omcc,2)
# [1] 0.00 0.06 0.02 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.06 0.00 0.00 0.00 0.00 0.00 0.10
# [22] 0.00 0.00 0.82 0.68
round(imcc,2)
# [1] 0.00 0.02 0.03 0.09 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.03 0.00 0.00 0.00 0.00 0.00 0.01
# [22] 0.00 0.00 0.85 0.71
write.csv(cbind(imcc, omcc), "MCC_results.csv")
###########################################
## check subpop effect on CC (local)
###########################################
for (i in 1:nsim){
## read in latent grouping variable
  grp = read.csv(paste0("/Users/wzheng/Documents/method/subpop_analysis/GUIDE_testing/sim_data_dscr/Test_Sets", i, "/subgroup_identifier_training_data.csv"))
  id1 = which(grp == 0.8)
## read in  trainig data 
  train = read.csv(paste0("/Users/wzheng/Documents/method/subpop_analysis/GUIDE_testing/sim_data_dscr/Test_Sets", i, "/causal_network_training_data.csv"))
  pdf(paste0("test_", i, "_boxplot_by_groups.pdf"))
  par(mfrow=c(1, 2))
  boxplot(train$CC[id1] ~ train$D0[id1], xlab="D0", ylab="CC", main="group1")
  boxplot(train$CC[-id1] ~ train$D0[-id1], xlab="D0", ylab="CC", main="group2")
  dev.off()
  
  pdf(paste0("test_", i, "_density_plot.pdf"))
  ggplot(train[id1,], aes(x=CC, fill= as.factor(D0))) + geom_density(alpha=0.25)
  ggplot(train[-id1,], aes(x=CC, fill= as.factor(D0))) + geom_density(alpha=0.25)
  dev.off()
  
  
  pdf(paste0("test_", i, "_residual_plot.pdf"), width=12, height=6)
  lm1 = lm(CC ~ D0, data = train)
  train$RES = lm1$residuals
  id1 = (train$D0==1)
  par(mfrow=c(1, 2))
  plot(train$G1_L1[id1], train$RES[id1], xlab="G1_L1", ylim=range(train$RES), ylab="Residuals", main = "D0=1", col= train$D0[id1])
  lines(lowess(train$G1_L1[id1], train$RES[id1]), col= "blue")
  plot(train$G1_L1[!id1], train$RES[!id1], xlab="G1_L1", ylim=range(train$RES), ylab="Residuals", main = "D0 = -1", col= 1-train$D0[!id1])
  lines(lowess(train$G1_L1[!id1], train$RES[!id1]), col= "blue")
  dev.off()
  
}#for i

nsim = 25
pp = matrix(0, nsim, 4)
simtype = c("dscr", "cntn")
interm = c("D0", "G0_L1")
for (k in 1:2){
  for (j in 1:2){
    for (i in 1:nsim){
      ## read in  trainig data 
      train = read.csv(paste0("/Users/wzheng/Documents/method/subpop_analysis/GUIDE_testing/sim_data_", simtype[j], "/Test_Sets", i, "/causal_network_training_data.csv"))
      
      m = lm(as.formula(paste("CC ~ D0 + G1_L1 + ", interm[k], ":G1_L1")), data=train)
      if (j==1) col = k else col = k+2
      pp[i, col] = summary(m)$coefficients[4,4]
      
    }#i
  }#j
}#k
colnames(pp) = paste(rep(simtype, each=2), rep(interm,2), sep="_")
pp

pp<0.05
