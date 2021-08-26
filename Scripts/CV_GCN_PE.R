#
library(RasperGade)
# functions to be added
source("RasperGade_hidden_state_prediction.R")
source("RasperGade_node_tracking.R")
# read in command arguments
arg = commandArgs(TRUE)
# setup parallel running
if(Sys.info()["sysname"]=='Windows'){
  numCores = 1
}else{
  numCores = min(detectCores()-1,as.numeric(arg[1]))
  if(is.na(numCores)) numCores = 1
}
#
numCV = as.numeric(arg[2])
if(is.na(numCV)) numCV = 1
#
by.CV = numCV>10
#
rescaled.data = readRDS("Reference/rescaled_data_model.RDS")
rescaled.tree = makeNodeLabel(rescaled.data$phy)
tree = readRDS("Reference/homogeneous_data.RDS")$phy
#
NSTD.cutoff = 10^seq(-3,0,length.out = 10)[1:9]
true.NSTD.cutoff = c(0,NSTD.cutoff[-1])
#
trait = rescaled.data$dat
#
if(!is.na(arg[3])) set.seed(seed = as.numeric(arg[3]))
test.set.tips = lapply(1:length(true.NSTD.cutoff),function(j){
  tip.random = tree$tip.label[randperm(Ntip(tree),Ntip(tree))]
  tip.partition = unname(split(tip.random,mod(1:Ntip(tree),50)))
  return(tip.partition)
})
saveRDS(test.set.tips,"CV/GCN.testset.RDS")
#
#
PE.CV = mclapply(1:length(true.NSTD.cutoff),function(j){
  this.batch.CV = mclapply(test.set.tips[[j]][1:min(numCV,50)],function(query.tip){
    test.data = partition_test_set_with_minimum_distance(
      phy = tree,trait = trait[tree$tip.label],
      test.set = query.tip,distance = true.NSTD.cutoff[j],location = FALSE)
    this.ref.tree = drop.tip(rescaled.tree,setdiff(rescaled.tree$tip.label,
                                                   test.data$phy$tip.label))
    test.data = partition_test_set_with_minimum_distance(
      phy = this.ref.tree,trait = trait[this.ref.tree$tip.label],
      test.set = query.tip,distance = 0)
    CV.FMR = fullMarginalReconstructionWithPE(
      phy = test.data$ref$phy,x = test.data$ref$dat,
      params = rescaled.data$model,laplace = FALSE,approximate = 1)
    CV.FMR$key = make_traceable_edge_label(test.data$ref$phy)
    this.tree.match = match(c(test.data$ref$phy$tip.label,test.data$ref$phy$node.label),
                            c(rescaled.tree$tip.label,rescaled.tree$node.label))
    CV.FMR$scale = rescaled.data$scale.branch[this.tree.match]*0+1
    CV.FMR$epsilon.branch = rescaled.data$epsilon.branch[this.tree.match]*0
    CV.res = predictHiddenStateWithPE(FMR = CV.FMR,query.keys = test.data$hash)
    return(CV.res)
  },mc.cores = max(numCores*as.numeric(by.CV),1))
  return(this.batch.CV)
},mc.cores = max(numCores*as.numeric(!by.CV),1))
saveRDS(PE.CV,"CV/GCN.PE.CV.RDS")
#