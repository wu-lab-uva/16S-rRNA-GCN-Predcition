#
library(RasperGade)
library(castor)
#
source("RasperGade_node_tracking.R")
#
arg = commandArgs(TRUE)
numCores = as.numeric(arg[1])
if(is.na(numCores)) numCores = 1
# load reference data
cat("Loading reference data\n")
tree = readRDS("Reference/homogeneous_data.RDS")$phy
trait = readRDS("Reference/homogeneous_data.RDS")$dat
# load tips for test sets
test.set.tips = readRDS("CV/GCN.testset.RDS")
#
NSTD.cutoff = 10^seq(-3,0,length.out = 10)[1:9]
true.NSTD.cutoff = c(0,NSTD.cutoff[-1])
# get reference trees for each CV
CV.ref.trees = lapply(1:length(test.set.tips),function(i){
  this.NSTD = test.set.tips[[i]]
  lapply(this.NSTD,function(this.tips){
    this.ref = partition_test_set_with_minimum_distance(phy = tree,trait = trait,test.set = this.tips,
                                                        distance = true.NSTD.cutoff[i],location=FALSE)
    return(c(this.ref$ref,list(insert=this.ref$phy)))
  })
})
# CV with castor
reroot_tree_at_tip = function(phy,tip){
  if(is.character(tip)){
    tip.idx = which(phy$tip.label==tip)
  }else{
    tip.idx = tip
  }
  phy = makeNodeLabel(phy)
  pseudo.tree = rtree(2,br = 1)
  pseudo.tree  = makeNodeLabel(phy = pseudo.tree,prefix = "PseudoNode")
  pseudo.tree$tip.label = paste0("Pseudo_",pseudo.tree$tip.label)
  combined.tree = bind.tree(x = phy,y = pseudo.tree,where = tip.idx,position = 0)
  combined.rt.tree = root.phylo(combined.tree,
                                node = get_mrca_of_set(tree = combined.tree,
                                                       descendants = pseudo.tree$tip.label))
  rt.tree = drop.tip(combined.rt.tree,
                     tip = pseudo.tree$tip.label,collapse.singles = FALSE,trim.internal = FALSE)
  return(rt.tree)
}
#
CV.castor.rt = list()
for(i in 1:length(CV.ref.trees)){
  print(i)
  this.res = mclapply(CV.ref.trees[[i]],function(this.CV){
    query.tip = setdiff(this.CV$insert$tip.label,names(this.CV$dat[!is.na(this.CV$dat)]))
    CV.res = do.call(rbind,lapply(1:length(query.tip),function(j){
      this.tree = drop.tip(phy = this.CV$insert,tip = query.tip[-j])
      this.tree = reroot_tree_at_tip(this.tree,query.tip[j])
      this.asr = asr_independent_contrasts(tree = this.tree,tip_states = this.CV$dat[this.tree$tip.label],
                                           weighted = TRUE,include_CI = TRUE)
      return(data.frame(node=-1,label=query.tip[j],x=this.asr$ancestral_states[1],
                        var=this.asr$standard_errors[1]^2,stringsAsFactors = FALSE))
    }))
    return(CV.res)
  },mc.cores = numCores)
  CV.castor.rt[[i]] = this.res
  saveRDS(CV.castor.rt,"CV/GCN.BM.CV.RDS")
}