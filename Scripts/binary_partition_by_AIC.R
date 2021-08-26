#
library(RasperGade)
source("RasperGade16S_partition.R")
# read in key arguments
arg=commandArgs(TRUE)
#
data.file = as.character(arg[1])
if(is.na(data.file)) data.file = "Reference/homogeneous_data.RDS"
model.file = as.character(arg[2])
if(is.na(model.file)) model.file = "Reference/homogeneous_models.RDS"
save.name = as.character(arg[3])
if(is.na(save.name)) save.name = "Reference/binary_partition.RDS"
#
res.name = as.character(arg[4])
if(is.na(res.name)) res.name = "Reference/rescaled_data_model.RDS"
#
numCores = as.numeric(arg[5])
if(is.na(numCores)) numCores = 1

repRange = 10
ini.multiplier = expand.grid(lambda=c(1,repRange,1/repRange),
                             size=c(1,repRange,1/repRange),
                             KEEP.OUT.ATTRS = TRUE)
#
tree = readRDS(data.file)$phy
trait = readRDS(data.file)$dat
#
homogeneous.model = readRDS(model.file)$params$pee
#
cat("Starting binary partition\n")
#
iteration = 1
res = list()
res[[iteration]] = list(model = list(regular = homogeneous.model,
                                     slow = homogeneous.model*c(0.01,1,1,1)),
                        partition=list())
#
reg.ancs = 
  reconstructAncestralStates(phy = tree,x = trait[tree$tip.label],
                            rate = total.process.variance(res[[iteration ]]$model$regular),
                            epsilon = res[[iteration ]]$model$regular["epsilon"])
slow.ancs = 
  reconstructAncestralStates(phy = tree,x = trait[tree$tip.label],
                             rate = total.process.variance(res[[iteration ]]$model$slow),
                             epsilon = res[[iteration ]]$model$slow["epsilon"])
reg.AIC = sapply(1:Nnode(tree),function(i){
  -2*log(dPEpoisnorm(x = reg.ancs$contrast[i],t = reg.ancs$l[i],
                  lambda=unname(res[[iteration ]]$model$regular[1]),
                  size=unname(res[[iteration ]]$model$regular[2]),
                  sigma=unname(res[[iteration ]]$model$regular[3]),
                  epsilon=unname(reg.ancs$epsilon[i])))
})
slow.AIC = sapply(1:Nnode(tree),function(i){
  -2*log(dPEpoisnorm(x = slow.ancs$contrast[i],t = slow.ancs$l[i],
                  lambda=unname(res[[iteration ]]$model$slow[1]),
                  size=unname(res[[iteration ]]$model$slow[2]),
                  sigma=unname(res[[iteration ]]$model$slow[3]),
                  epsilon=unname(slow.ancs$epsilon[i])))
})
#
dAIC = reg.AIC-slow.AIC
total.reg.AIC = 0
total.slow.AIC = -sum(dAIC)
part.flip = partition_phylogeny_by_AIC_single(tree = tree,dAIC = dAIC,flip = TRUE)
refined.part = refine_partition(tree = tree,dAIC = dAIC,
                                       group = part.flip$group)
res[[iteration]]$partition = refined.part$group
saveRDS(res,save.name)
warnings()
# Fit model and partition tree with the newly fitted model
new.reg.model = mclapply(1:dim(ini.multiplier)[1],function(i){
  fitPE.phy(phy = tree,x = trait[tree$tip.label],
            start.value = initialize.parameters.SP(c(lambda=ini.multiplier$lambda[i],
                                                     size=ini.multiplier$size[i],
                                                     sigma=1e-6,epsilon=1)),
            fixed=list(sigma=log(1e-6)),
            min.value = c(0,0,0,1e-6),
            laplace = FALSE,bootstrap = 
              which(refined.part$group[Ntip(tree)+(1:Nnode(tree))]>0))
},mc.cores = numCores)
new.reg.model = new.reg.model[[which.min(sapply(new.reg.model,get.AIC))]]
#
new.slow.model = mclapply(1:dim(ini.multiplier)[1],function(i){
  fitPE.phy(phy = tree,x = trait[tree$tip.label],
            start.value = initialize.parameters.SP(c(lambda=ini.multiplier$lambda[i],
                                                     size=ini.multiplier$size[i],
                                                     sigma=1e-6,epsilon=1)),
            fixed=list(size=unname(log(new.reg.model$params[2])),sigma=log(1e-6)),
            min.value = c(0,0,0,1e-6),
            laplace = FALSE,bootstrap = 
              which(refined.part$group[Ntip(tree)+(1:Nnode(tree))]<0))
},mc.cores = numCores)
new.slow.model = new.slow.model[[which.min(sapply(new.slow.model,get.AIC))]]
#
warnings()
iteration = iteration + 1
res[[iteration]] = list(model = list(regular = new.reg.model$params,
                                     slow = new.slow.model$params),
                        partition=list())
print(res[[iteration]]$model)
#
#
reg.ancs = 
  reconstructAncestralStates(phy = tree,x = trait[tree$tip.label],
                             rate = total.process.variance(res[[iteration ]]$model$regular),
                             epsilon = res[[iteration ]]$model$regular["epsilon"])
slow.ancs = 
  reconstructAncestralStates(phy = tree,x = trait[tree$tip.label],
                             rate = total.process.variance(res[[iteration ]]$model$slow),
                             epsilon = res[[iteration ]]$model$slow["epsilon"])
reg.AIC = sapply(1:Nnode(tree),function(i){
  -2*log(dPEpoisnorm(x = reg.ancs$contrast[i],t = reg.ancs$l[i],
                     lambda=unname(res[[iteration ]]$model$regular[1]),
                     size=unname(res[[iteration ]]$model$regular[2]),
                     sigma=unname(res[[iteration ]]$model$regular[3]),
                     epsilon=unname(reg.ancs$epsilon[i])))
})
slow.AIC = sapply(1:Nnode(tree),function(i){
  -2*log(dPEpoisnorm(x = slow.ancs$contrast[i],t = slow.ancs$l[i],
                     lambda=unname(res[[iteration ]]$model$slow[1]),
                     size=unname(res[[iteration ]]$model$slow[2]),
                     sigma=unname(res[[iteration ]]$model$slow[3]),
                     epsilon=unname(slow.ancs$epsilon[i])))
})
#
dAIC = reg.AIC-slow.AIC
refined.part = refine_partition(tree = tree,dAIC = dAIC,
                                group = res[[iteration-1]]$partition)
res[[iteration]]$partition = refined.part$group
saveRDS(res,save.name)
warnings()
#
rescaled.tree = 
  rescale_tree_by_partition_and_model(tree = tree,
    group = as.numeric(refined.part$group<0)+1,
    models = list(new.reg.model$params,new.slow.model$params),epsilon.as.PE = TRUE)
rescaled.tree$dat = trait[rescaled.tree$phy$tip.label]
rescaled.tree$original = tree
rescaled.tree$parttion = refined.part$group
saveRDS(rescaled.tree,res.name)
#






