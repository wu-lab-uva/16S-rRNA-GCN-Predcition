library(randomForest)
# load newly added functions
source("added_functions.R")
# read in command arguments
arg = commandArgs(TRUE)
# setup parallel running
if(Sys.info()["sysname"]=='Windows'){
  numCores = 1
}else{
  numCores = min(detectCores()-1,as.numeric(arg[1]))
  if(is.na(numCores)) numCores = 1
}
# load data
cat("Loading reference data\n")
all.trait = readRDS("Reference/homogeneous_data.RDS")$dat
CV.res.summarized = readRDS("CV/GCN.PE.CV.RDS")
#
#
if(is.na(arg[3])){
  sim.meta = readRDS("Sim/sim_2env_30otu_2f_beta_meta.RDS")
}else{
  sim.meta = readRDS(arg[3])
}

if(is.na(arg[2])){
  CV.community.sim = readRDS("Sim/sim_2env_30otu_2f_beta_data.RDS")
}else{
  CV.community.sim = readRDS(arg[2])
}
#
cat("Data loaded\n")
#
CV.sim.table = lapply(1:length(CV.community.sim),function(i){
  mclapply(CV.community.sim[[i]],function(this.batch){
    gene.table=as.data.frame(data.table::rbindlist(lapply(this.batch,function(x){
      data.frame(t(x$gene/sum(x$gene)))
    }),fill=TRUE))
    cell.table=as.data.frame(data.table::rbindlist(lapply(this.batch,function(x){
      data.frame(t(x$cell/sum(x$cell)))
    }),fill=TRUE))
    gene.table[is.na(gene.table)] = 0
    cell.table[is.na(cell.table)] = 0
    return(list(gene=gene.table,cell=cell.table))
  },mc.cores = numCores)
})
#
cat("Table generated\n")
#
savename = arg[4]
if(is.na(savename)) savename = "Sim/Sim_30_2"
#
CV.sim.rf = list()
for(i in 1:length(CV.community.sim[1])){
  res.summarized = CV.res.summarized[[i]]
  this.res = mclapply(1:length(CV.community.sim[[i]]),function(j){
    this.batch = CV.community.sim[[i]][[j]]
    this.table = CV.sim.table[[i]][[j]]
    this.pred.idx = match(colnames(this.table$gene),res.summarized[[this.batch[[1]]$CV.idx]]$hsp$label)
    pred.error = res.summarized[[this.batch[[1]]$CV.idx]]$error[this.pred.idx]
    pred.GCN = res.summarized[[this.batch[[1]]$CV.idx]]$hsp$x[this.pred.idx]
    gene.rf = randomforest_cross_validation(predictor = as.matrix(this.table$gene),
                                            response = sim.meta$V1)
    cell.rf = randomforest_cross_validation(predictor = as.matrix(this.table$cell),
                                            response = sim.meta$V1)
    correct.rf = correct_randomforest_with_confidence(abundance.table = this.table$gene,error = pred.error,
                                                      metadata = sim.meta,n = 100)
    return(list(gene=gene.rf,cell=cell.rf,correct = correct.rf))
  },mc.cores = numCores)
  CV.sim.rf[[i]] = this.res
  saveRDS(CV.sim.rf,paste0(savename,"rf_test.RDS"))
}
