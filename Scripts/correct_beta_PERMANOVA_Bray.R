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
CV.sim.beta = list()
CV.sim.PERMANOVA = list()
for(i in 1:length(CV.community.sim[1])){
  print(i)
  res.summarized = CV.res.summarized[[i]]
  this.res = mclapply(1:length(CV.community.sim[[i]]),function(j){
    this.batch = CV.community.sim[[i]][[j]]
    this.table = CV.sim.table[[i]][[j]]
    this.pred.idx = match(colnames(this.table$gene),res.summarized[[this.batch[[1]]$CV.idx]]$hsp$label)
    pred.error = res.summarized[[this.batch[[1]]$CV.idx]]$error[this.pred.idx]
    pred.GCN = res.summarized[[this.batch[[1]]$CV.idx]]$hsp$x[this.pred.idx]
    gene.dist = vegdist(this.table$gene)
    gene.PERMANOVA = adonis(formula = gene.dist~V1,data = sim.meta,permutations = 999)
    cell.dist = vegdist(this.table$cell)
    cell.PERMANOVA = adonis(formula = cell.dist~V1,data = sim.meta,permutations = 999)
    correct.PERMANOVA = correct_PERMANOVA_with_confidence(abundance.table = this.table$gene,error = pred.error,
                                                          method="bray",
                                                      metadata = sim.meta,n = 1000)
    return(list(gene=gene.dist,cell=cell.dist,correct=correct.PERMANOVA$dist,CI=correct.PERMANOVA$CI,
                gene.PERMANOVA=gene.PERMANOVA$aov.tab$R2,cell.PERMANOVA = cell.PERMANOVA$aov.tab$R2,
                R2 = correct.PERMANOVA$correct$aov.tab$R2,R2.CI = correct.PERMANOVA$R2.CI))
  },mc.cores = numCores)
  CV.sim.beta[[i]] = lapply(this.res,function(x){x[1:4]})
  CV.sim.PERMANOVA[[i]] = lapply(this.res,function(x){x[5:8]})
  saveRDS(CV.sim.beta,paste0(savename,"beta.RDS"))
  saveRDS(CV.sim.PERMANOVA,paste0(savename,"PERMANOVA.RDS"))
}
