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
  sim.meta = readRDS("Sim/sim_2env_5otu_2f_Beta_meta.RDS")
}else{
  sim.meta = readRDS(arg[3])
}

if(is.na(arg[2])){
  CV.community.sim = readRDS("Sim/sim_2env_5otu_2f_beta_data.RDS")
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
if(is.na(savename)) savename = "Sim/Sim_5_2"
#
t_and_wilcox_with_GCN.2sample = function(abundance.table,error,metadata,focal,
                                         method="t.test",n=1000,numCores=1){
  sample.idx = lapply(unique(metadata),function(x){
    which(metadata==x)
  })
  sample.count = sapply(sample.idx,length)
  max.comb = choose(n = length(metadata),k = min(sample.count))
  if(max.comb<n) n = max.comb
  if(max.comb>1e5){
    perm.sample = lapply(1:n,function(i){
      sample(1:length(metadata),size = min(sample.count),replace = FALSE)
    })
  }else{
    perm.sample = sample(combn(1:length(metadata),min(sample.count),simplify = FALSE),
                         size = n,replace = FALSE)
  }
  perm.GCN = sapply(1:n,function(nn){
    this.GCN = sapply(error,function(x){
      rMixNormal(n = 1,mean = x$x,sd = sqrt(x$var),probs = x$probs)
      })
    this.GCN[this.GCN<=1] = 1
    return(this.GCN)
  })
  GCN = round(apply(perm.GCN,1,median))
  correct.table = correct_abundance_table(abundance.table,GCN)
  correct.test = lapply(focal,function(i){
    do.call(method,list(correct.table[sample.idx[[1]],i],
                        correct.table[sample.idx[[2]],i]))
    })
  perm.test = do.call(cbind,mclapply(1:dim(perm.GCN)[2],function(k){
    this.GCN = round(perm.GCN[,k])
    this.table = correct_abundance_table(abundance.table,this.GCN)
    this.test = sapply(focal,function(i){
      do.call(method,list(correct.table[perm.sample[[k]],i],
                          correct.table[-perm.sample[[k]],i]))$statistic})
  },mc.cores = numCores))
  correct.p.values = sapply(1:length(focal),function(i){
    correct.stat = correct.test[[i]]$statistic
    baseline.stat = perm.test[i,]
    this.p = (sum(correct.stat>baseline.stat)+
      0.5*sum(correct.stat==baseline.stat))/length(baseline.stat)
    return(2*min(this.p,1-this.p))
  })
  return(list(test=correct.test,correct = list(p.value=correct.p.values)))
}
#
sim.diff = list()
#
for(i in 1:length(CV.community.sim[1])){
  print(i)
  res.summarized = CV.res.summarized[[i]]
  this.res = mclapply(1:length(CV.community.sim[[i]]),function(j){
    this.batch = CV.community.sim[[i]][[j]]
    this.table = CV.sim.table[[i]][[j]]
    this.pred.idx = match(colnames(this.table$gene),res.summarized[[this.batch[[1]]$CV.idx]]$hsp$label)
    pred.error = res.summarized[[this.batch[[1]]$CV.idx]]$error[this.pred.idx]
    pred.GCN = res.summarized[[this.batch[[1]]$CV.idx]]$hsp$x[this.pred.idx]
    cell.test = sapply(1:dim(this.table$cell)[2],function(i){
      wilcox.test(x = this.table$cell[1:10,i],y = this.table$cell[11:20,i])$p.value
    })
    gene.test = sapply(1:dim(this.table$gene)[2],function(i){
      wilcox.test(x = this.table$gene[1:10,i],y = this.table$gene[11:20,i])$p.value
    })
    cell.fold.diff = apply(this.table$cell[1:10,],2,mean)/
      apply(this.table$cell[11:20,],2,mean)
    gene.fold.diff = apply(this.table$gene[1:10,],2,mean)/
      apply(this.table$gene[11:20,],2,mean)
    correct.test = list()
    #  t_and_wilcox_with_GCN.2sample(abundance.table = this.table$gene,error = pred.error,
    #                                metadata = sim.meta$V1,focal = 1:dim(this.table$gene)[2],
    #                                method = "wilcox.test",n = 1000,numCores = 1)
    return(list(cell=list(diff = cell.fold.diff, test = cell.test),
                gene=list(diff = gene.fold.diff, test = gene.test),correct=correct.test))
  },mc.cores = numCores)
  sim.diff[[i]] = this.res
  saveRDS(sim.diff,paste0(savename,"_diff.RDS"))
}
