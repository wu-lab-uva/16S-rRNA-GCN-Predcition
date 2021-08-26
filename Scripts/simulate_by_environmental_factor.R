#
library(sads)
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
CV.PE = readRDS("CV/GCN.PE.CV.RDS")
#
cat("Data loaded\n")
# simulate communities for analyses
Sn = 100
Nn = 2000
numCategory = as.numeric(arg[2])
if(is.na(numCategory)|numCategory<1) numCategory = 1
effect.OTU = as.numeric(arg[3])
if(is.na(effect.OTU)) effect.OTU = 30
effect.size = as.numeric(arg[4])
if(is.na(effect.size)) effect.size = 2
numRep = as.numeric(arg[5])
if(is.na(numRep)) numRep = 10
#
OTU.pool.idx = 1:length(CV.PE[[1]])
sim.category = rep(1:numCategory,each=numRep)
sim.meta = data.frame(V1=sim.category)
CV.community.sim = lapply(CV.PE,function(res.summarized){
  mclapply(1:length(OTU.pool.idx),function(i){
    this.CV = all.trait[res.summarized[[OTU.pool.idx[i]]]$hsp$label]
    this.CV.group = numeric(length(this.CV))*0
    for(this.category in unique(sim.category)){
      this.CV.group[this.CV.group<=0][1:effect.OTU] = this.category
    }
    this.CV.group = this.CV.group[randperm(a = length(this.CV),k = length(this.CV))]
    names(this.CV.group) = names(this.CV)
    lapply(1:(numRep*numCategory),function(j){
      this.category = sim.category[j]
      OTUs = sample(x = this.CV,size = min(Sn,length(this.CV)),replace = FALSE)
      OTUs.group = this.CV.group[names(OTUs)]
      OTUs.max = exp(OTUs*0)
      OTUs.max[names(which(OTUs.group==this.category))] = effect.size
      OTUs.weight = runif(n = length(OTUs),min = 0,max = OTUs.max)
      names(OTUs.weight) = names(OTUs)
      cell.abundance = sort(rls(n = min(Sn,length(this.CV)),N = Nn,alpha = 20),decreasing = TRUE)
      names(cell.abundance) = names(sort(OTUs.weight,decreasing = TRUE))
      gene.abundance = OTUs[names(cell.abundance)]*cell.abundance
      return(list(cell=cell.abundance,gene=gene.abundance,GCN=OTUs,
                  category = this.category,group=OTUs.group[names(cell.abundance)],CV.idx=OTU.pool.idx[i]))
    })
  },mc.cores = numCores)
})
cat("Data simulated\n")
#
savename = arg[6]
if(is.na(savename)) savename = "CV_sim_by_env_factor"
saveRDS(CV.community.sim,paste0(savename,"_data.RDS"))
saveRDS(sim.meta,paste0(savename,"_meta.RDS"))
