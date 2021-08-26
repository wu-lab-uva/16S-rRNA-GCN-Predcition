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
CV.PE = readRDS("CV/GCN.PE.CV.RDS")
#
if(is.na(arg[2])){
  CV.community.sim = readRDS("Sim/sim_1env_20otu_1f_RA_data.RDS")
}else{
  CV.community.sim = readRDS(arg[2])
}
#
cat("Data loaded\n")
#
savename = arg[3]
if(is.na(savename)) savename = "Sim/Sim_1_1RA.RDS"
#
CV.correct.RA = list()
for(i in 1:length(CV.community.sim)){
  print(i)
  res.summarized = CV.PE[[i]]
  ra.res = mclapply(CV.community.sim[[i]],function(this.batch){
    lapply(this.batch,function(x){
      this.pred.idx = match(names(x$gene),res.summarized[[x$CV.idx]]$hsp$label)
      pred.error = res.summarized[[x$CV.idx]]$error[this.pred.idx]
      pred.GCN = res.summarized[[x$CV.idx]]$hsp$x[this.pred.idx]
      correct.RA = correct_RAD_with_confidence(abundance = x$gene,
                                               error = pred.error,
                                               n = 1000)
      return(c(list(true = relative_abundance(x$cell),raw = relative_abundance(x$gene)),correct.RA))
    })
  },mc.cores = numCores)
  CV.correct.RA[[i]] = ra.res
  saveRDS(CV.correct.RA,savename)
}
#

