#
library(sads)
# functions to be added
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
# Read in HQ metadata
HMP.meta = readRDS("HMP_v13_batch00/HMP_v13_HQ.meta.RDS")
# Read in predicted GCN
pred.error = readRDS("HMP_v13_batch00/HMP_v13.pred.GCN.RDS.RDS")
pred.GCN = pred.error$discrete$x
names(pred.GCN) = pred.error$discrete$label
pred.continuous = pred.error$contiguous$hsp$x
names(pred.continuous) = pred.error$contiguous$hsp$label
pred.error = pred.error$contiguous$error
#
HMP.abundance = readRDS("HMP_v13_batch00/HMP_v13_HQ.abundance.RDS")
#
HMP.RAD = mclapply(HMP.abundance,function(this.abundance){
  this.error = pred.error[match(names(this.abundance),names(pred.GCN))]
  correct_RAD_with_confidence(abundance = this.abundance,error = this.error)
},mc.cores = numCores)
#
saveRDS(HMP.RAD,"HMP/HMP_v13_HQ.RAD.RDS")

