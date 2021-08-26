#
library(RasperGade)
source("RasperGade_hidden_state_prediction.R")
source("RasperGade_node_tracking.R")
#
arg = commandArgs(TRUE)
#
jplace.dir = arg[1]
if(is.na(jplace.dir)) stop("EPA-NG output missing!")
#
ref.file = arg[2]
if(is.na(ref.file)) ref.file = "Reference/prepared_reference.RDS"
#
ref.data = readRDS(ref.file)
ref.tree = ref.data$original
#
numCores = as.numeric(arg[3])
if(is.na(numCores)) numCores = 1
cat(sprintf("Working on %d cores.\n",numCores))
#
save.name = arg[4]
if(is.na(save.name)) save.name = paste0(tree.dir,"GCN_prediction")
#
insert.locations = parse_jplace(jplace.dir,split = numCores)
saveRDS(insert.locations,paste0(save.name,".locations.RDS"))
cat(sprintf("Predicting GCN for %d inserted sequences\n",
            sum(sapply(insert.locations,function(x){length(x$query)}))))
#
insert.prediction = mclapply(insert.locations,function(this.insert){
  this.res = predictHiddenStateWithPE(FMR = ref.data$FMR,query.keys = this.insert$hash,laplace = FALSE)
  return(this.res)
},mc.cores = numCores)
saveRDS(insert.prediction,paste0(save.name,".RDS"))
#
cat("Prediction completed!\n")
#
insert.res = list(hsp=do.call(rbind,lapply(insert.prediction,function(x){x$hsp})),
                  error=do.call(c,lapply(insert.prediction,function(x){x$error})))
insert.discrete.res = discretizeResult(res = insert.res$hsp,
                                       error = insert.res$error,laplace = FALSE,numCores=numCores)
#
saveRDS(list(contiguous=insert.res,discrete=insert.discrete.res),paste0(save.name,".RDS"))