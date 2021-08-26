#
library(RasperGade)
source("RasperGade_hidden_state_prediction.R")
source("RasperGade_node_tracking.R")
#
arg = commandArgs(TRUE)
#
ref.file = arg[1]
if(is.na(ref.file)) ref.file = "Reference/rescaled_data_model.RDS"
ref.data = readRDS(ref.file)
#
ref.FMR = 
  fullMarginalReconstructionWithPE(
    phy = ref.data$phy,x = ref.data$dat[ref.data$phy$tip.label],
    params = ref.data$model,laplace = FALSE)
ref.FMR$scale = ref.data$scale.branch
ref.FMR$epsilon.branch = ref.data$epsilon.branch
ref.FMR$key = make_traceable_edge_label(ref.data$phy)
ref.data$FMR = ref.FMR
#
save.file = arg[2]
if(is.na(save.file)) save.file = "Reference/prepared_reference.RDS"
saveRDS(ref.data,save.file)
