#
library(RasperGade)
library(castor)
#
source("RasperGade_node_tracking.R")
source("RasperGade_hidden_state_prediction.R")
#
arg = commandArgs(TRUE)
location.file = arg[1]
if(is.na(location.file)) stop("Location file must be specified.")
save.name = arg[2]
if(is.na(save.name)) save.name = "GCN.NSTD.RDS"
#
ref.data = readRDS("Reference/prepared_reference.RDS")
nearest.ref = find_nearest_tips(tree = ref.data$original,only_descending_tips = FALSE)
nearest.tip = c(nearest.ref$nearest_tip_per_tip,nearest.ref$nearest_tip_per_node)
nearest.adj.ref = find_nearest_tips(tree = ref.data$phy,only_descending_tips = FALSE)
nearest.adj.tip = c(nearest.adj.ref$nearest_tip_per_tip,nearest.adj.ref$nearest_tip_per_node)
#
insert.locations = unname(do.call(c,lapply(readRDS(location.file),function(x){x$hash})))
#
cat("Location loaded\nCalculating NSTD\n")
#
insert.NSTD = sapply(insert.locations,function(x){
  des.nodes = ref.data$FMR$key[[x$hash[1]]]
  des.l = x$l[1:2]
  if(is.null(des.nodes)){
    des.nodes = ref.data$FMR$key[[x$hash[2]]]
    des.l = rev(des.l)
    if(is.null(des.nodes)) return(c(NA,NA))
  }
  parent.node = intersect(des.nodes,sapply(des.nodes,getAncestor,phy=ref.data$original))
  if(length(parent.node)<1) parent.node = getAncestor(phy = ref.data$original,des.nodes[1])
  neighbor.raw.NSTD = get_pairwise_distances(tree = ref.data$original,
                                    A = des.nodes,
                                    B = nearest.tip[des.nodes])
  neighbor.adjusted.NSTD = get_pairwise_distances(tree = ref.data$phy,
                                         A = des.nodes,
                                         B = nearest.adj.tip[des.nodes])
  min.node = unname(min(neighbor.raw.NSTD+des.l))
  min.adj.node = unname(min(neighbor.adjusted.NSTD+des.l/ref.data$scale.branch[parent.node]))
  return(c(raw = min.node+x$l[3],
           adjusted=min.adj.node+x$l[3]/ref.data$scale.branch[parent.node]+
             ref.data$epsilon.branch[parent.node]))
})
colnames(insert.NSTD) = sapply(insert.locations,function(x){x$label})
rownames(insert.NSTD) = c("raw","adjusted")
cat("Calculation finished\nSaving...\n")
saveRDS(insert.NSTD,save.name)
#