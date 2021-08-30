# Load required packages
library(RasperGade)
# Read in data
reference.data = readRDS("Reference/prepared_reference.RDS")
# Link nodes and edges
node2edge = lapply(names(reference.data$FMR$key),function(x){
  edge.nodes = reference.data$FMR$key[[x]]
  edge.index = intersect(c(which(reference.data$phy$edge[,2]==edge.nodes[1]),
                           which(reference.data$phy$edge[,1]==edge.nodes[1])),
                         c(which(reference.data$phy$edge[,2]==edge.nodes[2]),
                           which(reference.data$phy$edge[,1]==edge.nodes[2])))
  if(length(edge.index)<1) edge.index = which(reference.data$phy$edge[,1]==getRoot(reference.data$phy))
  return(edge.index)
})
names(node2edge) = names(reference.data$FMR$key)
# Read in SILVA results
SILVA.insertion = do.call(c,lapply(readRDS("SILVA/SILVA.pred.GCN.RDS.locations.RDS"),
                                   function(x){x$hash}))
SILVA.res = readRDS("SILVA/SILVA.pred.GCN.RDS.RDS")$discrete
SILVA.NSTD = readRDS("SILVA/SILVA.NSTD.RDS")
# Read in SILVA taxonomy
SILVA.taxa = read.table(file = "SILVA/SILVA132NR99.sorted.id.3level.txt",
                        sep = ";",stringsAsFactors = FALSE)[,1:2]
SILVA.taxa$V1 = sapply(SILVA.taxa$V1,function(x){strsplit(x = x,split = " ",fixed = TRUE)[[1]][1]})
SILVA.taxa$acc = sapply(SILVA.taxa$V1,function(x){strsplit(x = x,split = ".",fixed = TRUE)[[1]][1]})
SILVA.phyla = table(SILVA.taxa$V2)
# Find top phyla in OTU count
top.phyla = names(sort(SILVA.phyla,decreasing = TRUE))[1:6]
# Check rate group for SILVA prediction
SILVA.group = sapply(SILVA.insertion,function(x){
  parent.node = unique(reference.data$original$edge[node2edge[[x$hash[1]]],1])
  2*as.numeric(reference.data$scale.branch[parent.node]<2)-1
})
names(SILVA.group) = sapply(SILVA.insertion,function(x){x$label})
# Empty some memory space
rm(SILVA.insertion)
# Make Table 2
summary.by.phylum = do.call(rbind,lapply(top.phyla,function(x){
  this.pred = SILVA.res$probs[SILVA.res$label%in%SILVA.taxa$V1[SILVA.taxa$V2==x]]
  this.otus = as.character(SILVA.res$label)[SILVA.res$label%in%SILVA.taxa$V1[SILVA.taxa$V2==x]]
  this.NSTD = SILVA.NSTD[2,this.otus]
  this.summary = data.frame(phylum=x,total=length(this.pred),
                            NSTD = median(this.NSTD),
                            slow = sum(SILVA.group[this.otus]<0)/length(this.otus),
                            confident=sum(this.pred>=0.95),moderate=sum(sum(this.pred>=0.5)))
  this.summary$confident.ratio = this.summary$confident/this.summary$total
  this.summary$moderate.ratio = this.summary$moderate/this.summary$total
  return(this.summary)
}))
summary.by.superkingdom = data.frame(phylum="Bacteria",total = length(SILVA.res$probs),
                                     NSTD = median(SILVA.NSTD[2,]),
                                     slow = sum(SILVA.group<0)/length(SILVA.group),
                                     confident=sum(SILVA.res$probs>=0.95),
                                     moderate =sum(SILVA.res$probs>=0.5),
                                     confident.ratio = sum(SILVA.res$probs>=0.95)/length(SILVA.res$probs),
                                     moderate.ratio = sum(SILVA.res$probs>=0.5)/length(SILVA.res$probs))
# Print Table 2
write.table(rbind(summary.by.superkingdom,summary.by.phylum),file = "Table_2.txt",
            col.names = TRUE,row.names = FALSE,quote=FALSE,sep = "\t")




