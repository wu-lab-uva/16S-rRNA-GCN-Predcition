# Load required packages
library(RasperGade)
library(RasperGade16S)
# Read in data
homo.data = readRDS("Reference/homogeneous_data.RDS")
homo.model = readRDS("Reference/homogeneous_models.RDS")
# Fit BM model
homo.BM.model = fitPE.phy(phy = homo.data$phy,x = homo.data$dat[homo.data$phy$tip.label],
                          start.value = 
                            initialize.parameters.SP(c(lambda=0,size=0,sigma=1e5,epsilon=0)),
                          fixed = list(lambda=log(0),size=log(0),epsilon=log(0)),
                          laplace = FALSE,ignore.zero = TRUE,eval.only = FALSE)
# Read in partition results
partition.data = readRDS("Reference/binary_partition.RDS")
numTransitions = 
  get_partition_transition(tree = homo.data$phy,
                           group = partition.data[[2]]$partition)
# Fit models in two groups
reg.BM.model = 
  fitPE.phy(phy = homo.data$phy,
     x = homo.data$dat[homo.data$phy$tip.label],
     start.value = 
       initialize.parameters.SP(c(lambda=0,size=0,sigma=1e5,epsilon=0)),
     fixed = list(lambda=log(0),size=log(0),epsilon=log(0)),
     laplace = FALSE,ignore.zero = TRUE,
     bootstrap = which(partition.data[[2]]$partition[(1:Nnode(homo.data$phy))+Ntip(homo.data$phy)]>0))
#
reg.BME.model = 
  fitPE.phy(phy = homo.data$phy,
            x = homo.data$dat[homo.data$phy$tip.label],
            start.value = 
              initialize.parameters.SP(c(lambda=0,size=0,sigma=20,epsilon=1)),
            fixed = list(lambda=log(0),size=log(0)),
            laplace = FALSE,ignore.zero = FALSE,
            bootstrap = which(partition.data[[2]]$partition[(1:Nnode(homo.data$phy))+Ntip(homo.data$phy)]>0))
reg.PE.model = 
  fitPE.phy(phy = homo.data$phy,
            x = homo.data$dat[homo.data$phy$tip.label],
            start.value = 
              initialize.parameters.SP(partition.data[[2]]$model$regular),
            fixed = list(sigma=log(partition.data[[2]]$model$regular[3])),
            laplace = FALSE,eval.only = TRUE,ignore.zero = FALSE,
            bootstrap = which(partition.data[[2]]$partition[(1:Nnode(homo.data$phy))+Ntip(homo.data$phy)]>0))
#
slow.BM.model =
  fitPE.phy(phy = homo.data$phy,
            x = homo.data$dat[homo.data$phy$tip.label],
            start.value = 
              initialize.parameters.SP(c(lambda=0,size=0,sigma=1e5,epsilon=0)),
            fixed = list(lambda=log(0),size=log(0),epsilon=log(0)),
            laplace = FALSE,ignore.zero = TRUE,
            bootstrap = which(partition.data[[2]]$partition[(1:Nnode(homo.data$phy))+Ntip(homo.data$phy)]<0))
slow.BME.model =
  fitPE.phy(phy = homo.data$phy,
            x = homo.data$dat[homo.data$phy$tip.label],
            start.value = 
              initialize.parameters.SP(c(lambda=0,size=0,sigma=20,epsilon=1)),
            fixed = list(lambda=log(0),size=log(0)),
            laplace = FALSE,ignore.zero = FALSE,
            bootstrap = which(partition.data[[2]]$partition[(1:Nnode(homo.data$phy))+Ntip(homo.data$phy)]<0))
slow.PE.model = 
  fitPE.phy(phy = homo.data$phy,
            x = homo.data$dat[homo.data$phy$tip.label],
            start.value = 
              initialize.parameters.SP(partition.data[[2]]$model$slow),
            fixed = list(size = log(partition.data[[2]]$model$slow[2]),
                         sigma=log(partition.data[[2]]$model$slow[3])),
            laplace = FALSE,eval.only = TRUE,ignore.zero = FALSE,
            bootstrap = which(partition.data[[2]]$partition[(1:Nnode(homo.data$phy))+Ntip(homo.data$phy)]<0))
# Get AIC
hetero.BM.AIC = round(slow.BM.model$AIC) + round(reg.BM.model$AIC)+2*numTransitions
hetero.BME.AIC = round(slow.BME.model$AIC) + round(reg.BME.model$AIC)+2*numTransitions
hetero.PE.AIC = round(slow.PE.model$AIC) + round(reg.PE.model$AIC)+2*numTransitions
# Make Table 1
table1.df = data.frame(label=c("Homogeneous","Hetero-regular","Hetero-slow","Hetero-full"),
                       BM = round(c(homo.BM.model$AIC,reg.BM.model$AIC,
                                    slow.BM.model$AIC,hetero.BM.AIC)),
                       BME = round(c(homo.model$stat[1],reg.BME.model$AIC,
                                    slow.BME.model$AIC,hetero.BME.AIC)),
                       PE = round(c(homo.model$stat[2],reg.PE.model$AIC,
                                    slow.PE.model$AIC,hetero.PE.AIC)))
# Print Table 1
write.table(table1.df,file = "Table_1.txt",
            col.names = TRUE,row.names = FALSE,quote=FALSE,sep = "\t")
# Some additional statistics
detailed.partition = check_partition(tree = homo.data$phy,
                                     dAIC = rep(0,Nnode(homo.data$phy)),
                                     group = partition.data[[2]]$partition)
#
group.count = table(sapply(detailed.partition,function(x){x$group}))
group.tip.size = table(partition.data[[2]]$partition[1:Ntip(homo.data$phy)])
group.node.size = table(partition.data[[2]]$partition[Ntip(homo.data$phy)+
                                                        (1:Nnode(homo.data$phy))])
# Print them
cat("Neighborhood count per group\n")
print(group.count)
cat("Tip count per group\n")
print(group.tip.size)
cat("Node count per group\n")
print(group.node.size)


