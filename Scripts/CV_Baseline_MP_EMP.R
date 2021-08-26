# using PICRUST2 to predict GCN
library(RasperGade)
library(castor)
#
source("RasperGade_node_tracking.R")
#
arg = commandArgs(TRUE)
numCores = as.numeric(arg[1])
if(is.na(numCores)) numCores = 1
# load reference data
cat("Loading reference data\n")
tree = readRDS("Reference/homogeneous_data.RDS")$phy
trait = readRDS("CV/Baseline.trait.RDS")
trait = trait-min(trait)
# load tips for test sets
test.set.tips = readRDS("CV/Baseline.testset.RDS")
#
NSTD.cutoff = 10^seq(-3,0,length.out = 10)[1:9]
true.NSTD.cutoff = c(0,NSTD.cutoff[-1])
# get reference trees for each CV
CV.ref.trees = lapply(1:length(test.set.tips),function(i){
  this.NSTD = test.set.tips[[i]]
  lapply(this.NSTD,function(this.tips){
    this.ref = partition_test_set_with_minimum_distance(phy = tree,trait = trait,test.set = this.tips,
                                                        distance = true.NSTD.cutoff[i],location=FALSE)
    return(c(this.ref$ref,list(insert=this.ref$phy)))
  })
})
# CV with PICRUST
system("mkdir PICRUST2")
run.out = mclapply(1:length(CV.ref.trees),function(i){
  print(i)
  this.res = lapply(1:length(CV.ref.trees[[i]]),function(j){
    this.CV = CV.ref.trees[[i]][[j]]
    tree.path = paste0("PICRUST2/CV_",i,"_",j,".tre")
    trait.path = paste0("PICRUST2/CV_",i,"_",j,".txt")
    pred.path = paste0("PICRUST2/CV_",i,"_",j,".pred")
    CI.path = paste0("PICRUST2/CV_",i,"_",j,".CI")
    write.tree(phy = this.CV$insert,file = tree.path)
    write.table(x = data.frame(sequence=names(this.CV$dat),x=this.CV$dat),file = trait.path,
                quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
    cmd1 = sprintf("Rscript castor_hsp_95CI_with_prob.R %s %s mp 0.5 TRUE FALSE %s.mp %s.mp None",
                   tree.path,trait.path,pred.path,CI.path)
    cmd2 = sprintf("Rscript castor_hsp_95CI_with_prob.R %s %s emp_prob 0 TRUE FALSE %s.emp_prob %s.emp_prob None",
                   tree.path,trait.path,pred.path,CI.path)
    cmd.out = system(command = cmd1,intern = TRUE)
    cmd.out = system(command = cmd2,intern = TRUE)
    return()
  })
},mc.cores = numCores)
#
cat("Prediction finished.\nCompiling results.\n")
# compile results
CV.mp_emp.summarized = lapply(1:length(CV.ref.trees),function(i){
  lapply(1:length(CV.ref.trees[[i]]),function(j){
    MP.pred.file = paste0("PICRUST2/CV_",i,"_",j,".pred.mp")
    MP.CI.file = paste0("PICRUST2/CV_",i,"_",j,".pred.mp.prob")
    EMP.pred.file = paste0("PICRUST2/CV_",i,"_",j,".pred.emp_prob")
    EMP.CI.file = paste0("PICRUST2/CV_",i,"_",j,".pred.emp_prob.prob")
    mp.pred = read.table(MP.pred.file,header = TRUE,sep = "\t",stringsAsFactors = FALSE)
    mp.CI = read.table(MP.CI.file,header = TRUE,sep = "\t",stringsAsFactors = FALSE)
    emp.pred = read.table(EMP.pred.file,header = TRUE,sep = "\t",stringsAsFactors = FALSE)
    emp.CI = read.table(EMP.CI.file,header = TRUE,sep = "\t",stringsAsFactors = FALSE)
    this.res = list(HomoMP = data.frame(node=-1,label=mp.pred$sequence,x=mp.pred[,2],
                                        probs=mp.CI[match(mp.pred$sequence,mp.CI$sequence),2],
                                        stringsAsFactors = FALSE),
                    HomoEMP = data.frame(node=-1,label=emp.pred$sequence,x=emp.pred[,2],
                                         probs=emp.CI[match(emp.pred$sequence,emp.CI$sequence),2],
                                         stringsAsFactors = FALSE))
    return(this.res)
  })
})
# save compiled CV results
saveRDS(CV.mp_emp.summarized,"CV/Baseline.MP_EMP.CV.RDS")