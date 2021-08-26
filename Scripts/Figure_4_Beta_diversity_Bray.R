#
library(ggplot2)
library(ggpubr)
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
# load data
cat("Loading data\n")
#
NSTD.cutoff = 10^seq(-3,0,length.out = 10)
true.NSTD.cutoff = c(0,NSTD.cutoff[-1])
#
trait = readRDS("Reference/homogeneous_data.RDS")$dat
#
CV.adjusted.NSTD = readRDS("CV/CV.NSTD.RDS")
#
CV.res.summarized = readRDS("CV/GCN.PE.CV.RDS")
#
cat("Data loaded\n")
# simulate communities for RAD and beta diversity
sim.meta = readRDS("Sim/sim_2env_5otu_2f_Beta_meta.RDS")
CV.community.sim = readRDS("Sim/sim_2env_5otu_2f_Beta_data.RDS")
cat("Data simulated\n")
#
CV.sim.table = lapply(1:length(CV.community.sim),function(i){
  mclapply(CV.community.sim[[i]],function(this.batch){
    gene.table=as.data.frame(data.table::rbindlist(lapply(this.batch,function(x){
      data.frame(t(x$gene)/sum(x$gene))
    }),fill=TRUE))
    cell.table=as.data.frame(data.table::rbindlist(lapply(this.batch,function(x){
      data.frame(t(x$cell)/sum(x$cell))
    }),fill=TRUE))
    gene.table[is.na(gene.table)] = 0
    cell.table[is.na(cell.table)] = 0
    return(list(gene=gene.table,cell=cell.table))
  },mc.cores = numCores)
})
#
CV.sim.rf = readRDS("Sim/Sim_5_2rf_test.RDS")
CV.sim.beta = readRDS("Sim/Sim_5_2beta.RDS")
CV.sim.PERMANOVA = readRDS("Sim/Sim_5_2PERMANOVA.RDS")
CV.sim.diff = readRDS("Sim/Sim_5_2_diff.RDS")
# plot example PCA
i=1
j=1
adj.sim.meta = data.frame(sampleID = 1:dim(sim.meta)[1],environment = sim.meta$V1)
example.gene.physeq = phyloseq(sample_data(adj.sim.meta),
                               otu_table(CV.sim.table[[i]][[j]]$gene,taxa_are_rows = FALSE))
example.cell.physeq = phyloseq(sample_data(adj.sim.meta),
                               otu_table(CV.sim.table[[i]][[j]]$cell,taxa_are_rows = FALSE))
sample_names(example.gene.physeq) = paste0("Gene",sample_names(example.gene.physeq))
sample_names(example.cell.physeq) = paste0("Cell",sample_names(example.cell.physeq))
sample_data(example.gene.physeq)$type="Gene"
sample_data(example.cell.physeq)$type="Cell"
example.correct.physeq = example.gene.physeq
pred.GCN = round(CV.res.summarized[[i]][[j]]$hsp$x)
names(pred.GCN) = CV.res.summarized[[i]][[j]]$hsp$label
true.GCN = trait[CV.res.summarized[[i]][[j]]$hsp$label]
if(attr(example.gene.physeq@otu_table,which = "taxa_are_rows")){
  correct.example.table = correct_abundance_table(abundance = t(example.gene.physeq@otu_table),
                                                  GCN = pred.GCN[rownames(example.gene.physeq@otu_table)])
}else{
  correct.example.table = correct_abundance_table(abundance = example.gene.physeq@otu_table,
                                                  GCN = pred.GCN[colnames(example.gene.physeq@otu_table)])
}
example.correct.physeq@otu_table = otu_table(correct.example.table,taxa_are_rows = FALSE)
sample_names(example.correct.physeq) = paste0("Correct",sample_names(example.correct.physeq))
sample_data(example.correct.physeq)$type = "Correct"
example.combined.physeq = merge_phyloseq(example.cell.physeq,
                                         example.gene.physeq,
                                         example.correct.physeq)
example.combined.ord = ordinate(example.combined.physeq,"PCoA",distance = "bray")
correct.plot = plot_ordination(physeq = example.combined.physeq,
                               ordination = example.combined.ord,color = "environment")
#
example.plot.data = correct.plot$data
example.plot.data$environment = as.character(example.plot.data$environment)
example.arrow.data = data.frame(gene.PC1=example.plot.data$Axis.1[example.plot.data$type=="Gene"],
                                gene.PC2=example.plot.data$Axis.2[example.plot.data$type=="Gene"],
                                correct.PC1=example.plot.data$Axis.1[example.plot.data$type=="Correct"],
                                correct.PC2=example.plot.data$Axis.2[example.plot.data$type=="Correct"],
                                cell.PC1=example.plot.data$Axis.1[example.plot.data$type=="Cell"],
                                cell.PC2=example.plot.data$Axis.2[example.plot.data$type=="Cell"],
                                environment=example.plot.data$environment)
example.plot = ggplot()+
  geom_point(mapping = aes(x=Axis.1,y=Axis.2,group=type,shape=type,color=environment),data = example.plot.data,
             alpha=0.5)+
  geom_segment(mapping = aes(x=cell.PC1,y=cell.PC2,xend=correct.PC1,yend=correct.PC2,color=environment),
               linetype="dotted",alpha=0.5,data = example.arrow.data)+
  geom_segment(mapping = aes(x=cell.PC1,y=cell.PC2,xend=gene.PC1,yend=gene.PC2,color=environment),
               linetype="solid",alpha=0.5,data = example.arrow.data)+
  coord_cartesian(xlim=c(-0.35,0.45),ylim = c(-0.4,0.4))+
  scale_x_continuous(breaks = seq(-1,1,0.25))+
  scale_shape_manual(values = c(Correct=16,Gene=17,Cell=15),
                     labels=c(Correct="Corrected abundance",
                              Gene="Gene abundance",
                              Cell="True cell abundance"))+
  guides(shape=guide_legend(title="Method",order = 1),color=guide_legend(title = "Environment",order = 2))+
  xlab("PC1")+ylab("PC2")+
  theme(legend.position = c(0.9,0.15),legend.box="vertical",
        panel.background = element_blank(),
        axis.line = element_line())
#
beta.data = lapply(1:length(CV.sim.beta),function(i){
  good.sample = which(sapply(CV.sim.beta[[i]],function(x){!is.null(x$gene)}))
  this.cutoff = CV.sim.beta[[i]][good.sample]
  print(length(this.cutoff))
  do.call(rbind,lapply(1:length(this.cutoff),function(j){
    this.batch = this.cutoff[[j]]
    total.shift = c(gene=total_shift(this.batch$gene,this.batch$cell,squared=FALSE),
                    correct=total_shift(this.batch$correct,this.batch$cell,squared=FALSE))
    mean.shift = total.shift/dim(as.matrix(this.batch$cell))
    mean.dist = mean(as.matrix(vegdist(as.matrix(this.batch$cell),
                                  method = "euclidean"))[1:10,11:20])
    mean.within.dist = mean(c(vegdist(as.matrix(this.batch$cell)[1:10,1:10],
                                      method = "euclidean"),
                              vegdist(as.matrix(this.batch$cell)[11:20,11:20],
                                      method = "euclidean")))
    element.cp = c(gene=elementwise_coverage_dist(dist1 = this.batch$gene,CI = this.batch$CI,detail = FALSE),
                   cell=elementwise_coverage_dist(dist1 = this.batch$cell,CI = this.batch$CI,detail = FALSE))
    this.NSTD = sapply(CV.community.sim[[i]][[good.sample[j]]],function(x){
      exp(sum(log(CV.adjusted.NSTD[[i]][2,names(x$cell)])*relative_abundance(x$cell)))
    })
    data.frame(shift=t(mean.shift),cp=t(element.cp),
               NSTD = mean(this.NSTD),improvement=unname(1-total.shift[2]/total.shift[1]),
               mean.pairwise.dist = mean.dist,mean.within.dist = mean.within.dist,
               relative.shift = t(mean.shift)/mean.dist)
  }))
})
#
PERMANOVA.data = lapply(1:length(CV.sim.PERMANOVA),function(i){
  good.sample = which(sapply(CV.sim.PERMANOVA[[i]],function(x){!is.null(x$gene)}))
  this.cutoff = CV.sim.PERMANOVA[[i]][good.sample]
  print(length(this.cutoff))
  do.call(rbind,lapply(1:length(this.cutoff),function(j){
    this.batch = this.cutoff[[j]]
    R2 = c(gene= this.batch$gene.PERMANOVA[1],
           cell = this.batch$cell.PERMANOVA[1],
           correct=this.batch$R2[1])
    dR2 = c(gene= abs(this.batch$gene.PERMANOVA[1]-this.batch$cell.PERMANOVA[1]),
            correct=abs(this.batch$R2[1]-this.batch$cell.PERMANOVA[1]))
    cp = c(gene=as.numeric((this.batch$gene.PERMANOVA[1]<this.batch$R2.CI[2])&
                             (this.batch$gene.PERMANOVA[1]>this.batch$R2.CI[1])),
           cell=as.numeric((this.batch$cell.PERMANOVA[1]<this.batch$R2.CI[2])&
                             (this.batch$cell.PERMANOVA[1]>this.batch$R2.CI[1])))
    this.NSTD = sapply(CV.community.sim[[i]][[good.sample[j]]],function(x){
      exp(sum(log(CV.adjusted.NSTD[[i]][2,names(x$cell)])*relative_abundance(x$cell)))
    })
    data.frame(dR2=t(dR2),cp=t(cp),R2=t(R2),
               NSTD = mean(this.NSTD),improvement=unname(1-dR2[2]/dR2[1]))
  }))
})
rf.data = lapply(1:length(CV.sim.rf),function(i){
  this.cutoff = CV.sim.rf[[i]]
  do.call(rbind,lapply(1:length(this.cutoff),function(j){
    this.batch = this.cutoff[[j]]
    this.sim = CV.community.sim[[i]][[j]]
    this.group = do.call(c,lapply(this.sim,function(x){x$group}))
    this.group = this.group[unique(names(this.group))]
    this.feature = names(which(this.group>0))
    gene.recall = sum(which(names(this.batch$gene$rank)%in%this.feature)<=10)
    cell.recall = sum(which(names(this.batch$cell$rank)%in%this.feature)<=10)
    correct.recall = sum(which(names(this.batch$correct$correct$rank)%in%this.feature)<=10)
    this.NSTD = sapply(CV.community.sim[[i]][[j]],function(x){
      exp(sum(log(CV.adjusted.NSTD[[i]][2,names(x$cell)])*relative_abundance(x$cell)))
    })
    return(data.frame(gene=gene.recall,cell=cell.recall,correct=correct.recall,NSTD=mean(this.NSTD)))
  }))
})
#
diff.data = lapply(1:length(CV.sim.diff),function(i){
  this.cutoff = CV.sim.diff[[i]]
  do.call(rbind,lapply(this.cutoff,function(this.batch){
    data.frame(cell.diff = exp(log(this.batch$cell$diff)),
               cell.p = this.batch$cell$test,
               gene.diff = exp(log(this.batch$gene$diff)),
               gene.p = this.batch$gene$test,
               ddiff = exp(abs(log(this.batch$cell$diff/this.batch$gene$diff))))
  }))
}) 
#
PERMANOVA.plot = ggplot()+
  geom_point(mapping=aes(x=R2.cell,y=R2.gene),data = PERMANOVA.data[[1]],alpha=0.5)+
  geom_abline(slope = 1,intercept = 0,linetype="dashed",color="red")+
  xlab(expression(R^2*" using true cell abundance"))+
  ylab(expression(R^2*" using gene abundance"))+
  coord_cartesian(xlim=c(0.05,0.2),ylim = c(0.05,0.2))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line())
#
diff.plot = ggplot()+
  geom_point(mapping=aes(x=cell.diff,y=gene.diff),data = diff.data[[1]],alpha=0.25)+
  geom_abline(slope = 1,intercept = 0,linetype="dashed",color="red")+
  scale_x_continuous(trans="log10",breaks = 10^seq(-2,2,1),
                     labels = c(expression("10"[-2]),expression("10"[-1]),expression("10"[0]),
                                expression("10"[1]),expression("10"[2])))+
  scale_y_continuous(trans="log10",breaks = 10^seq(-2,2,1),
                     labels = c(expression("10"[-2]),expression("10"[-1]),expression("10"[0]),
                                expression("10"[1]),expression("10"[2])))+
  coord_cartesian(xlim=c(0.01,100),ylim=c(0.01,100))+
  xlab("Fold change in true cell abundance")+
  ylab("Fold change in gene abundance")+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line())
#
pvalue.plot = ggplot()+
  geom_point(mapping=aes(x=cell.p,y=gene.p),data = diff.data[[1]],alpha=0.25)+
  geom_abline(slope = 1,intercept = 0,linetype="dashed",color="red")+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  xlab("Fold change in true cell abundance")+
  ylab("Fold change in gene abundance")+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line())
#
right.panel = ggarrange(plotlist = list(PERMANOVA.plot,diff.plot),
                        labels = c("B","C"),ncol = 1,nrow = 2,
                        align = "hv",common.legend = FALSE)
#
combined.plot = ggarrange(plotlist = list(example.plot,right.panel),
                          labels = c("A",""),ncol = 2,nrow = 1,
                          align = "hv",common.legend = FALSE)
if(Sys.info()["sysname"]=='Windows'){
  ggsave(filename = "Fig_4.png",plot = combined.plot,device = "png",
         width = 9,height = 6,units = "in",
         dpi = "print",scale = 1.5,type = "cairo")
}else{
  ggsave(filename = "Fig_4.png",plot = combined.plot,device = "png",
         width = 9,height = 6,units = "in",
         dpi = "print",scale = 1.5)
}
ggsave(filename = "Fig_4.pdf",plot = combined.plot,device = "pdf",
       width = 9,height = 6,units = "in",scale = 1.5)


