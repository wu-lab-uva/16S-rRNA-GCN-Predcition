# functions to be added
library(ggplot2)
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
# Read in data
HMP.physeq = readRDS("../20210601_Manuscript/HMP/phyloseq.rarefy.1000.random.100.site.RDS")
HMP.physeq@otu_table = otu_table(apply(HMP.physeq@otu_table,2,function(x){x/sum(x)}),
                                 taxa_are_rows = TRUE)
# Read in predicted GCN
pred.error = readRDS("HMP_v13_batch00/HMP_v13.pred.GCN.RDS.RDS")
pred.GCN = pred.error$discrete$x
names(pred.GCN) = pred.error$discrete$label
pred.error = pred.error$contiguous$error
#
if(attr(HMP.physeq@otu_table,which = "taxa_are_rows")){
  correct.table = correct_abundance_table(abundance = t(HMP.physeq@otu_table),GCN = pred.GCN[rownames(HMP.physeq@otu_table)])
}else{
  correct.table = correct_abundance_table(abundance = HMP.physeq@otu_table,GCN = pred.GCN[colnames(HMP.physeq@otu_table)])
}
#
sample_data(HMP.physeq)$type = "Gene"
HMP.correct.physeq = phyloseq(otu_table(correct.table,taxa_are_rows = FALSE),
                            sample_data(HMP.physeq),phy_tree(HMP.physeq))
sample_data(HMP.correct.physeq)$type = "Correct"
sample_names(HMP.correct.physeq) = 
  paste0("correct_",sample_names(HMP.correct.physeq))
HMP.combined.physeq = merge_phyloseq(HMP.physeq,HMP.correct.physeq)
#
HMP.physeq.dist = distance(HMP.physeq,method = "bray")
HMP.correct.physeq.dist = distance(HMP.correct.physeq,method = "bray")

#
HMP.combined.ord = ordinate(HMP.combined.physeq,method = "PCoA",distance = "bray")
HMP.combined.plot = plot_ordination(physeq = HMP.combined.physeq,ordination = HMP.combined.ord,
                                       type="samples", color="HMPBodySite",shape="type")
#
combined.ord.data = HMP.combined.plot$data
combined.arrow.data = data.frame(gene.PC1=combined.ord.data$Axis.1[combined.ord.data$type=="Gene"],
                                 gene.PC2=combined.ord.data$Axis.2[combined.ord.data$type=="Gene"],
                                 correct.PC1=combined.ord.data$Axis.1[combined.ord.data$type=="Correct"],
                                 correct.PC2=combined.ord.data$Axis.2[combined.ord.data$type=="Correct"],
                                 HMPBodySite=combined.ord.data$HMPBodySite)
#
HMP.overlay.plot = ggplot()+
  geom_point(mapping = aes(x=Axis.1,y=Axis.2,group=type,shape=type,color=HMPBodySite),data = combined.ord.data,
             alpha=0.5)+
  geom_segment(mapping = aes(x=gene.PC1,y=gene.PC2,xend=correct.PC1,yend=correct.PC2,color=HMPBodySite),
               linetype="dotted",alpha=0.5,data = combined.arrow.data)+
  scale_shape_discrete(labels = c(Gene="Gene abundance",Correct="Corrected abundance"))+
  guides(color=guide_legend(nrow=2,byrow = TRUE),
         shape=guide_legend(title = "Method"))+
  theme(legend.position = "bottom",legend.box="vertical",
        panel.background = element_blank(),
        axis.line = element_line())
#
saveRDS(HMP.overlay.plot,"HMP_v13_batch00/HMP.beta.RDS")
#
HMP.PERMANOVA = adonis(formula = HMP.physeq.dist~HMPBodySite,data = unclass(sample_data(HMP.physeq)))
gene.PERMANOVA = adonis(formula = HMP.correct.physeq.dist~HMPBodySite,data = unclass(sample_data(HMP.correct.physeq)))
#
HMP.bodysites = unique(sample_data(HMP.physeq)$HMPBodySite)
#
bodysite.combs = combn(x = HMP.bodysites,m = 2)
#
pairwise.rf = list()
pairwise.PERMANOVA = list()
pairwise.fold = list()
for(i in 1:dim(bodysite.combs)[2]){
  print(i)
  this.sites = bodysite.combs[,i]
  this.gene.physeq = subset_samples(HMP.physeq,HMPBodySite%in%this.sites)
  this.correct.physeq = subset_samples(HMP.correct.physeq,HMPBodySite%in%this.sites)
  #this.gene.physeq.dist = distance(physeq = this.gene.physeq,method = "bray")
  #this.correct.physeq.dist = distance(physeq = this.correct.physeq,method = "bray")
  #this.correct.PERMANOVA = adonis(formula =this.correct.physeq.dist~HMPBodySite,
  #                            data = unclass(sample_data(this.correct.physeq)))
  #this.gene.PERMANOVA = adonis(formula = this.gene.physeq.dist~HMPBodySite,
  #                             data = unclass(sample_data(this.gene.physeq)))
  #gene.rf = randomforest_cross_validation(predictor = t(otu_table(this.gene.physeq)),
  #                                        response = sample_data(this.gene.physeq)$HMPBodySite)
  #correct.rf = randomforest_cross_validation(predictor = otu_table(this.correct.physeq),
  #                                           response = sample_data(this.correct.physeq)$HMPBodySite)
  correct.idx = split(x = 1:dim(sample_data(this.correct.physeq))[1],
              f = sample_data(this.correct.physeq)$HMPBodySite)
  gene.idx = split(x = 1:dim(sample_data(this.gene.physeq))[1],
                   f = sample_data(this.gene.physeq)$HMPBodySite)
  correct.fold.diff = apply(otu_table(this.correct.physeq)[correct.idx[[1]],],2,mean)/
    apply(otu_table(this.correct.physeq)[correct.idx[[2]],],2,mean)
  gene.fold.diff = apply(otu_table(this.gene.physeq)[,gene.idx[[1]]],1,mean)/
    apply(otu_table(this.gene.physeq)[,gene.idx[[2]]],1,mean)
  correct.test = sapply(which(is.finite(correct.fold.diff)),function(i){
    wilcox.test(x = as.vector(otu_table(this.correct.physeq)[correct.idx[[1]],i]),
                y = as.vector(otu_table(this.correct.physeq)[correct.idx[[2]],i]))$p.value
  })
  gene.test = sapply(which(is.finite(gene.fold.diff)),function(i){
    wilcox.test(x = as.vector(otu_table(this.gene.physeq)[i,gene.idx[[1]]]),
                y = as.vector(otu_table(this.gene.physeq)[i,gene.idx[[2]]]))$p.value
  })
  pairwise.fold[[i]] = list(correct=list(diff = correct.fold.diff, test = correct.test),
              gene=list(diff = gene.fold.diff, test = gene.test))
  #pairwise.rf[[i]] = list(gene=gene.rf,correct=correct.rf)
  #pairwise.PERMANOVA[[i]] = list(gene=this.gene.PERMANOVA,correct=this.correct.PERMANOVA)
}
#
pairwise.rf.stat = sapply(pairwise.rf,function(x){
  c(matchness=rank_square_difference(rank1 = names(x$gene$rank),
                         rank2 = names(x$correct$rank),normalize = TRUE),
    top10=sum(names(x$gene$rank[1:10])%in%names(x$correct$rank[1:10])),
    top100=sum(names(x$gene$rank[1:100])%in%names(x$correct$rank[1:100])))
})
pairwise.PERMANOVA.stat = sapply(pairwise.PERMANOVA,function(x){
  c(gene = x$gene$aov.tab$R2[1],
    correct = x$correct$aov.tab$R2[1])
})
pairwise.fold.stat = sapply(pairwise.fold,function(x){
  calculate_R2cv(trait = x$correct$diff[is.finite(x$correct$diff)],
                 pred = x$gene$diff[is.finite(x$gene$diff)])
})
saveRDS(pairwise.rf,"HMP_v13_batch00/HMP.rf.RDS")
saveRDS(pairwise.PERMANOVA,"HMP_v13_batch00/HMP.PERMANOVA.RDS")
saveRDS(pairwise.fold,"HMP_v13_batch00/HMP.diff.RDS")
