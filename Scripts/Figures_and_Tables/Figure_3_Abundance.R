#
library(RasperGade16S)
library(ggplot2)
library(ggpubr)
# load data
cat("Loading data\n")
#
NSTD.cutoff = 10^seq(-3,0,length.out = 10)[1:9]
true.NSTD.cutoff = c(0,NSTD.cutoff[-1])
#
CV.res.summarized = readRDS("CV/GCN.PE.CV.RDS")
CV.NSTD = readRDS("CV/CV.NSTD.RDS")
cat("Data loaded\n")
# simulate communities for RAD and beta diversity
CV.community.sim = readRDS("Sim/sim_1env_20otu_1f_RA_data.RDS")
cat("Data simulated\n")
#
CV.sim.RAD = readRDS("Sim/Sim_1_1RA.RDS")
cat("RAD calculated\n")
# Theoretical fold-change by abundance
theoretical_fold_change = function(abundance,ACN.fold){
  this.fold = 1/(ACN.fold+(1-ACN.fold)*abundance)
  return(this.fold)
}
#
theoretical.data = expand.grid(abundance=seq(0,1,0.05),ACN.fold=c(0.2,0.5,1,2,5),KEEP.OUT.ATTRS = TRUE)
theoretical.data$gene = sapply(1:dim(theoretical.data)[1],function(i){
  theoretical_fold_change(abundance = theoretical.data[i,2],
                          ACN.fold = theoretical.data[i,1])
})
theoretical.GCN.data = expand.grid(abundance=c(0.01,0.2,0.4,0.6,0.8),
                                   ACN.fold=exp(seq(-log(5),log(5),length.out = 21)),KEEP.OUT.ATTRS = TRUE)
theoretical.GCN.data$gene = sapply(1:dim(theoretical.GCN.data)[1],function(i){
  theoretical_fold_change(abundance = theoretical.GCN.data[i,2],
                          ACN.fold = theoretical.GCN.data[i,1])
})
#
theoretical.abundance.plot = ggplot(data=theoretical.data)+
  geom_line(mapping = aes(x=abundance,y=gene,group=ACN.fold,color=as.character(ACN.fold)))+
  geom_point(mapping = aes(x=abundance,y=gene,group=ACN.fold,color=as.character(ACN.fold)))+
  scale_y_continuous(trans="log10")+
  scale_color_discrete(name="GCN to average GCN ratio")+
  xlab("True relative cell abundance")+ylab("Fold change in\nrelative cell abundance")+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())

theoretical.GCN.plot = ggplot(data=theoretical.GCN.data)+
  geom_line(mapping = aes(x=ACN.fold,y=gene,group=abundance,color=as.character(abundance)))+
  geom_point(mapping = aes(x=ACN.fold,y=gene,group=abundance,color=as.character(abundance)))+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log10")+
  scale_color_discrete(name="True relative abundance")+
  xlab("Average GCN")+ylab("Fold change in\nrelative cell abundance")+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())
# Relative abundance results
CV.RAD.data = lapply(1:(length(CV.sim.RAD)),function(i){
  this.cutoff.data = do.call(rbind,lapply(1:length(CV.sim.RAD[[i]]),function(j){
    this.batch = CV.sim.RAD[[i]][[j]]
    as.data.frame(t(apply(sapply(this.batch,function(x){
      c(true.cover = coverage_length(rad = x$true,CI = x$CI[,names(x$true)]),
      gene.cover = coverage_length(rad = x$raw,CI = x$CI[,names(x$raw)]),
      rad.cover = coverage_length(rad = x$rad,CI = x$CI[,names(x$rad)]),
      rad.diff = exp(mean(abs(log(x$rad[names(x$true)]/x$true)))),
      raw.diff = exp(mean(abs(log(x$raw[names(x$true)]/x$true)))),
      top.support = x$support[1,1],
      raw.support = x$gene.support[1,1],
      top.hit = as.numeric(names(x$true)[1]==names(x$rad)[1]),
      gene.hit = as.numeric(names(x$true)[1]==names(sort(x$raw,decreasing = TRUE))[1]),
      top.agree = as.numeric(names(x$rad)[1]==names(sort(x$raw,decreasing = TRUE))[1]),
      med.NSTD = exp(sum(log(CV.NSTD[[i]][2,names(x$rad)])*relative_abundance(x$rad)))) 
    }),1,mean)))
  }))
})
#
CV.RAD.line.data = do.call(rbind,lapply(CV.RAD.data,function(x){
  do.call(cbind,apply(x,2,function(xx){
    data.frame(mean=mean(xx),se=std_err(xx))
  }))
}))
#
CV.RAD.data = do.call(rbind,CV.RAD.data)
#
CV.diff.plot.data = data.frame(abundance = c(CV.RAD.data$rad.diff,CV.RAD.data$raw.diff),
                               type=rep(c("With correction","Without correction"),each=dim(CV.RAD.data)[1]),
                               distance = rep(CV.RAD.data$med.NSTD,2))
CV.diff.line.plot.data = data.frame(abundance = c(CV.RAD.line.data$rad.diff.mean,CV.RAD.line.data$raw.diff.mean),
                                    abundance.se = c(CV.RAD.line.data$rad.diff.se,CV.RAD.line.data$raw.diff.se),
                               type=rep(c("With correction","Without correction"),each=dim(CV.RAD.line.data)[1]),
                               distance = rep(CV.RAD.line.data$med.NSTD.mean,2))
#
cat(sprintf("Average fold-change without GCN correction: %.1f\n",mean(CV.RAD.line.data$raw.diff.mean)))
cat(sprintf("Average hit-rate of most abundance OTU without GCN correction: %.0f%%\n",
            mean(CV.RAD.line.data$gene.hit.mean)*100))
cat(sprintf("Best average fold-change with GCN correction: %.1f\n",CV.RAD.line.data$rad.diff.mean[1]))
cat(sprintf("Best average hit-rate of most abundance OTU with GCN correction: %.0f%%\n",
            CV.RAD.line.data$top.hit.mean[1]*100))
#
abun.diff.plot = ggplot(mapping = aes(group=type,color=type))+
  geom_line(data=CV.diff.line.plot.data,
            mapping = aes(x=distance,y=abundance),linetype="solid",size=1)+
  geom_errorbar(data=CV.diff.line.plot.data,
                mapping = aes(x=distance,ymin=abundance-abundance.se,ymax=abundance+abundance.se),
                linetype="dashed",width=0.2)+
  xlab("Adjusted NSTI")+ylab("Average fold change in\nrelative gene abundance")+
  scale_x_continuous(trans="log10")+
  scale_y_continuous(breaks = seq(0,2,0.2))+
  scale_color_manual(name = "",
                     values = c("With correction"="#619CFF",
                                "Without correction"="#F8766D"),
                     labels=c("With correction"="Fold-change with correction",
                              "Without correction"="Fold-change without correction"))+
  guides(color=guide_legend(nrow=2,byrow = FALSE))+
  coord_cartesian(ylim=c(1,2))+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())
#
CV.conf.plot.data = data.frame(abundance = CV.RAD.data$true.cover,
                               distance = CV.RAD.data$med.NSTD)
CV.conf.line.plot.data = data.frame(abundance = CV.RAD.line.data$true.cover.mean,
                                    abundance.se = CV.RAD.line.data$true.cover.se,
                               distance = CV.RAD.line.data$med.NSTD.mean)
abun.cover.plot = ggplot()+
  geom_line(data=CV.conf.line.plot.data,
            mapping = aes(x=distance,y=abundance),linetype="solid",size=1)+
  geom_errorbar(data=CV.conf.line.plot.data,
                mapping = aes(x=distance,ymin=abundance-abundance.se,ymax=abundance+abundance.se),
                linetype="dashed",width=0.2)+
  geom_hline(yintercept = 0.95,linetype="dashed",color="red")+
  xlab("Adjusted NSTI")+ylab("Coverage probability")+
  scale_x_continuous(trans="log10")+
  scale_y_continuous(breaks = seq(0,1,0.1),labels = paste0(round(seq(0,1,0.1)*100),"%"))+
  coord_cartesian(ylim=c(0,1))+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())
#
cutoff.plot.data = data.frame(diff = CV.RAD.data$raw.diff/CV.RAD.data$rad.diff-1,
                              cover = CV.RAD.data$gene.cover)
#
cutoff.plot = ggplot()+
  geom_point(mapping = aes(x=cover,y=diff,),alpha=0.25,
             data = cutoff.plot.data)+
  geom_hline(yintercept = 0,color="red",linetype = "dashed")+
  geom_vline(xintercept = 0.95,color="red",linetype = "dashed")+
  scale_x_continuous(trans = "exp",breaks=seq(0,1,0.1),
                     labels = paste0(as.character(round(100*seq(0,1,0.1))),"%"))+
  scale_y_continuous(breaks=seq(-1,1,0.1),labels = paste0(as.character(round(100*seq(-1,1,0.1))),"%"))+
  coord_cartesian(xlim=c(0,1),ylim = c(-0.2,0.8))+
  xlab("Coverage probability to gene abundance")+ylab("Relative improvement")+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())
# Top identity plot
CV.top.line.plot.data = data.frame(abundance = c(CV.RAD.line.data$top.hit.mean,
                                                 CV.RAD.line.data$gene.hit.mean,
                                                 CV.RAD.line.data$top.support.local.mean),
                                    abundance.se = c(CV.RAD.line.data$top.hit.se,
                                                     CV.RAD.line.data$gene.hit.se,
                                                     CV.RAD.line.data$top.support.local.se),
                                   method = rep(c("With correction",
                                                  "Without correction",
                                                  "Support value"),
                                              each=dim(CV.RAD.line.data)[1]),
                                    distance = rep(CV.RAD.line.data$med.NSTD.mean,3))
#
top.plot = ggplot(mapping = aes(group=method,color=method))+
  geom_line(data=CV.top.line.plot.data,
            mapping = aes(x=distance,y=abundance),linetype="solid",size=1)+
  geom_errorbar(data=CV.top.line.plot.data,
                mapping = aes(x=distance,ymin=abundance-abundance.se,ymax=abundance+abundance.se),
                linetype="dashed",width=0.2)+
  xlab("Adjusted NSTI")+ylab("Empirical probability")+
  scale_x_continuous(trans="log10")+
  scale_y_continuous(breaks = seq(0,1,0.1),labels = paste0(round(seq(0,1,0.1)*100),"%"),
                     sec.axis = sec_axis(trans = ~ . , breaks = seq(0,1,0.1),
                                         labels = seq(0,1,0.1),name = "Support value"))+
  coord_cartesian(ylim=c(0,1))+
  scale_color_manual(name = "",
                     values = c("With correction"="#619CFF",
                                "Without correction"="#F8766D",
                                "Support value"="#00BA38"),
                     labels=c("With correction"="Empirical probability with correction",
                              "Without correction"="Empirical probability without correction",
                              "Support value"="Support value"))+
  guides(color=guide_legend(nrow=2,byrow = FALSE))+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())
#
HMP.abundance = readRDS("HMP/HMP_v13_HQ.abundance.RDS")
HMP.RAD = readRDS("HMP/HMP_v13_HQ.RAD.RDS")
#
HMP.RAD.data = do.call(rbind,lapply(1:length(HMP.abundance),function(i){
  dx = HMP.abundance[[i]]/HMP.RAD[[i]]$rad[names(HMP.abundance[[i]])]
  mean.dx = exp(mean(abs(log(dx))))
  gene.CI = coverage_length(rad = HMP.abundance[[i]][names(HMP.RAD[[i]]$rad)],CI = HMP.RAD[[i]]$CI)
  top.agree = as.numeric(names(HMP.RAD[[i]]$rad)[1] == names(which.max(HMP.abundance[[i]])))
  support = HMP.RAD[[i]]$support[1,1]
  return(data.frame(dx = mean.dx,cp=gene.CI,agree=top.agree,support = support))
}))
#
HMP.ecdf.data = data.frame(cp=seq(0,1,0.001),
                           ecdf = sapply(seq(0,1,0.001),function(x){
                             sum(HMP.RAD.data$cp<=x)/dim(HMP.RAD.data)[1]}))
HMP.ecdf.line.data = data.frame(x=c(0,0.95),xend=c(0.95,0.95),
                                y=c(HMP.ecdf.data$ecdf[951],0),
                                yend=rep(HMP.ecdf.data$ecdf[951],2),
                                idx = 1:2)
HMP.ecdf.plot = ggplot()+
  geom_line(mapping = aes(x=cp,y=ecdf),data = HMP.ecdf.data)+
  geom_segment(mapping = aes(x=x,xend=xend,y=y,yend=yend,group=idx),
               data=HMP.ecdf.line.data,
               color="red",linetype="dashed")+
  xlab("Coverage probability to gene abundance")+ylab("Cumulative probability")+
  scale_x_continuous(breaks = seq(0,1,0.1),
                     labels = paste0(as.character(round(100*seq(0,1,0.1))),"%"))+
  scale_y_continuous(breaks = seq(0,1,0.1),
                     labels = paste0(as.character(round(100*seq(0,1,0.1))),"%"))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line())
#
top.panels = ggarrange(plotlist = list(theoretical.abundance.plot,abun.diff.plot,top.plot),
                       ncol = 3,nrow = 1,labels = c("A","B","C"),align = "hv",common.legend = FALSE)
bottom.panels = ggarrange(plotlist = list(abun.cover.plot,cutoff.plot,HMP.ecdf.plot),
                          ncol = 3,nrow = 1,labels = c("D","E","F"),align = "hv",common.legend = FALSE)
#
combined.plot = ggarrange(plotlist = list(top.panels,bottom.panels),
                          ncol = 1,nrow = 2,align = "hv",common.legend = FALSE)
if(Sys.info()["sysname"]=='Windows'){
  ggsave(filename = "Fig_3.png",plot = combined.plot,device = "png",
         width = 10,height = 6,units = "in",
         dpi = "print",scale = 1.5,type = "cairo")
}else{
  ggsave(filename = "Fig_3.png",plot = combined.plot,device = "png",
         width = 10,height = 6,units = "in",
         dpi = "print",scale = 1.5)
}
ggsave(filename = "Fig_3.pdf",plot = combined.plot,device = "pdf",
       width = 10,height = 6,units = "in",scale = 1.5)