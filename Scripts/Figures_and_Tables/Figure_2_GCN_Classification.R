# Load required packages
library(RasperGade)
library(ggplot2)
library(ggpubr)
# Read in trait values
pe.trait = readRDS("Reference/homogeneous_data.RDS")$dat
min.trait = 0
# Set NSTD thresholds
NSTD.cutoff = 10^seq(-3,0,length.out = 10)[1:9]
true.NSTD.cutoff = c(0,NSTD.cutoff[-1])
# Read in CV results
CV.BM = readRDS("CV/GCN.BM.CV.RDS")
CV.PE = readRDS("CV/GCN.PE.CV.RDS")
CV.MP_EMP = readRDS("CV/GCN.MP_EMP.CV.RDS")
# Read in NSTD at different threshold
CV.NSTD = readRDS("CV/CV.NSTD.RDS")
# Estimate uncertainty
CV.uncertainty = lapply(1:length(CV.BM),function(i){
  sapply(1:length(CV.PE[[i]]),function(j){
    sapply(c(list(HomoPE=CV.PE[[i]][[j]],HomoBM=CV.BM[[i]][[j]]),CV.MP_EMP[[i]][[j]]),function(x){
      if(is.data.frame(x)){
        within.range = which(CV.NSTD[[i]][1,x$label]<(NSTD.cutoff[i]*(10^(1/3))))
        if(colnames(x)[4]=="var"){
          mean(1-discretizeResult(res = x[within.range,],error = NULL,epsilon = 0,laplace = FALSE)$probs)
        }else{
          mean(1-x$probs)
        }
      }else{
        within.range = which(CV.NSTD[[i]][1,x$hsp$label]<(NSTD.cutoff[i]*(10^(1/3))))
        mean(1-discretizeResult(res = x$hsp[within.range,],error = x$error[within.range],
                                epsilon = 0,laplace = FALSE)$probs)
      }
    })
  })
})
# Get classification profile
this.alpha = 0.05
CV.error = lapply(1:length(CV.BM),function(i){
  lapply(1:length(CV.PE[[i]]),function(j){
    lapply(c(list(HomoPE=CV.PE[[i]][[j]],HomoBM=CV.BM[[i]][[j]]),CV.MP_EMP[[i]][[j]]),function(x){
      if(is.data.frame(x)){
        within.range = which(CV.NSTD[[i]][1,x$label]<(NSTD.cutoff[i]*(10^(1/3))))
        if(colnames(x)[4]=="var"){
          calculateExclusionErrorRate(obs = pe.trait[x$label[within.range]],
                                      pred = x$x[within.range],error = x$var[within.range],
                                      tolerance = 0,alpha = this.alpha,epsilon = 0,laplace = FALSE,
                                      discrete = TRUE)
        }else{
          list(obs=(pe.trait[x$label[within.range]]-min.trait)!=(x$x[within.range]),
               pred=x$probs[within.range]<(1-this.alpha))
        }
      }else{
        within.range = which(CV.NSTD[[i]][1,x$hsp$label]<(NSTD.cutoff[i]*(10^(1/3))))
        calculateExclusionErrorRate(obs = pe.trait[x$hsp$label[within.range]],
                                    pred = x$hsp$x[within.range],error = x$error[within.range],
                                    tolerance = 0,alpha = this.alpha,epsilon = 0,laplace = FALSE,
                                    discrete = TRUE)
      }
    })
  })
})
# Calculate precision
CV.precision = lapply(CV.error,function(res){
  sapply(res,function(xx){
    sapply(xx,function(xxx){
      unname(precision_recall(obs = !xxx$obs,pred = !xxx$pred,as.na = TRUE)[1])
    })
  })
})
# Calculate recall
CV.recall = lapply(CV.error,function(res){
  sapply(res,function(xx){
    sapply(xx,function(xxx){
      unname(precision_recall(obs = !xxx$obs,pred = !xxx$pred,as.na = TRUE)[2])
    })
  })
})
# Calculate R2
CV.R2CV = lapply(1:length(CV.BM),function(i){
  sapply(1:length(CV.PE[[i]]),function(j){
    sapply(c(list(HomoPE=CV.PE[[i]][[j]],HomoBM=CV.BM[[i]][[j]]),CV.MP_EMP[[i]][[j]]),function(x){
      if(is.data.frame(x)){
        within.range = which(CV.NSTD[[i]][1,x$label]<(NSTD.cutoff[i]*(10^(1/3))))
        if(colnames(x)[4]=="var"){
          calculate_R2cv(trait = pe.trait[x$label[within.range]],
                         pred = round(x$x[within.range]))
        }else{
          calculate_R2cv(trait = pe.trait[x$label[within.range]]-min.trait,
                         pred = round(x$x[within.range]))
        }
      }else{
        within.range = which(CV.NSTD[[i]][1,x$hsp$label]<(NSTD.cutoff[i]*(10^(1/3))))
        calculate_R2cv(trait = pe.trait[x$hsp$label[within.range]],
                       pred = round(x$hsp$x[within.range]))
      }
    })
  })
})
# Make Figure 2
method2show = c("HomoPE","HomoBM","HomoMP","HomoEMP")
label2show = c(HomoPE = "RasperGade16S (PE)",
           HomoBM = "PICRUST2 (pic)",
           HomoMP = "PICRUST2 (mp)",
           HomoEMP = "PICRUST2 (emp)")
all.line.data = do.call(rbind,lapply(method2show,function(x){
  this.method.res = do.call(rbind,lapply(1:length(CV.BM),function(i){
    this.R2.mean = mean(CV.R2CV[[i]][x,],na.rm = TRUE)
    this.R2.CI = unname(std_err(CV.R2CV[[i]][x,]))
    this.precision = CV.precision[[i]][x,]
    this.recall = CV.recall[[i]][x,]
    precision.mean = mean(this.precision,na.rm=TRUE)
    if(is.nan(precision.mean)) precision.mean = NA 
    recall.mean = mean(this.recall,na.rm = TRUE)
    if(is.nan(recall.mean)) recall.mean = NA 
    precision.CI = unname(std_err(this.precision[!is.na(this.precision)]))
    if(is.null(precision.CI)) precision.CI = NA
    recall.CI = unname(std_err(this.recall[!is.na(this.recall)]))
    if(is.null(recall.CI)) recall.CI = NA
    uncertainty.mean = mean(CV.uncertainty[[i]][x,])
    uncertainty.CI = std_err(CV.uncertainty[[i]][x,])
    this.res = data.frame(NSTD = mean(CV.NSTD[[i]][1,do.call(c,lapply(CV.error[[i]],function(res){
      names(res[[x]]$obs)
    }))]),
    precision = precision.mean,
    precision.CI = precision.CI,
    recall = recall.mean,
    recall.CI = recall.CI,
    R2 = this.R2.mean,
    R2.CI = this.R2.CI,
    probs = uncertainty.mean,
    probs.CI = uncertainty.CI,
    model=x,
    stringsAsFactors = FALSE)
    return(this.res)
  }))
  return(this.method.res)
}))
# Panel D
fprecision.plot = ggplot(mapping = aes(group=model,color=model),
                         data=all.line.data)+
  geom_line(mapping = aes(x=NSTD,y=precision),
            size=1,linetype="solid")+
  geom_point(mapping = aes(x=NSTD,y=precision),
            size=1.25)+
  geom_errorbar(mapping = aes(x=NSTD,ymin=precision-precision.CI,ymax=precision+precision.CI),
                linetype="dashed",width=0.2)+
  coord_cartesian(xlim = c(1e-3,1),ylim = c(0,1))+
  scale_x_continuous(trans="log10",breaks=c(0.001,0.01,0.1,1),labels = c(0.001,0.01,0.1,1))+
  scale_y_continuous(breaks=seq(0,1,0.2),labels = seq(0,1,0.2))+
  scale_color_discrete(name="Method",breaks=method2show,
                       labels=label2show[method2show])+
  guides(color=guide_legend(nrow=2,byrow = TRUE))+
  xlab("Mean NSTD")+ylab("Precision")+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())
# Panel C
frecall.plot =ggplot(mapping = aes(group=model,color=model),
                     data=all.line.data)+
  geom_line(mapping = aes(x=NSTD,y=recall),
            size=1,linetype="solid")+
  geom_point(mapping = aes(x=NSTD,y=recall),
            size=1.25)+
  geom_errorbar(mapping = aes(x=NSTD,ymin=recall-recall.CI,ymax=recall+recall.CI),
                linetype="dashed",width=0.2)+
  coord_cartesian(xlim = c(1e-3,1),ylim = c(0,1))+
  scale_x_continuous(trans="log10",breaks=c(0.001,0.01,0.1,1),labels = c(0.001,0.01,0.1,1))+
  scale_y_continuous(breaks=seq(0,1,0.2),labels = seq(0,1,0.2))+
  scale_color_discrete(name="Method",breaks=method2show,
                       labels=label2show[method2show])+
  guides(color=guide_legend(nrow=2,byrow = TRUE))+
  xlab("Mean NSTD")+ylab("Recall")+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())
# Panel B
R2.plot =ggplot(mapping = aes(group=model,color=model),
                data=all.line.data)+
  geom_line(mapping = aes(x=NSTD,y=R2),
            size=1,linetype="solid")+
  geom_point(mapping = aes(x=NSTD,y=R2),
            size=1.25)+
  geom_errorbar(mapping = aes(x=NSTD,ymin=R2-R2.CI,ymax=R2+R2.CI),
                linetype="dashed",width=0.2)+
  coord_cartesian(xlim = c(1e-3,1),ylim = c(0,1))+
  scale_x_continuous(trans="log10",breaks=c(0.001,0.01,0.1,1),labels = c(0.001,0.01,0.1,1))+
  scale_y_continuous(breaks=seq(-1,1,0.2),labels = seq(-1,1,0.2))+
  scale_color_discrete(name="Method",breaks=method2show,
                       labels=label2show[method2show])+
  guides(color=guide_legend(nrow=2,byrow = TRUE))+
  xlab("Mean NSTD")+ylab(expression(R^2))+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())
# Panel A
probs.plot = ggplot(mapping = aes(group=model,color=model),
                    data=all.line.data)+
  geom_line(mapping = aes(x=NSTD,y=probs),
            size=1,linetype="solid")+
  geom_point(mapping = aes(x=NSTD,y=probs),
            size=1.25)+
  geom_errorbar(mapping = aes(x=NSTD,ymin=probs-probs.CI,ymax=probs+probs.CI),
                linetype="dashed",width=0.2)+
  coord_cartesian(xlim = c(1e-3,1),ylim = c(0,1))+
  scale_x_continuous(trans="log10",breaks=c(0.001,0.01,0.1,1),labels = c(0.001,0.01,0.1,1))+
  scale_color_discrete(name="Method",breaks=method2show,
                       labels=label2show[method2show])+
  scale_y_continuous(name="Predicted uncertainty",
                     sec.axis = sec_axis(name = "Confidence",trans = ~1-.)
  )+
  guides(color=guide_legend(nrow=2,byrow = TRUE))+
  xlab("Mean NSTD")+#ylab("Predicted uncertainty")+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())
#
combined.plot = ggarrange(plotlist = list(probs.plot,R2.plot,frecall.plot,fprecision.plot),
                          ncol = 2,nrow = 2,labels = "AUTO",align = "hv",
                          common.legend = TRUE,legend = "bottom")
if(Sys.info()["sysname"]=='Windows'){
  ggsave(filename = "Fig_2.png",plot = combined.plot,device = "png",
         width = 8,height = 6,units = "in",
         dpi = "print",scale = 1.5,type = "cairo")
}else{
  ggsave(filename = "Fig_2.png",plot = combined.plot,device = "png",
         width = 8,height = 6,units = "in",
         dpi = "print",scale = 1.5)
}
ggsave(filename = "Fig_2.pdf",plot = combined.plot,device = "pdf",
       width = 8,height = 6,units = "in",scale = 1.5)
