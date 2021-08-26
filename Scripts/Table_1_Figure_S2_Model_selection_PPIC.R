#
library(ggplot2)
library(ggpubr)
library(RasperGade)
source("RasperGade16S_partition.R")
#
homo.data = readRDS("Reference/homogeneous_data.RDS")
homo.model = readRDS("Reference/homogeneous_models.RDS")
#
homo.BM.model = fitPE.phy(phy = homo.data$phy,x = homo.data$dat[homo.data$phy$tip.label],
                          start.value = 
                            initialize.parameters.SP(c(lambda=0,size=0,sigma=1e5,epsilon=0)),
                          fixed = list(lambda=log(0),size=log(0),epsilon=log(0)),
                          laplace = FALSE,ignore.zero = TRUE,eval.only = FALSE)
#
partition.data = readRDS("Reference/binary_partition.RDS")
numTransitions = 
  get_partition_transition(tree = homo.data$phy,
                           group = partition.data[[2]]$partition)
#
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
#
hetero.BM.AIC = round(slow.BM.model$AIC) + round(reg.BM.model$AIC)+2*numTransitions
hetero.BME.AIC = round(slow.BME.model$AIC) + round(reg.BME.model$AIC)+2*numTransitions
hetero.PE.AIC = round(slow.PE.model$AIC) + round(reg.PE.model$AIC)+2*numTransitions
#
table1.df = data.frame(label=c("Homogeneous","Hetero-regular","Hetero-slow","Hetero-full"),
                       BM = round(c(homo.BM.model$AIC,reg.BM.model$AIC,
                                    slow.BM.model$AIC,hetero.BM.AIC)),
                       BME = round(c(homo.model$stat[1],reg.BME.model$AIC,
                                    slow.BME.model$AIC,hetero.BME.AIC)),
                       PE = round(c(homo.model$stat[2],reg.PE.model$AIC,
                                    slow.PE.model$AIC,hetero.PE.AIC)))
#
print(table1.df)
#
detailed.partition = check_partition(tree = homo.data$phy,
                                     dAIC = rep(0,Nnode(homo.data$phy)),
                                     group = partition.data[[2]]$partition)
#
group.count = table(sapply(detailed.partition,function(x){x$group}))
group.tip.size = table(partition.data[[2]]$partition[1:Ntip(homo.data$phy)])
group.node.size = table(partition.data[[2]]$partition[Ntip(homo.data$phy)+
                                                        (1:Nnode(homo.data$phy))])
#
homo.BM.pic = pic(x = homo.data$dat[homo.data$phy$tip.label],phy = homo.data$phy)
homo.BM.pic = homo.BM.pic/sqrt(var(homo.BM.pic))
homo.BME.ppic = pseudo.pic(phy = homo.data$phy,x = homo.data$dat[homo.data$phy$tip.label],
                           pfunc = "pPEpoisnorm",params = homo.model$params$bme)
homo.PEE.ppic = pseudo.pic(phy = homo.data$phy,x = homo.data$dat[homo.data$phy$tip.label],
                           pfunc = "pPEpoisnorm",params = homo.model$params$pee)
std.norm.data = data.frame(x=seq(-4.95,4.95,0.1),
                           density = sapply(seq(-5,4.9,0.1),function(x){
                             pnorm(q = x+0.1)-pnorm(q = x)
                           })/0.1)
#
reg.nodes = which(partition.data[[2]]$partition[-(1:Ntip(homo.data$phy))]>0)
slow.nodes = which(partition.data[[2]]$partition[-(1:Ntip(homo.data$phy))]<0)
#
hetero.BM.pic = pic(x = homo.data$dat[homo.data$phy$tip.label],phy = homo.data$phy)
hetero.BM.pic[reg.nodes] = homo.BM.pic[reg.nodes]/sqrt(var(homo.BM.pic[reg.nodes]))
hetero.BM.pic[slow.nodes] = homo.BM.pic[slow.nodes]/sqrt(var(homo.BM.pic[slow.nodes]))
#
hetero.BME.ppic = hetero.PEE.ppic = 1:Nnode(homo.data$phy)
hetero.BME.ppic[reg.nodes] = pseudo.pic(phy = homo.data$phy,x = homo.data$dat[homo.data$phy$tip.label],
                                        pfunc = "pPEpoisnorm",params = reg.BME.model$params)[reg.nodes]
hetero.BME.ppic[slow.nodes] = pseudo.pic(phy = homo.data$phy,x = homo.data$dat[homo.data$phy$tip.label],
                                         pfunc = "pPEpoisnorm",params = slow.BME.model$params)[slow.nodes]
#
hetero.PEE.ppic[reg.nodes] = pseudo.pic(phy = homo.data$phy,x = homo.data$dat[homo.data$phy$tip.label],
                                        pfunc = "pPEpoisnorm",params = partition.data[[2]]$model$regular)[reg.nodes]
hetero.PEE.ppic[slow.nodes] = pseudo.pic(phy = homo.data$phy,x = homo.data$dat[homo.data$phy$tip.label],
                                         pfunc = "pPEpoisnorm",params = partition.data[[2]]$model$slow)[slow.nodes]
#
homo.plot.data = data.frame(ppic=c(sort(homo.BM.pic),sort(homo.BME.ppic),sort(homo.PEE.ppic)),
                            emp_p = rep(seq(0,1,length.out=length(homo.BM.pic)+1)[-1]-0.5/length(homo.BM.pic),3),
                            model=rep(c("BM","BME","PEE"),each=length(homo.BM.pic)))
homo.plot.data$theo_p = pnorm(homo.plot.data$ppic)
hetero.plot.data = data.frame(ppic=c(sort(hetero.BM.pic),sort(hetero.BME.ppic),sort(hetero.PEE.ppic)),
                              emp_p = rep(seq(0,1,length.out=length(homo.BM.pic)+1)[-1]-0.5/length(homo.BM.pic),3),
                              model=rep(c("BM","BME","PEE"),each=length(hetero.BM.pic)))
hetero.plot.data$theo_p = pnorm(hetero.plot.data$ppic)
#
homo.ppic.plot = ggplot()+
  geom_histogram(mapping = aes(x=ppic,y = ..density..),binwidth = 0.1,data = homo.plot.data)+
  geom_line(mapping = aes(x=x,y=density),color="red",linetype="dashed",data = std.norm.data)+
  facet_wrap(facets = vars(model),nrow = 1,ncol = 3)+
  coord_cartesian(xlim=c(-5,5),ylim = c(0,0.5))+
  xlab("Pseudo-PIC")+ylab("Density")+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line())
homo.pp.plot = ggplot()+
  geom_line(mapping = aes(x=emp_p,y = theo_p),data = homo.plot.data)+
  geom_abline(slope = 1,intercept = 0,color="red",linetype="dashed")+
  facet_wrap(facets = vars(model),nrow = 1,ncol = 3)+
  coord_cartesian(xlim=c(0,1),ylim = c(0,1))+
  xlab("Empirical quantile")+ylab("Theoretical quantile")+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line())
#
hetero.ppic.plot =   ggplot()+
  geom_histogram(mapping = aes(x=ppic,y = ..density..),binwidth = 0.1,data = hetero.plot.data)+
  geom_line(mapping = aes(x=x,y=density),color="red",linetype="dashed",data = std.norm.data)+
  facet_wrap(facets = vars(model),nrow = 1,ncol = 3)+
  coord_cartesian(xlim=c(-5,5),ylim = c(0,0.5))+
  xlab("Pseudo-PIC")+ylab("Density")+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line())
hetero.pp.plot = ggplot()+
  geom_line(mapping = aes(x=emp_p,y = theo_p),data = hetero.plot.data)+
  geom_abline(slope = 1,intercept = 0,color="red",linetype="dashed")+
  facet_wrap(facets = vars(model),nrow = 1,ncol = 3)+
  coord_cartesian(xlim=c(0,1),ylim = c(0,1))+
  xlab("Empirical quantile")+ylab("Theoretical quantile")+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line())
#
combined.plot = ggarrange(plotlist = list(homo.ppic.plot,homo.pp.plot,
                                          hetero.ppic.plot,hetero.pp.plot),
                          ncol = 1,nrow = 4,labels = "AUTO",align = "hv")
#
if(Sys.info()["sysname"]=='Windows'){
  ggsave(filename = "Fig_S2.png",plot = combined.plot,device = "png",
         width = 6,height = 8,units = "in",
         dpi = "print",scale = 1.5,type = "cairo")
}else{
  ggsave(filename = "Fig_S2.png",plot = combined.plot,device = "png",
         width = 6,height = 8,units = "in",
         dpi = "print",scale = 1.5)
}
#