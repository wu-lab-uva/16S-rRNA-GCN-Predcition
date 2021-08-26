# Plotting reference phylogeny and SILVA insertion
library(ggplot2)
library(ggtree)
library(RasperGade)
source("RasperGade_node_tracking.R")
#
reference.data = readRDS("Reference/prepared_reference.RDS")
#
taxid =readRDS("Reference/taxids.RDS")[reference.data$original$tip.label]
taxlin =readRDS("Reference/lineage_table.RDS")
tip.tax = sapply(reference.data$original$tip.label,function(x){
  sapply(taxlin[taxlin$taxid==as.character(taxid[x]),],as.character)
})
#
tip.tax[is.na(tip.tax)] = "Unassigned"
#
tips.per.node = get_descendant_tips_for_each_node(reference.data$original)
family.per.node = lapply(tips.per.node,function(x){
  unique(tip.tax["family",x])
})
monophyletic.status = sapply(family.per.node,length)
maximal.mono = sapply(1:(Ntip(reference.data$original)+
                           Nnode(reference.data$original)),function(i){
                             if(i==getRoot(reference.data$original)) return(FALSE)
                             return((monophyletic.status[i]==1)&
                               (monophyletic.status[getAncestor(reference.data$original,i)]>1))
                           })
maximal.mono.node = which(maximal.mono)
maximal.mono.size = sapply(maximal.mono.node,function(i){
  length(tips.per.node[[i]])
})
maximal.mono.family = sapply(family.per.node[maximal.mono.node],identity)
#
reference.tree = reference.data$original
reference.tree$edge.length = sqrt(reference.data$original$edge.length)
attr(reference.tree,"rate") =
  sapply(1:(Ntip(reference.data$original)+
              Nnode(reference.data$original)),function(x){
                this.parent = getAncestor(phy = reference.data$original,x = x)
                if(is.na(this.parent)) return("0")
                this.rate = as.character(as.numeric(reference.data$scale.branch[this.parent]>2))
              })
#
reference.SILVA.plot = ggtree(tr = reference.tree,
                        mapping = aes(color=rate),layout = "rectangular")+
  scale_color_manual(values = c("0"="black","1"="red"),
                        labels=c("0"="Regular","1"="Slowly-evolving"))+
  scale_x_continuous(expand=expansion(mult=c(0.05,0.1),add=c(0,0.2)))+
  guides(color=guide_legend(title = "Rate group"))+
  theme(legend.key.size = unit(0.75, 'in'), #change legend key size
        legend.key.height = unit(0.75, 'in'), #change legend key height
        legend.key.width = unit(0.75, 'in'), #change legend key width
        legend.title = element_text(size=28), #change legend title font size
        legend.text = element_text(size=24)) #change legend text font size
#
for(i in which(maximal.mono.size>=10)){
  reference.SILVA.plot = reference.SILVA.plot + geom_cladelabel(node=maximal.mono.node[i],label = maximal.mono.family[i])
}
#
ggplot2::ggsave(filename = "Fig_S2.pdf",plot = reference.SILVA.plot,device = "pdf",
       width = 12,height = 96,units = "in",scale = 1,limitsize = FALSE)
#
