# Load required packages
library(ape)
library(castor)
library(ggplot2)
library(ggpubr)
# Read in data
tree = read.tree("Reference/reference.tre")
taxid =readRDS("Reference/taxids.RDS")[tree$tip.label]
taxlin =readRDS("Reference/lineage_table.RDS")
# Compile tip taxonomy
tip.tax = sapply(tree$tip.label,function(x){
  sapply(taxlin[taxlin$taxid==as.character(taxid[x]),],as.character)
})
# Set up taxonomic levels
tax.level = c("phylum","class","order","family","genus","species")
higher.tax.level = c("superkingdom","phylum","class","order","family","genus")
names(higher.tax.level) = tax.level
# Get unique taxonomic group in each level
tax.groups = lapply(tax.level,function(x){
  tax = unique(tip.tax[x,])
  tax[!is.na(tax)]
})
names(tax.groups) = tax.level
# Get NSTD by taxonomic level
all.NSTD.by.group = lapply(tax.level,function(x){
  res = lapply(tax.groups[[x]],function(xx){
    g1 = which(tip.tax[x,] == xx)
    hightax = unique(tip.tax[higher.tax.level[x],g1])
    g2 = which(tip.tax[higher.tax.level[x],] == hightax)
    g3 = setdiff(g2,g1)
    if(length(g3)<1) return(numeric(0))
    this.res = find_nearest_tips(
      tree = tree,target_tips = g3,
      only_descending_tips = FALSE)
    return(this.res$nearest_distance_per_tip[g1])
  })
  names(res) = tax.groups[[x]]
  return(list(level=x,NSTD=res))
})
# Calculate the mean NSTD
mean.NSTD.by.group = lapply(all.NSTD.by.group,function(x){
  res = do.call(c,lapply(x$NSTD,function(xx){
    if(length(xx)<1) return(numeric(0))
    exp(mean(log(xx)))
  }))
  return(list(level=x$level,NSTD=res))
})
# Compile NSTD data
mean.NSTD.by.group.data = do.call(rbind,lapply(mean.NSTD.by.group,function(x){
  data.frame(NSTD=x$NSTD,level=x$level,stringsAsFactors = FALSE)
}))
mean.NSTD.by.group.data$level = factor(mean.NSTD.by.group.data$level,levels = tax.level)
# Make Figure S3
NSTD.plot = ggplot()+
  geom_violin(mapping = aes(x=level,y=NSTD),
              draw_quantiles = c(0.25,0.5,0.75),
                data=mean.NSTD.by.group.data)+
  geom_hline(yintercept = 0.5,color="red",linetype="dashed")+
  scale_y_continuous(trans="sqrt",
                     breaks = c(0.05,0.1,0.5,1))+
  scale_x_discrete(breaks=tax.level,labels=sapply(tax.level,stringr::str_to_sentence))+
  xlab("Taxonomic level")+ylab("NSTD")+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())
#
if(Sys.info()["sysname"]=='Windows'){
  ggsave(filename = "Fig_S3.png",plot = NSTD.plot,device = "png",
         width = 4,height = 3,units = "in",
         dpi = "print",scale = 1.5,type = "cairo")
}else{
  ggsave(filename = "Fig_S3.png",plot = NSTD.plot,device = "png",
         width = 4,height = 3,units = "in",
         dpi = "print",scale = 1.5)
}
ggsave(filename = "Fig_S3.pdf",plot = NSTD.plot,device = "pdf",
       width = 4,height = 3,units = "in",scale = 1.5)

