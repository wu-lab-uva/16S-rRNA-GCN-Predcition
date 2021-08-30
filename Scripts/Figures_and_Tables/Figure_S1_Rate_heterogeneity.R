# Load required packages
library(RasperGade)
library(ggplot2)
# Read in data
tree = readRDS("Reference/homogeneous_data.RDS")$phy
trait = readRDS("Reference/homogeneous_data.RDS")$dat
taxid =readRDS("Reference/taxids.RDS")[tree$tip.label]
taxlin =readRDS("Reference/lineage_table.RDS")
# Get taxonomy information for each tip
tip.tax = sapply(tree$tip.label,function(x){
  sapply(taxlin[taxlin$taxid==as.character(taxid[x]),],as.character)
})
# Set up taxonomic levels
tax.level = c("phylum","class","order","family","genus","species")
higher.tax.level = c("superkingdom","phylum","class","order","family","genus")
names(higher.tax.level) = tax.level
# Unique taxonomic groups in each level
tax.groups = lapply(tax.level,function(x){
  tax = unique(tip.tax[x,])
  tax[!is.na(tax)]
})
names(tax.groups) = tax.level
# Get node index for each genus' MRCA
genus.node = sapply(tax.groups$genus,function(x){
  this.tip = which(tip.tax["genus",]==x)
  if(length(this.tip)<10) return(NA)
  getMRCA(phy = tree,tip = this.tip)
})
names(genus.node) = tax.groups$genus
genus.node = genus.node[!is.na(genus.node)]
# Calculate average rate as the variance of PIC
ave.rate = sapply(names(genus.node),function(x){
  subtree = extract.clade(tree,genus.node[x])
  subtree = drop.tip(subtree,subtree$tip.label[tip.tax["genus",subtree$tip.label]!=x])
  var(pic(phy = subtree,x = trait[subtree$tip.label]))
})
# Compile rate data 
rate.data = data.frame(rate=ave.rate,log10rate=log10(ave.rate+0.1),level="genus")
# Make Figure S1
rate.plot = ggplot()+
  geom_histogram(mapping = aes(x=log10rate),binwidth = 0.25,data=rate.data)+
  scale_x_continuous(breaks = seq(-1,6,1),
                     labels = sapply(seq(-1,6,1),function(i){
                       parse(text=sprintf("10^%d",i))
                     }))+
  xlab(expression("Log"[10]*" rate"))+ylab("Count")+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line())
#
if(Sys.info()["sysname"]=='Windows'){
  ggsave(filename = "Fig_S1.png",plot = rate.plot,device = "png",
         width = 4,height = 3,units = "in",
         dpi = "print",scale = 1.5,type = "cairo")
}else{
  ggsave(filename = "Fig_S1.png",plot = rate.plot,device = "png",
         width = 4,height = 3,units = "in",
         dpi = "print",scale = 1.5)
}
ggsave(filename = "Fig_S1.pdf",plot = rate.plot,device = "pdf",
       width = 4,height = 3,units = "in",scale = 1.5)
