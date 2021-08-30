# Load required packages
library(ggplot2)
library(ggpubr)
# Read in data
NSTD.data = read.table("EBI/EBI.adj.NSTI.txt",sep = "\t",stringsAsFactors = FALSE,header=TRUE)
# Parse biome information
NSTD.data$Lfull = sapply(NSTD.data$biomes,function(x){
  xx = gsub(x,pattern = "root:",replacement = "")
})
NSTD.data$L1 = sapply(NSTD.data$biomes,function(x){
  xx = strsplit(x = x,split = ":",fixed = TRUE)[[1]]
  if(length(xx)<2) xx = c(xx,"Unassigned")
  xx[2]
})
NSTD.data$L2 = sapply(NSTD.data$biomes,function(x){
  xx = strsplit(x = x,split = ":",fixed = TRUE)[[1]]
  if(length(xx)<3) xx = c(xx,"Unassigned")
  paste0(xx[2:min(3,length(xx))],collapse = ":")
})
NSTD.data$L3 = sapply(NSTD.data$biomes,function(x){
  xx = strsplit(x = x,split = ":",fixed = TRUE)[[1]]
  if(length(xx)<4) xx = c(xx,"Unassigned")
  paste0(xx[2:min(4,length(xx))],collapse = ":")
})
# count L3 biomes
count.L3 = sort(table(NSTD.data$L3),decreasing = TRUE)
# Show top 35 biomes
L3toShow = setdiff(names(count.L3)[count.L3>=85],"Unassigned")
# Make Figure 5
NSTI.L3.plot = ggplot(data = NSTD.data[NSTD.data$L3 %in% L3toShow,])+
  geom_histogram(mapping = aes(x=adj.NSTI,fill=L1),bins = 50)+
  geom_vline(xintercept = 0.3,color="red",linetype="dashed")+
  coord_cartesian(xlim=c(1e-3,1))+
  scale_x_continuous(trans="log10")+
  scale_fill_discrete(name="")+
  xlab("Adjusted NSTI")+ylab("Count")+
  facet_wrap(vars(L3),ncol = 5,scales = "free_y")+
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.line = element_line())
#
if(Sys.info()["sysname"]=='Windows'){
  ggsave(filename = "Fig_5.png",plot = NSTI.L3.plot,device = "png",
         width = 10,height = 6,units = "in",
         dpi = "print",scale = 2,type = "cairo")
}else{
  ggsave(filename = "Fig_5.png",plot = NSTI.L3.plot,device = "png",
         width = 10,height = 6,units = "in",
         dpi = "print",scale = 2)
}
#
ggsave(filename = "Fig_5.pdf",plot = NSTI.L3.plot,device = "pdf",
       width = 10,height = 6,units = "in",scale = 2)