# Calculate NSTD for cross-validation
library(RasperGade)
library(castor)
# load rescaled tree
rescaled.data = readRDS("Reference/rescaled_data_model.RDS")
# set up cutoff
NSTD.cutoff = 10^seq(-3,0,length.out = 10)
true.NSTD.cutoff = c(0,NSTD.cutoff[-1])
#
CV.NSTD = lapply(true.NSTD.cutoff,function(x){
  this.NSTD = sapply(rescaled.data$phy$tip.label,function(this.tip){
      res = find_nearest_tips(tree = rescaled.data$original,
                              only_descending_tips = FALSE,target_tips = this.tip)
      adj.res = find_nearest_tips(tree = rescaled.data$phy,
                                  only_descending_tips = FALSE,target_tips = this.tip)
      all.unadjusted.NSTD = 
        res$nearest_distance_per_tip[-which(rescaled.data$phy$tip.label==this.tip)]
      all.adjusted.NSTD = 
        adj.res$nearest_distance_per_tip[-which(rescaled.data$phy$tip.label==this.tip)]
      raw.NSTD = min(all.unadjusted.NSTD[all.unadjusted.NSTD>x])
      adj.NSTD = min(all.adjusted.NSTD[all.unadjusted.NSTD>x])
      return(c(raw=raw.NSTD,adj=adj.NSTD))
  })
  return(this.NSTD)
})
saveRDS(CV.NSTD,"CV/CV.NSTD.RDS")
#