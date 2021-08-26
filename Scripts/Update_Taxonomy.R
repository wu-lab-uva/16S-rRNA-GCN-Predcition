#
library(taxize)
library(myTAI)
library(plyr)
#
taxids <- readRDS("Reference/taxids.RDS")
#
taxtable = readRDS("Reference/lineage_table.RDS")
#
taxon_summary <- ncbi_get_taxon_summary(id = unique(taxids))
#
match.idx= match(taxon_summary$uid,taxtable$taxid)
missing.idx = which(is.na(match.idx))
if(length(missing.idx)>0){
  
}else{
  diff.idx = which(taxon_summary$name!=as.character(taxtable$species[match.idx]))
}
#
df_list <- list()
for (i in diff.idx){
  tax  <- taxonomy(organism = taxon_summary[i,]$name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list[[i]] <- df
  Sys.sleep(0.5)
}
df_list = df_list[diff.idx]
# for those with brackets in the name, run again
for(i in diff.idx[which(sapply(df_list,dim)[1,]<1)]){
  ncbi.name = taxon_summary[i,]$name
  new.name=gsub("\\s*\\([^\\)]+\\)","",ncbi.name)
  new.name = gsub(pattern = "[",replacement = "",x = new.name,fixed = TRUE)
  new.name = gsub(pattern = "]",replacement = "",x = new.name,fixed = TRUE)
  tax  <- taxonomy(organism = new.name, db = "ncbi",output = "classification")
  df <- data.frame(lapply(tax$name, function(x) data.frame(x)))
  colnames(df) <- tax$rank
  df_list[[which(diff.idx==i)]] <- df
  Sys.sleep(0.5)
}
#
taxa.level = c("superkingdom","phylum","class","order","family","genus","species")
combined_df <- do.call(rbind.fill, df_list)[,taxa.level]
combined_df$taxid = taxon_summary$uid[diff.idx]
#
taxtable = rbind(taxtable[-match.idx[diff.idx],],combined_df)
#
#saveRDS(taxtable,"Reference/lineage_table.RDS")