#


# setup the environment for model fitting
library(RasperGade)
# read in key arguments for the fitting
arg=commandArgs(TRUE)
# interpret the arguments
if(arg[1]==arg[2]){
  # read in phylogeny and trait from 1 RDS file
  tree = readRDS(arg[1])$phy
  trait = readRDS(arg[1])$dat
}else{
  # read in the phylogeny
  tree = read.tree(arg[1])
  # read in the trait value
  trait = read.table(file = arg[2],stringsAsFactors = FALSE)[,2]
  names(trait) = read.table(file = arg[2],stringsAsFactors = FALSE)[,1]
}
# set up a name for saving the result of this analysis
job.name = arg[3]
if(is.na(job.name)) stop("Job name is required!")
# how many cores to use in the fitting
numCores=ceiling(as.numeric(arg[4]))
if(is.na(numCores)) numCores=1
#
repRange = 10
#
ignore.zero = FALSE
#
AIC.cutoff=2
#
reltol = 1e-4
# set scaling coefficient
scale.coef = 1
# organize data for easy transfer to other analysis
this.data = list(phy=tree,dat=trait[tree$tip.label])
if(any(is.na(this.data$dat))) stop("NAs are not supported in tip values!")
saveRDS(this.data,sprintf("%s_data.RDS",job.name))
# estimating scaling coefficient, epsilon and sigma term so the optimizer have a reasonable starting point
est.sigma = var(pic(x = trait[tree$tip.label],phy=tree))*scale.coef
median.l = median(tree$edge.length)
#
terminal.pic = findTerminalPIC(phy = tree,x = sqrt(scale.coef)*trait[tree$tip.label],
                               remove.zero = FALSE)
term.pee = fitPE(x = sapply(terminal.pic$x,function(x){x[1]-x[2]}),
                 l = sapply(terminal.pic$l,sum),numCores = 1,laplace = FALSE,
                 start.value = list(lambda=log(1/median.l),size=log(est.sigma*median.l),epsilon=log(1)),
                 min.value = c(0,0,0,1e-4),
                 fixed = list(sigma=log(0)))
term.pebme = fitPE(x = sapply(terminal.pic$x,function(x){x[1]-x[2]}),
                   l = sapply(terminal.pic$l,sum),numCores = 1,laplace = FALSE,
                   min.value = c(0,0,0,1e-4),
                   start.value = list(lambda=log(1/median.l),size=log(0.9*est.sigma*median.l),sigma=log(0.1*est.sigma),epsilon=log(1)))
# fitting the models
model.functions = c(bme = "fitPE.phy",
                    pee = "fitPE.phy", pebme ="fitPE.phy")
fixed.params = list(bme=list(lambda=log(0),size=log(0)),
                    pee = list(sigma=log(0)),pebme =list())
initial.parameters = list(bme = c(lambda=0,size=0,sigma=est.sigma,epsilon=1),
                          pee = term.pee$params,
                          pebme = term.pebme$params)
min.values = list(bme = c(rep(0,3),1e-6),
                  pee = c(rep(0,3),1e-6), pebme = c(rep(0,3),1e-6))
ini.multiplier = expand.grid(lambda=c(1,repRange,1/repRange),
                             size=c(1,repRange,1/repRange),
                             KEEP.OUT.ATTRS = TRUE)
#
res = list(bme = list(), pee = list(), pebme = list())
for(i in rev(1:3)){
  print(names(res)[i])
  ini.func = "initialize.parameters.SP"
  res[[i]] = mclapply(1:dim(ini.multiplier)[1],function(j){
    if(i>1){
      start.params = initialize.parameters.SP(initial.parameters[[i]]*
        c(ini.multiplier$lambda[j],ini.multiplier$size[j],1,1))
    }else{
      start.params = initialize.parameters.SP(initial.parameters[[i]]*
        c(1,1,ini.multiplier$lambda[j],ini.multiplier$size[j]))
    }
    mdl = NULL
    if(any(sapply(start.params,is.nan))) return(mdl)
    try({
      mdl= fit.model.byNM.2step(fit.func = model.functions[[i]],
                                start.params = list(start.value=start.params),
                                params1 =list(fixed=fixed.params[[i]],min.value = min.values[[i]],
                                              ignore.zero=ignore.zero,laplace = FALSE,
                                              phy=tree,x=sqrt(scale.coef)*trait[tree$tip.label],
                                              numCores=1,reltol=reltol),
                                params2 =list(fixed=fixed.params[[i]],min.value = min.values[[i]],
                                              ignore.zero=ignore.zero,laplace = FALSE,
                                              phy=tree,x=sqrt(scale.coef)*trait[tree$tip.label],
                                              numCores=1),
                                AIC.func = "get.AIC",
                                params.func = "get.params",
                                ini.func = ini.func,
                                AIC.cutoff = AIC.cutoff)
    })
    return(mdl)
  },mc.cores = numCores)
  #
  saveRDS(list(scale=scale.coef,res = res,jitter=0),
          sprintf("%s_models.RDS",job.name))
}
#
model.summary = lapply(res,function(x){
  x = x[!is.null(x)]
  x = x[[which.min(sapply(x,get.AIC))]]
  return(list(AIC=get.AIC(x),params=get.params(x)))
})
#
saveRDS(list(scale=scale.coef,jitter=0,res = res,stat=sapply(model.summary,get.AIC),params=lapply(model.summary,get.params)),
        sprintf("%s_models.RDS",job.name))
#
warnings()
