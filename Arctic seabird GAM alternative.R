### Arctic Seabird Counts and adding a GAM smooth
# 
# #setwd("M:/My Documents/State of Birds/SOCB/seabird data 2017")
 library(dplyr)
 library(tibble)
 library(tidyr)
 library(mgcv)
 library(rjags)
 
 
 
 ## reading in the arctic colony data used in teh 2012 SOCB report
 ## this will have to suffice, given the lack of surveys over the last decade.
 
 all = read.csv("data/ArcticSeabirds 2012 socb data.csv",stringsAsFactors = F)
 
 
 
 
 
 
weighted = T #change to FALSE if unweighted gam desired.















sptomod <- unique(all$species)
#spmissed <- names(ncountsbysp)[which(ncountsbysp <= 19)]

source("Functions/transparency function.r")

j = 0
for(ss in sptomod){
#  for(ss in re.run){
    
  j = j+1
  spdat <- all[which(all$species == ss & !is.na(all$estimate)),]
 
  
  

spdat$col.f <- as.integer(factor(spdat$colony))


colonylist = unique(spdat[,c("colony","col.f")])
colony = spdat$col.f

wtlst1 <- tapply(spdat$estimate,spdat$colony,mean,na.rm = T)
wtlst <- wtlst1/sum(wtlst1)
wts <- wtlst[colony]
count <- round(spdat$estimate)
  ncounts = length(count)
  ncolony <- max(colony)
  if(ncolony == 1) {
    print(paste("only one colony for",ss))
  }
  year <- spdat$Year - (min(spdat$Year)-1)
  ymax <- max(year)
  ymin <- min(year)
  colonies <- levels(factor(spdat$colony))
  
 

  
  nknots <- 4
if(ss == "THGU"){nknots = 2}
  

    
###### building the GAM basis function

# Following Crainiceanu, C. M., Ruppert, D. & Wand, M. P. (2005). Bayesian Analysis for Penalized Spline Regression Using WinBUGS. Journal of Statistical Softare, 14 (14), 1-24.

knotsX<- round(seq(ymin,ymax,length=(nknots+2))[-c(1,nknots+2)])
X_K<-(abs(outer(year,knotsX,"-")))^3
X_OMEGA_all<-(abs(outer(knotsX,knotsX,"-")))^3
X_svd.OMEGA_all<-svd(X_OMEGA_all)
X_sqrt.OMEGA_all<-t(X_svd.OMEGA_all$v  %*% (t(X_svd.OMEGA_all$u)*sqrt(X_svd.OMEGA_all$d)))
X.basis<-t(solve(X_sqrt.OMEGA_all,t(X_K)))

nyears = length(ymin:ymax)
ymaxpred = 2009-(min(spdat$Year)-1)
nyearspred = length(ymin:ymaxpred)

X_Kp<-(abs(outer(ymin:ymaxpred,knotsX,"-")))^3
X_OMEGA_allp<-(abs(outer(knotsX,knotsX,"-")))^3
X_svd.OMEGA_allp<-svd(X_OMEGA_allp)
X_sqrt.OMEGA_allp<-t(X_svd.OMEGA_allp$v  %*% (t(X_svd.OMEGA_allp$u)*sqrt(X_svd.OMEGA_allp$d)))
X.pred <- t(solve(X_sqrt.OMEGA_allp,t(X_Kp)))

if(weighted){
dat = list(X.basis = X.basis,
           X.pred = X.pred,
           colony = colony,
           ncounts = ncounts,
           ncolony = ncolony,
           count = count,
           nknots = nknots,
           nyearspred = nyearspred,
 wts = wts,
 wtlst = wtlst)

if(ss == "NOGA"){ #use the Negative binomial model for NOGA
  gamfitqp <- jags.model(data = dat,
                       file = paste0("models/GAM model ","NOGA"," NBaltW.txt"),
                       n.chains = 3)
}else{ #use the overdispersed Poisson for the remaining species
gamfitqp <- jags.model(data = dat,
                       file = paste0("models/GAM model ","NOGA"," quasialtW.txt"),
                       n.chains = 3)
}
}else{
dat = list(X.basis = X.basis,
           X.pred = X.pred,
           colony = colony,
           ncounts = ncounts,
           ncolony = ncolony,
           count = count,
           nknots = nknots,
           nyearspred = nyearspred)#,
           # wts = wts,
           # wtlst = wtlst)

gamfitqp <- jags.model(data = dat,
                       file = paste0("models/GAM model ","NOGA"," quasialt.txt"),
                       n.chains = 3)

}
adaptest <- adapt(object = gamfitqp,
                  n.iter = 10)

while(adaptest == F){
  adaptest <- adapt(object = gamfitqp,
                    n.iter = 1000)
  
}

if(ss == "NOGA"){nburn = 40000}else{nburn = 10000}

update(gamfitqp,n.iter = nburn)

poster <- coda.samples(gamfitqp,
                                c(#"sdnoise",
                                  "etapred",
                                  "beta.X",
                                  "x.gampred",
                                  "pop",
                                  "sdX"),
                                n.iter = 10000,
                                #inits = test$jags.ini,
                                n.burnin = nburn,
                                thin = 10)
#below supports the hierarchical parameterization - not implemented 
# poster <- coda.samples(gamfitqp,
#                        c("sdnoise",
#                          "colpredhyp",
#                          "colpred",
#                          "poppredhyp",
#                          "pop",
#                          "pophyp",
#                          "beta.X",
#                          "sdX",
#                          "sdbeta",
#                          "B.X"),
#                        n.iter = 20000,
#                        #inits = test$jags.ini,
#                        n.burnin = 10000,
#                        thin = 5)

out <- summary(poster)

if(ss == "NOGA"){
pdf(paste0("time series models/gam/",ss,"Arctic NB alt weighted mcmc diagnostic.pdf"))
plot(poster)
dev.off()
}else{

pdf(paste0("time series models/gam/",ss,"Arctic quasi alt weighted mcmc diagnostic.pdf"))
plot(poster)
dev.off()
}

qout <- out$quantiles


netapred = paste0("etapred[",rep(1:ncolony,each = nyearspred),",",rep(1:nyearspred,times = ncolony),"]")


indicescol = data.frame(qout[netapred,])
indicescol[,"node"] <- netapred
indicescol[,"year"] <- rep(ymin:ymaxpred,times = ncolony)+min(spdat$Year) 
indicescol[,"colony"] <-  rep(1:ncolony,each = nyearspred)
indicescol[,"species"] <- ss 
indicescol[,"colonyname"] <- rep(colonies,each = nyearspred) 
indicescol[,"nknots"] <- nknots 

popi <- paste0("pop[",1:nyearspred,"]")
indicesreg <- data.frame(qout[popi,])
indicesreg[,"node"] <- popi
indicesreg[,"year"] <- (ymin:ymaxpred)+min(spdat$Year)
indicesreg[,"species"] <- ss 
indicesreg[,"nknots"] <- nknots 


names(indicescol)[which(names(indicescol) == "X50.")] <- "index"
names(indicescol)[which(names(indicescol) == "X2.5.")] <- "index.lci"
names(indicescol)[which(names(indicescol) == "X97.5.")] <- "index.uci"


names(indicesreg)[which(names(indicesreg) == "X50.")] <- "index"
names(indicesreg)[which(names(indicesreg) == "X2.5.")] <- "index.lci"
names(indicesreg)[which(names(indicesreg) == "X97.5.")] <- "index.uci"


trendcol = rbind(colonylist,colonylist)
trendcol$time = rep(c("long-term","short-term"),each = ncolony)
trendcol[2*ncolony+1,"col.f"] <- NA
trendcol[2*ncolony+1,"colony"] <- "regional"
trendcol[2*ncolony+1,"time"] <- "long-term"
trendcol[2*ncolony+2,"col.f"] <- NA
trendcol[2*ncolony+2,"colony"] <- "regional"
trendcol[2*ncolony+2,"time"] <- "short-term"

trendcol[2*ncolony+3,"col.f"] <- NA
trendcol[2*ncolony+3,"colony"] <- "regional.sum"
trendcol[2*ncolony+3,"time"] <- "long-term"
trendcol[2*ncolony+4,"col.f"] <- NA
trendcol[2*ncolony+4,"colony"] <- "regional.sum"
trendcol[2*ncolony+4,"time"] <- "short-term"

for (cc in 1:ncolony){
  
  ccn = colonylist[which(colonylist$col.f == cc),"colony"]
  netapredi = paste0("etapred[",cc,",",ymin,"]")
  
  base = unlist(poster[,netapredi, drop=FALSE]) 
  
  netapredi = paste0("etapred[",cc,",",nyearspred-10,"]")
  
  base10 = unlist(poster[,netapredi, drop=FALSE]) 
  
  
  netapredi = paste0("etapred[",cc,",",nyearspred,"]")
  
  post = unlist(poster[,netapredi, drop=FALSE]) 
  
  for(tt in c("long-term","short-term")){
    if(tt == "long-term"){
      baset = base
      ymint = ymin
    }else{
      baset = base10
      ymint = nyearspred-10}
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"start.year"] <- ymint+min(spdat$Year-1)
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"end.year"] <- nyearspred+min(spdat$Year-1)
    
    
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"trend"] <- round(median(((post/baset)^(1/(nyearspred-ymint))-1)*100),3)
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"trend.lci"] <- round(quantile(((post/baset)^(1/(nyearspred-ymint))-1)*100,0.025),3)
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"trend.uci"] <- round(quantile(((post/baset)^(1/(nyearspred-ymint))-1)*100,0.975),3)
    
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"pop.change"] <- (round(quantile(post/baset,0.5),3)-1)*100
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"pop.change.lci"] <- (round(quantile(post/baset,0.025),3)-1)*100
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"pop.change.uci"] <- (round(quantile(post/baset,0.975),3)-1)*100
    
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"weight.p"] <- round(wtlst[which(names(wtlst) == ccn)],2)
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"individuals"] <- round(wtlst1[which(names(wtlst1) == ccn)],0)
    
    
    if(tt == "long-term"){
      tmpind = indicescol[which(indicescol$colony == cc),]
      minyy = tmpind[which.min(tmpind$index),"year"]
      miny = minyy-(min(tmpind$year)-1)
      maxyy = tmpind[which.max(tmpind$index),"year"]
      maxy = maxyy-(min(tmpind$year)-1)
      
      netapredi = paste0("etapred[",cc,",",miny,"]")
      
      basetmin = unlist(poster[,netapredi, drop=FALSE]) 
      netapredi = paste0("etapred[",cc,",",maxy,"]")
      
      basetmax = unlist(poster[,netapredi, drop=FALSE]) 
      
      
      trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"change.from.min"] <- (round(quantile(post/basetmin,0.5),3)-1)*100
      trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"change.from.min.lci"] <- (round(quantile(post/basetmin,0.025),3)-1)*100
      trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"change.from.min.uci"] <- (round(quantile(post/basetmin,0.975),3)-1)*100
      
      trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"change.from.max"] <- (round(quantile(post/basetmax,0.5),3)-1)*100
      trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"change.from.max.lci"] <- (round(quantile(post/basetmax,0.025),3)-1)*100
      trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"change.from.max.uci"] <- (round(quantile(post/basetmax,0.975),3)-1)*100
      
    }
    
  }#tt
  
  
  
}#cc

for (cc in c("regional")){
  
 # if(cc == "regional"){
    
    popi = paste0("pop[",ymin,"]")
    
    base = unlist(poster[,popi, drop=FALSE]) 
    
    popi = paste0("pop[",nyearspred-10,"]")
    
    base10 = unlist(poster[,popi, drop=FALSE]) 
    
    popi = paste0("pop[",nyearspred,"]")
    
    post = unlist(poster[,popi, drop=FALSE]) 
    
  
  for(tt in c("long-term","short-term")){
    if(tt == "long-term"){
      baset = base
      ymint = ymin
    }else{
      baset = base10
      ymint = nyearspred-10}
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"start.year"] <- ymint+min(spdat$Year-1)
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"end.year"] <- nyearspred+min(spdat$Year-1)
    
    
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"trend"] <- round(median(((post/baset)^(1/(nyearspred-ymint))-1)*100),3)
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"trend.lci"] <- round(quantile(((post/baset)^(1/(nyearspred-ymint))-1)*100,0.025),3)
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"trend.uci"] <- round(quantile(((post/baset)^(1/(nyearspred-ymint))-1)*100,0.975),3)
    
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"pop.change"] <- (round(quantile(post/baset,0.5),3)-1)*100
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"pop.change.lci"] <- (round(quantile(post/baset,0.025),3)-1)*100
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"pop.change.uci"] <- (round(quantile(post/baset,0.975),3)-1)*100
    
    if(tt == "long-term"){
      if(cc == "regional.sum"){
        tmpind = indicesregr
        minyy = tmpind[which.min(tmpind$index),"year"]
        miny = minyy-(min(tmpind$year)-1)
        maxyy = tmpind[which.max(tmpind$index),"year"]
        maxy = maxyy-(min(tmpind$year)-1)
        
        popir <- paste0("popr[",miny,"]")
        
        basetmin = unlist(poster[,popir, drop=FALSE]) 
        popir <- paste0("popr[",maxy,"]")
        
        basetmax = unlist(poster[,popir, drop=FALSE]) 
      }else{
        tmpind = indicesreg
        minyy = tmpind[which.min(tmpind$index),"year"]
        miny = minyy-(min(tmpind$year)-1)
        maxyy = tmpind[which.max(tmpind$index),"year"]
        maxy = maxyy-(min(tmpind$year)-1)
        
        popi <- paste0("pop[",miny,"]")
        
        basetmin = unlist(poster[,popi, drop=FALSE]) 
        popi <- paste0("pop[",maxy,"]")
        
        basetmax = unlist(poster[,popi, drop=FALSE]) 
      }
      
      trendcol[which(trendcol$colony == cc & trendcol$time == tt),"change.from.min"] <- (round(quantile(post/basetmin,0.5),3)-1)*100
      trendcol[which(trendcol$colony == cc & trendcol$time == tt),"change.from.min.lci"] <- (round(quantile(post/basetmin,0.025),3)-1)*100
      trendcol[which(trendcol$colony == cc & trendcol$time == tt),"change.from.min.uci"] <- (round(quantile(post/basetmin,0.975),3)-1)*100
      
      trendcol[which(trendcol$colony == cc & trendcol$time == tt),"change.from.max"] <- (round(quantile(post/basetmax,0.5),3)-1)*100
      trendcol[which(trendcol$colony == cc & trendcol$time == tt),"change.from.max.lci"] <- (round(quantile(post/basetmax,0.025),3)-1)*100
      trendcol[which(trendcol$colony == cc & trendcol$time == tt),"change.from.max.uci"] <- (round(quantile(post/basetmax,0.975),3)-1)*100
      
    }
    
  }#tt
  
}#cc

trendcol$species = ss

xplot = indicesreg$year
if(ss == "NOGA"){
  png(filename = paste0("time series models/gam/",ss," Arctic gam trajectory overplot NB weighted.png"),
    res = 300,
    height = 8,
    width = 12,
    units = "in") 
}else{
  png(filename = paste0("time series models/gam/",ss," Arctic gam trajectory overplot quasi poisson weighted.png"),
    res = 300,
    height = 8,
    width = 12,
    units = "in") 
}
par(mar = c(3,3,3,12))

plot(1,1,
     xlim = c(1960,2040),
     ylim = c(0,max(indicesreg$index*1.6)),
     type = "l",
     main = ss,
     bty = "l",
     ylab = "population by colony and total",
     xlab = "")

polygon(x = c(xplot,rev(xplot)),
        y = c(indicesreg$index.uci,rev(indicesreg$index.lci)),
        col = transp.func(grey(0.8),0.6),
        border = NA)
lines(y = indicesreg$index,
      x = xplot,
      col = "black")


for(cc in 1:ncolony){
  tmp = indicescol[which(indicescol$colony == cc),"index"]
  tmplci = indicescol[which(indicescol$colony == cc),"index.lci"]
  tmpuci = indicescol[which(indicescol$colony == cc),"index.uci"]
  #coli = unique(indicescol[which(indicescol$colony == j),"colonyname"])
  polygon(x = c(xplot,rev(xplot)),
          y = c(tmplci,rev(tmpuci)),
          col = transp.func(rainbow(length(unique(colony)))[cc],0.1),
          border = NA)
  lines(y = tmp,
        x = xplot,
        col = rainbow(length(unique(colony)))[cc])
  points(col = rainbow(length(unique(colony)))[cc],
         x = min(spdat$Year)+year[which(colony == cc)],
         y = count[which(colony == cc)])
  text(colonies[cc],
       x = 2018,
       y = tmp[length(tmp)],
       col = rainbow(length(unique(colony)))[cc],
       pos = 4)
}


dev.off()

if(j == 1){
  indicespop <- indicesreg
  
  indicescolony <- indicescol
  trendscolony = trendcol
}else{
  indicespop = rbind(indicespop,indicesreg)
  indicescolony = rbind(indicescolony,indicescol)
  trendscolony = rbind(trendscolony,trendcol)
}                      


write.csv(indicescolony,"Arctic colony level indices weighted qpoisson.csv")

write.csv(indicespop,"Arctic population level indices weighted qpoisson.csv")

write.csv(trendscolony,"Arctic population level trends weighted qpoisson.csv")


}#ss






