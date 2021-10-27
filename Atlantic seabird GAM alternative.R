### Updating Atlantic Seabird Counts
# 
# #setwd("M:/My Documents/State of Birds/SOCB/seabird data 2017")
 setwd("c:/SOCB/seabird data 2017")
 library(dplyr)
 library(tibble)
 library(tidyr)
 library(mgcv)
# 
# old <- read.csv("AtlanticSeabirds 2012.csv",
#                 stringsAsFactors = F)
# 
# new <- read.csv("AtlanticSeabirds 2017.csv",
#                 stringsAsFactors = F)
# 
# names(new) <- c("colony",
#                 "species",
#                 "Year",
#                 "estimate",
#                 "estimate.comment")
# 
# old <- old[,c("colony",
#               "species",
#               "Year",
#               "estimate",
#               "estimate.comment")]
# #above step removes column "year.comment", "se", and "X"
# # columns se and X were blank
# # column year.comment only had two entries, both of which related to the
# # mixed Common and Arctic Tern counts from NF, which will be dropped in the next line
# old <- old[-grep(old$species, pattern = "/",fixed = T),]
# 
# all <- as_tibble(rbind(old,new))
# 
# 
# tmp <- count(all,colony)
# write.csv(tmp,"colonies w counts.csv")
# 
# 
# 
# all[which(all$colony == "Gannet Islands Newfoundland"),"colony"] <- "Gannet Clusters"
# all[which(all$colony == "Pee Pee Island (Witless Bay)"),"colony"] <- "Pee Pee Island"
# all[which(all$colony == "Baccalieu Island, NL"),"colony"] <- "Baccalieu Island"
# all[which(all$colony == "Grand Colombier, St Pierre"),"colony"] <- "Grand Colombier"
# 
# 
# all <- arrange(all,colony,species,Year)
# 
# 
# ######### quebec data
# 
# qcold = read.csv("QuebecSeabirds.csv",stringsAsFactors = F)
# qcnew = read.csv("Quebec seabird data 2017.csv",stringsAsFactors = F)
# 
# qc <- rbind(qcold,qcnew)
# 
# 
# qc <- as_tibble(qc[,c("colony",
#             "Year",
#             "species",
#             "estimate",
#             "estimate.comment")])
# 
# table(qc[,c("colony","species")])
# table(all[,c("colony","species")])
# 
# all <- bind_rows(all,qc)
# 
# all[which(all$species == "HEGU"),"species"] = "HERG"
# all[which(all$species == "LHSP"),"species"] = "LESP"
# 
# 
# 
# allmat <- unique(all[,c("colony","species")])
# 
# for(y in min(all$Year,na.rm = T):max(all$Year,na.rm = T)){
#   allmat[,paste0("y",y)] <- NA
#   tmp <- filter(all,Year == y)
#   if(nrow(tmp) > 0){
#     for(cc in unique(tmp$colony)){
#       tmp2 <- filter(tmp,colony == cc)
#       for(ss in unique(tmp2$species)){
#         if(length(which(tmp2$species == ss)) > 1){
#           print(paste(ss,cc,y))
#           allmat[which(allmat$colony == cc &
#                          allmat$species == ss),paste0("y",y)] <- max(tmp2[which(tmp2$species == ss),"estimate"],na.rm = T)
# 
#         }else{
#         allmat[which(allmat$colony == cc &
#                        allmat$species == ss),paste0("y",y)] <- tmp2[which(tmp2$species == ss),"estimate"]
#         }
#         }
#     }
#   }
# }
# 
# 
# allmat <- as_tibble(allmat)
# allmatout <- select(allmat,colony,species,y1960:y2016)
# allmatout <- arrange(allmatout,species,colony)
# 
# for(i in 1:nrow(allmatout)){
#   allmatout[i,"nyears"] <- length(1960:2016)-length(which(is.na(allmatout[i,paste0("y",1960:2016)])))
# yrvc =  names(allmatout)[2+which(!is.na(allmatout[i,paste0("y",1960:2016)]))]
# yrrng = range(as.integer(substr(yrvc,2,5)))
# allmatout[i,"span"] <- diff(yrrng)
# rm("yrrng")
# allmatout[i,"mean count"] <- mean(as.integer(allmatout[i,paste0("y",1960:2016)]),na.rm = T)
#   }
# 
# drop <- which(rowSums(allmatout[,paste0("y",1960:2016)],na.rm = T) == 0)
# write.csv(allmatout[-drop,], "all atlantic colony summary 2017.csv",row.names = F)
# 
# 
# write.csv(all,"tempall.csv")
# # ###############
 # #5-year summaries
 # ###############
# allmat5 <- unique(all[,c("colony","species")])
# ys5 <- seq(from = 1960, to = 2015,
#            by = 5)
# for(y in ys5){
#   allmat5[,paste0("y",y)] <- NA
#   tmp <- filter(all,Year %in% seq(y-2,y+2,1))
#   if(nrow(tmp) > 0){
#     for(cc in unique(tmp$colony)){
#       tmp2 <- filter(tmp,colony == cc)
#       for(ss in unique(tmp2$species)){
#         #if(length(which(tmp2$species == ss)) > 1){
#           #print(paste(ss,cc,y))
#         #}
#         allmat5[which(allmat5$colony == cc &
#                        allmat5$species == ss),paste0("y",y)] <- length(which(tmp2$species == ss))
# 
# 
#         }
#     }
#   }
# }
# 
# for(i in 1:nrow(allmat5)){
# allmat5[i,"n 5years"] <- length(which(!is.na(allmat5[i,paste0("y",ys5)])))
# allmat5[i,"span"] <- diff(range(ys5[which(!is.na(allmat5[i,paste0("y",ys5)]))]))
# }
# 
# 
# write.csv(allmat5,"allmat5.csv")
# 
# ###############
# #decadal summaries
# ###############
# allmat10 <- unique(all[,c("colony","species")])
# ys10 <- seq(from = 1965, to = 2015,
#            by = 10)
# for(y in ys10){
#   allmat10[,paste0("y",y)] <- NA
#   tmp <- filter(all,Year %in% seq(y-5,y+4,1))
#   if(nrow(tmp) > 0){
#     for(cc in unique(tmp$colony)){
#       tmp2 <- filter(tmp,colony == cc)
#       for(ss in unique(tmp2$species)){
#         #if(length(which(tmp2$species == ss)) > 1){
#         #print(paste(ss,cc,y))
#         #}
#         allmat10[which(allmat10$colony == cc &
#                         allmat10$species == ss),paste0("y",y)] <- length(which(tmp2$species == ss))
# 
# 
#       }
#     }
#   }
# }
# 
# for(i in 1:nrow(allmat10)){
#   allmat10[i,"n 10years"] <- length(which(!is.na(allmat10[i,paste0("y",ys10)])))
#   allmat10[i,"span"] <- diff(range(ys10[which(!is.na(allmat10[i,paste0("y",ys10)]))]))
# }
# 
# 
# write.csv(allmat10,"allmat10.csv")
# 
# 
# colsw10 <- allmat10[which(allmat10$span > 20),]
# 
# table(colsw10$species)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### time series for each colony with a first difference model?
# 

# save.image("seabird trend prep.RData")

# load(("seabird trend prep.RData"))

# 
# waic <- function(counts = 1:ncounts,
#                  lambda = lambdal,
#                  count = spdat$estimate){
#   lppds <- vector(length = length(counts))
#   pwaic2s <- vector(length = length(counts))
#   j = 0  
#   for(k in counts){
#     j = j+1
#     logtest <- mean(dpois(count[k],unlist(lambda[,k,drop = F])))
#     lppds[j] <- log(max(logtest,exp(-100)))
#     pwaic2s[j] <-  var(dpois(count[k],unlist(lambda[,k,drop = F]),log = T))
#   }
#   
#   lppd <- mean(lppds)
#   
#   pwaic2 <- sum(pwaic2s)
#   
#   WAIC <- -2*(lppd-pwaic2)
#   return(WAIC)
# }#end waic function
# 
# 
# library(rjags)
# 
# allpost60 <- all[which(all$Year > 1959 & !is.na(all$estimate)),]
# ncountsbysp <- rev(sort(table(allpost60$species)))
# sptomod <- names(ncountsbysp)[which(ncountsbysp > 9)]
# #spmissed <- names(ncountsbysp)[which(ncountsbysp <= 19)]
# 
# source("c:/Functions/transparency function.r")
# 
# # gammodel <- c("GAM model.txt")
# # 
# # waics <- list()
# # length(waics) <- length(sptomod)
# # names(waics) <- sptomod
# 
# cvresults <- data.frame(species = sptomod)
# for(ss in sptomod){
# cvr <- which(cvresults$species == ss)
# spdat <- allpost60[which(allpost60$species == ss),]
# colony <- as.integer(factor(spdat$colony))
# count <- spdat$estimate
# ncounts = length(count)
# ncolony <- max(colony)
# if(ncolony == 1) {
#   print(paste("only one colony for",ss))
# }
# year <- spdat$Year - 1957
# ymax <- max(year)
# ymin <- min(year)
# colonies <- levels(factor(spdat$colony))
# 
# countbycol <- table(colony)
# mncount <- mean(countbycol)
# maxcount <- max(countbycol)
# nknots <- floor(mncount/2)
# nknotso <- min(4,nknots)
# cvresults[cvr,"nknots mncount"] <- nknotso
# cvresults[cvr,"ncolony"] <- ncolony
# cvresults[cvr,"ncounts"] <- ncounts
# cvresults[cvr,"mncount"] <- mncount
# cvresults[cvr,"maxcount"] <- maxcount
# # setting up the initial jags model code with jagam() mgcv package
# #knots <- quantile(year,probs = seq(0,1,length = nknots+1))
# 
# 
# # wts <- tapply(count,colony,mean,na.rm = T)
# # wts <- wts/sum(wts)
# # wtsvec <- wts[colony]
# 
# if (maxcount < 8){next}
# 
# if(ncolony > 1){
#    loocvp <- gam(count~ 1+factor(colony) + s(year),
#                 family = "poisson",
#                 method = "REML")
#    loocvqp <- gam(count~ 1+factor(colony) + s(year),
#                  family = "quasipoisson",
#                  method = "REML")
#  
#    
#    loocvnb <- gam(count~ 1+factor(colony) + s(year),
#                 family = nb(),
#                 method = "REML")
#    cvresults[cvr,"REML Poisson"] <- summary(loocvp)$sp.criterion
#    cvresults[cvr,"REML quasiPoisson"] <- summary(loocvqp)$sp.criterion
#    cvresults[cvr,"REML nb"] <- summary(loocvnb)$sp.criterion
#    cvresults[cvr,"AIC Poisson"] <-loocvp$aic
#    cvresults[cvr,"AIC quasiPoisson"] <-loocvqp$aic
#    cvresults[cvr,"AIC nb"] <-loocvnb$aic
#    cvresults[cvr,"df Poisson"] <- summary(loocvp)$edf
#    cvresults[cvr,"df quasiPoisson"] <- summary(loocvqp)$edf
#    cvresults[cvr,"df nb"] <- summary(loocvnb)$edf
# 
#    loocvnbo <- gam(count~ 1+factor(colony) + s(year,k = nknotso+1),
#                    family = nb(),
#                    method = "REML")
#    loocvnb9 <- gam(count~ 1+factor(colony) + s(year,k = 9),
#                    family = nb(),
#                    method = "REML")
#    
#    cvresults[cvr,"AIC nb nknoto"] <-loocvnbo$aic
#    cvresults[cvr,"AIC nb 9"] <-loocvnb9$aic
#    
#    remls <- c(summary(loocvp)$sp.criterion,
#               summary(loocvnb)$sp.criterion,
#               summary(loocvqp)$sp.criterion)
#    knotvec <- c(summary(loocvp)$edf,
#                 summary(loocvnb)$edf,
#                 summary(loocvqp)$edf)
#    
# 
#      nknots <- ceiling(knotvec[which.min(remls)])
# 
#    
# test <- jagam(count~ 1+factor(colony) + s(year,k = nknots+1),
#               #knots = list(knots = knots),
#               family = "poisson",
#               file = paste("models/GAM model",ss,"quasiblank.txt"))#,
# 
# test2 <- jagam(count~ 1+factor(colony) + s(year,k = nknotso+1),
#               #knots = list(knots = knots),
#               family = "poisson",
#               file = paste("models/GAM model",ss,"nknotso.txt"))#,
# 
# test3 <- jagam(count~ 1+factor(colony) + s(year,k = 9),
#                #knots = list(knots = knots),
#                family = "poisson",
#                file = paste("models/GAM model",ss,"nknots9.txt"))#,
# 
# 
# }else{
#   loocvp <- gam(count~ s(year),
#                 family = "poisson",
#                 method = "REML")
#   loocvqp <- gam(count~ s(year),
#                 family = "quasipoisson",
#                 method = "REML")
#   loocvnb <- gam(count~ s(year),
#                  family = nb(),
#                  method = "REML")
#   
#   cvresults[cvr,"REML Poisson"] <- summary(loocvp)$sp.criterion
#   cvresults[cvr,"REML nb"] <- summary(loocvnb)$sp.criterion
#   cvresults[cvr,"AIC Poisson"] <-loocvp$aic
#   cvresults[cvr,"AIC nb"] <-loocvnb$aic
#   cvresults[cvr,"df Poisson"] <- summary(loocvp)$ed
#   cvresults[cvr,"df nb"] <- summary(loocvnb)$edf
#   cvresults[cvr,"df quasiPoisson"] <- summary(loocvqp)$edf
#   cvresults[cvr,"AIC quasiPoisson"] <-loocvqp$aic
#   cvresults[cvr,"REML quasiPoisson"] <- summary(loocvqp)$sp.criterion
#   
#   loocvnbo <- gam(count~ s(year,k = nknotso+1),
#                   family = nb(),
#                   method = "REML")
#   loocvnb9 <- gam(count~ s(year,k = 9),
#                   family = nb(),
#                   method = "REML")
#   
#   cvresults[cvr,"AIC nb nknoto"] <-loocvnbo$aic
#   cvresults[cvr,"AIC nb 9"] <-loocvnb9$aic
#   
#   
#   remls <- c(summary(loocvp)$sp.criterion,
#              summary(loocvnb)$sp.criterion,
#              summary(loocvqp)$sp.criterion)
#   knotvec <- c(summary(loocvp)$edf,
#                summary(loocvnb)$edf,
#                summary(loocvqp)$edf)
#   
#   nknots = ceiling(knotvec[which.min(remls)])
#   
#   test <- jagam(count~ s(year,k = nknots+1),
#                 #knots = list(knots = knots),
#                 family = "poisson",
#                 file = paste("models/GAM model",ss,"quasiblank.txt"))#,
#   test2 <- jagam(count~ s(year,k = nknotso+1),
#                  #knots = list(knots = knots),
#                  family = "poisson",
#                  file = paste("models/GAM model",ss,"nknotso.txt"))#,
#   
#   test3 <- jagam(count~ s(year,k = 9),
#                  #knots = list(knots = knots),
#                  family = "poisson",
#                  file = paste("models/GAM model",ss,"nknots9.txt"))#,
#   
#   
#               }
# 
# }#ss
#  
# 
# write.csv(cvresults,"cvresults.csv",row.names = F)
# 



#weights = wtsvec)## must add in the site-level intercept
# test2 <- gam(count~ 1+factor(colony) + s(year,k = length(knots)),
#               knots = list(knots = knots),
#               family = "poisson")#,
#               #file = "GAM model.txt",
#               weights = wtsvec)
# gamfit <- jags(data = test$jags.data,
#                model.file = "GAM model.txt",
#                parameters.to.save = "mu",
#                n.chains = 3,
#                n.iter = 200000,
#                #inits = test$jags.ini,
#                n.burnin = 100000) 











##################### running jags GAMs based on cross validation results above
load(("data/seabird trend prep.RData"))
 pairestimates = grep(all$estimate.comment,pattern = "airs",fixed = T)
 all[pairestimates,"estimate"] = all[pairestimates,"estimate"]*2 
 
weighted = T #change to FALSE if unweighted gam desired.
library(rjags)

allpost60 <- all[which(all$Year > 1959 & !is.na(all$estimate)),]
allpost60 = allpost60[-which(allpost60$colony == "Baccalieu Island" & allpost60$species == "LESP" & allpost60$estimate == 6720000),]
ncountsbysp <- rev(sort(table(allpost60$species)))
sptomod <- names(ncountsbysp)[which(ncountsbysp > 9)]
#spmissed <- names(ncountsbysp)[which(ncountsbysp <= 19)]

source("c:/Functions/transparency function.r")

gamparam <- read.csv("cvresults.csv",
                     stringsAsFactors = F)
j = 0
run.order <- c(2:length(sptomod),1)
#re.run <- c("BLGU","NOGA","GRCO","BLKI","LESP")

for(ss in sptomod[run.order]){
#  for(ss in re.run){
    
  j = j+1
  parc <- which(gamparam$species == ss)
  spdat <- allpost60[which(allpost60$species == ss),]
  if(ss == "BLKI"){
    #removing all colonies with < 5% of the maximum colony size because there are so many that they swamp the calculation
    mcountsbycol = tapply(spdat$estimate,spdat$colony,mean,na.rm = T)
    mcountsbycol = sort(mcountsbycol)
    coltokeep = names(mcountsbycol)[which(mcountsbycol > max(mcountsbycol)*0.05)]
    spdat = spdat[which(spdat$colony %in% coltokeep),]
  }
  
  if(ss == "BLGU"){
    #removing all counts from pre 1977
    spdat = spdat[which(spdat$Year > 1976),]
    ncountsnew = tapply(spdat$Year,spdat$colony,length)
    spdat = spdat[which(spdat$colony %in% names(ncountsnew)[which(ncountsnew > 1)]),]
  }
  if(ss == "LESP"){
    # spdat = spdat[which(spdat$Year > 1983),]
    # ncountsnew = tapply(spdat$Year,spdat$colony,length)
    # spdat = spdat[which(spdat$colony %in% names(ncountsnew)[which(ncountsnew > 1)]),]
  ### removing all colonies with <1% of the largest colony
    mcountsbycol = tapply(spdat$estimate,spdat$colony,mean,na.rm = T)
    mcountsbycol = sort(mcountsbycol)
    coltokeep = names(mcountsbycol)[which(mcountsbycol > max(mcountsbycol)*0.01)]
    spdat = spdat[which(spdat$colony %in% coltokeep),]
  }
  
  if(ss == "NOGA"){
    # spdat = spdat[which(spdat$Year > 1983),]
    # ncountsnew = tapply(spdat$Year,spdat$colony,length)
    # spdat = spdat[which(spdat$colony %in% names(ncountsnew)[which(ncountsnew > 1)]),]
    ### removing all colonies with <10% of the largest colony
    
    mcountsbycol = tapply(spdat$estimate,spdat$colony,mean,na.rm = T)
    mcountsbycol = sort(mcountsbycol)
    coltokeep = names(mcountsbycol)[which(mcountsbycol > max(mcountsbycol)*0.1)]
   #above line removes 1 colony that accounts for about 5% of the largest colony but shows a different pattern in recent years
    spdat = spdat[which(spdat$colony %in% coltokeep),]
  }
  
  
  if(ss == "GRCO"){
    #removing all counts from pre 1975
    
    spdat = spdat[which(spdat$Year > 1975),]
    ncountsnew = tapply(spdat$Year,spdat$colony,length)
    spdat = spdat[which(spdat$colony %in% names(ncountsnew)[which(ncountsnew > 1)]),]
  }
  
  if(ss == "ARTE"){
    #removing all counts from pre 1975 because only 1 count at 1 colony and it's 3-5X larger than any other count
    
    spdat = spdat[which(spdat$Year > 1975),]
    ncountsnew = tapply(spdat$Year,spdat$colony,length)
    spdat = spdat[which(spdat$colony %in% names(ncountsnew)[which(ncountsnew > 1)]),]
  }
  
# tmp <- spdat[which(spdat$colony == names(which.max(wts))),]
# for(i in 1:10){
#   spdat <- rbind(spdat,tmp)
# }
spdat$col.f <- as.integer(factor(spdat$colony))
colonylist = unique(spdat[,c("colony","col.f")])
colony = spdat$col.f
wtlst1 <- tapply(spdat$estimate,spdat$colony,mean)
wtlst <- wtlst1/sum(wtlst1)
wts <- wtlst[colony]
count <- spdat$estimate
  ncounts = length(count)
  ncolony <- max(colony)
  if(ncolony == 1) {
    print(paste("only one colony for",ss))
  }
  year <- spdat$Year - 1957
  ymax <- max(year)
  ymin <- min(year)
  colonies <- levels(factor(spdat$colony))
  
 
 

 
  remls <- c(gamparam[parc,"REML.Poisson"],
             gamparam[parc,"REML.nb"],
             gamparam[parc,"REML.quasiPoisson"])
  knotvec <- c(gamparam[parc,"df.Poisson"],
               gamparam[parc,"df.nb"],
               gamparam[parc,"df.quasiPoisson"])
  
  
  nknots <- ceiling(knotvec[which.min(remls)])
  if(ss == "BLGU"){nknots = 8}
  if(ss == "NOGA"){nknots = 8}
  if(ss == "GRCO"){nknots = 5}
  if(ss == "BLKI"){nknots = 3}
  if(ss == "LESP"){nknots = 2}
  if(ss == "TBMU"){nknots = 2}
  if(ss == "ATPU"){nknots = 7}
  
  
#  
# if(ss == "TBMU"){
#   nknots <- 1 ### this is specific to TBMU with had too few observations per colony to run the cross validation
# }
# nknotso <- gamparam[parc,"nknots.mncount"]
# #  if(ncolony != 1) {
    
###### building the GAM basis function

# Following Crainiceanu, C. M., Ruppert, D. & Wand, M. P. (2005). Bayesian Analysis for Penalized Spline Regression Using WinBUGS. Journal of Statistical Softare, 14 (14), 1-24.

knotsX<- seq(ymin,ymax,length=(nknots+2))[-c(1,nknots+2)]
X_K<-(abs(outer(year,knotsX,"-")))^3
X_OMEGA_all<-(abs(outer(knotsX,knotsX,"-")))^3
X_svd.OMEGA_all<-svd(X_OMEGA_all)
X_sqrt.OMEGA_all<-t(X_svd.OMEGA_all$v  %*% (t(X_svd.OMEGA_all$u)*sqrt(X_svd.OMEGA_all$d)))
X.basis<-t(solve(X_sqrt.OMEGA_all,t(X_K)))

nyears = length(ymin:ymax)
ymaxpred = 2016-1957
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

if(ss == "NOGA"){
  gamfitqp <- jags.model(data = dat,
                       file = paste0("models/GAM model ","NOGA"," NBaltW.txt"),
                       n.chains = 3)
}else{
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

if(ss == "NOGA"){nburn = 20000}else{nburn = 10000}

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
                                #n.burnin = nburn,
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
pdf(paste0("time series models/gam/",ss,"NB alt weighted mcmc diagnostic.pdf"))
plot(poster)
dev.off()
}else{

pdf(paste0("time series models/gam/",ss,"quasi alt weighted mcmc diagnostic.pdf"))
plot(poster)
dev.off()
}
# pops <- qout[popi,"50%"]
# popslci <- qout[popi,"2.5%"]
# popsuci <- qout[popi,"97.5%"]
# x11()
# plot(y = pops,
#      x = xplot,
#      type = "l",
#      ylim = c(0,max(popsuci)))
# polygon(x = c(xplot,rev(xplot)),
#         y = c(popslci,rev(popsuci)),
#         col = grey(0.8),
#         border = NA)
# lines(y = pops,
#       x = xplot)
#
# 
# Below is code ot exract predicted values from the hierarchical GAM which 
# could not be fit well with these seabird data
# qout <- out$quantiles
# xplot <- ymin:ymax+1957
# popi <- paste0("pop[",1:length(ymin:ymax),"]")
# 
# pops <- qout[popi,"50%"]
# popslci <- qout[popi,"2.5%"]
# popsuci <- qout[popi,"97.5%"]
# x11()
# plot(y = pops,
#      x = xplot,
#      type = "l",
#      ylim = c(0,max(pops)*3))
# polygon(x = c(xplot,rev(xplot)),
#         y = c(popslci,rev(popsuci)),
#         col = grey(0.8),
#         border = NA)
# lines(y = pops,
#       x = xplot)
# 
# popi <- paste0("pophyp[",1:length(ymin:ymax),"]")
# 
# pops <- qout[popi,"50%"]
# popslci <- qout[popi,"2.5%"]
# popsuci <- qout[popi,"97.5%"]
# 
# 
# polygon(x = c(xplot,rev(xplot)),
#         y = c(popslci,rev(popsuci)),
#         col = grey(0.8),
#         border = NA)
# lines(y = pops,
#       x = xplot)
# 

# trendcol = data.frame(colony = 1:ncolony,
#                       colonyname = as.character(colonies),
#                       stringsAsFactors = F)
# trendcol[ncolony+1,"colony"] <- NA
# trendcol[ncolony+1,"colonyname"] <- "regional"


qout <- out$quantiles


netapred = paste0("etapred[",rep(1:ncolony,each = nyearspred),",",rep(1:nyearspred,times = ncolony),"]")


indicescol = data.frame(qout[netapred,])
indicescol[,"node"] <- netapred
indicescol[,"year"] <- rep(ymin:ymaxpred,times = ncolony)+1957 
indicescol[,"colony"] <-  rep(1:ncolony,each = nyearspred)
indicescol[,"species"] <- ss 
indicescol[,"colonyname"] <- rep(colonies,each = nyearspred) 
indicescol[,"nknots"] <- nknots 

popi <- paste0("pop[",1:nyearspred,"]")
indicesreg <- data.frame(qout[popi,])
indicesreg[,"node"] <- popi
indicesreg[,"year"] <- (ymin:ymaxpred)+1957
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

# trendcol[2*ncolony+3,"col.f"] <- NA
# trendcol[2*ncolony+3,"colony"] <- "regional.sum"
# trendcol[2*ncolony+3,"time"] <- "long-term"
# trendcol[2*ncolony+4,"col.f"] <- NA
# trendcol[2*ncolony+4,"colony"] <- "regional.sum"
# trendcol[2*ncolony+4,"time"] <- "short-term"

for (cc in 1:ncolony){
  
  ccn = as.character(colonylist[which(colonylist$col.f == cc),"colony"])
  netapredi = paste0("etapred[",cc,",",max(1,nyearspred-46),"]")
  
  base = unlist(poster[,netapredi, drop=FALSE]) 
  
  netapredi = paste0("etapred[",cc,",",nyearspred-10,"]")
  
  base10 = unlist(poster[,netapredi, drop=FALSE]) 
  
  
  netapredi = paste0("etapred[",cc,",",nyearspred,"]")
  
  post = unlist(poster[,netapredi, drop=FALSE]) 
  
  for(tt in c("long-term","short-term")){
    if(tt == "long-term"){
      baset = base
      ymint = indicesreg[max(1,nyearspred-46),"year"]
      nyrtrend = min(nyearspred,46) 
    }else{
      baset = base10
      ymint = 2006
      nyrtrend = 10
      }
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"start.year"] <- ymint
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"end.year"] <- 2016
    
    
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"trend"] <- round(median(((post/baset)^(1/(nyrtrend))-1)*100),3)
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"trend.lci"] <- round(quantile(((post/baset)^(1/(nyrtrend))-1)*100,0.025),3)
    trendcol[which(trendcol$col.f == cc & trendcol$time == tt),"trend.uci"] <- round(quantile(((post/baset)^(1/(nyrtrend))-1)*100,0.975),3)
    
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
  ymint = max(1,nrow(indicesreg)-46)
    popi = paste0("pop[",ymint,"]")
    
    base = unlist(poster[,popi, drop=FALSE]) 
    
    popi = paste0("pop[",nyearspred-10,"]")
    
    base10 = unlist(poster[,popi, drop=FALSE]) 
    
    popi = paste0("pop[",nyearspred,"]")
    
    post = unlist(poster[,popi, drop=FALSE]) 
    
  # }else{
  #   
  #   popi = paste0("popr[",ymin,"]")
  #   
  #   base = unlist(poster[,popi, drop=FALSE]) 
  #   
  #   popi = paste0("popr[",nyearspred-10,"]")
  #   
  #   base10 = unlist(poster[,popi, drop=FALSE]) 
  #   
  #   popi = paste0("popr[",nyearspred,"]")
  #   
  #   post = unlist(poster[,popi, drop=FALSE]) 
  #   
  # }
  # 
  # trendcol[which(trendcol$colony == cc),"long-term"] <- round(median(((post/base)^(1/(nyearspred-ymin))-1)*100),3)
  # trendcol[which(trendcol$colony == cc),"long-term.lci"] <- round(quantile(((post/base)^(1/(nyearspred-ymin))-1)*100,0.025),3)
  # trendcol[which(trendcol$colony == cc),"long-term.uci"] <- round(quantile(((post/base)^(1/(nyearspred-ymin))-1)*100,0.975),3)
  # 
  # trendcol[which(trendcol$colony == cc),"short-term"] <- round(median(((post/base10)^(1/(10))-1)*100),3)
  # trendcol[which(trendcol$colony == cc),"short-term.lci"] <- round(quantile(((post/base10)^(1/(10))-1)*100,0.025),3)
  # trendcol[which(trendcol$colony == cc),"short-term.uci"] <- round(quantile(((post/base10)^(1/(10))-1)*100,0.975),3)
  # 
    for(tt in c("long-term","short-term")){
      if(tt == "long-term"){
        baset = base
        ymint = indicesreg[max(1,nyearspred-46),"year"]
        nyrtrend = min(nyearspred,46) 
      }else{
        baset = base10
        ymint = 2006
        nyrtrend = 10
      }
      
      trendcol[which(trendcol$colony == cc & trendcol$time == tt),"start.year"] <- ymint
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"end.year"] <- nyearspred+min(spdat$Year-1)
    
    
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"trend"] <- round(median(((post/baset)^(1/(nyrtrend))-1)*100),3)
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"trend.lci"] <- round(quantile(((post/baset)^(1/(nyrtrend))-1)*100,0.025),3)
    trendcol[which(trendcol$colony == cc & trendcol$time == tt),"trend.uci"] <- round(quantile(((post/baset)^(1/(nyrtrend))-1)*100,0.975),3)
    
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



xplot = indicesreg$year
if(ss == "NOGA"){
  png(filename = paste0("time series models/gam/",ss," gam trajectory overplot NB weighted.png"),
    res = 300,
    height = 8,
    width = 12,
    units = "in") 
}else{
  png(filename = paste0("time series models/gam/",ss," gam trajectory overplot quasi poisson weighted.png"),
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
         x = 1957+year[which(colony == cc)],
         y = count[which(colony == cc)])
  text(colonies[cc],
       x = 2018,
       y = tmp[length(tmp)],
       col = rainbow(length(unique(colony)))[cc],
       pos = 4)
}


dev.off()
trendcol$species = ss
if(j == 1){
  indicespop <- indicesreg
  
  indicescolony <- indicescol
  trendscolony = trendcol
}else{
  indicespop = rbind(indicespop,indicesreg)
  indicescolony = rbind(indicescolony,indicescol)
  trendscolony = rbind(trendscolony,trendcol)
}                      


write.csv(indicescolony,"Atlantic colony level indices weighted qpoisson.csv")

write.csv(indicespop,"Atlantic population level indices weighted qpoisson.csv")

write.csv(trendscolony,"Atlantic population level trends weighted qpoisson.csv")


}#ss


# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#     testml <- gam(count~ 1+factor(colony) + s(year),
#                   #knots = list(knots = knots),
#                   family = "quasipoisson",
#                   method = "REML")
#                   
#   test <- jagam(count~ 1+factor(colony) + s(year,k = nknots+1),
#                 #knots = list(knots = knots),
#                 family = "poisson",
#                 file = paste("models/GAM model",ss,"nknotsblank.txt"))#,
#  
#   
#   predmat <- expand.grid(predyear = c(ymin:ymax),
#                          predcolony = 1:ncolony)
#   predmat$predcount = rpois(nrow(predmat),5)
#   testpred <- jagam(predcount ~ 1+factor(predcolony) + s(predyear,k = nknots+1),
#                     data = predmat,
#                     knots = test$pregam$smooth$qr,
#                     family = "poisson",
#                     file = paste("models/GAM model",ss,"nknotspredict.txt"))
#   names(predmat) = gsub(pattern = "pred",
#                         replacement = "",
#                         names(predmat))
#   
#   #  testo <- jagam(count~ 1+factor(colony) + s(year,k = nknotso+1),
#   #                #knots = list(knots = knots),
#   #                family = "poisson",
#   #                file = paste("models/GAM model",ss,"nknotsoblank.txt"))#,
#   # test9 <- jagam(count~ 1+factor(colony) + s(year,k = 9),
#   #                #knots = list(knots = knots),
#   #                family = "poisson",
#   #                file = paste("models/GAM model",ss,"nknots9blank.txt"))#,
#   # 
#   # 
#   # }else{
#   #   test <- jagam(count~ s(year,k = nknots+1),
#   #                 #knots = list(knots = knots),
#   #                 family = "poisson",
#   #                 file = paste("models/GAM model",ss,"nknotsblank.txt"))#,
#   #   testo <- jagam(count~ s(year,k = nknotso+1),
#   #                  #knots = list(knots = knots),
#   #                  family = "poisson",
#   #                  file = paste("models/GAM model",ss,"nknotsoblank.txt"))#,
#   #   test9 <- jagam(count~ s(year,k = 9),
#   #                  #knots = list(knots = knots),
#   #                  family = "poisson",
#   #                  file = paste("models/GAM model",ss,"nknots9blank.txt"))#,
#   #   
#   # }
#   
#   # 
#   #  gamfit2 <- jags.model(data = test$jags.data,
#   #                       file = paste0("models/GAM model ",ss,".txt"),
#   #                       n.chains = 3)
#   # 
#   # 
#   # adaptest <- adapt(object = gamfit2,
#   #                   n.iter = 1000)
#   # 
#   # while(adaptest == F){
#   #   adaptest <- adapt(object = gamfit2,
#   #                     n.iter = 1000)
#   #   
#   # }
#   # 
#   # gamfito <- jags.model(data = testo$jags.data,
#   #                       file = paste0("models/GAM model ",ss," nknotso nb.txt"),
#   #                       n.chains = 3)
#   # 
#   # 
#   # adaptest <- adapt(object = gamfito,
#   #                   n.iter = 1000)
#   # 
#   # while(adaptest == F){
#   #   adaptest <- adapt(object = gamfito,
#   #                     n.iter = 1000)
#   #   
#   # }
#   # 
#   # 
#   # gamfit9 <- jags.model(data = test9$jags.data,
#   #                       file = paste0("models/GAM model ",ss," nknots9 nb.txt"),
#   #                       n.chains = 3)
#   # 
#   # 
#   # adaptest <- adapt(object = gamfit9,
#   #                   n.iter = 10)
#   # 
#   # while(adaptest == F){
#   #   adaptest <- adapt(object = gamfit9,
#   #                     n.iter = 1000)
#   #   
#   # }
#   # 
#   # 
#   
#   
#   
#   dat = test$jags.data
#   
#   
#   
#   dat[["Xp"]] <- testpred$jags.data$X
#   dat[["npred"]] <- nrow(testpred$jags.data$X)
#   
#   gamfitqp <- jags.model(data = dat,
#                         file = paste0("models/GAM model ",ss," quasi2.txt"),
#                         n.chains = 3)
#   
#   
#   adaptest <- adapt(object = gamfitqp,
#                     n.iter = 10)
#   
#   while(adaptest == F){
#     adaptest <- adapt(object = gamfitqp,
#                       n.iter = 1000)
#     
#   }
#   
#     # 
#     # gamfitsamcoda <- coda.samples(gamfit2,
#     #                           c("b","rho","mu","r"),
#     #                           n.iter = 10000,
#     #                           #inits = test$jags.ini,
#     #                           n.burnin = 5000,
#     #                           thin = 5)
#     # 
#     # gamfitsamocoda <- coda.samples(gamfito,
#     #                            c("b","rho","mu","r"),
#     #                            n.iter = 10000,
#     #                            #inits = test$jags.ini,
#     #                            n.burnin = 5000,
#     #                            thin = 5)
#     # gamfitsam9coda <- coda.samples(gamfit9,
#     #                            c("b","rho","mu","r"),
#     #                            n.iter = 10000,
#     #                            #inits = test$jags.ini,
#     #                            n.burnin = 5000,
#     #                            thin = 5)
#     gamfitsamqpcoda <- coda.samples(gamfitqp,
#                                   c("b","rho","mu","sdnoise","etapred"),
#                                   n.iter = 20000,
#                                   #inits = test$jags.ini,
#                                   n.burnin = 10000,
#                                   thin = 5)
#     
# # gamfitsam <- jags.samples(gamfit2,
# #                           c("b","rho","mu"),
# #                           n.iter = 5000,
# #                           #inits = test$jags.ini,
# #                           thin = 5)
# # 
# # gamfitsamo <- jags.samples(gamfito,
# #                           c("b","rho","mu"),
# #                           n.iter = 5000,
# #                           #inits = test$jags.ini,
# #                           thin = 5)
# # gamfitsam9 <- jags.samples(gamfit9,
# #                            c("b","rho","mu"),
# #                            n.iter = 5000,
# #                            #inits = test$jags.ini,
# #                            thin = 5)
# gamfitsamqp <- jags.samples(gamfitqp,
#                                 c("b","rho","mu","sdnoise","etapred"),
#                                 n.iter = 5000,
#                                 thin = 5)
# 
# # }else{
# #   gamfitsamcoda <- coda.samples(gamfit2,
# #                             c("b","rho","mu","r"),
# #                             n.iter = 10000,
# #                             #inits = test$jags.ini,
# #                             n.burnin = 5000,
# #                             thin = 5)  
# #   gamfitsamocoda <- coda.samples(gamfito,
# #                             c("b","rho","mu","r"),
# #                             n.iter = 10000,
# #                             #inits = test$jags.ini,
# #                             n.burnin = 5000,
# #                             thin = 5)  
# #   gamfitsam9coda <- coda.samples(gamfit9,
# #                              c("b","rho","mu","r"),
# #                              n.iter = 10000,
# #                              #inits = test$jags.ini,
# #                              n.burnin = 5000,
# #                              thin = 5)
# #   gamfitsam <- jags.samples(gamfit2,
# #                             c("b","rho","mu","r"),
# #                             n.iter = 5000,
# #                             #inits = test$jags.ini,
# #                             thin = 5)  
# #   gamfitsamo <- jags.samples(gamfito,
# #                              c("b","rho","mu","r"),
# #                              n.iter = 5000,
# #                              #inits = test$jags.ini,
# #                              thin = 5)  
# #   gamfitsam9 <- jags.samples(gamfit9,
# #                              c("b","rho","mu","r"),
# #                              n.iter = 5000,
# #                              #inits = test$jags.ini,
# #                              thin = 5)
# #   
# # }
# # 
# # preds.mu <- as.vector(gamfit$BUGSoutput$summary[paste0("mu[",1:64,"]"),"50%"])
# # preds.lci <- as.vector(gamfit$BUGSoutput$summary[paste0("mu[",1:64,"]"),"2.5%"])
# # preds.uci <- as.vector(gamfit$BUGSoutput$summary[paste0("mu[",1:64,"]"),"97.5%"])
# 
#   # pdf(paste0("time series models/gam/",ss,"best mcmc diagnostic.pdf"))
#   # plot(gamfitsamcoda)
#   # dev.off()
#   # 
#   # pdf(paste0("time series models/gam/",ss,"nknotso mcmc diagnostic.pdf"))
#   # plot(gamfitsamocoda)
#   # dev.off()
#   # 
#   # pdf(paste0("time series models/gam/",ss,"nknots9 mcmc diagnostic.pdf"))
#   # plot(gamfitsam9coda)
#   # dev.off()
#   pdf(paste0("time series models/gam/",ss,"quasi mcmc diagnostic.pdf"))
#   plot(gamfitsamqpcoda)
#   dev.off()
# 
#   
#   
#   #gamfitsamqpcoda  
#   
# # pd <- data.frame(expand.grid(colony = unique(colony),
# #                              year = min(year):(2016-1957)))
# 
# netapred = paste0("etapred[",1:nrow(predmat),"]")
# 
# etapred = gamfitsamqpcoda[,netapred]
# 
# 
# trendcol = data.frame(colony = 1:ncolony,
#                       colonyname = as.character(colonies),
#                       stringsAsFactors = F)
# trendcol[ncolony+1,"colony"] <- NA
# trendcol[ncolony+1,"colonyname"] <- "regional"
# 
# indicescol <- predmat[,c("year","colony")]
# indicesreg = data.frame(year = ymin:ymax)
# 
# basel = list()
# length(basel) <- ncolony
# base10l = basel
# 
# for(y in ymin:ymax){
# for(cc in 1:ncolony){
# ip = which(indicescol$colony == cc & indicescol$year == y)
# 
# post =  unlist(gamfitsamqpcoda[,netapred[ip]])
# if(cc == 1){
#   postsave = data.frame(post)
# }else{
#   postsave[,cc] <- post
# }
# ## indices
# indicescol[ip,"index"] <- round(median(post))
# indicescol[ip,"index.lci"] <- round(quantile(post,0.025))
# indicescol[ip,"index.uci"] <- round(quantile(post,0.975))
# indicescol[ip,"CalYear"] <- y+1957
# indicescol[ip,"species"] <- ss 
# indicescol[ip,"colonyname"] <- colonies[cc] 
# 
# if(y == ymin){
#   basel[[cc]] = post
# }
# if (y == ymax-10){
#   base10l[[cc]] = post
# }
# if(y == ymax){
#   base = basel[[cc]]
#   base10 = base10l[[cc]]
#   trendcol[which(trendcol$colony == cc),"long-term"] <- round(median(((post/base)^(1/(ymax-ymin))-1)*100),3)
#   trendcol[which(trendcol$colony == cc),"long-term.lci"] <- round(quantile(((post/base)^(1/(ymax-ymin))-1)*100,0.025),3)
#   trendcol[which(trendcol$colony == cc),"long-term.uci"] <- round(quantile(((post/base)^(1/(ymax-ymin))-1)*100,0.975),3)
#  
#   trendcol[which(trendcol$colony == cc),"short-term"] <- round(median(((post/base10)^(1/(10))-1)*100),3)
#   trendcol[which(trendcol$colony == cc),"short-term.lci"] <- round(quantile(((post/base10)^(1/(10))-1)*100,0.025),3)
#   trendcol[which(trendcol$colony == cc),"short-term.uci"] <- round(quantile(((post/base10)^(1/(10))-1)*100,0.975),3)
#   
# }
#   
#   
# }#cc
#   
#   postr = rowSums(postsave)
#   ipr = which(indicesreg$year == y)
#   
#   indicesreg[ipr,"index"] <- round(median(postr))
#   indicesreg[ipr,"index.lci"] <- round(quantile(postr,0.025))
#   indicesreg[ipr,"index.uci"] <- round(quantile(postr,0.975))
#   indicesreg[ipr,"CalYear"] <- y+1957
#   indicesreg[ipr,"species"] <- ss 
#   
#   
#   if(y == ymin){
#     baser = rowSums(postsave)
#   }
#   if (y == ymax-10){
#     base10r = rowSums(postsave)
#   }
#   if(y == ymax){
#     trendcol[which(trendcol$colonyname == "regional"),"long-term"] <- round(median(((postr/baser)^(1/(ymax-ymin))-1)*100),3)
#     trendcol[which(trendcol$colonyname == "regional"),"long-term.lci"] <- round(quantile(((postr/baser)^(1/(ymax-ymin))-1)*100,0.025),3)
#     trendcol[which(trendcol$colonyname == "regional"),"long-term.uci"] <- round(quantile(((postr/baser)^(1/(ymax-ymin))-1)*100,0.975),3)
#     
#     trendcol[which(trendcol$colonyname == "regional"),"short-term"] <- round(median(((postr/base10r)^(1/(10))-1)*100),3)
#     trendcol[which(trendcol$colonyname == "regional"),"short-term.lci"] <- round(quantile(((postr/base10r)^(1/(10))-1)*100,0.025),3)
#     trendcol[which(trendcol$colonyname == "regional"),"short-term.uci"] <- round(quantile(((postr/base10r)^(1/(10))-1)*100,0.975),3)
#     
#   }
#   
#   
# }#y
# 
# xplot = indicesreg$year+1957
# 
# plot(1,1,
#      xlim = c(1960,2040),
#      ylim = c(0,max(indicesreg$index*1.6)),
#      type = "l",
#      main = ss,
#      bty = "l",
#      ylab = "population by colony and total",
#      xlab = "")
# 
# polygon(x = c(xplot,rev(xplot)),
#         y = c(indicesreg$index.uci,rev(indicesreg$index.lci)),
#         col = transp.func(grey(0.8),0.6),
#         border = NA)
# lines(y = indicesreg$index,
#       x = xplot,
#       col = "black")
# 
# }
# for(cc in 1:ncolony){
#   tmp = indicescol[which(indicescol$colony == cc),"index"]
#   tmplci = indicescol[which(indicescol$colony == cc),"index.lci"]
#   tmpuci = indicescol[which(indicescol$colony == cc),"index.uci"]
#   #coli = unique(indicescol[which(indicescol$colony == j),"colonyname"])
# polygon(x = c(xplot,rev(xplot)),
#         y = c(tmplci,rev(tmpuci)),
#         col = transp.func(rainbow(length(unique(colony)))[cc],0.1),
#         border = NA)
# lines(y = tmp,
#       x = xplot,
#       col = rainbow(length(unique(colony)))[cc])
# points(col = rainbow(length(unique(colony)))[cc],
#        x = 1957+year[which(colony == cc)],
#        y = count[which(colony == cc)])
# text(colonies[cc],
#      x = 2018,
#      y = tmp[length(tmp)],
#      col = rainbow(length(unique(colony)))[cc],
#      pos = 4)
# }
# 
# 
# 
# 
# 
# 
# 
# # gamjam <- sim2jam(sam = gamfitsam,pregam = test$pregam)
# # predplot <- predict(gamjam,newdata = pd,
# #                     type = "response",
# #                     se.fit = T)
# # 
# # gamjamo <- sim2jam(sam = gamfitsamo,pregam = testo$pregam)
# # predploto <- predict(gamjamo,newdata = pd,
# #                     type = "response",
# #                     se.fit = T)
# # gamjam9 <- sim2jam(sam = gamfitsam9,pregam = test9$pregam)
# # predplot9 <- predict(gamjam9,newdata = pd,
# #                      type = "response",
# #                      se.fit = T)
# 
# gamjamqp <- sim2jam(sam = gamfitsamqp,pregam = test$pregam)
# plot(gamjamqp,
#      trans = poisson()$linkinv,
#      shift = 1)#sum(coef(gamjamqp)[1]))#,
#      seWithMean = T)
# 
# predplotqp <- predict(gamjamqp,newdata = pd,
#                     type = "link",
#                     se.fit = T)
# # 
# # 
# # save(file = paste0("time series models/gam/",ss," gam jags save quasi.RData"),
# #      list = ls()[which( ls() %in% c("predplot",
# #               "pd",
# #               "gamfitsam",
# #               "test",
# #               "gamjam",
# #               "predploto",
# #               "gamfitsamo",
# #               "testo",
# #               "gamjamo",
# #               "predplot9",
# #               "gamfitsam9",
# #               "test9",
# #               "gamjam9",
# #               "predplotqp",
# #               "gamfitsamqp",
# #               "testqp",
# #               "gamjamqp"))])
# 
# 
# png(filename = paste0("time series models/gam/",ss," gam trajectory overplot quasi poisson cor.png"),
#     res = 300,
#     height = 8,
#     width = 12,
#     units = "in") 
# par(mar = c(3,3,3,12))
# 
# 
# predplot <- pd
# predplot[,"fit"] <- predplotqp$fit
# predplot[,"se.fit"] <- predplotqp$se.fit
# 
# predplot[,"fit.uci"] <- predplotqp$fit+(2*predplotqp$se.fit)
# predplot[,"fit.lci"] <- predplotqp$fit-(2*predplotqp$se.fit)
# 
# predplot[,"fit.resp"] <- exp(predplot$fit)
# predplot[,"fit.lci.resp"] <- exp(predplot$fit.lci)
# predplot[,"fit.uci.resp"] <- exp(predplot$fit.uci)
# 
# 
# tmpall <- tapply(exp(predplot$fit),pd$year,sum)
# tmpalllog <- log(tmpall)
# sqsum <- function(x){
#   y = sqrt(sum(x^2))
# }
# tmpseall <- tapply(predplot$se.fit,pd$year,sqsum)
# 
# tmpalllci <- vector(length = length(tmpall))
# tmpalluci <- vector(length = length(tmpall))
# for(y in 1:length(tmpseall)){
#   tmpalluci[y] <- exp(tmpalllog[y]+(1.96*tmpseall[y]))
#   tmpalllci[y] <- exp(tmpalllog[y]-(1.96*tmpseall[y]))
#   }
# 
# for(cc in sort(unique(colony))){
#   tmp <- predplot$fit.resp[which(pd$colony == cc)]
#   tmplci <- predplot$fit.lci.resp[which(pd$colony == cc)]
#   tmpuci <- predplot$fit.uci.resp[which(pd$colony == cc)]
#   xplot <- pd[which(pd$colony == cc),"year"]+1957
#   if(cc == 1){
#     plot(1,1,
#          xlim = c(1960,2040),
#          ylim = c(0,max(tmpall*1.6)),
#          type = "l",
#          main = ss,
#          bty = "l",
#          ylab = "population by colony and total",
#          xlab = "")
#     
#     polygon(x = c(xplot,rev(xplot)),
#             y = c(tmpalllci,rev(tmpalluci)),
#             col = transp.func(grey(0.8),0.6),
#             border = NA)
#     lines(y = tmpall,
#           x = xplot,
#           col = "black")
#     
#   }
#   polygon(x = c(xplot,rev(xplot)),
#           y = c(tmplci,rev(tmpuci)),
#           col = transp.func(rainbow(length(unique(colony)))[cc],0.1),
#           border = NA)
#   lines(y = tmp,
#         x = xplot,
#         col = rainbow(length(unique(colony)))[cc])
#   points(col = rainbow(length(unique(colony)))[cc],
#          x = 1957+year[which(colony == cc)],
#          y = count[which(colony == cc)])
#   text(colonies[cc],
#        x = 2018,
#        y = tmp[length(tmp)],
#        col = rainbow(length(unique(colony)))[cc],
#        pos = 4)
# }
# 
# 
# 
# dev.off()
# 
# 
# 
# 
# 
# inds <- data.frame(species = ss,
#                    years = min(year):(2016-1957),
#                    index = tmpall,
#                    indexlci = tmpalllci,
#                    indexuci = tmpalluci,
#                    nknots = nknots+1,
#                    counts = ncounts,
#                    colony = ncolony)
# if(j == 1){
#   indices <- inds
# }else{
#   indices = rbind(indices,inds)
# }                      
# 
# 
# # 
# # 
# # if(paic > nbaic){
# #   png(filename = paste0("time series models/gam/",ss," gam trajectory overplot nb nknotso.png"),
# #       res = 300,
# #       height = 8,
# #       width = 12,
# #       units = "in") 
# # }else{    
# #   
# #   png(filename = paste0("time series models/gam/",ss," gam trajectory overplot Poisson nknotso.png"),
# #       res = 300,
# #       height = 8,
# #       width = 12,
# #       units = "in") 
# # }
# # par(mar = c(3,3,3,12))
# # tmpall <- tapply(predploto$fit,pd$year,sum)
# # sqsum <- function(x){
# #   y = sqrt(sum(x^2))
# # }
# # tmpseall <- tapply(predploto$se.fit,pd$year,sqsum)
# # 
# # for(cc in sort(unique(colony))){
# #   tmp <- predploto$fit[which(pd$colony == cc)]
# #   tmpse <- predploto$se.fit[which(pd$colony == cc)]
# #   xplot <- pd[which(pd$colony == cc),"year"]+1957
# #   if(cc == 1){
# #     plot(1,1,
# #          xlim = c(1960,2040),
# #          ylim = c(0,(max(tmpall) + max(tmpseall*2))*1.1),
# #          type = "l",
# #          main = ss,
# #          bty = "l",
# #          ylab = "population by colony and total",
# #          xlab = "")
# #     
# #     polygon(x = c(xplot,rev(xplot)),
# #             y = c((tmpall-(2*tmpseall)),rev((tmpall+(2*tmpseall)))),
# #             col = transp.func(grey(0.8),0.6),
# #             border = NA)
# #     lines(y = tmpall,
# #           x = xplot,
# #           col = "black")
# #     
# #   }
# #   polygon(x = c(xplot,rev(xplot)),
# #           y = c((tmp-(2*tmpse)),rev((tmp+(2*tmpse)))),
# #           col = transp.func(rainbow(length(unique(colony)))[cc],0.1),
# #           border = NA)
# #   lines(y = tmp,
# #         x = xplot,
# #         col = rainbow(length(unique(colony)))[cc])
# #   points(col = rainbow(length(unique(colony)))[cc],
# #          x = 1957+year[which(colony == cc)],
# #          y = count[which(colony == cc)])
# #   text(colonies[cc],
# #        x = 2018,
# #        y = tmp[length(tmp)],
# #        col = rainbow(length(unique(colony)))[cc],
# #        pos = 4)
# # }
# # 
# # 
# # 
# # dev.off()
# # 
# # 
# # 
# # 
# # indso <- data.frame(species = ss,
# #                    years = min(year):(2016-1957),
# #                    index = tmpall,
# #                    indexse = tmpseall,
# #                    nknots = nknots+1,
# #                    counts = ncounts,
# #                    colony = ncolony)
# # if(j == 1){
# #   indiceso <- indso
# # }else{
# #   indiceso = rbind(indiceso,indso)
# # }                      
# # 
# # 
# # 
# # 
# # if(paic > nbaic){
# #   png(filename = paste0("time series models/gam/",ss," gam trajectory overplot nb nknots9.png"),
# #       res = 300,
# #       height = 8,
# #       width = 12,
# #       units = "in") 
# # }else{    
# #   
# #   png(filename = paste0("time series models/gam/",ss," gam trajectory overplot Poisson nknots9.png"),
# #       res = 300,
# #       height = 8,
# #       width = 12,
# #       units = "in") 
# # }
# # par(mar = c(3,3,3,12))
# # tmpall <- tapply(predplot9$fit,pd$year,sum)
# # sqsum <- function(x){
# #   y = sqrt(sum(x^2))
# # }
# # tmpseall <- tapply(predplot9$se.fit,pd$year,sqsum)
# # 
# # for(cc in sort(unique(colony))){
# #   tmp <- predplot9$fit[which(pd$colony == cc)]
# #   tmpse <- predplot9$se.fit[which(pd$colony == cc)]
# #   xplot <- pd[which(pd$colony == cc),"year"]+1957
# #   if(cc == 1){
# #     plot(1,1,
# #          xlim = c(1960,2040),
# #          ylim = c(0,(max(tmpall) + max(tmpseall*2))*1.1),
# #          type = "l",
# #          main = ss,
# #          bty = "l",
# #          ylab = "population by colony and total",
# #          xlab = "")
# #     
# #     polygon(x = c(xplot,rev(xplot)),
# #             y = c((tmpall-(2*tmpseall)),rev((tmpall+(2*tmpseall)))),
# #             col = transp.func(grey(0.8),0.6),
# #             border = NA)
# #     lines(y = tmpall,
# #           x = xplot,
# #           col = "black")
# #     
# #   }
# #   polygon(x = c(xplot,rev(xplot)),
# #           y = c((tmp-(2*tmpse)),rev((tmp+(2*tmpse)))),
# #           col = transp.func(rainbow(length(unique(colony)))[cc],0.1),
# #           border = NA)
# #   lines(y = tmp,
# #         x = xplot,
# #         col = rainbow(length(unique(colony)))[cc])
# #   points(col = rainbow(length(unique(colony)))[cc],
# #          x = 1957+year[which(colony == cc)],
# #          y = count[which(colony == cc)])
# #   text(colonies[cc],
# #        x = 2018,
# #        y = tmp[length(tmp)],
# #        col = rainbow(length(unique(colony)))[cc],
# #        pos = 4)
# # }
# # 
# # 
# # 
# # dev.off()
# # 
# # 
# # inds9 <- data.frame(species = ss,
# #                    years = min(year):(2016-1957),
# #                    index = tmpall,
# #                    indexse = tmpseall,
# #                    nknots = nknots+1,
# #                    counts = ncounts,
# #                    colony = ncolony)
# # if(j == 1){
# #   indices9 <- inds9
# # }else{
# #   indices9 = rbind(indices9,inds9)
# # }                      
# # 
# 
# 
# 
# write.csv(indices,"GAM indices quasi.csv",row.names = F)
# # write.csv(indiceso,"GAM indices nknotso.csv",row.names = F)
# # write.csv(indices9,"GAM indices nknots9.csv",row.names = F)
# 
# }#ss
# 
# ####################### on Monday am
# 
# ### add in the quasi poisson models
# 
# 
# 
# 
# 
# 
# # 
# # 
# # 
# # 
# # 
# # year.s <- 1
# # year.e <- 2018-1957## setting the outer bounds for predicted value
# # ymiss.s <- (ymin-1):year.s
# # ymiss.e <- (ymax+1):year.e
# # 
# # 
# # ### estimates of the log-scale means and precisions
# # ### of colony averages, assuming the colony means
# # ### have only Poisson-distribution variance (var/mean = 1)
# # ### and so ignoring the variance among years at each colony
# # colmn <- tapply(count,colony,mean,na.rm = T)
# # colcvar <- colmn/(colmn^2) ## cv of a poisson count variance/mean^2
# # colhat <- log((colmn^2)/sqrt(colmn+colmn^2)) # mean = log(mean^2/sqrt(var + mean^2))
# # taucolhat <- 1/(log(1+colcvar)) # var = log(1+ var/mean^2)
# # ## above equations are the standard equations for the mean and
# # ## variance of a lognormal distribution, derived from the mean
# # ## and variance of the normally distributed variate
# # ### Johnson, Norman L.; Kotz, Samuel; Balakrishnan, N. (1994), "14: Lognormal Distributions", Continuous univariate distributions. Vol. 1, Wiley Series in Probability and Mathematical Statistics: Applied Probability and Statistics (2nd ed.), New York: John Wiley & Sons, ISBN 978-0-471-58495-7, MR 1299979
# # ### 
# # 
# # upperchange <- 25#% change that has a 0.1% prior probability
# # upperchange01 <- upperchange/(qnorm(0.9995)) 
# # 
# # sdyear <- log((upperchange01/100)+1)
# # 
# # # x = rnorm(10000,0,upperchange01l)
# # # hist(x)
# # # (exp(quantile(x,probs = c(0.001,0.999)))-1)*100
# # 
# # 
# # 
# # jagsdat <- list(colony = colony,
# #                 count = count,
# #                 ncounts = ncounts,
# #                 ncolony = ncolony,
# #                 colhat = colhat,
# #                 taucolhat = taucolhat,
# #                 ymax = ymax,
# #                 ymin = ymin,
# #                 year = year,
# #                 ymiss.s = ymiss.s,
# #                 ymiss.e = ymiss.e,
# #                 year.e = year.e,
# #                 sdyear = sdyear)
# # 
# # params = c("yeareffect",
# #            "n",
# #            "sdyear",
# #            "lambda",
# #            "col")
# # 
# # adaptSteps = 500              # Number of steps to "tune" the samplers.
# # burnInSteps = 100000            # Number of steps to "burn-in" the samplers.
# # nChains = 3                   # Number of chains to run.
# # numSavedSteps=2000           # Total number of steps to save.
# # thinSteps=100                   # Number of steps to "thin" (1=keep every step).
# # nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# # 
# # t1 = Sys.time()
# # 
# # 
# # jagsMod = jags.model( difmodel, 
# #                       data= jagsdat ,   
# #                       n.chains= nChains , 
# #                       n.adapt= adaptSteps )
# # 
# # # Burn-in:
# # cat( "Burning in the MCMC chain...\n" )
# # update( jagsMod , n.iter=burnInSteps )
# # # The saved MCMC chain:
# # cd13 = coda.samples( jagsMod , variable.names=params ,
# #                      n.iter=nIter , thin=thinSteps )
# # # 
# # # 
# # # 
# # # cdlambda = coda.samples( jagsMod , variable.names=list("lambda") ,
# # #                      n.iter=nIter , thin=thinSteps )
# # 
# # 
# # t2 = Sys.time()
# # print(paste(round(t2-t1,2),ss,"finished"))
# # pdf(paste0("time series models/fixed prior sdyear/difference/",ss," difference diagnostic plot.pdf"))
# # plot(cd13)
# # dev.off()
# # 
# # lambdanms <- paste0("lambda[",1:ncounts,"]")
# # lambdal <- cd13[,lambdanms,drop = F]
# #  
# # waics[[ss]]["difference"] <- waic()
# # 
# # sum <- summary(cd13)
# # 
# # qsum <- sum$quantiles
# # 
# # wind <- paste0("n[",c(1:year.e),"]")
# # inds <- as.data.frame(qsum[wind,])
# # inds[,"species"] <- ss
# # inds[,"node"] <- as.character(row.names(inds))
# # inds[,"model"] <- "difference"
# # if(ss == sptomod[1]){indicesout = inds}else{indicesout = rbind(indicesout,inds)}
# # write.csv(inds,paste0("time series models/fixed prior sdyear/difference/",ss," difference annual indices.csv"))
# # write.csv(qsum,paste0("time series models/fixed prior sdyear/difference/",ss," difference posterior quantiles.csv"))
# # 
# # xplot <- c((1+1957):(year.e+1957))
# # 
# # png(filename = paste0("time series models/fixed prior sdyear/difference/",ss," difference trajectory overplot.png"),
# #     res = 300,
# #     height = 8,
# #     width = 8,
# #     units = "in") 
# # 
# # plot(x = xplot,
# #      y = inds[,"50%"],
# #      ylim = c(0,min(max(inds[,"97.5%"]),max(inds[,"50%"])*3)),
# #      type = "l")
# # polygon(x = c(xplot,rev(xplot)),
# #         y = c(inds[,"97.5%"],rev(inds[,"2.5%"])),
# #         col = grey(0.6),
# #         border = NA)
# # lines(x = xplot,
# #      y = inds[,"50%"])
# # colcols <- rainbow(ncolony)
# # plotcols <- colcols[colony]
# # 
# # points(x = year+1957,
# #        y = count,
# #        col = plotcols)
# # dev.off()
# # 
# # save(file = paste0("time series models/fixed prior sdyear/difference/",ss," difference post jags.RData"),
# #      list = c("cd13",
# #               "spdat",
# #               "jagsdat"))
# # 
# # }#ss
# 
# 
# # 
# # seconddifmodel <- c("seconddifmodel constrained year.txt")
# # 
# # for(ss in sptomod){
# #   
# #   
# #   
# #   spdat <- allpost60[which(allpost60$species == ss),]
# #   colony <- as.integer(factor(spdat$colony))
# #   count <- spdat$estimate
# #   ncounts = length(count)
# #   ncolony <- max(colony)
# #   year <- spdat$Year - 1957
# #   ymax <- max(year)
# #   ymin <- min(year)
# #  
# #   
# #    
# #   
# #   year.s <- 1
# #   year.e <- 2018-1957## setting the outer bounds for predicted values
# #   ymiss.s <- (ymin-1):year.s
# #   ymiss.e <- (ymax+1):year.e
# #   
# #   
# #   ### estimates of the log-scale means and precisions
# #   ### of colony averages, assuming the colony means
# #   ### have only Poisson-distribution variance (var/mean = 1)
# #   ### and so ignoring the variance among years at each colony
# #   colmn <- tapply(count,colony,mean,na.rm = T)
# #   colcvar <- colmn/(colmn^2) ## cv of a poisson count variance/mean^2
# #   colhat <- log((colmn^2)/sqrt(colmn+colmn^2)) # mean = log(mean^2/sqrt(var + mean^2))
# #   taucolhat <- 1/(log(1+colcvar)) # var = log(1+ var/mean^2)
# #   ## above equations are the standard equations for the mean and
# #   ## variance of a lognormal distribution, derived from the mean
# #   ## and variance of the normally distributed variate
# #   ### Johnson, Norman L.; Kotz, Samuel; Balakrishnan, N. (1994), "14: Lognormal Distributions", Continuous univariate distributions. Vol. 1, Wiley Series in Probability and Mathematical Statistics: Applied Probability and Statistics (2nd ed.), New York: John Wiley & Sons, ISBN 978-0-471-58495-7, MR 1299979
# #   ### 
# #   
# #   upperchange <- 25#% change that has a 0.1% prior probability
# #   upperchange01 <- upperchange/(qnorm(0.9995)) 
# #   
# #   sdyear <- log((upperchange01/100)+1)
# #   
# #   #
# #   jagsdat <- list(colony = colony,
# #                   count = count,
# #                   ncounts = ncounts,
# #                   ncolony = ncolony,
# #                   colhat = colhat,
# #                   taucolhat = taucolhat,
# #                   ymax = ymax,
# #                   ymin = ymin,
# #                   year = year,
# #                   ymiss.s = ymiss.s,
# #                   ymiss.e = ymiss.e,
# #                   year.e = year.e,
# #                   sdyear = sdyear)
# #   
# #   params = c("yeareffect",
# #              "n",
# #              "sdyear",
# #              "lambda",
# #              "col")
# #   
# #   adaptSteps = 500              # Number of steps to "tune" the samplers.
# #   burnInSteps = 100000            # Number of steps to "burn-in" the samplers.
# #   nChains = 3                   # Number of chains to run.
# #   numSavedSteps=2000           # Total number of steps to save.
# #   thinSteps=100                   # Number of steps to "thin" (1=keep every step).
# #   nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# #   
# #   t1 = Sys.time()
# #   
# #   
# #   jagsMod = jags.model( seconddifmodel, 
# #                         data= jagsdat ,   
# #                         n.chains= nChains , 
# #                         n.adapt= adaptSteps )
# #   
# #   # Burn-in:
# #   cat( "Burning in the MCMC chain...\n" )
# #   update( jagsMod , n.iter=burnInSteps )
# #   # The saved MCMC chain:
# #   cd13 = coda.samples( jagsMod , variable.names=params ,
# #                        n.iter=nIter , thin=thinSteps )
# #   
# #   t2 = Sys.time()
# #   pdf(paste0("time series models/fixed prior sdyear/second difference/",ss," second difference diagnostic plot.pdf"))
# #   plot(cd13)
# #   dev.off()
# #   
# #   lambdanms <- paste0("lambda[",1:ncounts,"]")
# #   lambdal <- cd13[,lambdanms,drop = F]
# #   
# #   waics[[ss]]["second difference"] <- waic()
# #   
# #   sum <- summary(cd13)
# #   
# #   qsum <- sum$quantiles
# #   
# #   wind <- paste0("n[",c(1:year.e),"]")
# #   
# #   xplot <- c((1+1957):(year.e+1957))
# #   inds <- as.data.frame(qsum[wind,])
# #   inds[,"species"] <- ss
# #   inds[,"node"] <- as.character(row.names(inds))
# #   inds[,"model"] <- "second difference"
# #   
# # indicesout = rbind(indicesout,inds)
# # 
# #   
# #   write.csv(inds,paste0("time series models/fixed prior sdyear/second difference/",ss," second difference annual indices.csv"))
# #   write.csv(qsum,paste0("time series models/fixed prior sdyear/second difference/",ss," second difference posterior quantiles.csv"))
# #   
# # 
# #   png(filename = paste0("time series models/fixed prior sdyear/second difference/",ss," second difference trajectory overplot.png"),
# #       res = 300,
# #       height = 8,
# #       width = 8,
# #       units = "in") 
# #   plot(x = xplot,
# #        y = inds[,"50%"],
# #        ylim = c(0,min(max(inds[,"97.5%"]),max(inds[,"50%"])*3)),
# #        type = "l")
# #   polygon(x = c(xplot,rev(xplot)),
# #           y = c(inds[,"97.5%"],rev(inds[,"2.5%"])),
# #           col = grey(0.6),
# #           border = NA)
# #   lines(x = xplot,
# #         y = inds[,"50%"])
# #   colcols <- rainbow(ncolony)
# #   plotcols <- colcols[colony]
# #   
# #   points(x = year+1957,
# #          y = count,
# #          col = plotcols)
# #   dev.off()
# #   
# #   save(file = paste0("time series models/fixed prior sdyear/second difference/",ss," second difference post jags.RData"),
# #        list = c("cd13",
# #                 "spdat",
# #                 "jagsdat"))
# #   print(paste(round(t2-t1,2),ss,"finished"))
# #   
# # }#ss
# # 
# # 
# # 
# # 
# # 
# # slopemodel <- c("slopemodel constrained year.txt")
# # 
# # for(ss in sptomod){
# #   
# #   
# #   
# #   spdat <- allpost60[which(allpost60$species == ss),]
# #   colony <- as.integer(factor(spdat$colony))
# #   count <- spdat$estimate
# #   ncounts = length(count)
# #   ncolony <- max(colony)
# #   year <- spdat$Year - 1957
# #   ymax <- max(year)
# #   ymin <- min(year)
# #   fixedyear <- floor(ymax-(diff(c(ymin,ymax))/2))
# #   
# #   year.s <- 1
# #   year.e <- 2018-1957## setting the outer bounds for predicted values
# #   ymiss.s <- (ymin-1):year.s
# #   ymiss.e <- (ymax+1):year.e
# #   
# #   
# #   ### estimates of the log-scale means and precisions
# #   ### of colony averages, assuming the colony means
# #   ### have only Poisson-distribution variance (var/mean = 1)
# #   ### and so ignoring the variance among years at each colony
# #   colmn <- tapply(count,colony,mean,na.rm = T)
# #   colcvar <- colmn/(colmn^2) ## cv of a poisson count variance/mean^2
# #   colhat <- log((colmn^2)/sqrt(colmn+colmn^2)) # mean = log(mean^2/sqrt(var + mean^2))
# #   taucolhat <- 1/(log(1+colcvar)) # var = log(1+ var/mean^2)
# #   ## above equations are the standard equations for the mean and
# #   ## variance of a lognormal distribution, derived from the mean
# #   ## and variance of the normally distributed variate
# #   ### Johnson, Norman L.; Kotz, Samuel; Balakrishnan, N. (1994), "14: Lognormal Distributions", Continuous univariate distributions. Vol. 1, Wiley Series in Probability and Mathematical Statistics: Applied Probability and Statistics (2nd ed.), New York: John Wiley & Sons, ISBN 978-0-471-58495-7, MR 1299979
# #   ### 
# #   
# #   upperchange <- 25#% change that has a 0.1% prior probability
# #   upperchange01 <- upperchange/(qnorm(0.9995)) 
# #   
# #   sdyear <- log((upperchange01/100)+1)
# #   
# #   #
# #   jagsdat <- list(colony = colony,
# #                   count = count,
# #                   ncounts = ncounts,
# #                   ncolony = ncolony,
# #                   colhat = colhat,
# #                   taucolhat = taucolhat,
# #                   # ymax = ymax,
# #                   # ymin = ymin,
# #                   year = year,
# #                   # ymiss.s = ymiss.s,
# #                   # ymiss.e = ymiss.e,
# #                   year.e = year.e,
# #                   fixedyear = fixedyear,
# #                   sdyear = sdyear)
# #   
# #   
# #   # jagsdat <- list(colony = colony,
# #   #                 count = count,
# #   #                 ncounts = ncounts,
# #   #                 ncolony = ncolony,
# #   #                 ymax = ymax,
# #   #                 ymin = ymin,
# #   #                 year = year,
# #   #                 fixedyear = fixedyear)
# #   # 
# #   params = c("yeareffect",
# #              "n",
# #              "sdbeta",
# #              "col",
# #              "sdyear",
# #              "B",
# #              "beta",
# #              "lambda")
# #   
# #   adaptSteps = 1500              # Number of steps to "tune" the samplers.
# #   burnInSteps = 200000            # Number of steps to "burn-in" the samplers.
# #   nChains = 3                   # Number of chains to run.
# #   numSavedSteps=2000           # Total number of steps to save.
# #   thinSteps=100                   # Number of steps to "thin" (1=keep every step).
# #   nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
# #   
# #   t1 = Sys.time()
# #   
# #   
# #   jagsMod = jags.model(slopemodel, 
# #                         data= jagsdat ,   
# #                         n.chains= nChains , 
# #                         n.adapt= adaptSteps )
# #   
# #   # Burn-in:
# #   cat( "Burning in the MCMC chain...\n" )
# #   update( jagsMod , n.iter=burnInSteps )
# #   # The saved MCMC chain:
# #   cd13 = coda.samples( jagsMod , variable.names=params ,
# #                        n.iter=nIter , thin=thinSteps )
# #   
# #   t2 = Sys.time()
# #   print(paste(round(t2-t1,2),ss,"finished"))
# #   pdf(paste0("time series models/fixed prior sdyear/slope/",ss," slope diagnostic plot.pdf"))
# #   plot(cd13)
# #   dev.off()
# #   
# #   lambdanms <- paste0("lambda[",1:ncounts,"]")
# #   lambdal <- cd13[,lambdanms,drop = F]
# #   
# #   waics[[ss]]["slope"] <- waic()
# #   
# #   sum <- summary(cd13)
# #   
# #   qsum <- sum$quantiles
# #   
# #   
# #   xplot <- c((1+1957):(year.e+1957))
# #   
# #   wind <- paste0("n[",c(1:year.e),"]")
# #   inds <- as.data.frame(qsum[wind,])
# #   inds[,"species"] <- ss
# #   inds[,"node"] <- as.character(row.names(inds))
# #   inds[,"model"] <- "slope"
# #   
# #   indicesout = rbind(indicesout,inds)
# #   
# #   write.csv(inds,paste0("time series models/fixed prior sdyear/slope/",ss," slope annual indices.csv"))
# #   write.csv(qsum,paste0("time series models/fixed prior sdyear/slope/",ss," slope posterior quantiles.csv"))
# #   
# #  
# #   png(filename = paste0("time series models/fixed prior sdyear/slope/",ss," slope trajectory overplot.png"),
# #       res = 300,
# #       height = 8,
# #       width = 8,
# #       units = "in")
# #   plot(x = xplot,
# #        y = inds[,"50%"],
# #        ylim = c(0,min(max(inds[,"97.5%"]),max(inds[,"50%"])*3)),
# #        type = "l")
# #   polygon(x = c(xplot,rev(xplot)),
# #           y = c(inds[,"97.5%"],rev(inds[,"2.5%"])),
# #           col = grey(0.6),
# #           border = NA)
# #   lines(x = xplot,
# #         y = inds[,"50%"])
# #   colcols <- rainbow(ncolony)
# #   plotcols <- colcols[colony]
# #   
# #   points(x = year+1957,
# #          y = count,
# #          col = plotcols)
# #   dev.off()
# #   
# #   save(file = paste0("time series models/fixed prior sdyear/slope/",ss," slope post jags.RData"),
# #        list = c("cd13",
# #                 "spdat",
# #                 "jagsdat"))
# #   
# # }#ss
# # 
# # write.csv(indicesout,"time series models/fixed prior sdyear/all species indices.csv",row.names = F)
# # save(file = "time series models/fixed prior sdyear/waic summary.RData",
# #      list = "waics")
# # 
# # # then sum all time-series for a given species and year
# # 
# # ## consider the "standard-year" approach of the living planet index (move everything to an 5-year schedule)
# # 
# # 
# 
# 
# 
# 
# 
# 
# 
# 





