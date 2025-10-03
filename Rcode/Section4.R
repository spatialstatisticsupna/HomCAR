################################################################################################
## R code for the reproduce the analyses of the paper entitled "A proposal for homoskedastic  ##
## modelling with conditional auto-regressive distributions" (Martínez-Beneito et al., 2025). ##
################################################################################################
rm(list=ls())
library(INLA)  
library(LearnBayes)
library(MASS)
library(RColorBrewer)
library(sf)
library(spdep)
library(tmap)
library(units)

# Option for the management of cartographies
sf_use_s2(FALSE)


##################################################################
## SECTION 4: A PRACTICAL ASSESSMENT OF THE HOMCAR DISTRIBUTION ##
##################################################################
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

##################################################################################
## INLA-based function for computing marginal variances from a precision matrix ##
## i.e, diagonal elements of the Moore-Penrose inverse                          ##
##################################################################################
inla.ginv.diag <- function(Q, constr = NULL, eps = sqrt(.Machine$double.eps)) {
  
  marg.var <- rep(0, nrow(Q))
  Q <- inla.as.sparse(Q)
  g <- inla.read.graph(Q)
  
  if(is.null(constr)) constr <- list(A=matrix(1,1,nrow(Q)), e=0)
  
  for (k in seq_len(g$cc$n)) {
    i <- g$cc$nodes[[k]]
    n <- length(i)
    QQ <- Q[i, i, drop = FALSE]
    if (n == 1) {
      QQ[1, 1] <- 1
      marg.var[i] <- 1
    }
    else {
      cconstr <- constr
      if (!is.null(constr)) {
        cconstr$A <- constr$A[, i, drop = FALSE]
        eeps <- eps
      }
      else {
        eeps <- 0
      }
      idx.zero <- which(rowSums(abs(cconstr$A)) == 0)
      if (length(idx.zero) > 0) {
        cconstr$A <- cconstr$A[-idx.zero, , drop = FALSE]
        cconstr$e <- cconstr$e[-idx.zero]
      }
      res <- inla.qinv(QQ + Diagonal(n) * max(diag(QQ)) * 
                         eeps, constr = cconstr)
      
      # fac <- exp(mean(log(diag(res))))
      # QQ <- fac * QQ
      # marg.var[i] <- diag(res)/fac
      
      marg.var[i] <- diag(res)
    }
    Q[i,i] <- QQ
  }
  return(marg.var)
}


##############################################################
## Compute weighted precision matrices for the HomCAR prior ##
##############################################################
load("./Carto_files.Rdata")

prec.matrices <- lapply(list(VR=CartoVR, CL=CartoCL, AR=CartoAR, CM=CartoCM), function(x){
  nb <- poly2nb(st_sf(x))
  
  W <- nb2mat(nb, style="B")
  dimnames(W) <- NULL
  attr(W,"call") <- NULL
  
  Q <- diag(colSums(W))-W
  vars.marginal <- inla.ginv.diag(Q)
  
  Q.weighted <- t(sqrt(vars.marginal)*t(sqrt(vars.marginal)*Q))
  
  return(list(Q=Q, Q.weighted=Q.weighted))
})

## Summary statistics of the marginal variances for the HomCAR precision matrices ##
lapply(prec.matrices, function(x) summary(inla.ginv.diag(x$Q.weighted)))


###########################################
## Generation of the simulated data sets ##
###########################################
carto.list <- list('Valencian Region'=CartoVR,
                   'Castile and Leon'=CartoCL,
                   'Aragon'=CartoAR,
                   'Castile-La Mancha'=CartoCM)

n.sim <- 100
E <- 5

simu.data <- lapply(names(carto.list), function(x){
  
  carto <- carto.list[[x]]
  n.data <- nrow(carto)
  cat(x," (n=",n.data,")",sep="")
  
  cat("\n + generating data sets...")
  t1 <- system.time({
    
    ## Define the covariance matrix ##
    carto <- st_transform(carto, 25830)
    centroids <- suppressWarnings(st_centroid(carto))
    distMatrix <- drop_units(st_distance(centroids))
    covMatrix <- (0.2^2)*exp(-3*distMatrix/quantile(c(distMatrix),0.25))
    
    ## Generation of simulated data sets ##
    set.seed(1)
    log.risks <- t(rmnorm(n=n.sim, mean=rep(0,n.data), varcov=covMatrix))
    Obs <- apply(log.risks, 2, function(x) rpois(n.data, E*exp(x)))
    
  })
  cat(" elapsed time:",t1[3],"sec\n\n")

  return(list(Obs=Obs, log.risks=log.risks))
})
names(simu.data) <- c("VR","CL","AR","CM")


################################################################
## INLA calls for the traditional BYM and homoscedastic model ##
################################################################

## Hyperprior distribution for the standard deviations ##
sdunif="expression:
  logdens=log(0.5)-log_precision/2;
  return(logdens)"


## Valencian Region
####################
carto <- CartoVR
Q <- prec.matrices$VR$Q
Q.weighted <- prec.matrices$VR$Q.weighted

## Objects to save ##
post.mean.VR <- list(`CAR`=matrix(nrow=nrow(carto), ncol=n.sim),
                     `HomCAR`=matrix(nrow=nrow(carto), ncol=n.sim))

DIC.VR <- list(`CAR`=data.frame(DIC=numeric(n.sim), WAIC=numeric(n.sim), LS=numeric(n.sim)),
               `HomCAR`=data.frame(DIC=numeric(n.sim), WAIC=numeric(n.sim), LS=numeric(n.sim)))

## Fit models with INLA ##
for(i in 1:n.sim){
  print(sprintf("Simulation %d", i))
  
  data.INLA <- data.frame(O=simu.data$VR$Obs[,i], E=E,
                          ID.area=1:nrow(carto), ID.area2=1:nrow(carto))
  
  ## Heterocedastic BYM model ##
  f.CAR <- O ~ f(ID.area, model="bym", graph=Q, constr=TRUE,
                 hyper=list(prec.spatial=list(prior=sdunif),
                            prec.unstruct=list(prior=sdunif)))
  
  M1 <- inla(f.CAR, family="poisson", data=data.INLA, E=E,
             control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=FALSE),
             inla.mode="compact", num.threads="8:1")

  ## Homoscedastic BYM model ##
  f.HomCAR <- O ~ f(ID.area, model="generic0", Cmatrix=Q.weighted, constr=TRUE,
                    hyper=list(prec=list(prior=sdunif))) + 
    f(ID.area2, model="iid", constr=TRUE, hyper=list(prec=list(prior=sdunif)))
  
  M2 <- inla(f.HomCAR, family="poisson", data=data.INLA, E=E,
             control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=FALSE),
             inla.mode="compact", num.threads="8:1")
  
  ## Save results ##
  post.mean.VR$CAR[,i] <- M1$summary.linear.predictor$mean
  post.mean.VR$HomCAR[,i] <- M2$summary.linear.predictor$mean

  DIC.VR$CAR[i,] <- c(M1$dic$dic, M1$waic$waic, -sum(log(M1$cpo$cpo)))
  DIC.VR$HomCAR[i,] <- c(M2$dic$dic, M2$waic$waic, -sum(log(M2$cpo$cpo)))
}


## Castile and Leon
####################
carto <- CartoCL
Q <- prec.matrices$CL$Q
Q.weighted <- prec.matrices$CL$Q.weighted

## Objects to save ##
post.mean.CL <- list(`CAR`=matrix(nrow=nrow(carto), ncol=n.sim),
                     `HomCAR`=matrix(nrow=nrow(carto), ncol=n.sim))

DIC.CL <- list(`CAR`=data.frame(DIC=numeric(n.sim), WAIC=numeric(n.sim), LS=numeric(n.sim)),
               `HomCAR`=data.frame(DIC=numeric(n.sim), WAIC=numeric(n.sim), LS=numeric(n.sim)))


## Fit models with INLA ##
for(i in 1:n.sim){
  print(sprintf("Simulation %d", i))
  
  data.INLA <- data.frame(O=simu.data$CL$Obs[,i], E=E,
                          ID.area=1:nrow(carto), ID.area2=1:nrow(carto))
  
  ## Heterocedastic BYM model ##
  f.CAR <- O ~ f(ID.area, model="bym", graph=Q, constr=TRUE,
                 hyper=list(prec.spatial=list(prior=sdunif),
                            prec.unstruct=list(prior=sdunif)))
  
  M1 <- inla(f.CAR, family="poisson", data=data.INLA, E=E,
             control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=FALSE),
             inla.mode="compact", num.threads="8:1")
  
  ## Homoscedastic BYM model ##
  f.HomCAR <- O ~ f(ID.area, model="generic0", Cmatrix=Q.weighted, constr=TRUE,
                    hyper=list(prec=list(prior=sdunif))) + 
    f(ID.area2, model="iid", constr=TRUE, hyper=list(prec=list(prior=sdunif)))
  
  M2 <- inla(f.HomCAR, family="poisson", data=data.INLA, E=E,
             control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=FALSE),
             inla.mode="compact", num.threads="8:1")
  
  ## Save results ##
  post.mean.CL$CAR[,i] <- M1$summary.linear.predictor$mean
  post.mean.CL$HomCAR[,i] <- M2$summary.linear.predictor$mean
  
  DIC.CL$CAR[i,] <- c(M1$dic$dic, M1$waic$waic, -sum(log(M1$cpo$cpo)))
  DIC.CL$HomCAR[i,] <- c(M2$dic$dic, M2$waic$waic, -sum(log(M2$cpo$cpo)))
}


## Aragon
##########
carto <- CartoAR
Q <- prec.matrices$AR$Q
Q.weighted <- prec.matrices$AR$Q.weighted

## Objects to save ##
post.mean.AR <- list(`CAR`=matrix(nrow=nrow(carto), ncol=n.sim),
                     `HomCAR`=matrix(nrow=nrow(carto), ncol=n.sim))

DIC.AR <- list(`CAR`=data.frame(DIC=numeric(n.sim), WAIC=numeric(n.sim), LS=numeric(n.sim)),
               `HomCAR`=data.frame(DIC=numeric(n.sim), WAIC=numeric(n.sim), LS=numeric(n.sim)))


## Fit models with INLA ##
for(i in 1:n.sim){
  print(sprintf("Simulation %d", i))
  
  data.INLA <- data.frame(O=simu.data$AR$Obs[,i], E=E,
                          ID.area=1:nrow(carto), ID.area2=1:nrow(carto))
  
  ## Heterocedastic BYM model ##
  f.CAR <- O ~ f(ID.area, model="bym", graph=Q, constr=TRUE,
                 hyper=list(prec.spatial=list(prior=sdunif),
                            prec.unstruct=list(prior=sdunif)))
  
  M1 <- inla(f.CAR, family="poisson", data=data.INLA, E=E,
             control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=FALSE),
             inla.mode="compact", num.threads="8:1")
  
  ## Homoscedastic BYM model ##
  f.HomCAR <- O ~ f(ID.area, model="generic0", Cmatrix=Q.weighted, constr=TRUE,
                    hyper=list(prec=list(prior=sdunif))) + 
    f(ID.area2, model="iid", constr=TRUE, hyper=list(prec=list(prior=sdunif)))
  
  M2 <- inla(f.HomCAR, family="poisson", data=data.INLA, E=E,
             control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=FALSE),
             inla.mode="compact", num.threads="8:1")
  
  ## Save results ##
  post.mean.AR$CAR[,i] <- M1$summary.linear.predictor$mean
  post.mean.AR$HomCAR[,i] <- M2$summary.linear.predictor$mean

  DIC.AR$CAR[i,] <- c(M1$dic$dic, M1$waic$waic, -sum(log(M1$cpo$cpo)))
  DIC.AR$HomCAR[i,] <- c(M2$dic$dic, M2$waic$waic, -sum(log(M2$cpo$cpo)))
}


## Castile-La Mancha
#####################
carto <- CartoCM
Q <- prec.matrices$CM$Q
Q.weighted <- prec.matrices$CM$Q.weighted

## Objects to save ##
post.mean.CM <- list(`CAR`=matrix(nrow=nrow(carto), ncol=n.sim),
                     `HomCAR`=matrix(nrow=nrow(carto), ncol=n.sim))

DIC.CM <- list(`CAR`=data.frame(DIC=numeric(n.sim), WAIC=numeric(n.sim), LS=numeric(n.sim)),
               `HomCAR`=data.frame(DIC=numeric(n.sim), WAIC=numeric(n.sim), LS=numeric(n.sim)))


## Fit models with INLA ##
for(i in 1:n.sim){
  print(sprintf("Simulation %d", i))
  
  data.INLA <- data.frame(O=simu.data$CM$Obs[,i], E=E,
                          ID.area=1:nrow(carto), ID.area2=1:nrow(carto))
  
  ## Heterocedastic BYM model ##
  f.CAR <- O ~ f(ID.area, model="bym", graph=Q, constr=TRUE,
                 hyper=list(prec.spatial=list(prior=sdunif),
                            prec.unstruct=list(prior=sdunif)))
  
  M1 <- inla(f.CAR, family="poisson", data=data.INLA, E=E,
             control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=FALSE),
             inla.mode="compact", num.threads="8:1")
  
  ## Homoscedastic BYM model ##
  f.HomCAR <- O ~ f(ID.area, model="generic0", Cmatrix=Q.weighted, constr=TRUE,
                    hyper=list(prec=list(prior=sdunif))) + 
    f(ID.area2, model="iid", constr=TRUE, hyper=list(prec=list(prior=sdunif)))
  
  M2 <- inla(f.HomCAR, family="poisson", data=data.INLA, E=E,
             control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, return.marginals.predictor=FALSE),
             inla.mode="compact", num.threads="8:1")
  
  ## Save results ##
  post.mean.CM$CAR[,i] <- M1$summary.linear.predictor$mean
  post.mean.CM$HomCAR[,i] <- M2$summary.linear.predictor$mean

  DIC.CM$CAR[i,] <- c(M1$dic$dic, M1$waic$waic, -sum(log(M1$cpo$cpo)))
  DIC.CM$HomCAR[i,] <- c(M2$dic$dic, M2$waic$waic, -sum(log(M2$cpo$cpo)))
}


##############################################################################
## FIGURE 3:                                                                ##
## Empirical variances of posterior mean estimates for the log-risks across ##
## 100 simulated data sets. Left, variances for the traditional BYM model,  ##
## and right, variances for that including the HomCAR distribution.         ##
## Darker municipalities stand for higher variances.                        ##
##############################################################################

## Valencian Region ##
carto <- CartoVR
aux <- lapply(post.mean.VR, function(x) apply(sweep(x,2,apply(x,2,mean)),1,var))

carto$Var.CAR <- aux$CAR
carto$Var.HomCAR <- aux$HomCAR

breaks <- quantile(carto$Var.CAR,(0:7)/7)
labels <- sprintf("%.3f - %.3f", head(breaks,-1), tail(breaks,-1))

Fig3.VR <- tm_shape(carto) +
  tm_polygons(fill=c("Var.CAR","Var.HomCAR"),
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend(" ", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3, fill.free=FALSE) + 
  tm_layout(panel.labels=c("BYM model","HomCAR model"),
            panel.label.frame=FALSE,
            panel.label.fontface="bold",
            panel.label.size=1.2)


## Castile and Leon ##
carto <- CartoCL
aux <- lapply(post.mean.CL, function(x) apply(sweep(x,2,apply(x,2,mean)),1,var))

carto$Var.CAR <- aux$CAR
carto$Var.HomCAR <- aux$HomCAR

breaks <- quantile(carto$Var.CAR,(0:7)/7)
labels <- sprintf("%.3f - %.3f", head(breaks,-1), tail(breaks,-1))

Fig3.CL <- tm_shape(carto) +
  tm_polygons(fill=c("Var.CAR","Var.HomCAR"),
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend(" ", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3, fill.free=FALSE) + 
  tm_layout(panel.labels=c("BYM model","HomCAR model"),
            panel.label.frame=FALSE,
            panel.label.fontface="bold",
            panel.label.size=1.2)


## Aragon ##
carto <- CartoAR
aux <- lapply(post.mean.AR, function(x) apply(sweep(x,2,apply(x,2,mean)),1,var))

carto$Var.CAR <- aux$CAR
carto$Var.HomCAR <- aux$HomCAR

breaks <- quantile(carto$Var.CAR,(0:7)/7)
labels <- sprintf("%.3f - %.3f", head(breaks,-1), tail(breaks,-1))

Fig3.AR <- tm_shape(carto) +
  tm_polygons(fill=c("Var.CAR","Var.HomCAR"),
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend(" ", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3, fill.free=FALSE) + 
  tm_layout(panel.labels=c("BYM model","HomCAR model"),
            panel.label.frame=FALSE,
            panel.label.fontface="bold",
            panel.label.size=1.2)


## Castile-La Mancha ##
carto <- CartoCM
aux <- lapply(post.mean.CM, function(x) apply(sweep(x,2,apply(x,2,mean)),1,var))

carto$Var.CAR <- aux$CAR
carto$Var.HomCAR <- aux$HomCAR

breaks <- quantile(carto$Var.CAR,(0:7)/7)
labels <- sprintf("%.3f - %.3f", head(breaks,-1), tail(breaks,-1))

Fig3.CM <- tm_shape(carto) +
  tm_polygons(fill=c("Var.CAR","Var.HomCAR"),
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend(" ", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3, fill.free=FALSE) + 
  tm_layout(panel.labels=c("BYM model","HomCAR model"),
            panel.label.frame=FALSE,
            panel.label.fontface="bold",
            panel.label.size=1.2)

## Save the maps ##
Fig3 <- tmap_arrange(Fig3.VR, Fig3.CL, Fig3.AR, Fig3.CM, nrow=4, ncol=1)
tmap_save(Fig3, filename="Figure3.pdf", width=10, height=12)


######################################################################################################################
## TABLE 1: Standard deviations for the municipal empirical variances of posterior mean estimates for the log-risks ##
######################################################################################################################
Table1 <- data.frame('Valencian Region'=unlist(lapply(post.mean.VR, function(x) sd(apply(x,1,var)))),
                     'Castile Leon'=unlist(lapply(post.mean.CL, function(x) sd(apply(x,1,var)))),
                     'Aragon'=unlist(lapply(post.mean.AR, function(x) sd(apply(x,1,var)))),
                     'Castile-La Mancha'=unlist(lapply(post.mean.CM, function(x) sd(apply(x,1,var)))))

round(Table1,4)


#########################################################################
## FIGURE 5:                                                           ##
## Differences in variance of the relative risks between HomCAR-based  ##
## and ICAR-based BYM models. Brown areas correspond to locations with ##
## higher variance for the HomCAR model while green areas show lower   ##
## variance for this model.                                            ##
#########################################################################
rm(list=ls())

## Load cartography files ##
load("./Carto_files.Rdata")

## Load results from the traditional ICAR prior ##
load("./Results_ICARmodels.Rdata")
CartoVR["Var.ICAR"] <- apply(sweep(Res.VR,2, apply(Res.VR,2,mean)),1,var)
CartoCL["Var.ICAR"] <- apply(sweep(Res.CL,2, apply(Res.CL,2,mean)),1,var)
CartoAR["Var.ICAR"] <- apply(sweep(Res.AR,2, apply(Res.AR,2,mean)),1,var)
CartoCM["Var.ICAR"] <- apply(sweep(Res.CM,2, apply(Res.CM,2,mean)),1,var)

## Load results from the proposed HomCAR prior ##
load("./Results_HomCARmodels.Rdata")
CartoVR["Var.HomCAR"] <- apply(sweep(Res.VR.HomCAR,2, apply(Res.VR.HomCAR,2,mean)),1,var)
CartoCL["Var.HomCAR"] <- apply(sweep(Res.CL.HomCAR,2, apply(Res.CL.HomCAR,2,mean)),1,var)
CartoAR["Var.HomCAR"] <- apply(sweep(Res.AR.HomCAR,2, apply(Res.AR.HomCAR,2,mean)),1,var)
CartoCM["Var.HomCAR"] <- apply(sweep(Res.CM.HomCAR,2, apply(Res.CM.HomCAR,2,mean)),1,var)

## Compute relative differences of the municipal variances ##
CartoVR["dif.var"] <- (CartoVR$Var.HomCAR-CartoVR$Var.ICAR)/CartoVR$Var.ICAR
CartoCL["dif.var"] <- (CartoCL$Var.HomCAR-CartoCL$Var.ICAR)/CartoCL$Var.ICAR
CartoAR["dif.var"] <- (CartoAR$Var.HomCAR-CartoAR$Var.ICAR)/CartoAR$Var.ICAR
CartoCM["dif.var"] <- (CartoCM$Var.HomCAR-CartoCM$Var.ICAR)/CartoCM$Var.ICAR

Fig4.VR <- tm_shape(CartoVR) +
  tm_polygons(fill="dif.var",
              fill.scale=tm_scale(values=brewer.pal(7,"BrBG")[7:1],
                                  midpoint=0,
                                  breaks=c(-Inf,-0.5,-0.25,-0.1,0.1,0.25,0.5,Inf),
                                  labels=c("less than -50%","-50% to 10%","-25% to -10%","-10% to 10%","10% to 25%","25% to 50%","50% or more")),
              fill.legend=tm_legend("Changes in variance", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("Valencian Region", fontface="bold") + 
  tm_options(component.autoscale = FALSE)

Fig4.CL <- tm_shape(CartoCL) +
  tm_polygons(fill="dif.var",
              fill.scale=tm_scale(values=brewer.pal(7,"BrBG")[7:1],
                                  midpoint=0,
                                  breaks=c(-Inf,-0.5,-0.25,-0.1,0.1,0.25,0.5,Inf),
                                  labels=c("less than -50%","-50% to 10%","-25% to -10%","-10% to 10%","10% to 25%","25% to 50%","50% or more")),
              fill.legend=tm_legend("Changes in variance", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("Castile and Leon", fontface="bold") + 
  tm_options(component.autoscale = FALSE)

Fig4.AR <- tm_shape(CartoAR) +
  tm_polygons(fill="dif.var",
              fill.scale=tm_scale(values=brewer.pal(7,"BrBG")[7:1],
                                  midpoint=0,
                                  breaks=c(-Inf,-0.5,-0.25,-0.1,0.1,0.25,0.5,Inf),
                                  labels=c("less than -50%","-50% to 10%","-25% to -10%","-10% to 10%","10% to 25%","25% to 50%","50% or more")),
              fill.legend=tm_legend("Changes in variance", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("Aragon", fontface="bold") + 
  tm_options(component.autoscale = FALSE)

Fig4.CM <- tm_shape(CartoCM) +
  tm_polygons(fill="dif.var",
              fill.scale=tm_scale(values=brewer.pal(7,"BrBG")[7:1],
                                  midpoint=0,
                                  breaks=c(-Inf,-0.5,-0.25,-0.1,0.1,0.25,0.5,Inf),
                                  labels=c("less than -50%","-50% to 10%","-25% to -10%","-10% to 10%","10% to 25%","25% to 50%","50% or more")),
              fill.legend=tm_legend("Changes in variance", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("Castile-La Mancha", fontface="bold") + 
  tm_options(component.autoscale = FALSE)

Fig4 <- tmap_arrange(Fig4.VR, Fig4.CL, Fig4.AR, Fig4.CM, nrow=2, ncol=2)
tmap_save(Fig4, filename="Figure4.pdf", width=12, height=10)
