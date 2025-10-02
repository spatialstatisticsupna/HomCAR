################################################################################################
## R code for the reproduce the analyses of the paper entitled "A proposal for homoskedastic  ##
## modelling with conditional auto-regressive distributions" (Martínez-Beneito et al., 2025). ##
################################################################################################
rm(list=ls())
library(INLA)  
library(RColorBrewer)
library(sf)
library(spdep)
library(tmap)

# Option for the management of cartographies
sf_use_s2(FALSE)


#################################################################
## SECTION 2: THE ICAR DISTRIBUTION AND ITS MARGINAL VARIANCES ##
#################################################################
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


####################################################################################
## FIGURE 1:                                                                      ##
## Marginal variance for the municipalities of 4 autonomous regions of Spain with ##
## markedly different shapes. Variances correspond to an ICAR distribution over   ##
## the corresponding region. Darker municipalities stand for higher variances.    ##
####################################################################################

## Load cartography files and compute the marginal variances of their corresponding spatial precision matrices ##
load("./Carto_files.Rdata")

vars.marginal <- lapply(list(VR=CartoVR, CL=CartoCL, AR=CartoAR, CM=CartoCM), function(x){
  nb <- poly2nb(st_sf(x))
  
  W <- nb2mat(nb, style="B")
  dimnames(W) <- NULL
  attr(W,"call") <- NULL
  
  Q <- diag(colSums(W))-W
  vars.marginal <- inla.ginv.diag(Q)
  
  return(vars.marginal)
})


## Valencian Region ##
CartoVR["PriorVar"] <- vars.marginal$VR

breaks <- quantile(CartoVR$PriorVar,(0:7)/7)
labels <- sprintf("%.2f - %.2f", head(breaks,-1), tail(breaks,-1))

Fig1.VR <- tm_shape(CartoVR) +
  tm_polygons(fill="PriorVar",
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend("", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("ICAR Variances, Valencian Region", fontface="bold") + 
  tm_options(component.autoscale = FALSE)


## Castile and Leon ##
CartoCL["PriorVar"] <- vars.marginal$CL

breaks <- quantile(CartoCL$PriorVar,(0:7)/7)
labels <- sprintf("%.2f - %.2f", head(breaks,-1), tail(breaks,-1))

Fig1.CL <- tm_shape(CartoCL) +
  tm_polygons(fill="PriorVar",
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend("", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("ICAR Variances, Castile and Leon", fontface="bold") + 
  tm_options(component.autoscale = FALSE)


## Aragon ##
CartoAR["PriorVar"] <- vars.marginal$AR

breaks <- quantile(CartoAR$PriorVar,(0:7)/7)
labels <- sprintf("%.2f - %.2f", head(breaks,-1), tail(breaks,-1))

Fig1.AR <- tm_shape(CartoAR) +
  tm_polygons(fill="PriorVar",
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend("", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("ICAR Variances, Aragon", fontface="bold") + 
  tm_options(component.autoscale = FALSE)


## Castile-La Mancha ##
CartoCM["PriorVar"] <- vars.marginal$CM

breaks <- quantile(CartoCM$PriorVar,(0:7)/7)
labels <- sprintf("%.2f - %.2f", head(breaks,-1), tail(breaks,-1))

Fig1.CM <- tm_shape(CartoCM) +
  tm_polygons(fill="PriorVar",
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend("", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("ICAR Variances, Castile-La Mancha", fontface="bold") + 
  tm_options(component.autoscale = FALSE)


## Save the maps ##
Fig1 <- tmap_arrange(Fig1.VR, Fig1.CL, Fig1.AR, Fig1.CM, nrow=2, ncol=2)
tmap_save(Fig1, filename="Figure1.pdf", width=12, height=10)


###################################################################################
## FIGURE 2:                                                                     ##
## Empirical variances of the posterior mean estimates of the log-risks for each ##
## municipality across the 100 causes of death analyzed in each region.          ##
## Darker municipalities stand for higher variances.                             ##
###################################################################################

## Load the results /posterior mean estimates) from previously fitted BYM models ##
## for 100 causes of death at the municipality level                             ##
load("./Results_models.Rdata")


## Valencian Region ##
CartoVR$Var.Q <- apply(sweep(Res.VR,2,apply(Res.VR,2,mean)),1,var)

breaks <- quantile(CartoVR$Var.Q,(0:7)/7)
labels <- sprintf("%.3f - %.3f", head(breaks,-1), tail(breaks,-1))

Fig2.VR <- tm_shape(CartoVR) +
  tm_polygons(fill="Var.Q",
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend("", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("Variance of posterior estimates, Valencian Region", fontface="bold", size=1.2) + 
  tm_options(component.autoscale = FALSE)


## Castile and Leon ##
CartoCL$Var.Q <- apply(sweep(Res.CL,2,apply(Res.CL,2,mean)),1,var)

breaks <- quantile(CartoCL$Var.Q,(0:7)/7)
labels <- sprintf("%.3f - %.3f", head(breaks,-1), tail(breaks,-1))

Fig2.CL <- tm_shape(CartoCL) +
  tm_polygons(fill="Var.Q",
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend("", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("Variance of posterior estimates, Castile and Leon", fontface="bold", size=1.2) + 
  tm_options(component.autoscale = FALSE)


## Aragon ##
CartoAR$Var.Q <- apply(sweep(Res.AR,2,apply(Res.AR,2,mean)),1,var)

breaks <- quantile(CartoAR$Var.Q,(0:7)/7)
labels <- sprintf("%.3f - %.3f", head(breaks,-1), tail(breaks,-1))

Fig2.AR <- tm_shape(CartoAR) +
  tm_polygons(fill="Var.Q",
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend("", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("Variance of posterior estimates, Aragon", fontface="bold", size=1.2) + 
  tm_options(component.autoscale = FALSE)


## Castile-La Mancha ##
CartoCM$Var.Q <- apply(sweep(Res.CM,2,apply(Res.CM,2,mean)),1,var)

breaks <- quantile(CartoCM$Var.Q,(0:7)/7)
labels <- sprintf("%.3f - %.3f", head(breaks,-1), tail(breaks,-1))

Fig2.CM <- tm_shape(CartoCM) +
  tm_polygons(fill="Var.Q",
              fill.scale=tm_scale(values=brewer.pal(7,"YlOrBr"),
                                  breaks=breaks, labels=labels),
              fill.legend=tm_legend("", show=TRUE, reverse=TRUE,
                                    position=tm_pos_out("right","center"),
                                    frame=FALSE),
              col_alpha=0.3) + 
  tm_title("Variance of posterior estimates, Castile-La Mancha", fontface="bold", size=1.2) + 
  tm_options(component.autoscale = FALSE)


## Save the maps ##
Fig2 <- tmap_arrange(Fig2.VR, Fig2.CL, Fig2.AR, Fig2.CM, nrow=2, ncol=2)
tmap_save(Fig2, filename="Figure2.pdf", width=12, height=10)


################################################################################################
## Correlation between the prior variances (Fig 1) and the variances of the estimates (Fig 2) ##
################################################################################################
lapply(list('Valencian Region'=CartoVR,
            'Castile and Leon'=CartoCL,
            'Aragon'=CartoAR,
            'Castile-La Mancha'=CartoCM),
       function(x){
         round(cor(x$PriorVar, x$Var.Q),2)
})
