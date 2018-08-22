###############################################################################################
## Diagnostics graphics pour les modèles linéaires mixtes de type "mfl"                      ##
###############################################################################################                                                                                          ##
## Input:                                                                                    ##
##   mfl, un modèle linéaire mixte généré par le module lmixedm de JF Martin (fonction       ##
##        lmRepeated2FF), i.e. : un modèle linéaire mixte de type lmer créé par la formule   ##
##        mfl <- lmer( vd ~ time + fixfact + time:fixfact + (1| subject), ids)               ##
##        => noms de colonnes importants                                                     ##
##        => 1 seul effet aléatoire sur l'intercept                                          ##
###############################################################################################
## Output :                                                                                  ##
##   Les graphics, les calculs associés et les notations* utilisées dans le script suivent   ##
##   l'article de Singer et al (2016) Graphical Tools for detedcting departures from linear  ##
##   mixed model assumptions and some remedial measures, International Statistical Review    ##
##   (doi:10.1111/insr.12178)                                                                ##
##   * ajout d'une ou 2 lettres de typage de la variables                                    ##
##
##   Script adapté de http://www.ime.unicamp.br/~cnaber/residdiag_nlme_v22.R pour fonction-  ##
##   ner avec un modèle lmer (et non lme), des sujets avec des identifiants non numériques,  ##
##   et des observations non ordonnées sujet par sujet (dernier point à vérifier.)           ##
## ############################################################################################
## Remarques sur les calculs numériques                                                      ##
##    - l'inverse d'une matrice est calculée à partir de la fonction ginv du package MASS    ##
##       (Moore- Penrose generalized inverse) au lieu de la fonction "solve"                 ##
##    - la racine carrée des matrices sont calculées grâce à une déomposition SVD (car       ##
##      s'applique normalement ici uniquement à des matrices symétriques positives). A       ##
##      remplacer par la fonction sqrtm{expm} si des erreurs apparaissent ???                ##
###############################################################################################


library(ggplot2)
library(gridExtra)
library(grid)


##-------------------------------------------------------------------------------------------------##
## Helpers
##-------------------------------------------------------------------------------------------------##

## square root of a matrix
## http://www.cs.toronto.edu/~jepson/csc420/notes/introSVD.pdf (page 6)
## (for matMN a n x n matrix that symetric and non-negative definite)

sqrtmF<-function(matMN){
   matMN <- as.matrix(matMN)
   # ## check that matMN is symetric
   # if(!all(t(matMN==matMN)))
   #   stop("matMN must be symetric.")
   svd_dec <- svd(matMN)
   invisible(svd_dec$u%*%sqrt(diag(svd_dec$d))%*%t(svd_dec$v))
}


## qqplotF
## adapted from https://gist.github.com/rentrop/d39a8406ad8af2a1066c


qqplotF <- function(x, 
                    distribution = "norm", ..., 
                    line.estimate = NULL, 
                    conf = 0.95,
                    labels = names(x)){
   q.function <- eval(parse(text = paste0("q", distribution)))
   d.function <- eval(parse(text = paste0("d", distribution)))
   x <- na.omit(x)
   ord <- order(x)
   n <- length(x)
   P <- ppoints(length(x))
   daf <- data.frame(ord.x = x[ord], z = q.function(P, ...))
   
   if(is.null(line.estimate)){
      Q.x <- quantile(daf$ord.x, c(0.25, 0.75))
      Q.z <- q.function(c(0.25, 0.75), ...)
      b <- diff(Q.x)/diff(Q.z)
      coef <- c(Q.x[1] - b * Q.z[1], b)
   } else {
      coef <- coef(line.estimate(ord.x ~ z))
   }
   
   zz <- qnorm(1 - (1 - conf)/2)
   SE <- (coef[2]/d.function(daf$z,...)) * sqrt(P * (1 - P)/n)
   fit.value <- coef[1] + coef[2] * daf$z
   daf$upper <- fit.value + zz * SE
   daf$lower <- fit.value - zz * SE
   
   if(!is.null(labels)){ 
      daf$label <- ifelse(daf$ord.x > daf$upper | daf$ord.x < daf$lower, labels[ord],"")
   }
   
   p <- ggplot(daf, aes(x=z, y=ord.x)) +
      geom_point() + 
      geom_abline(intercept = coef[1], slope = coef[2], col = "red") +
      geom_line(aes(x=z, y = lower),daf,  col = "red", linetype = "dashed") +
      geom_line(aes(x=z, y = upper),daf,  col = "red", linetype = "dashed") +
      #geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2)+
      xlab("")+ylab("")
   if(!is.null(labels)) p <- p + geom_text( aes(label = label))
   
   return(p)
   #print(p)
   #coef
}


## histogramm
histF <- function(x, sd_x = NULL, breaks = "scott"){
   
   if(is.null(sd_x)){ sd_x <- sd(x)}
   
   ## Bandwith estimation (default is Scott)
   if(!breaks %in% c("sqrt", "sturges", "rice", "scott", "fd"))
      breaks <- "scott"
   
   if(breaks %in% c("sqrt", "sturges", "rice")){
      k <- switch(breaks,
                  sqrt = sqrt(length(x)),
                  sturges = floor(log2(x))+1,
                  rice = floor(2*length(x)^(1/3))
      ) 
      bw <- diff(range(x))/k
   }else{
      bw <- switch(breaks,
                   scott = 3.5*sd_x/length(x)^(1/3),
                   fd = diff(range(x))/(2*IQR(x)/length(x)^(1/3))
      )
   }
   
   
   daf <- data.frame(x=x)
   ## graph
   return(ggplot(data=daf, aes(x)) + 
             geom_histogram(aes(y = ..density..),
                            col="black", fill = "grey", binwidth = bw)+
             geom_density(size = 1.2,
                          col = "blue",
                          linetype = "blank",
                          fill = rgb(0,0,1, 0.1))+
             stat_function(fun = dnorm,
                           args = list(mean = 0, sd = sd_x),
                           col = "blue", size = 1.2)+
             theme(legend.position="none")+
             xlab(""))
}



##-------------------------------------------------------------------------------------------------##
## Main function
##-------------------------------------------------------------------------------------------------##


diagmflF <- function(mfl, title = "",
                     outlier.limit = 3, 
                     pvalCutof = 0.05,
                     least.confounded = FALSE){
   
   ## initialisation -------------------------------------------------------------------------------------
   responseC<- "vd"
   unitC <- "subject"
   timeC<- "time" 
   fixfactC<-"fixfact"
   hlimitN <- outlier.limit
   
   
   
   ## extracting information from mfl models -------------------------------------------------------------
   df <- mfl@frame
   yVn <- df[, responseC]
   nobsN <- length(yVn)
   idunitVc <- unique(df[, unitC])
   nunitN <- length(unique(idunitVc))
   xMN <- mfl@pp$X
   pN <- ncol(xMN)    
   zMN <- as.matrix(t(mfl@pp$Zt))
   gMN <- as.matrix(as.data.frame(VarCorr(mfl))[1, "vcov"]) ## Valable seulement pour 1 seul effet aléatoire
   gammaMN <- as.matrix(kronecker(diag(nunitN), gMN))  
   sigsqN <- as.data.frame(VarCorr(mfl))[2, "vcov"]
   rMN <- sigsqN*diag(nobsN)
   vMN <- (zMN%*%gammaMN%*%t(zMN)) + rMN
   invvMN <- MASS::ginv(vMN)
   hMN <- MASS::ginv(t(xMN)%*%invvMN%*%xMN)
   qMN <- invvMN - invvMN%*%xMN%*%(hMN)%*%t(xMN)%*%invvMN
   eblueVn<-mfl@beta 
   eblupVn<-gammaMN%*%t(zMN)%*%invvMN%*%(yVn-xMN%*%eblueVn) ## equivalent de ranef(mfl)
   rownames(eblupVn) <- colnames(zMN)
   
   
   ##  Calculs of matrices and vectors used in graph diagnosics ---------------------------------------------
   
   ## Marginal and individual predictions, residuals and variances
   marpredVn <- xMN%*%eblueVn
   marresVn <- yVn - marpredVn
   marvarMN <- vMN - xMN%*%hMN%*%t(xMN)
   
   condpredVn <- marpredVn + zMN%*%eblupVn
   condresVn <- yVn - condpredVn
   condvarMN <- rMN%*%qMN%*%rMN
   
   
   ## Analysis of marginal and conditional residuals
   stmarresVn <-stcondresVn <- rep(0,nobsN)
   lesverVn <- rep(0, nunitN)
   names(lesverVn )<- idunitVc
   
   for(i in 1:nunitN){
      
      idxiVn <- which(df[, unitC] == idunitVc[i]) ## position des observations du sujet i
      miN <- length(idxiVn)
      
      ## standardization of marginal residual
      stmarresVn[idxiVn] <- as.vector(solve(sqrtmF(marvarMN[idxiVn,idxiVn]))%*%marresVn[idxiVn])
      
      ##Standardized Lessafre and Verbeke's measure
      auxMN <- diag(1, ncol = miN, nrow =miN)- stmarresVn[idxiVn]%*%t(stmarresVn[idxiVn])
      lesverVn[i] <- sum(diag(auxMN%*%t(auxMN)))
      
      ## standardization of conditional residual
      stcondresVn[idxiVn] <- as.vector(solve(sqrtmF(condvarMN[idxiVn,idxiVn]))%*%condresVn[idxiVn])
   }
   lesverVn <- lesverVn/sum(lesverVn)
   
   
   ##  Least confounded residuals (à valider !)
   ## résultats différents des valeurs estimées via le script de Singer, car 
   ## non unicité des vecteurs propres de la décomposition spectrale ?
   if(least.confounded){
      sqrMN <- sqrtmF(rMN)
      specDec <- eigen((sqrMN%*%qMN%*%sqrMN), symmetric = TRUE, only.values = FALSE)
      
      cMN <- t(sqrt(solve(diag((specDec$values[1:(nobsN -pN)])))) %*% t(specDec$vectors[1:(nobsN -pN),1:(nobsN -pN)]) %*%
                  solve(sqrtmF(rMN[1:(nobsN -pN),1:(nobsN -pN)])) )
      
      lccondresVn <- (cMN%*%condresVn[1:(nobsN -pN)])/sqrt(diag(cMN%*%condvarMN[1:(nobsN -pN),1:(nobsN -pN)]%*%t(cMN)))
   }
   
   
   ##  EBLUP analysis (Mahalanobis' distance)
   varbMN <- gammaMN%*%t(zMN)%*%qMN%*%zMN%*%gammaMN
   mdistVn <- rep(0, nunitN)
   for(i in 1:nunitN){
      mdistVn[i] <- eblupVn[i]^2/varbMN[i, i]
   }
   mdistVn <-  mdistVn/sum(mdistVn)
   
   
   ## Combine data and results in 2 data frames for easy plotting with ggplot2 -----------------------------
   
   ## long data frame (all observations)

   df <- data.frame(df,
                    mar.pred = marpredVn,
                    mar.res = marresVn,
                    st.mar.res = stmarresVn,
                    cond.pred = condpredVn,
                    cond.res = condresVn,
                    st.cond.res = stcondresVn)
   
   if(!("fixfact" %in% colnames(df)))
      df$fixfact <- rep(1, nrow(df))
   df$numTime <- as.numeric(levels(as.factor(df$time)))
   
   
   ## short data frame (1 row per unit)
   
   unitDf <- data.frame(unit = idunitVc,
                        eblup = eblupVn,
                        lvm = lesverVn,
                        mal = mdistVn)
   
   unitDf$fixfact <- sapply(1:nrow(unitDf),
                            function(i){
                               unique(df[which(df[, unitC] == unitDf$unit[i]),
                                         fixfactC])
                            })
   
   unitDf$se <- rep(NA, nrow(unitDf))
   re <- ranef(mfl, condVar =TRUE)
   for(i in 1:nrow(unitDf))
      unitDf$se[i] <- sqrt(attr(re[[unitC]], "postVar")[1,1,i])
   unitDf$upper <- unitDf$eblup+unitDf$se*1.96
   unitDf$lower <- unitDf$eblup-unitDf$se*1.96
   
   
   ## Outliers "annotations"
   df$marres.out <- rep("", nrow(df))
   df$marres.out[abs(df$st.mar.res)>hlimitN] <- 
      paste(df[abs(df$st.mar.res)>hlimitN, unitC], 
            df[abs(df$st.mar.res)>hlimitN, timeC],
            sep = ".")
   df$marres.out <- paste(" ", df$marres.out, sep ="")
   
   df$condres.out <- rep("", nrow(df))
   df$condres.out[abs(df$st.cond.res)>hlimitN] <- 
      paste(df[abs(df$st.cond.res)>hlimitN, unitC], 
            df[abs(df$st.cond.res)>hlimitN, timeC],
            sep = ".")
   df$condres.out <- paste(" ", df$condres.out, sep ="")
   
   ## Diagnostic Plots -------------------------------------------------------------------------------------
  
   ## Linearity of effect and outlying observations
   p1 <- ggplot(data = df, 
                aes(x=mar.pred, 
                    y=st.mar.res, 
                    colour=fixfact)) +
      geom_point(size =2) + 
      geom_hline(yintercept = 0, col = "grey")+
      geom_smooth(aes(x=mar.pred, y=st.mar.res), data = df,  se = FALSE, col = "blue", method = "loess")+
      ggtitle("Linearity of effects/outlying obervations")+
      xlab("Marginal predictions")+
      ylab("Standardized marginal residuals")+
      theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))+
      geom_hline(yintercept = c(-1,1)*hlimitN, linetype = "dashed")+
      geom_text(aes(label = marres.out, col = fixfact), hjust=0, vjust=0)
   
   # p1hist <- histF(df$st.mar.res, sd_x =1)+
   #   xlab("Standardized marginal residuals")
   
   
   ## Presence of outlying observations and homoscedacity of residuals
   p2 <- ggplot(data = df, 
                aes(x=cond.pred, 
                    y=st.cond.res, 
                    colour=fixfact)) +
      geom_point(size =2) + 
      geom_hline(yintercept = 0, col = "grey")+
      geom_smooth(aes(x=cond.pred,  y=st.cond.res), data = df,  se = FALSE, col = "blue", method = "loess")+
      ggtitle("Homosedasticity of conditional residuals/outlying observations")+
      xlab("Individual predictions")+
      ylab("Standardized conditional residuals")+
      theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))+
      geom_hline(yintercept = c(-1,1)*hlimitN, linetype = "dashed")+
      geom_text(aes(label = condres.out, col = fixfact), hjust=0, vjust=0)
   
   # p2hist <- histF(df$st.cond.res, sd_x =1)+
   #   xlab("Standardized conditional residuals")
   
   
   ## Normality of residuals
   if(least.confounded){
      p3 <- qqplotF(x = lccondresVn, 
                    distribution = "norm", 
                    line.estimate = NULL,
                    conf = 0.95)+
         xlab("Standard normal quantiles")+
         ylab("Least confounded conditional residual quantiles")+
         ggtitle("Normality of conditional error (least confounded)")+
         theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))
      
      # p3hist <- histF(lccondresVn, sd_x =1)+
      #     xlab("Least confounded conditional residuals")
   }else{
      p3 <-qqplotF(x = df$st.cond.res, 
                   distribution = "norm", 
                   line.estimate = NULL,
                   conf = 0.95)+
         xlab("Standard normal quantiles")+
         ylab("Standardized conditional residual quantiles")+
         ggtitle("Normality of conditional error")+
         theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))
      # p3hist <-p2hist
   }

   ## Within-units covariance structure
   p4 <- ggplot(data = unitDf, 
                aes(x=unit, 
                    y=lvm, 
                    colour=fixfact)) +
      geom_point(size =2) +
      theme(legend.position="none")+
      xlab(unitC)+
      ylab("Standardized Lesaffre-Verbeke measure")+
      geom_hline(yintercept = 2*mean(unitDf$lvm), linetype = "dashed")+
      geom_text(aes(label = unit, col = fixfact), 
                data = unitDf[unitDf$lvm>2*mean(unitDf$lvm), ],
                hjust=0, vjust=0)+
      ggtitle("Within-units covariance matrice")+
      theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))
   
   ## EBLUP modif lmerTest v3 mais plante si pas d'effet aléatoire
   #pvl <- ranova(model=mfl, reduce.terms = TRUE)
   #pvrnd <- pvl[[6]][2]
   #ggtitle(paste("Random effect on intercept (LRT p-value = ",round(pvrnd,digits = 5), ")", sep = ""))+
   
   p5 <-
      ggplot(aes(x = eblup, y = unit, col = fixfact), data = unitDf)+
      geom_point(size =3)+
      geom_segment(aes(xend = lower, yend = unit)) + 
      geom_segment(aes(xend = upper, yend = unit))+
      ggtitle("Random effect on intercept")+
      theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))+
      geom_vline(xintercept = 0, linetype = "dashed")+
      ylab(unitC)
   
   ## p6
   p6 <-qqplotF(x = unitDf$mal, 
                distribution = "chisq", 
                df= 1,
                line.estimate = NULL,
                conf = 0.95)+
      xlab("Chi-squared quantiles")+
      ylab("Standadized Mahalanobis distance")+
      ggtitle("Normality of random effect")+
      theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))

   ## Outlying subjects
   p7 <-
      ggplot(aes(y = mal, x= unit, col = fixfact), data = unitDf)+
      geom_point(size =3)+
      ylab("Standardized Mahalanobis distance")+
      geom_vline(xintercept = 0, linetype = "dashed")+
      theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))+
      geom_hline(yintercept = 2*mean(unitDf$mal), linetype = "dashed")+
      geom_text(aes(label = unit, col = fixfact), 
                data = unitDf[unitDf$mal>2*mean(unitDf$mal), ],
                hjust=1, vjust=0)+
      ggtitle("Outlying subjects")+
      xlab(unitC)

   ## "Data" and "modeling" Plots --------------------------------------------------------------------------
   
   ## Individual time-course
   rawPlot <- ggplot(data = df, 
                     aes(x=time, y=vd, colour=fixfact, group=subject)) +
      geom_line() +  ggtitle("Individual time-courses ")+
      theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))

   ## Post-hoc estimates (modification due to lmerTest v3)
   ddlsm1  <- data.frame(difflsmeans(mfl,test.effs=NULL))
   colnames(ddlsm1)[9] <- "pvalue"
   # ddlsm1$name <- sapply(rownames(ddlsm1),
   #                       function(nam){
   #                          strsplit(nam, split = " ", fixed =TRUE)[[1]][1]
   #                       }) 
   # ddlsm1$detail <- sapply(rownames(ddlsm1),
   #                         function(nam){
   #                            paste(strsplit(nam, split = " ", fixed =TRUE)[[1]][-1],
   #                                  collapse= "")
   #                         })
   # 
   # colnames(ddlsm1)<- make.names(colnames(ddlsm1))
   ddlsm1$Significance <- rep("NS", nrow(ddlsm1))
   
   ## modif JF pour tenir compte du seuil de pvalues defini par le user  
   ddlsm1$Significance[which(ddlsm1$pvalue <pvalCutof & ddlsm1$pvalue >=0.01)] <- "p-value < threshold"
   ddlsm1$Significance[which(ddlsm1$pvalue <0.01 & ddlsm1$pvalue >=0.005)] <- "p-value < 0.01"
   ddlsm1$Significance[which(ddlsm1$pvalue <0.005)] <- "p-value < 0.005"

   phPlot <- ggplot(ddlsm1, aes(x = levels, y = Estimate))+
      facet_grid(facets = ~term, ddlsm1,scales = "free", space = "free")+
      geom_bar( aes(fill = Significance), stat="identity")+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      scale_fill_manual(
         values = c("NS" = "grey",
                    "p-value < threshold"  = "yellow",
                    "p-value < 0.01"  = "orange",
                    "p-value < 0.005" = "red"))+
      geom_errorbar(aes(ymin = lower, ymax =upper ), width=0.25)+
      ggtitle("Post-hoc estimates")+xlab("")+
      theme(plot.title = element_text(size = rel(1.2), face = "bold"))

   ## Final plotting
   
   blank<-rectGrob(gp=gpar(col="white"))
   rectspacer<-rectGrob(height = unit(0.1, "npc"), gp=gpar(col="grey")) 
   
   grid.arrange(blank,
                rawPlot,  blank, phPlot,
                rectspacer,
                p1,blank, p2, blank, p3, blank, p4,
                blank,
                p5,blank, p6, blank,  p7, blank,
                blank,blank,
                top = textGrob(title,gp=gpar(fontsize=40,font=4)),
                layout_matrix = matrix(c(rep(1,7),
                                         2, 3, rep(4,3), 20,21,
                                         rep(5,7),
                                         6:12, 
                                         rep(13,7),
                                         14:18, rep(19,2)),
                                       ncol=7, nrow=6, byrow=TRUE),
                heights= c(0.1/3, 0.3, 0.1/3, 0.3, 0.1/3, 0.3 ),
                widths = c(0.22, 0.04, 0.22,0.04 , 0.22, 0.04, 0.22))
   
   invisible(NULL)

}

#diagmflF(mfl, title = "",  outlier.limit = 3, least.confounded =TRUE)

plot.res.Lmixed <- function(mfl, df, title = "", pvalCutof = 0.05) {

   nameVar <- colnames(df)[4]
   fflab <- colnames(df)[1]
   ## Individual time-course
   rawPlot <- 
      ggplot(data = df, aes(x=df[[2]], y=df[[4]], colour=df[[1]], group=df[[3]])) +
      geom_line() +  ggtitle("Individual time-courses (raw data)")+ 
      ylab(nameVar) +
      xlab(label = colnames(df)[2])+
      theme(legend.title = element_blank() , legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))
   
   ## Boxplot of fixed factor
   bPlot <- 
      ggplot(data = df, aes(y=df[[4]], x=df[[1]], color=df[[1]]))+
      geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+
      ggtitle(paste("Boxplot by ",fflab,sep=""))+xlab("")+ylab("")+
      theme(legend.title = element_blank(), plot.title = element_text(size = rel(1.2), face = "bold"))
   
   ## Post-hoc estimates
   
   ddlsm1  <- mfl
   ddlsm1$name <- rownames(ddlsm1)
   # ddlsm1$name <- sapply(rownames(ddlsm1),
   #                       function(nam){
   #                          strsplit(nam, split = " ", fixed =TRUE)[[1]][1]
   #                       })
   # ddlsm1$detail <- sapply(rownames(ddlsm1),
   #                         function(nam){
   #                            paste(strsplit(nam, split = " ", fixed =TRUE)[[1]][-1],
   #                                  collapse= "")
   #                         })
   # 
   #colnames(ddlsm1)<- make.names(colnames(ddlsm1))
   ddlsm1$Significance <- rep("NS", nrow(ddlsm1))
   
   ## modif JF pour tenir compte du seuil de pvalues defini par le user 
   options("scipen"=100, "digits"=5)
   pvalCutof <- as.numeric(pvalCutof)
   bs=0.05; bm=0.01; bi=0.005
   if (pvalCutof >bm) {bs <- pvalCutof} else
      if (pvalCutof <bm & pvalCutof >bi) {bm <- pvalCutof; bs <- pvalCutof} else
         if (pvalCutof < bi) {bi <- pvalCutof; bm <- pvalCutof; bs <- pvalCutof}
   lbs <- paste("p-value < ",bs, sep="")
   lbm <- paste("p-value < ",bm, sep="")
   lbi <- paste("p-value < ",bi, sep="")
   
   cols <- paste("p-value < ",bs, sep="")
   colm <- paste("p-value < ",bm, sep="")
   coli <- paste("p-value < ",bi, sep="")
   valcol <- c("grey","yellow","orange","red")
   names (valcol) <- c("NS",lbs,lbm,lbi)
   ddlsm1$Significance[which(ddlsm1$p.value<= bs)] <- lbs    
   ddlsm1$Significance[which(ddlsm1$p.value<bs & ddlsm1$p.value>= bm)] <- lbm
   ddlsm1$Significance[which(ddlsm1$p.value<bi)] <- lbi

   phPlot <- 
      ggplot(ddlsm1, aes(x = levels, y = Estimate))+
      facet_grid(facets = ~term, ddlsm1,scales = "free", space = "free")+
      geom_bar( aes(fill = Significance), stat="identity")+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
      scale_fill_manual(
         # values = c("NS" = "grey", "p-value < threshold" = "yellow","p-value < 0.01" = "orange","p-value < 0.005" = "red"))+
         # values = c("NS" = 'grey', "pvalue < 0.05 "= 'yellow',"p-value < 0.01" = 'orange',"p-value < 0.005" = 'red'))+
         # values = c("NS = grey", "p-value < threshold = yellow","p-value < 0.01 = orange","p-value < 0.005 = red"))+
         values = valcol )+ 
      geom_errorbar(aes(ymin = Lower.CI, ymax =Upper.CI ), width=0.25)+
      # ggtitle(paste("Post-hoc estimates with p-value threshold = ",pvalCutof,sep=""))+xlab("")+
      ggtitle("Post-hoc estimates ")+xlab("")+
      theme(plot.title = element_text(size = rel(1.2), face = "bold"))
   
   ## Final plotting
   grid.arrange(arrangeGrob(rawPlot,bPlot,ncol=2),
                phPlot,nrow=2,
                top = textGrob(title,gp=gpar(fontsize=32,font=4))
   )
   
} 