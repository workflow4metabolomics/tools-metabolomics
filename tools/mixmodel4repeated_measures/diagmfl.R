#' Calcul des grandeurs "diagnostiques"
#'
#'  Script adapte de http://www.ime.unicamp.br/~cnaber/residdiag_nlme_v22.R pour fonctionner
#'  avec un modele lmer (et non lme), des sujets avec des identifiants non numeriques,
#'  et des observations non ordonnees sujet par sujet (dernier point a verifier.)
#'
#'  @detail Les graphiques, les calculs associés et les notations utilisees dans le script suivent
#'   l'article de Singer et al (2016) Graphical Tools for detedcting departures from linear
#'    mixed model assumptions and some remedial measures, International Statistical Review
#'       (doi:10.1111/insr.12178)
#'
#' @param mfl A linear mixed model fitted via lmer or a data frame containing data
#' @return A list
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export lmer.computeDiag



lmer.computeDiag <- function(mfl){
  
  ## Check arguments ---------------------------------------------------------
  
  if(length(mfl@flist)>1)
    stop("Several 'grouping level' for random effect not implemented yet.")
  
  
  ## extracting information from mfl models -------------------------------------------------------------
  
  # data
  df <- mfl@frame
  responseC <- names(df)[1]
  unitC <- names(mfl@flist)[1]
  
  # observations
  yVn <- df[, responseC]
  nobsN <- length(yVn)
  
  # units
  idunitVc <- levels(mfl@flist[[1]])#unique(df[, unitC])
  nunitN <- length(unique(idunitVc))
  
  #X
  xMN <- mfl@pp$X
  pN <- ncol(xMN)
  
  #Z
  zMN <- t(as.matrix(mfl@pp$Zt))
  
  # Estimated covariance matrix of random effects (Gam)
  aux <- VarCorr(mfl)[[1]] ## assuming only one level of grouping
  aux2 <- attr(aux, "stddev")
  gMN <- attr(aux, "correlation")*(aux2%*%t(aux2))
  gammaMN <- as.matrix(kronecker(diag(nunitN), gMN))
  q <- dim(gMN)[1]
  
  # Estimated covariance matrix of conditonal error (homoskedastic conditional independance model)
  sigsqN <- attr(VarCorr(mfl), "sc")^2
  rMN <- sigsqN*diag(nobsN)
  
  # Estimated covariance matrix of Y
  vMN <- (zMN%*%gammaMN%*%t(zMN)) + rMN
  invvMN <- MASS::ginv(vMN)
  
  # H and Q matrix
  hMN <- MASS::ginv(t(xMN)%*%invvMN%*%xMN)
  qMN <- invvMN - invvMN%*%xMN%*%(hMN)%*%t(xMN)%*%invvMN
  
  # eblue and eblup
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
  
  
  ## Least confounded conditional residuals
  
  
  
  
  ##  EBLUP analysis (Mahalanobis' distance)
  varbMN <- gammaMN%*%t(zMN)%*%qMN%*%zMN%*%gammaMN
  mdistVn <- rep(0, nunitN)
  
  # Initial coding: works only for 1 single random effect
  # for(i in 1:nunitN){
  #   mdistVn[i] <- eblupVn[i]^2/varbMN[i, i]
  # }
  # mdistVn <-  mdistVn/sum(mdistVn)
  
  qm <- q-1
  for(j in 1:nunitN){
    gbi <- varbMN[(q*j-qm):(q*j), (q*j-qm):(q*j)]
    eblupi <- eblupVn[(q*j-qm):(q*j)]
    mdistVn[j] <- t(eblupi)%*%ginv(gbi)%*%eblupi
  }
  #pmdistVn <-  mdistVn/sum(mdistVn)
  names(mdistVn) <- levels(mfl@flist[[1]])
  
  
  ## output ----------------------------------------------
  
  return(list(
    data = df,
    q = q,
    eblue = eblueVn,
    eblup = eblupVn,
    marginal.prediction = marpredVn,
    conditional.prediction = condpredVn,
    std.marginal.residuals = stmarresVn ,
    std.conditional.residuals = stcondresVn ,
    mahalanobis.distance = mdistVn,
    std.mahalanobis.distance = mdistVn/sum(mdistVn),
    std.lesaffreverbeke.measure = lesverVn
  ))
}


#' Wrapper function for diagnostic plots of 'lmer' linear mixed models
#'
#' (W4M mixmod)
#'
#' @param mfl A linear mixed model fitted via lmer or a data frame containing data
#' @param  title aa
#' @param  outlier.limit aa
#' @param  pvalCutof aa
#' @param  resC aa
#' @param  uniC aa
#' @param  fixC aa
#' @param  lest.confounded Not used yet.
#' @return NULL
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export diagmflF


diagmflF <- function(mfl,
                     title = "",
                     outlier.limit = 3,
                     pvalCutof = 0.05,
                     resC = "vd",
                     uniC = "subject",
                     timC = "time",
                     fixC = "fixfact",
                     least.confounded = FALSE){
  
  
  ## diagnostics
  
  diagLs <- lmer.computeDiag(mfl)
  
  ## plots
  blank<-rectGrob(gp=gpar(col="white"))
  rectspacer<-rectGrob(height = unit(0.1, "npc"), gp=gpar(col="grey"))
  
  grid.arrange(blank,
               plot_timeCourse(mfl,
                               responseC = resC,
                               timeC = timC,
                               subjectC = uniC,
                               fixfactC =fixC,
                               offset_subject = FALSE,
                               plotL = FALSE,
                               colorType = "FIXFACT",
                               shapeType = "none",
                               lineType = "FIXFACT"),
               blank,
               plot_posthoc(mfl,
                            pvalCutof = pvalCutof,
                            plotL =FALSE,
                            titC = "Post-hoc estimates (uncorrected p-value)"),
               rectspacer,
               plot_linearity(diagLs, hlimitN =outlier.limit, plotL = FALSE,
                              label_factor = c(uniC,fixC, timC) ),
               blank,
               plot_conditionalResiduals(diagLs, hlimitN =outlier.limit, plotL = FALSE,
                                         label_factor = c(uniC,fixC, timC)),
               blank,
               plot_condresQQplot(diagLs,  plotL = FALSE),
               blank,
               plot_lesaffreVeerbeke(diagLs,  plotL = FALSE),
               blank,
               plot_randomEffect(mfl, plotL = FALSE)[[1]],
               blank,
               plot_mahalanobisKhi2(diagLs,  plotL = FALSE),
               blank,
               plot_mahalanobis(diagLs,  plotL = FALSE),
               blank,
               blank,
               blank,
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
  
  
}

#######################################################################################################
## Raw data time courses
#######################################################################################################

#' Visualization of raw time course
#'
#' Une
#'
#' @param mfl A linear mixed model fitted via lmer or a data frame containing data
#' @param  responseC Name of the 'response' variable
#' @param  timeC   Name of the 'time' variable
#' @param  subjectC  Name of the 'subject' variable
#' @param  fixfactC  Name of the 'fixed factor' variable (e.g.treatment)
#' @param  offset_subject Boolean indicating if an offset value (subject's mean) should substracted to each data point. Default is FALSE
#' @param  plotL Boolean
#' @param  colorType One of NA, FIXFACT or SUBJECT
#' @param  shapeType One of NA, FIXFACT or SUBJECT
#' @param  lineType One of NA, FIXFACT or SUBJECT
#' @return A plot
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export plot_timeCourse



plot_timeCourse <- function (mfl,
                             responseC,
                             timeC,
                             subjectC,
                             fixfactC =NULL,
                             offset_subject = FALSE,
                             plotL = TRUE,
                             colorType = NA, ## subject, fixfact, none or NA
                             shapeType = NA, ## subject, fixfact, none or NA
                             lineType = NA ## subject, fixfact, none or NA
){
  
  ## Data -----
  
  if (class(mfl) %in% c("merModLmerTest", "lmerMod", "lmerModLmerTest")){
    DF <- mfl@frame
  } else if (class(mfl) =="data.frame") {
    DF <- mfl
  } else {
    stop ("'mfl' argument must be a linear mixed effect model or a data frame.")
  }
  
  
  
  ## Format data -----
  
  if(is.null(fixfactC)){
    DF <- DF[, c(responseC,  timeC, subjectC)]
    colnames(DF) <- c("DV", "TIME", "SUBJECT")
    meanDF <- aggregate(DF$DV,
                        by = list(SUBJECT = DF$SUBJECT,
                                  TIME = DF$TIME),
                        FUN = mean,
                        na.rm = TRUE)
    colnames(meanDF) <- c("SUBJECT", "TIME", "DV")
    meanDF$GROUP <- meanDF$SUBJECT
    
  } else{
    
    DF <- DF[, c(responseC, fixfactC, timeC, subjectC)]
    colnames(DF) <- c("DV","FIXFACT", "TIME", "SUBJECT")
    
    meanDF <- aggregate(DF$DV,
                        by = list(SUBJECT = DF$SUBJECT,
                                  TIME = DF$TIME,
                                  FIXFACT = DF$FIXFACT),
                        FUN = mean,
                        na.rm = TRUE)
    colnames(meanDF) <- c("SUBJECT", "TIME","FIXFACT", "DV")
    meanDF$GROUP <- paste(meanDF$SUBJECT, meanDF$FIXFACT, sep = "_")
    
  }
  
  
  ## Offset -----
  
  
  if(offset_subject){
    
    offsetMN <- aggregate(DF$DV, by = list(DF$SUBJECT), mean, na.rm = TRUE )
    offsetVn <- offsetMN[, 2]
    names(offsetVn) <- offsetMN[,1]
    rm(offsetMN)
    DF$DV <- DF$DV - offsetVn[DF$SUBJECT]
    meanDF$DV <- meanDF$DV - offsetVn[as.character(meanDF$SUBJECT)]
  }
  
  ## Graphical parameters -----
  
  xlabC <-  timeC
  ylabC <- responseC
  titC <- "Individual time-courses"
  
  if(offset_subject){
    ylabC <- paste(ylabC, "minus 'within-subject' empirical mean")
    titC <- paste(titC, "('within-subject' empirical mean offset)")
  }
  
  
  ## color
  
  if(is.na(colorType)){ ## automaticatical attribution
    
    if(is.null(fixfactC)){
      colorType <- "SUBJECT"
    } else {
      colorType <- "FIXFACT"
    }
    colTxt <- paste(", colour=", colorType)
    
  } else if (colorType == "none"){
    colTxt  <- ""
  } else {
    colTxt <- paste(", colour=", colorType)
  }
  
  
  ## lineType
  if(is.na(lineType)){ ## automaticatical attribution
    
    if(is.null(fixfactC)){
      linTxt  <- ""
    } else {
      linTxt <- paste(", linetype=",
                      ifelse(colorType == "SUBJECT", "FIXFACT", "SUBJECT"))
    }
    
  } else if (lineType == "none"){
    linTxt  <- ""
  } else {
    linTxt  <-  paste(", linetype=", lineType)
  }
  
  ## shapeType
  if(is.na(shapeType)){ ## automaticatical attribution
    
    if(is.null(fixfactC)){
      shaTxt  <- ""
    } else {
      shaTxt <- paste(", shape=",
                      ifelse(colorType == "SUBJECT", "FIXFACT", "SUBJECT"))
    }
    
  } else if (shapeType == "none"){
    shaTxt  <- ""
  } else {
    shaTxt  <-  paste(", shape=", shapeType)
  }
  
  
  ## aes mapping
  txtMap <- paste("aes(x = TIME, y = DV",
                  colTxt, shaTxt, ")", sep = "")
  
  txtLineMap <- paste("aes(x = TIME, y = DV, group = GROUP ",
                      colTxt, linTxt,  ")", sep = "")
  
  
  ## plot and output
  p <- ggplot(data = DF,mapping = eval(parse(text = txtMap))) +
    ggtitle(titC)+
    xlab(xlabC)+ylab(ylabC) +
    theme(legend.position="left",
          plot.title = element_text(size = rel(1.2), face = "bold"))+
    geom_point()+
    geom_line(eval(parse(text = txtLineMap)), data = meanDF)+
    theme_bw()+
    NULL
  
  if(plotL) plot(p)
  
  invisible(p)
}



#######################################################################################################
## Post-hoc estimate
#######################################################################################################

#' Visualization of fixed effects (post-hoc estimates)
#'
#' Description
#'
#' @param mfl A linear mixed model fitted via lmer or a data frame containing data
#' @param  pvalCutof User pvalue cut of
#' @param plotL Boolean
#' @param titC Title of the plot
#' @return A plot
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export plot_posthoc

plot_posthoc <- function(mfl, pvalCutof = 0.05, plotL =TRUE, titC = "Post-hoc estimates"){
  
  #ddlsm1  <- data.frame(difflsmeans(mfl,test.effs=NULL)$diffs.lsmeans.table) ## OLD versions
  ddlsm1  <- as.data.frame(difflsmeans(mfl,test.effs=NULL))
  
  colnames(ddlsm1)[ncol(ddlsm1)] <- "pvalue"
  ddlsm1$Significance <- rep("NS", nrow(ddlsm1))
  
  ## modif JF pour tenir compte du seuil de pvalues defini par le user
  ddlsm1$Significance[which(ddlsm1$pvalue <pvalCutof)] <- paste("p-value < ", pvalCutof, sep = "")
  ddlsm1$Significance[which(ddlsm1$pvalue <pvalCutof/5)] <- paste("p-value < ", pvalCutof/5, sep = "")
  ddlsm1$Significance[which(ddlsm1$pvalue <pvalCutof/10)] <- paste("p-value < ", pvalCutof/10, sep = "")
  
  ddlsm1$levels <- rownames(ddlsm1)
  ddlsm1$term <- sapply(rownames(ddlsm1),function(namC){
    strsplit(namC, split = " ", fixed = TRUE)[[1]][1]
  })
  
  
  
  colValue <- c("grey", "yellow","orange","red")
  names(colValue) <- c("NS",
                       paste("p-value < ", pvalCutof, sep = ""),
                       paste("p-value < ", pvalCutof/5, sep = ""),
                       paste("p-value < ", pvalCutof/10, sep = ""))
  
  p <- ggplot(ddlsm1, aes(x = levels, y = Estimate))+
    facet_grid(facets = ~term, ddlsm1,scales = "free", space = "free")+
    geom_bar( aes(fill = Significance), stat="identity")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_fill_manual(values = colValue)+
    #geom_errorbar(aes(ymin = Lower.CI, ymax =Upper.CI ), width=0.25)+ ## OLD versions
    geom_errorbar(aes(ymin = lower, ymax =upper ), width=0.25)+
    ggtitle(titC)+xlab("")+
    #theme(plot.title = element_text(size = rel(1.2), face = "bold"))+
    NULL
  
  if(plotL) plot(p)
  invisible(p)
  
}




#######################################################################################################
## Visualisation des effets aléatoires
#######################################################################################################

#' Visualization of random effects
#'
#' Equivalent of dotplot(ranef)
#'
#' @param mfl A linear mixed model fitted via lmer or a data frame containing data
#' @param  plotL Logical
#' @return A plot
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export plot_randomEffect


plot_randomEffect <- function(mfl, plotL = TRUE){
  
  ## Estimation et format des effets aléatoires
  randomEffect <- ranef(mfl, condVar = TRUE)
  
  DF <- data.frame(randomEffect = rep(names(randomEffect),
                                      times  = sapply(1:length(randomEffect),
                                                      function(lsi){return(length(unlist(randomEffect[[lsi]])))})))
  
  DF$condVar <- DF$estimate <- DF$x2 <- DF$x1 <- rep(NA,nrow(DF))
  
  for(rafC in names(randomEffect)){
    
    eff <-randomEffect[[rafC]]
    
    DF$x1[which(DF$randomEffect == rafC)] <- rep(colnames(eff), each = nrow(eff))
    DF$x2[which(DF$randomEffect == rafC)] <- rep(rownames(eff), ncol(eff))
    DF$estimate[which(DF$randomEffect == rafC)] <- unlist(eff)
    
    condvar <- attr(randomEffect[[rafC]], "postVar")
    
    se <- NULL
    for(coli in 1:ncol(eff)){
      se <- c(se,
              sapply(1:nrow(eff),
                     function(i){return(condvar[coli,coli,i])}))
    }
    DF$condVar[which(DF$randomEffect == rafC)] <- se
    
  }
  DF$se <- sqrt(DF$condVar)
  DF$lower <- DF$estimate-1.96*DF$se
  DF$upper <- DF$estimate+1.96*DF$se
  
  
  
  
  ## Plot
  plotLs <-vector("list", length(randomEffect))
  names(plotLs) <- names(randomEffect)
  
  for(pi in 1:length(plotLs)){
    
    subDF <- DF[DF$randomEffect == names(plotLs)[pi], ]
    subDF <- subDF[order(subDF$x1, subDF$estimate, decreasing = FALSE),]
    #subDF$x2 <-factor(subDF$x2, levels=subDF$x2)
    
    
    p <- ggplot(data = subDF,
                mapping = aes(x = estimate, y = reorder(x2, estimate))
    )+
      geom_point(size =3)+
      geom_segment(aes(xend = lower, yend = x2)) +
      geom_segment(aes(xend = upper, yend = x2))+
      facet_wrap(~x1, ncol = length(unique(subDF$x1)))+
      ylab("")+xlab("")+
      ggtitle(paste("Random effect - ",names(plotLs)[pi], sep = ""))+
      theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))+
      geom_vline(xintercept = 0, linetype = "dashed")+
      theme_bw()
    
    
    plotLs[[pi]] <- p
    
    if(plotL) plot(p)
  }
  
  
  
  invisible(plotLs)
  
}




#######################################################################################################
## Linearité des effets et outlying observations
#######################################################################################################


#' Linarity of the fixed effect with regard to the continuous time
#'
#' @param diagLs diagnostic list
#' @param hlimitN Limit value for outliers (e.g.2 or 3)
#' @param plotL Boolean
#' @param label_factor Column of observation names used to label outlying values
#' @return A plot
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export plot_linearity
#'


plot_linearity <- function(diagLs, hlimitN, plotL = TRUE, label_factor = NULL){
  
  df <- cbind.data.frame(diagLs$data,
                         marginal.prediction = diagLs$marginal.prediction,
                         standardized.marginal.residuals = diagLs$std.marginal.residuals)
  
  # outlier annotation
  df$outliers <- rep("", nrow(df))
  outidx <- which(abs(df$standardized.marginal.residuals)>hlimitN)
  df[outidx, "outliers"] <- (1:nrow(df))[outidx]
  
  if(length(label_factor) >=1){
    df[outidx, "outliers"] <- paste(df[outidx, "outliers"],
                                    df[outidx, label_factor[1]],
                                    sep = "_")
    if(length(label_factor) >1){
      for(li in 2:length(label_factor)){
        df[outidx, "outliers"] <- paste(df[outidx, "outliers"],
                                        df[outidx, label_factor[li]],
                                        sep = ".")
      }
    }
    
  }
  
  # if(!is.null(label_factor)){
  #   df[outidx, "outliers"] <- paste(df[outidx, "outliers"],
  #                                   ifelse(length(label_factor)==1,
  #                                          df[outidx, label_factor],
  #                                          apply(df[outidx, label_factor], 1, paste, collapse = ".")),
  #                                   sep = "_")
  # }
  
  
  
  p <- ggplot(data = df,
              aes(x=marginal.prediction,
                  y=standardized.marginal.residuals)) +
    geom_point(size =2) +
    geom_hline(yintercept = 0, col = "grey")+
    geom_smooth(aes(x=marginal.prediction,
                    y=standardized.marginal.residuals), data = df,  se = FALSE, col = "blue", method = "loess")+
    ggtitle("Linearity of effects/outlying obervations")+
    xlab("Marginal predictions")+
    ylab("Standardized marginal residuals")+
    theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))+
    geom_hline(yintercept = c(-1,1)*hlimitN, linetype = "dashed")+
    geom_text(aes(label = outliers), hjust=0, vjust=0)
  
  if(plotL) plot(p)
  
  invisible(p)
  
}








#######################################################################################################
## EBLUP
#######################################################################################################


#' Mahalanobis distance
#'
#' @param diagLs diagnostic list
#' @param plotL Boolean
#' @return A plot
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export plot_mahalanobis
#'


plot_mahalanobis <- function(diagLs,  plotL = TRUE){
  
  unitDf <- data.frame(unit = names(diagLs$std.mahalanobis.distance),
                       mal = diagLs$std.mahalanobis.distance)
  
  
  ## Outlying subjects
  p <-
    ggplot(aes(y = mal, x= unit), data = unitDf)+
    geom_point(size =3)+
    ylab("Standardized Mahalanobis distance")+
    geom_vline(xintercept = 0, linetype = "dashed")+
    theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))+
    geom_hline(yintercept = 2*mean(unitDf$mal), linetype = "dashed")+
    geom_text(aes(label = unit),
              data = unitDf[unitDf$mal>2*mean(unitDf$mal), ],
              hjust=1, vjust=0)+
    ggtitle("Outlying unit")+
    xlab("unit")
  
  if(plotL) plot(p)
  
  invisible(p)
  
}






#' Mahalanobis distance (Chi2)
#'
#' @param diagLs diagnostic list
#' @param plotL aa
#' @return A plot
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export plot_mahalanobisKhi2
#'


plot_mahalanobisKhi2 <- function(diagLs,  plotL = TRUE){
  
  unitDf <- data.frame(unit = names(diagLs$std.mahalanobis.distance),
                       mal = diagLs$mahalanobis.distance)
  
  p <-qqplotF(x = unitDf$mal,
              distribution = "chisq",
              df= diagLs$q,
              line.estimate = NULL,
              conf = 0.95)+
    xlab("Chi-squared quantiles")+
    ylab("Mahalanobis distance")+
    ggtitle("Normality of random effect")+
    theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))
  
  if(plotL) plot(p)
  invisible(p)
}







#######################################################################################################
## Residus conditionels
#######################################################################################################

## Presence of outlying observations and homoscedacity of residuals

#' Homoskedacity of conditionalresiduals
#'
#' @param diagLs diagnostic list
#' @param hlimitN Limit value for outliers (e.g.2 or 3)
#' @param plotL Boolean
#' @param label_factor Column of observation names used to label outlying values
#' @return A plot
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export plot_conditionalResiduals
#'


plot_conditionalResiduals <-  function(diagLs, hlimitN, plotL = TRUE, label_factor = NULL){
  
  df <- cbind.data.frame(diagLs$data,
                         conditional.prediction = diagLs$conditional.prediction,
                         standardized.conditional.residuals = diagLs$std.conditional.residuals)
  
  # outlier annotation
  df$outliers <- rep("", nrow(df))
  outidx <- which(abs(df$standardized.conditional.residuals)>hlimitN)
  df[outidx, "outliers"] <- (1:nrow(df))[outidx]
  
  if(length(label_factor) >=1){
    df[outidx, "outliers"] <- paste(df[outidx, "outliers"],
                                    df[outidx, label_factor[1]],
                                    sep = "_")
    if(length(label_factor) >1){
      for(li in 2:length(label_factor)){
        df[outidx, "outliers"] <- paste(df[outidx, "outliers"],
                                        df[outidx, label_factor[li]],
                                        sep = ".")
      }
    }
    
  }
  
  p <- ggplot(data = df,
              aes(x=conditional.prediction,
                  y=standardized.conditional.residuals)) +
    geom_point(size =2) +
    geom_hline(yintercept = 0, col = "grey")+
    geom_smooth(aes(x=conditional.prediction,
                    y=standardized.conditional.residuals),
                data = df,  se = FALSE, col = "blue", method = "loess")+
    ggtitle("Homoscedasticity of conditional residuals/outlying observations")+
    xlab("Individual predictions")+
    ylab("Standardized conditional residuals")+
    theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))+
    geom_hline(yintercept = c(-1,1)*hlimitN, linetype = "dashed")+
    geom_text(aes(label = outliers), hjust=0, vjust=0)
  
  if(plotL) plot(p)
  invisible(p)
}




#' Normality of conditionalresiduals
#'
#' @param diagLs diagnostic list
#' @param plotL aa
#' @return A plot
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export plot_condresQQplot
#'


plot_condresQQplot <-  function(diagLs, plotL = TRUE){
  
  df <- cbind.data.frame(diagLs$data,
                         conditional.prediction = diagLs$conditional.prediction,
                         standardized.conditional.residuals = diagLs$std.conditional.residuals)
  
  p <-qqplotF(x = df$standardized.conditional.residuals,
              distribution = "norm",
              line.estimate = NULL,
              conf = 0.95)+
    xlab("Standard normal quantiles")+
    ylab("Standardized conditional residual quantiles")+
    ggtitle("Normality of conditional error")+
    theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))
  
  if(plotL) plot(p)
  invisible(p)
}





#######################################################################################################
## Within-units covariance structure
#######################################################################################################


#' Lesaffre-Veerbeke measure
#'
#' @param diagLs diagnostic list
#' @param plotL aa
#' @return A plot
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export plot_lesaffreVeerbeke
#'


plot_lesaffreVeerbeke <- function(diagLs,  plotL = TRUE){
  
  unitDf <- data.frame(unit = names(diagLs$std.lesaffreverbeke.measure),
                       lvm = diagLs$std.lesaffreverbeke.measure)
  
  p <- ggplot(data = unitDf,
              aes(x=unit,
                  y=lvm)) +
    geom_point(size =2) +
    theme(legend.position="none")+
    xlab("units")+
    ylab("Standardized Lesaffre-Verbeke measure")+
    geom_hline(yintercept = 2*mean(unitDf$lvm), linetype = "dashed")+
    geom_text(aes(label = unit),
              data = unitDf[unitDf$lvm>2*mean(unitDf$lvm), ],
              hjust=0, vjust=0)+
    ggtitle("Within-units covariance matrice")+
    theme(legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))
  
  if(plotL) plot(p)
  invisible(p)
}

##-------------------------------------------------------------------------------------------------##
## Helpers
##-------------------------------------------------------------------------------------------------##


## square root of a matrix
## From Rocha, Singer and Nobre

#' square root of a matrix (Rocha)
#'
#' Description
#'
#' @param mat Matrix
#' @return A list
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export sqrt.matrix

sqrt.matrix <- function(mat) {
  mat <- as.matrix(mat)  # new line of code
  singular_dec<-svd(mat,LINPACK=F)
  U<-singular_dec$u
  V<-singular_dec$v
  D<-diag(singular_dec$d)
  sqrtmatrix<-U%*%sqrt(D)%*%t(V)
  #  return(list(sqrt=sqrtmatrix))
}


## square root of a matrix
## http://www.cs.toronto.edu/~jepson/csc420/notes/introSVD.pdf (page 6)
## (for matMN a n x n matrix that symetric and non-negative definite)

#' square root of a matrix (Rocha)
#'
#' @param mat Matrix
#' @return A list
#' @author Natacha Lenuzza
#' @examples
#' print("hello !")
#'
#' @export sqrtmF

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


plot.res.Lmixed <- function(mfl, df, title = "", pvalCutof = 0.05) {
  ## define subscript of the different columns depending if we have only time (ncol(df)=3) or not
  if (ncol(df) >3) {
    varidx <- 4
    ffidx <- 1
    timidx <- 2
    individx <- 3
    
  } else 
  {
    varidx <- 3
    ffidx <- 1
    timidx <- 1
    individx <-2
  }
  nameVar <- colnames(df)[varidx]
  fflab <- colnames(df)[ffidx]
  ## Individual time-course
  rawPlot <- 
    ggplot(data = df, aes(x=df[[timidx]], y=df[[varidx]], colour=df[[ffidx]], group=df[[individx]])) +
    geom_point()+
    geom_line() +  ggtitle("Individual time-courses (raw data)")+ 
    ylab(nameVar) +
    xlab(label = colnames(df)[2])+
    theme(legend.title = element_blank() , legend.position="none", plot.title = element_text(size = rel(1.2), face = "bold"))
  
  ## Boxplot of fixed factor
  bPlot <- 
    ggplot(data = df, aes(y=df[[varidx]], x=df[[ffidx]], color=df[[ffidx]]))+
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
  
  ddlsm1$levels <- rownames(ddlsm1)
  ddlsm1$term <- sapply(rownames(ddlsm1),function(namC){
    strsplit(namC, split = " ", fixed = TRUE)[[1]][1]
  })
  
  
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