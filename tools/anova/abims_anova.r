#!/usr/local/public/bin/Rscript
# version="1.1"

# date: 06-06-2012
# update: 18-02-2014
# **Authors** Gildas Le Corguille  ABiMS - UPMC/CNRS - Station Biologique de Roscoff - gildas.lecorguille|at|sb-roscoff.fr

# abims_anova.r version 20140218

library(batch)


# function avova
anova = function (file, sampleinfo, varinfo, mode="column", condition=1, interaction=F, method="BH", threshold=0.01, selection_method="intersection", sep=";", dec=".", outputdatapvalue="anova.data.output", outputdatasignif="anova.datasignif.output") {


	if (sep=="tabulation") sep="\t"
    	if (sep=="semicolon") sep=";"
    	if (sep=="comma") sep=","

	anova_formula_operator = "+"
	if (interaction) anova_formula_operator = "*"

  	# -- import --
	data=read.table(file, header = TRUE, row.names=1, sep = sep, quote="\"", dec = dec, fill = TRUE, comment.char="",na.strings = "NA", check.names=FALSE)

  	if (mode == "row") {data=t(data)} else {data=as.matrix(data)}

	sampleinfoTab=read.table(sampleinfo, header = TRUE, row.names=1, sep = sep, quote="\"")
	rownames(sampleinfoTab) = make.names(rownames(sampleinfoTab))

	varinfoTab=read.table(varinfo, header = TRUE, sep = sep, quote="\"")
	if(sum(colnames(data)!=varinfoTab[,1])!=0){ # if ID not exactly identical
		if(sum(colnames(data)[order(colnames(data))]!=varinfoTab[order(varinfoTab[,1]),1])==0){
			# reorder data matrix to match variable metadata order
			data = data[,match(varinfoTab[,1],colnames(data))]
		}else{
			stop(paste0("\nVariables' ID do not match between your data matrix and your variable",
			"metadata file. \nPlease check your data."))
		}
	}

	# -- group --
	match_data_sampleinfoTab = match(rownames(data),rownames(sampleinfoTab))
	if (sum(is.na(match_data_sampleinfoTab)) > 0) {
	  write("ERROR: There is a problem during to match sample names from the data matrix and from the sample info (presence of NA).", stderr())
	  write("You may need to use change the mode (column/row)", stderr())
	  write("10 first sample names in the data matrix:", stderr())
	  write(head(colnames(data)), stderr())
	  write("10 first sample names in the sample info:", stderr())
	  write(head(rownames(sampleinfoTab)), stderr())
	  quit("no",status=10)
	}


	# -- anova --

  	# formula
	grps=list()
	anova_formula_s = "data ~ "
	cat("\ncontrasts:\n")
	for (i in 1:length(condition)) {
	  grps[[i]] = factor(sampleinfoTab[,condition[i]][match_data_sampleinfoTab])
	  anova_formula_s = paste(anova_formula_s, "grps[[",i,"]]",anova_formula_operator, sep="")
	  cat(condition[i],"\t",levels(grps[[i]]),"\n")
	# write("Current groups: ", stderr())
	# write(grp[[i]], stderr())
	}
	anova_formula_s = substr(anova_formula_s, 1, nchar(anova_formula_s)-1)
	anova_formula = as.formula(anova_formula_s)



	# anova
	manovaObjectList = manova(anova_formula)
	manovaList = summary.aov(manovaObjectList)

  	# condition renaming
	manovaRownames = gsub(" ","",rownames(manovaList[[1]]))
	manovaNbrPvalue = length(manovaRownames)-1
	manovaRownames = manovaRownames[-(manovaNbrPvalue+1)]

	for (i in 1:length(condition)) {
	  manovaRownames = sub(paste("grps\\[\\[",i,"\\]\\]",sep=""),condition[i],manovaRownames)
  	  anova_formula_s = sub(paste("grps\\[\\[",i,"\\]\\]",sep=""),condition[i],anova_formula_s)
	}

  	# log
  	cat("\nanova_formula",anova_formula_s,"\n")

	# p-value
	aovPValue = sapply(manovaList,function(x){x[-(manovaNbrPvalue+1),5]})
	if(length(condition) == 1) aovPValue = t(aovPValue)
	rownames(aovPValue) = paste("pvalue_",manovaRownames,sep="")

	# p-value adjusted
	if(length(condition) == 1) {
		aovAdjPValue = t(p.adjust(aovPValue,method=method))
	} else {
		aovAdjPValue = t(apply(aovPValue,1,p.adjust, method=method))
	}
	rownames(aovAdjPValue) = paste("pval.",method,".",manovaRownames,sep="")

	# selection
	colSumThreshold = colSums(aovAdjPValue <= threshold)
	if (selection_method == "intersection") {
		datafiltered = data[,colSumThreshold == nrow(aovAdjPValue ), drop=FALSE]
	} else {
		datafiltered = data[,colSumThreshold != 0, drop=FALSE]
	}
	selected.var = rep("no",ncol(data))
	selected.var[colnames(data)%in%colnames(datafiltered)] = "yes"

	#data=rbind(data, aovPValue, aovAdjPValue)
	varinfoTab=cbind(varinfoTab, round(t(aovAdjPValue),10), selected.var)

	# group means
	for (i in 1:length(condition)) {
		for(j in levels(grps[[i]])){
			subgp = rownames(sampleinfoTab[which(sampleinfoTab[,condition[i]]==j),])
			modmean = colMeans(data[which(rownames(data)%in%subgp),],na.rm=TRUE)
			varinfoTab=cbind(varinfoTab, modmean)
			colnames(varinfoTab)[ncol(varinfoTab)] = paste0("Mean_",condition[i],".",j)
		}
	}
	colnames(varinfoTab) = make.unique(colnames(varinfoTab))

	# pdf for significant variables
	pdf(outputdatasignif)
	### Venn diagram
	if(nrow(aovAdjPValue)>5){
		pie(100,labels=NA,main=paste0("Venn diagram only available for maximum 5 terms\n",
		    "(your analysis concerns ",nrow(aovAdjPValue)," terms)\n\n",
			"Number of significant variables relatively to\nyour chosen threshold and ",
			"selection method: ",ncol(datafiltered)),cex.main=0.8)
	}else{
		vennlist = list(NULL)
		names(vennlist) = rownames(aovAdjPValue)[1]
		if(length(colnames(aovAdjPValue))==0){colnames(aovAdjPValue)=varinfoTab[,1]}
		for(i in 1:nrow(aovAdjPValue)){
			vennlist[[rownames(aovAdjPValue)[i]]]=colnames(aovAdjPValue[i,which(aovAdjPValue[i,]<=threshold),drop=FALSE])
		}
		if(length(vennlist)==0){ pie(100,labels=NA,main="No significant ions was found.")
		}else{ library(venn) ; venn(vennlist, zcolor="style", cexil=2, cexsn=1.5) }
	}
	### Boxplot
	par(mfrow=c(2,2),mai=rep(0.5,4))
	data <- as.data.frame(data)
	for(i in 1:nrow(aovAdjPValue)){
		factmain = gsub(paste0("pval.",method,"."),"",rownames(aovAdjPValue)[i])
		factsignif = aovAdjPValue[i,which(aovAdjPValue[i,]<=threshold),drop=FALSE]
		if((ncol(factsignif)!=0)&(factmain%in%colnames(sampleinfoTab))){
			for(j in 1:ncol(factsignif)){
				varsignif = gsub(" Response ","",colnames(factsignif)[j])
				boxplot(as.formula(paste0("data$",varsignif," ~ sampleinfoTab$",factmain)),
				        main=paste0(factmain,"\n",varsignif), col="grey", mai=7)
			}
		}
	}
	dev.off()

	# summary for significant variables
	cat("\n\n- - - - - - - number of significant variables - - - - - - -\n\n")
	for(i in 1:nrow(aovAdjPValue)){
		cat(rownames(aovAdjPValue)[i],"-",
		    sum(aovAdjPValue[i,]<=threshold),"significant variable(s)\n")
	}
	cat("\nIf some of your factors are missing here, this may be due to\neffects",
	    "not estimable; your design may not be balanced enough.\n")

	# -- output / return --
	write.table(varinfoTab, outputdatapvalue, sep=sep, quote=F, row.names=FALSE)

	# log
	cat("\nthreshold:",threshold,"\n")
	cat("result:",ncol(datafiltered),"/",ncol(data),"\n\n")

	quit("no",status=0)
}

# log
cat("ANOVA\n\n")
cat("Arguments\n")
args <- commandArgs(trailingOnly = TRUE)
print(args)

listArguments = parseCommandArgs(evaluate=FALSE)
do.call(anova, listArguments)
