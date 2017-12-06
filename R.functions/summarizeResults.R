#################################
# R function
#
#  summarize.RF and summarize.KF
#  November 19th, 2016
#
#  Author: Ariadna Morales (ariadna.biologia@gmail.com)
#
#  Description: This function estimates the proportion of loci supporting a given model
#				
#				Before using this function look at https://github.com/ariadnamorales,
#				there might be an updated version.
#
#  Input: tables of retained distasced calculated by 'calcDistances.reject'

###################
#### Arguments ####

#pathRetainedRF

summarize.RF<-function(pathRetainedRF){

	## Summarize RF results
	list.summaryRF<-list.files(path=pathRetainedRF, pattern="*.txt")
	summaryRF<-c("locus"=c(), "model"=c(), "prop.Retained"=c(), "thresholValueRF"=c())

	locusTypeRF<-c("locusSupporting.tmp"=c(), "modelBestSupported.tmp"=c(), "propRet.modelBestSupported.tmp"=c())
	for(file in 1:length(list.summaryRF)){
		summary.tmp<-read.table(file=paste(pathRetainedRF, "/", list.summaryRF[file],sep=""), header=TRUE)
		summaryRF<-rbind(summaryRF, summary.tmp)
	
		## create a table of loci and what model support 
		locusSupporting.tmp<-summary.tmp[which(summary.tmp$prop.Retained == max(summary.tmp$prop.Retained)),1]
		modelBestSupported.tmp<-summary.tmp[which(summary.tmp$prop.Retained == max(summary.tmp$prop.Retained)),2]
		propRet.modelBestSupported.tmp<-summary.tmp[which(summary.tmp$prop.Retained == max(summary.tmp$prop.Retained)),3]
	
		locusTypeRF.tmp<-as.data.frame(cbind(locusSupporting.tmp, modelBestSupported.tmp, propRet.modelBestSupported.tmp))
		locusTypeRF<-rbind(locusTypeRF, locusTypeRF.tmp)
	}
	uninformativeLoci.RF<-locusTypeRF[which(duplicated(locusTypeRF$locusSupporting.tmp) == TRUE),]		## loci that appear more than once with equal %% values
	informativeLoci.RF<-locusTypeRF[!(duplicated(locusTypeRF$locusSupporting.tmp) | duplicated(locusTypeRF$locusSupporting.tmp, fromLast = TRUE) == TRUE), ]
	
	summaryTable.RF<-list("uninformativeLoci.RF"=uninformativeLoci.RF, "informativeLoci.RF"=informativeLoci.RF)
	return(summaryTable.RF)
	
	
	
}


###################
#### Arguments ####

#pathRetainedKF

summarize.KF<-function(pathRetainedKF){
	## Summarize KF results
	list.summaryKF<-list.files(path=pathRetainedKF, pattern="*.txt")
	summaryKF<-c("locus"=c(), "model"=c(), "prop.Retained"=c(), "thresholValueKF"=c())

	locusTypeKF<-c("locusSupporting.tmp"=c(), "modelBestSupported.tmp"=c(), "propRet.modelBestSupported.tmp"=c())
	for(file in 1:length(list.summaryKF)){
		summary.tmp<-read.table(file=paste(pathRetainedKF, "/", list.summaryKF[file],sep=""), header=TRUE)
		summaryKF<-rbind(summaryKF, summary.tmp)
	
		## create a table of loci and what model support 
		locusSupporting.tmp<-summary.tmp[which(summary.tmp$prop.Retained == max(summary.tmp$prop.Retained)),1]
		modelBestSupported.tmp<-summary.tmp[which(summary.tmp$prop.Retained == max(summary.tmp$prop.Retained)),2]
		propRet.modelBestSupported.tmp<-summary.tmp[which(summary.tmp$prop.Retained == max(summary.tmp$prop.Retained)),3]
	
		locusTypeKF.tmp<-as.data.frame(cbind(locusSupporting.tmp, modelBestSupported.tmp, propRet.modelBestSupported.tmp))
		locusTypeKF<-rbind(locusTypeKF, locusTypeKF.tmp)
	}
	uninformativeLoci.KF<-locusTypeKF[which(duplicated(locusTypeKF$locusSupporting.tmp) == TRUE),]		## loci that appear more than once with equal %% values
	informativeLoci.KF<-locusTypeKF[!(duplicated(locusTypeKF$locusSupporting.tmp) | duplicated(locusTypeKF$locusSupporting.tmp, fromLast = TRUE) == TRUE), ]
	
	summaryTable.KF<-list("uninformativeLoci.KF"=uninformativeLoci.KF, "uninformativeLoci.KF"=informativeLoci.KF)
	return(summaryTable.KF)
}