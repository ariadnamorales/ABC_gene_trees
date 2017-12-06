#################################
# R function
#
#  calculate.RF.FK.distance.R
#  November 19th, 2016
#
#  Author: Ariadna Morales (ariadna.biologia@gmail.com)
#
#  Description: This function calculates pairwise distances between gene trees
#               Robinson-Foulds (raw and normalized)
#               Kuhner-Felsestein (branch scrore distance)
#				Before using this function look at https://github.com/ariadnamorales,
#				there might be an updated version.
#
#  Dependencies in R: ape, phangorn

###################
#### Arguments ####

#simulatedTrees for models specified in 'modelRange'
#observedTrees
#modelRange (from a migrationArray file created in PHRAPL)


###############################################################################################################
#### Function: Loop to calculate RF an KF distance between all simulatedTrees and all subsampled geneTrees ####
###############################################################################################################
 
loopRF.dist<-function(simulatedTrees, observedTrees, modelRange){
	## Load dependencies 
    require(ape)
	require(phangorn)
	
    ## Create empty objects where results will be stored
    rfDist<-c()
	rfDist.norm<-c()
	kfDist<-c()
    
    ## Start loop
	for(gene in 1:length(observedTrees)){
		for(nTree in 1:length(observedTrees[[gene]])){
			
			## Calculate Robinson-Foulds distance
			rfDist.tmp<-RF.dist(simulatedTrees, observedTrees[[gene]][[nTree]], normalize=FALSE)
			rfDist<-append(rfDist, c(rfDist.tmp))
			
			## Calculate normalized Robinson-Foulds distance
			rfDist.norm.tmp<-RF.dist(simulatedTrees, observedTrees[[gene]][[nTree]], normalize=TRUE)
			rfDist.norm<-append(rfDist.norm, c(rfDist.norm.tmp))
            
			## Calculate KF distance (Kuhner-Felsestein branch scrore distance)
			kfDist.tmp<-KF.dist(simulatedTrees, observedTrees[[gene]][[nTree]])
			kfDist<-append(kfDist, c(kfDist.tmp))
		}
		addModel<-rep(modelRange, length(simulatedTrees)*length(observedTrees[[gene]]))
	}
	result.distances<-as.data.frame(cbind(addModel, rfDist, rfDist.norm, kfDist))
	return(result.distances)
}


####### END OF FUNCTION  loopRF.dist


################################################################################################
# R function
#
#  calcDistances.reject.R
#  
#
#  Author: Ariadna Morales (ariadna.biologia@gmail.com)
#
#  Description: Calculate Robinson-Foulds and Kuhner-Felsestein distances bewteen gene Trees
#                it also runs a rejection step based on 2% threshold for the distribution of KF
#
#  Dependencies in R: ape, phangorn, gsubfn
#                    function loopRF.dist

###################
#### Arguments ####

#simulatedTrees for models specified in 'modelRange'
#observedTrees
#modelRange (from a migrationArray file created in PHRAPL)

calcDistances.reject<-function(locus, pathOutput=NULL, pathSubsampledTrees, pathSimulatedTrees){

    ##################
    ### Load libraries 
    #library(gsubfn)
    library(ape)
    library(phangorn)


    #####################
    #### Run details ####
		if(is.null(pathOutput)){
		pathOutput<-getwd()
	}    
    pathSubsampledTrees=pathSubsampledTrees

    ### Load simulated Trees (in ms) for each model
    list.simulatedTrees<-list.files(path=pathSimulatedTrees, pattern="*.tre", full.names=TRUE)
    
    for(nModel in 1:length(list.simulatedTrees)){
    	assign(paste0("simulatedTrees.model", nModel), simulatedTrees.model1<-read.tree(list.simulatedTrees[nModel]))
    	print(paste0("Simulated Trees for model ", nModel," -----> Loaded)"))
    }
    
    ### Load subsampled geneTrees (in phrapl) ---> an object named "observedTrees"[[#gene]][[#subsampledTree]]
    load(paste(pathSubsampledTees, "/subsampledTrees_", locus, "_locus.rda",sep=""))
        print("Subsampled Trees  -----> Loaded")

    ### Run function loopRF.dis to calculate distances between simulated and observed Trees
        print("Calculating distances between simulated and observed trees ...")				##This step can take a while
        
        for(nModel in 1:length(list.simulatedTrees)){
        	assign(paste0("rfDist.model", nModel), loopRF.dist(get(paste0("simulatedTrees.model", nModel)), observedTrees, modelRange=nModel))
        	print(paste0("Distances for model ", nModel," -----> Calculated"))
        	}
        
	

    allModelsData<-do.call("rbind",mget(paste0("rfDist.model", 1:length(list.simulatedTrees))))
    #allModelsData<-as.data.frame(allModelsData)
    
    ### Save results
        print("Saving results ...")
        ifelse(!dir.exists(paste0(pathOutput, "/tablesDistances")), dir.create(paste0(pathOutput, "/tablesDistances")), FALSE)
    write.table(allModelsData, file=paste(pathOutput, "/tablesDistances/distances_simData_locus",locus,".txt",sep=""), quote = FALSE, row.names = FALSE)
    
    retainedData.KF<-data.frame("addModel"=c(), "rfDist"=c(), "rfDist.norm"=c(), "kfDist"=c())
    retainedData.RF<-data.frame("addModel"=c(), "rfDist"=c(), "rfDist.norm"=c(), "kfDist"=c())


	### Rejection step based on a threshold
        print("Starting rejection ... ")
    rejectionThresholdKF<-min(allModelsData$kfDist)+2*(max(allModelsData$kfDist)-min(allModelsData$kfDist))/100 	##2%
         print(paste("KF threshold value for this locus: ", rejectionThresholdKF,sep=""))

	rejectionThresholdRF<-min(allModelsData$rfDist)+2*(max(allModelsData$rfDist)-min(allModelsData$rfDist))/100 	##2%
        print(paste("RF threshold value for this locus: ", rejectionThresholdRF,sep=""))


	##############################
    ### Rejection using KF dist
    for(nrow in 1:nrow(allModelsData)){
        if(allModelsData$kfDist[nrow] <= rejectionThresholdKF)
			retainedData.KF<-rbind(retainedData.KF, allModelsData[nrow,])
    }
        print("Rejection KF step  -----> Done")
        
        
    ### Calculate proporton of retained models
    for(nModel in 1:length(list.simulatedTrees)){
    	p.retained.tmp<-length(retainedData.KF$kfDist.norm[which(retainedData.KF$addModel == nModel)])/length(retainedData.KF$kfDist)
    	assign(paste0("p.retainedModel", nModel, ".KF") , p.retained.tmp)
    	
    	print(paste("Proportion of Trees retained for model ", nModel, ": ", get(paste0("p.retainedModel", nModel, ".KF")),sep=""))
    	
    	assign(paste0("summaryRetained.model", nModel, ".KF"), c(locus, nModel, get(paste0("p.retainedModel", nModel, ".KF")), rejectionThresholdKF))
    	    	
    }
    
    summaryRetained.KF<-do.call("rbind",(mget(paste0("summaryRetained.model", 1:length(list.simulatedTrees), ".KF"))))
    summaryRetained.KF<-as.data.frame(summaryRetained.KF)
    
    colnames(summaryRetained.KF)[1]<-"locus"
    colnames(summaryRetained.KF)[2]<-"model"
    colnames(summaryRetained.KF)[3]<-"prop.Retained"
    colnames(summaryRetained.KF)[4]<-"thresholValueKF"
    
    
    
    ############################
    ### Rejection using RF dist
    for(nrow in 1:nrow(allModelsData)){
        if(allModelsData$rfDist[nrow] <= rejectionThresholdRF)
			retainedData.RF<-rbind(retainedData.RF, allModelsData[nrow,])
    }
        print("Rejection RF step  -----> Done")
        
        
    ### Calculate proporton of retained models
    for(nModel in 1:length(list.simulatedTrees)){
    	p.retained.tmp<-length(retainedData.RF$rfDist.norm[which(retainedData.RF$addModel == nModel)])/length(retainedData.RF$rfDist)
    	assign(paste0("p.retainedModel", nModel, ".RF") , p.retained.tmp)
    	
    	print(paste("Proportion of Trees retained for model ", nModel, ": ", get(paste0("p.retainedModel", nModel, ".RF")),sep=""))
    	
    	assign(paste0("summaryRetained.model", nModel, ".RF"), c(locus, nModel, get(paste0("p.retainedModel", nModel, ".RF")), rejectionThresholdRF))
    	    	
    }
    
    summaryRetained.RF<-do.call("rbind", (mget(paste0("summaryRetained.model", 1:length(list.simulatedTrees), ".RF"))))
    summaryRetained.RF<-as.data.frame(summaryRetained.RF)
    
    colnames(summaryRetained.RF)[1]<-"locus"
    colnames(summaryRetained.RF)[2]<-"model"
    colnames(summaryRetained.RF)[3]<-"prop.Retained"
    colnames(summaryRetained.RF)[4]<-"thresholValueKF"
            
    
    ### Save results
        print("Saving results ...")
    ifelse(!dir.exists(paste0(pathOutput, "/tablesRetained.KF")), dir.create(paste0(pathOutput, "/tablesRetained.KF")), FALSE)
	ifelse(!dir.exists(paste0(pathOutput, "/tablesRetained.RF")), dir.create(paste0(pathOutput, "/tablesRetained.RF")), FALSE)
	ifelse(!dir.exists(paste0(pathOutput, "/summaryRetained.KF")), dir.create(paste0(pathOutput, "/summaryRetained.KF")), FALSE)
	ifelse(!dir.exists(paste0(pathOutput, "/summaryRetained.RF")), dir.create(paste0(pathOutput, "/summaryRetained.RF")), FALSE)
    
    
    write.table(retainedData.KF, file=paste(pathOutput, "/tablesRetained.KF/KF.retained_2percent_locus",locus,".txt",sep=""), quote = FALSE, row.names = FALSE)
    write.table(summaryRetained.KF, file=paste(pathOutput, "/summaryRetained.KF/KF.summRetained_2percent_locus",locus,".txt",sep=""), quote = FALSE, row.names = FALSE)
    
    write.table(retainedData.RF, file=paste(pathOutput, "/tablesRetained.RF/RF.retained_2percent_locus",locus,".txt",sep=""), quote = FALSE, row.names = FALSE)
    write.table(summaryRetained.RF, file=paste(pathOutput, "/summaryRetained.RF/RF.summRetained_2percent_locus",locus,".txt",sep=""), quote = FALSE, row.names = FALSE)
    
    
        print("DONE")
    }
### END FUNCTION
