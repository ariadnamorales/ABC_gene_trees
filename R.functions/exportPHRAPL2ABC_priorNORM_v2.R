#################################
# R function
#
#  exportPHRAPL2ABC_priorNORM.R
#
#
#  Author: Ariadna Morales (ariadna.biologia@gmail.com)
#
#  Description: This function simulates trees and calculates summary statistics in ms under specified demographic models (built in PHRAPL, in a migrationArray object) 
#					the main goal is to use the output gene trees to compare them with empirical gene trees
#					Prior distributions for theta and migration is uniform, and for divergence times is normal.
#
#  Dependencies: ms, PHRAPL, modified functions from PHRAPL (batchMS.R)

###################
#### Arguments ####

#migrationArray				#object created in PHRAPL with demographic models see: https://github.com/ariadnamorales/phrapl-manual/blob/master/3.Generate_set_of_models.Rmd 
#Priorsize 					#number of trees that will be simulated per model
#popVector					#vector with number of individuals/tips per population/species that will be simulated
#nTrees=1 					#number of trees that will be simulated in each draw, to avoid redundancy set it at 1
#modelRange 				#range of models from 'migrationArray' that will be used for simulations
#thetaRange 				#minimum and maximum values for theta (popSize) that will be used for simulations assuming a uniform distribution

#nDivEvents
#meanDivEvents
#sdDivEvents

#migrationRange 			#minimum and maximum values for migration rate that will be used for simulations assuming a uniform distribution
#pathOutput 				#path where simulated trees will be saved, if null they will be saved in current working directory
#msLocation="ms" 			#path to ms
#SampleStats=FALSE			#specify if summary stats will be calculated in ms and save in text file. IMPORTANT, if TRUE "sample_stats" must be installed in 'ms'
#print.ms.string=FALSE		#print an save string used for simulations in ms
#print.screenout.ABC=FALSE	#print an save string screenout of ms

#################################################################################################
#### Function export phrapl models to ms ----> simulate Trees and calculate Sample_Stats ABC ####
#################################################################################################

exportPHRAPL2ABC_normPrior<-function(migrationArray=migrationArray, Priorsize, popVector, nTrees=1, modelRange, thetaRange, nDivEvents, meanDivEvents, sdDivEvents, migrationRange, pathOutput=NULL, msLocation="ms", print.ms.string=FALSE, print.screenout.ABC=FALSE, SampleStats=FALSE){
    library(phrapl)
    
    ms.strings.ALL<-c()
    trees.sim.ms<-c()
    modelRange<-as.numeric(modelRange)
	
	if(length(popVector) != dim(migrationArray[[1]]$collapseMatrix)[1]){
		stop("error in popVector ---> wrong number of Pops")
	}
	
	if(is.null(pathOutput)){
		pathOutput<-getwd()
	}
	ifelse(!dir.exists(paste0(pathOutput, "/simulatedTrees")), dir.create(paste0(pathOutput, "/simulatedTrees")), FALSE)
	
	#load(paste(migrationArray, sep=""))
    for(model in modelRange){
        
        migrationIndividual<-migrationArray[[model]]
        parameterVector<-msIndividualParameters(migrationIndividual)
        #ms.string<-createMSstringSpecific(popVector,migrationIndividual,parameterVector,nTrees)
		
		
		screenout.ABC<-c()
    	for(priorsize in 1:Priorsize){
            	#
            	ms.string<-createMSstringSpecific(popVector,migrationIndividual,parameterVector,nTrees)
            	#Generate random numbers (NORMAL distribution for divergence times)
            		theta<-runif(1,min(thetaRange), max(thetaRange))
            		mig<-runif(1,min(migrationRange), max(migrationRange))
            		## Generate one variable per divergence event
            		for(divEventID in 1: nDivEvents){
            			assign(paste0("div", divEventID),rnorm(1, meanDivEvents[divEventID], sdDivEvents[divEventID]))	
            		}

            	#Replace paramvalues in ms string
            	ms.string$opts<-gsub("n0multiplier_1", theta, ms.string$opts)
            	ms.string$opts<-gsub("migration_1", mig, ms.string$opts)
            	for(divEventID in 1: nDivEvents){
            		ms.string$opts<-gsub(paste0("collapse_", divEventID), get(paste0("div", divEventID)), ms.string$opts)
            		
            	}
            	          
				#merge string parameters -----> command Line in a File
            	ms.string.cmd<-paste(msLocation, ms.string$nsam, ms.string$nreps, "-t ",runif(1,min(thetaRange), max(thetaRange)),ms.string$opts,sep=" ")
            
           		#merge string parameters ---> print Trees
            	ms.string.Trees<-paste(ms.string.cmd, " | sed '5!d' >> ",paste(pathOutput, "/simulatedTrees/ABC.output.Trees_model",model,".tre", sep=""),sep=" ")

				# write ms.string.cmd to file
				if(print.ms.string == TRUE){
					ifelse(!dir.exists(paste0(pathOutput, "/ms_cmdLines")), dir.create(paste0(pathOutput, "/ms_cmdLines")), FALSE)
            		fileName<-file(paste(pathOutput,"/ms_cmdLines/ms_cmdLine_model",model,".sh",sep=""))
            			writeLines(c(ms.string.cmd),fileName)
            		close(fileName)
            	}

           	#run ms
              	system(ms.string.Trees)
            	#print(ms.string)

            	# write screenout ABC file
            	if(print.screenout.ABC == TRUE){
            		ifelse(!dir.exists(paste0(pathOutput, "/ABC.screenout")), dir.create(paste0(pathOutput, "/ABC.screenout")), FALSE)
            		screenout.ABC.temp<-paste(model, theta, mig, as.character(unlist(mget(paste0("div", 1:nDivEvents)))))
            		screenout.ABC<-append(screenout.ABC, c(screenout.ABC.temp))
            		write(screenout.ABC, file=paste(pathOutput,"/ABC.screenout/ABC.screenout_model",model,".txt",sep=""), append=FALSE)
            	}
            	
            
            	if(SampleStats == TRUE){
            		#merge string parameters ---> SampleStats
            		ifelse(!dir.exists(paste0(pathOutput, "/sampleStats")), dir.create(paste0(pathOutput, "/sampleStats")), FALSE)
            		ms.string.SampleStats<-paste(ms.string.cmd, " | sample_stats >>",paste(pathOutput, "/sampleStats/ABC.output.SumStats_model",model,".txt", sep=""),sep=" ")
            	
            		## Run ms and Sample_Satats
            		system(ms.string.SampleStats)
            	
            		## Save priors and output (sumStats) in the same file
            		paste.string<-paste0("paste ", pathOutput,"/ABC.screenout_model",model,".txt ", pathOutput, "/sampleStats/ABC.output.SumStats_model",model,".txt > ", pathOutput, "/simpleABC.prior_model",model,".txt")
        			system(paste.string)
        			#print(paste.string)
           		}else{
           			print(paste0("Sample stats for rep ",priorsize," in model ",model," will not be calculated"))
           		}
           
     	}
     	print(paste("done model ", model, " with ", priorsize, "reps", sep=""))
    }
}
####### END OF FUNCTION
