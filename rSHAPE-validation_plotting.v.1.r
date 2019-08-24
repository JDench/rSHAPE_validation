# This script is to make the plots for my landscapeSimulation study
library(RSQLite) # In case we need more info
library(scales)
library(cowplot) # this for the grid function for ggdraw
library(reshape2)
library(Cairo)
library(lemon) # this allows coord_capped_cart to cap ends of plot lines
library(lme4)
library(lmerTest)
library(tidyverse)
library(rSHAPE)

rm(list=ls())

options(stringsAsFactors = FALSE)

############################## FUNCTIONS ################################

# This allows the opacity, saturation and hue of colours to be adjusted.  The purpose is to allow
# a range of basic named colours to be dynamically manipulted into altered vectors of the colour.
colTransparency <- function(func_cols, func_opacity=1, func_scaleSaturation = NA, func_scaleValue = NA){
     # In case the user has passed scaling values outside of the 0-1 range we simply return them to being values of 1
     funcParms = list("Saturation"= func_scaleSaturation,
                      "Value"= func_scaleValue,
                      "Opacity"= func_opacity)
     # We check that the opacity value is within the range of 0-1
     for(thisParm in names(funcParms)){
          # If the user has not defined anything, then let it pass....
          if(all(!is.na(funcParms[[thisParm]]))){
               if(any(funcParms[[thisParm]] > 1)){
                    funcParms[[thisParm]][which(funcParms[thisParm] > 1)] <- 1
               } else if (any(funcParms[[thisParm]] < 0)){
                    funcParms[[thisParm]][which(funcParms[thisParm] < 0)] <- 0
               }
          }
     }

     # After extracting the RGB balues of a colour we change the values based on the opacity.
     tmpReturn <- rgb2hsv(sapply(func_cols, col2rgb))
     return( sapply(1:ncol(tmpReturn),function(x){
          hsv(h = tmpReturn[1,x],
              s = if(any(is.na(funcParms[["Saturation"]]))){
                   tmpReturn[2,x]
              } else {
                   funcParms[["Saturation"]]
              },
              v = if(any(is.na(funcParms[["Value"]]))){
                   tmpReturn[3,x]
              } else {
                   funcParms[["Value"]]
              },
              alpha= funcParms[["Opacity"]])
     }) )
}

############################ PRE-DEFINITIONS ############################
# This is the post-analysis script for my rSHAPE-validation
workDir <- "R:/Research_Data/SimulationRuns/rSHAPE-validation/"
postDir <- paste(workDir,"postAnal/",sep="")
outDir <- postDir

# This is the base batchString
ref_batchString <- "valid-rSHAPE"
# These are the parameters of the experiment
load(paste(postDir,"jobParameters_from_valid-rSHAPE.RData",sep=""))
# These are the files from the experiment
load(paste(postDir,"allFiles_from_valid-rSHAPE.RData",sep=""))

# I find all the post analysis files
all_popDemo <- list.files(postDir, pattern = "popDemographics_from_popDemo")
all_repFiles <- list.files(postDir, pattern = "repeatabilityData_from_Repeat")

# This is the line of descent string
ref_lineString <- "_->_"

############################### BODY ####################################
# I collect the fitness from my experiments
fitnessValues <- NULL
# This is the convenience of the columns we'll be keeping for popDemo
cols_popDemo <- c("simJob","simSet","Step","minFit","meanFit","maxFit")
for(thisPop in all_popDemo){
     # We create a new environment to make loading easier, as well as tidy
     tmpEnv <- new.env()
     load(paste(postDir,thisPop,sep=""), envir = tmpEnv)
     tmpData <- eval(as.name(sub("popDemographics_from_","",sub("\\.RData","",thisPop))), envir = tmpEnv)
     # Now, we simply gather the information about time step, particular replicate,
     # mean, min, max fitness throughout that time
     tmpAdd <- map_dfr(names(tmpData), function(theseResults){
                         # From "theseResults" I want to know what are the batch and job ID values
                         # This for later ease of collection
                         tmp_jobInfo <- all_jobInfo[,paste(theseResults,"_1",sep="")]
                         # From the rownames we extract the step information
                         as_tibble(tmpData[[theseResults]]) %>%
                           mutate(., Step = as.numeric(sub("Step_","",rownames(tmpData[[theseResults]]))),
                                        simJob = tmp_jobInfo["jobID"],
                                        simSet = tmp_jobInfo["setID"]) %>%
                           select(., one_of(cols_popDemo))
                    })
     # Now we return this to our fitnessValues
     fitnessValues <- rbind(fitnessValues, tmpAdd)
     rm( tmpEnv )
}


# It takes time for this, so we save it as a processed file
fileName_mutantValues <- paste(outDir,"processed_mutantValues_forPlot.RData",sep="")
if(!file.exists(fileName_mutantValues)){
     # Now, the more tricky part will be collecting the mutational information
     mutantValues <- NULL
     # This is the convenience of the columns we'll be keeping for popDemo
     cols_repInfo <- c("transition","absRank", "transitionStep",
                       "progenitor_genotypeID", "progenitor_numMuts", "progenitor_fitness",
                       "offspring_genotypeID", "offspring_numMuts", "offspring_fitness")
     for(thisMutant in all_repFiles){
          # We create a new environment to make loading easier, as well as tidy
          tmpEnv <- new.env()
          load(paste(postDir,thisMutant,sep=""), envir = tmpEnv)
          tmpData <- eval(as.name(sub("repeatabilityData_from_","",sub("\\.RData","",thisMutant))), envir = tmpEnv)
          # This is a string needed when clearing out the prefix information for the file spefic strings
          tmp_removePrefix <- "processed_runData_from_"
          # As a convenience, we gather all the names of all the sub list elements of interest
          tmp_evalGrid <- map_dfr(names(tmpData), function(x){
                                   # We first gather all the information about names in the list
                                   tmp_jobInfo <- all_jobInfo[,paste(x,"_1",sep="")]
                                   tibble("groupNames" = x,
                                             "fileNames" = names(tmpData[[x]]),
                                             "Job" = tmp_jobInfo["jobID"],
                                             "Set" = tmp_jobInfo["setID"])
                              }) %>%
                         # Now we can add in the Rep information as being the last piece of the name
                         mutate(., Rep = sapply(strsplit(fileNames,"_"),function(x){ x[length(x)] }))
          # With the grid built we can run apply more efficiently to gather the information from each
          # of the sub-items.  We use map to ensure a list is returned and nothing is converted
          tmpAdd <- mapply(function(thisGroup, thisFile, thisJob, thisSet, thisRep){
                              # So, we first find the final dominant lineage after gathering the focal info
                              tmp_subData <- tmpData[[thisGroup]][[thisFile]]
                              tmp_finalDom <- tmp_subData$final_domLineage
                              # If this is not a vector of 4 elements, where we can find a single
                              # line of descent there are problems
                              if(length(tmp_finalDom) != 4 || length(tmp_finalDom["line_ofDescent"]) != 1){
                                   stop(paste(thisGroup, thisFile, thisJob, thisSet, thisRep, sep = ", "))
                              }
                              # Now, line of descent may be different lines but these will have been pasted
                              # together with the "__and__" string
                              tmp_finalLine <- unlist(strsplit(tmp_finalDom["line_ofDescent"],"__and__"))
                              # So, we can proceed by defining the pairwise steps of the line_ofDescent
                              # This double unlisting is a caution in case it's doubly nested.
                              tmpMutants <- lapply(tmp_finalLine,function(x){ unlist(strsplit(x, "split" = ref_lineString)) })
                              # I handle it this way for instances where there are multiple lines of descent
                              tmp_pairList <- lapply(tmpMutants,function(x){
                                                       paste(x[-length(x)],
                                                             x[-1],
                                                             sep = ref_lineString)
                                                  })
                              tmp_pairWise <- unlist(tmp_pairList, use.names = FALSE)
                              # Now we gather the transitions and find the rows related to the transitions listed
                              tmp_allTransitions <- tmp_subData$transitions
                              tmp_returnTransitions <- filter(tmp_allTransitions,
                                                                 is.element(transition, tmp_pairWise)) %>%
                                                       select(., one_of(cols_repInfo)) %>%
                                                       # Now, it happens sometimes that transitions repeat,
                                                       # in which case I'll record the first
                                                       mutate(., transitionStep = as.numeric(sapply(strsplit(as.character(transitionStep),"_"),
                                                                                                         function(x){x[1]})))
                              # If any of the transitions are missing we'll need to get their information praying something is in transitions
                              if(!all(is.element(tmp_pairWise, tmp_returnTransitions$transition))){
                                   # We can try to grab the information from the transition matrix
                                   tmpMissing <- tmp_pairWise[which(!is.element(tmp_pairWise, tmp_returnTransitions$transition))]
                                   tmp_missPairs <- strsplit(tmpMissing, ref_lineString)
                                   tmpSingles <- unique(unlist(tmp_missPairs))
                                   # Now, if any singles cannot be found in the transitions we'll need to grab those
                                   # from the landscape and so we build the secondary vector
                                   tmp_missSingles <- tmpSingles[which(!is.element(tmpSingles,unlist(strsplit(tmp_allTransitions$transition,ref_lineString))))]
                                   # So, if something cannot be found as pre or post step then we gather
                                   if(length(tmp_missSingles) > 0){
                                        # This printout was used for debugging
                                        #print(paste(thisGroup, thisFile, thisJob, thisSet, thisRep, sep = ", "))

                                        # This means we should update tmpSingles
                                        tmpSingles <- setdiff(tmpSingles, tmp_missSingles)
                                        # We make a connection to the landscape file of this job
                                        tmp_landscapeFile <- list.files(paste(workDir, thisGroup, "/",sep=""),
                                                                        pattern = "Landscape(.)+\\.sqlite",
                                                                        full.names=TRUE)
                                        # There should only be one, but in case....
                                        if(length(tmp_landscapeFile) != 1){
                                             stop(paste("Not one landscape",thisGroup, thisFile, thisJob, thisSet, thisRep, sep = ", "))
                                        }
                                        tmpCon <- reset_shapeDB(func_conName = tmp_landscapeFile,
                                                                func_type = "connect")
                                        # We find the tables
                                        tmp_conTables <- RSQLite::dbListTables(tmpCon)
                                        # We look for out genotypeID, binaryString, fitnesss info, then shape this properly
                                        tmp_missSingles <- RSQLite::dbGetQuery(tmpCon,
                                                                                paste("SELECT fitness,binaryString,genotypeID FROM ",
                                                                                      tmp_conTables,
                                                                                      " WHERE genotypeID IN (",
                                                                                      paste(tmp_missSingles,collapse=","),
                                                                                      ')',
                                                                                      sep="", collapse=" UNION ")) %>%
                                                            # I assume something will have been found, it BETTER!
                                                            distinct(.) %>%
                                                            # we get the numMuts from the length of the genotypeID
                                                            mutate(., this_numMuts = sapply(strsplit(binaryString,"_"),length),
                                                                      binaryString = NULL) %>%
                                                            # we rename
                                                            rename(., thisSingle = genotypeID,
                                                                      this_fitness = fitness) %>%
                                                            # we order the columns
                                                            select(., thisSingle, this_numMuts, this_fitness)
                                        # We disconnect the DB
                                        dbDisconnect(tmpCon)
                                        # We can now add this to the
                                   }
                                   # We find the rows in the tmp_allTransitions where we can find the singles
                                   # as firsts, and then as seconds.  We also add in anything we had to scavenge from out
                                   tmpFound <- rbind(map_dfr(tmpSingles, function(thisSingle){
                                                       # This is for debugging
                                                       #print(thisSingle)
                                                       # We look for it as a first or second in that order
                                                       tmp_rowFind <- list("first" = which(grepl(paste("^",thisSingle,"_",sep=""), tmp_allTransitions$transition)),
                                                                           "second" = which(grepl(paste("_",thisSingle,"$",sep=""), tmp_allTransitions$transition)))
                                                       if(length(tmp_rowFind[["first"]]) > 0){
                                                            return( data.frame("thisSingle" = thisSingle,
                                                                                "this_numMuts" = tmp_allTransitions$progenitor_numMuts[tmp_rowFind[["first"]][1]],
                                                                                "this_fitness" = tmp_allTransitions$progenitor_fitness[tmp_rowFind[["first"]][1]]) )
                                                       } else if(length(tmp_rowFind[["second"]]) > 0){
                                                            return( data.frame("thisSingle" = thisSingle,
                                                                                "this_numMuts" = tmp_allTransitions$offspring_numMuts[tmp_rowFind[["second"]][1]],
                                                                                "this_fitness" = tmp_allTransitions$offspring_fitness[tmp_rowFind[["second"]][1]]) )
                                                       } else {
                                                            stop(paste(thisGroup, thisFile, thisJob, thisSet, thisRep, "singles:",thisSingle, sep = ", "))
                                                       }
                                                  }),
                                                  tmp_missSingles)
                                   # We now build the required elements using the appropriate rows
                                   tmpReturn <- data.frame("progenitor_genotypeID" = sapply(tmp_missPairs,function(x) x[1]) ,
                                                            "offspring_genotypeID" = sapply(tmp_missPairs,function(x) x[2]),
                                                            "transition" = tmpMissing) %>%
                                                  # By joining on the progenitor we have that info, then we do the same for offspring
                                                  left_join(.,
                                                            tmpFound,
                                                            by = c("progenitor_genotypeID" = "thisSingle")) %>%
                                                  rename(., progenitor_fitness = this_fitness,
                                                            progenitor_numMuts = this_numMuts) %>%
                                                  left_join(.,
                                                            tmpFound,
                                                            by = c("offspring_genotypeID" = "thisSingle")) %>%
                                                  rename(., offspring_fitness = this_fitness,
                                                            offspring_numMuts = this_numMuts) %>%
                                                  # We'll be missing on columns which we add here as null values
                                                  mutate(., transitionStep = as.numeric(NA),
                                                            absRank = NA)
                                   # So, we just need to weave this information into the tmp_allTransitions
                                   tmp_returnTransitions <- rbind(tmp_returnTransitions,
                                                                  tmpReturn[, colnames(tmp_returnTransitions)]) %>%
                                                            # We order this by joining in the tmp_pairWise
                                                            left_join(.,
                                                                      data.frame("transition" = tmp_pairWise,
                                                                                 "order" = unlist(lapply(tmp_pairList,function(x){ 1:length(x) }),use.names=FALSE)),
                                                                      by = "transition") %>%
                                                            arrange(., order) %>%
                                                            select(., -order) %>%
                                                            distinct(.)
                                                            #mutate(., transition = factor(transition, levels = unique(tmp_pairWise))) %>%
                                                            #arrange(., transition) %>%
                                                            #mutate(., transition = as.character(transition))


                                   # Now, for each of the NA values in transitionStep we take the mean of the above and below
                                   # this is admitedly an estimate, but it's good enough.
                                   tmpUpdate <- which(is.na(tmp_returnTransitions$transitionStep))
                                   # If we're updating the first or last row then we handle this as being mean diff
                                   # from first or last steps, we do this iteratively in case of sequentialy NA's
                                   if(length(tmpUpdate) > 0){
                                        # This existed for debugging
                                        #print(thisFile)
                                        for(thisUpdate in tmpUpdate){
                                             # As a convenience we extract the vector of indexes from next step to last step
                                             tmp_considerIndexes <- (thisUpdate+1):nrow(tmp_returnTransitions)
                                             # This is the current potential next number for mean
                                             tmp_putativeMean <- suppressWarnings(min(tmp_returnTransitions$transitionStep[tmp_considerIndexes],na.rm=TRUE))
                                             # We now update and replace first and last row elements by first and last step values
                                             tmp_returnTransitions$transitionStep[thisUpdate] <- mean(c(ifelse(thisUpdate > 1,
                                                                                                                   tmp_returnTransitions$transitionStep[thisUpdate-1],
                                                                                                                   1),
                                                                                                            ifelse(c(thisUpdate < nrow(tmp_returnTransitions) &&
                                                                                                                        tmp_putativeMean != Inf),
                                                                                                                   tmp_putativeMean,
                                                                                                                   all_parmInfo$shape_numGenerations[which(all_parmInfo$jobName == thisGroup)])))
                                        }
                                   }
                                   # Now we update to be only distinct aspects and order by transition step
                                   tmp_returnTransitions <- arrange(tmp_returnTransitions, transitionStep) %>%
                                                            distinct(.)

                              } # This closes the loop for having missing transitions
                              # Now to actually return the information
                              return( mutate(tmp_returnTransitions,
                                             Rep = thisRep,
                                             Set = thisSet,
                                             Job = thisJob,
                                             orderedTransition = 1:nrow(tmp_returnTransitions)) )

                         },
                         thisGroup = tmp_evalGrid$groupNames,
                         thisFile = tmp_evalGrid$fileNames,
                         thisJob = tmp_evalGrid$Job,
                         thisSet = tmp_evalGrid$Set,
                         thisRep = tmp_evalGrid$Rep,
                         SIMPLIFY = FALSE)
          # We join the df's into a single one and add it
          mutantValues <- rbind(mutantValues,
                                do.call(rbind,tmpAdd))
          # Tidy our environment
          rm( tmpEnv )
     }
     # Sanity checking if there are any NA of Inf values in the mutant frame
     if(any(apply(mutantValues,2,function(x){ any(x == Inf) })) ||
        any(apply(mutantValues,2,function(x){ any(is.na(unlist(x))) }))){
          stop("There are either NA or Inf values in the mutantValues frame")
     }
     # We now save the processed information as this takes time
     save(mutantValues, file = fileName_mutantValues)
} else {
     load( fileName_mutantValues )
}

# Ok, for plotting we'll gather the bottleneck size as an element since this is something
# I want to compare across the different plots
# We create a grouping based on sequence so that we effectively "smooth" our lines
stepSmooth = 20
plotFitness <- mutate(fitnessValues, 
                      jobName = name_batchString(ref_batchString,simSet,simJob, func_sepString="_"),
                      stepGroup = floor(Step/stepSmooth) * stepSmooth) %>%
                # then we join in the extra information
                left_join(.,
                          mutate(all_parmInfo[,c("jobName","distFactor","shape_simModel")], 
                                 colString = paste(distFactor, shape_simModel, sep="_")),
                          by = "jobName") %>%
                # This gets the mean value for the elements of the stepGroup
                group_by(., stepGroup, colString) %>%
                        summarise_at(., .vars = vars(c("meanFit","minFit","maxFit")), .funs = mean) %>%
                        ungroup(.) %>%
                mutate(., distFactor = sapply(strsplit(colString,"_"),function(x){ x[1]}),
                          shape_simModel = sapply(strsplit(colString,"_"),function(x){ x[2]}))
                
# We'll define colours by the shape_simModel and bottleneck for shade
plotCols <- expand.grid(distFactor = unique(plotFitness$distFactor),
                        shape_simModel = unique(plotFitness$shape_simModel)) %>%
                # Then we define our colours and plotting group strings
                mutate(., col = colTransparency(c("blue","red"),
                                                func_scaleSaturation = seq(0.2,0.9,length.out = n_distinct(distFactor)),
                                                func_scaleValue = seq(0.8,0.3,length.out = n_distinct(distFactor))) %>%
                                as.vector(.),
                        colString = paste(distFactor,shape_simModel,sep="_"))
# named plotting colours 
plot_namedCols <- setNames(plotCols$col,
                           paste(plotCols$distFactor,
                                 plotCols$shape_simModel,
                                 sep="_"))

# We now create the sub mean and polygon data frames
plot_meanFrame <- distinct(plotFitness[,c("stepGroup","meanFit","colString")])
                        
# We join the inverse of the second values
plot_polygonFrame <- rbind(rename(plotFitness[,c("stepGroup","minFit","colString")], polyFit = minFit),
                           rename(plotFitness[nrow(plotFitness) : 1,c("stepGroup","maxFit","colString")], polyFit = maxFit))
plotCex <- c("baseLeg" = 22,
             "baseMain" = 36,
             "lwdMain" = 2,
             "lwdAxis" = 1)
plotLegend <- ggplot(plotCols, 
                     aes(x = as.numeric(as.factor(distFactor)), 
                         y = as.numeric(as.factor(shape_simModel)))) +
                theme_classic(base_size = plotCex["baseLeg"]) +
                geom_tile(aes(fill = factor(col, levels = unique(col)),
                              colour = factor(1)),
                          width = 0.925, height = 0.925) +
                scale_fill_manual(name = "", values = plotCols$col) +
                scale_colour_manual(name = "", values = "black") +
                theme(legend.position = "none",
                      axis.text.y = element_text(margin = margin(r = plotCex["baseLeg"]/2)),
                      axis.text.x = element_text(margin = margin(b = plotCex["baseLeg"]/2)),
                      axis.line = element_blank(),
                      axis.ticks = element_blank(),
                      axis.title = element_blank(),
                      plot.margin = margin(1,1,1,1),
                      line = element_line(lineend = "round")) + 
                scale_x_continuous(breaks = 1:n_distinct(plotCols$distFactor),
                                 labels = unique(plotCols$distFactor),
                                 position = "top",
                                 expand = c(0,0)) +
                scale_y_continuous(breaks = 1:n_distinct(plotCols$shape_simModel),
                                 labels = levels(plotCols$shape_simModel),
                                 position = "left",
                                 expand = c(0,0))

# We make the main plot here
# I need to work on the polygon but this is a good start.
#text = element_text(family = windowsFonts()$serif),   ## THIS didn't help as expected
plotMain <- ggplot(plot_meanFrame, aes(x = stepGroup, y = meanFit, 
                                 group = colString,
                                 colour = colString)) +
                theme_classic(base_size = plotCex["baseMain"]) +
                theme(legend.position = "none",
                      axis.line = element_line(size = plotCex["lwdAxis"],
                                               lineend = "square"),
                      axis.ticks = element_line(size = plotCex["lwdAxis"],
                                                lineend = "butt"),
                      axis.text.y.left = element_text(margin = margin(r = plotCex["baseMain"]/4,
                                                                      l = plotCex["baseMain"]/2)),
                      axis.text.x.bottom = element_text(margin = margin(b = plotCex["baseMain"]/2,
                                                                        t = plotCex["baseMain"]/4)),
                      plot.margin = margin(plotCex["baseMain"]/8,
                                           plotCex["baseMain"]/8,
                                           plotCex["baseMain"]/8,
                                           plotCex["baseMain"]/8)) +
                geom_polygon(aes(y = polyFit, fill = colString, colour = NA), 
                             data = plot_polygonFrame) +
                geom_line(size = plotCex["lwdMain"]) +
                scale_x_continuous(name = "Generation Number") +
                scale_y_continuous(name = "Mean community fitness") +
                scale_fill_manual(values = setNames(alpha(plot_namedCols, 0.2),names(plot_namedCols))) +
                scale_colour_manual(values = plot_namedCols)
                
    

fullPlot <- ggdraw() + 
                draw_plot(plotMain, x = 0 , y = 0,
                          width = 1 , height = 01) +
                draw_plot(plotLegend, x = 0.2, y = 0.95,
                          width = 0.25, height = 0.25,
                          hjust = 0.15, vjust = 1)


ggsave(paste(outDir,"plot_rSHAPE_validation_fitnessIncrease.png",sep=""),
       plot = fullPlot, dpi = 120, width = 320, height = 180, units = "mm")

#
# Check that the rate of difference depends on the shape_simModel (this because the diff b/ lines
# seems affected by the model type interacting with distFactor.


############################### PLOT OF NUMBER OF STEPS AND STEP SIZE ###################################
# We explicitly define the step size and a string for grouping work, this will involve
# joining in the information about distFactor and shape_simModel
plot_stepsData <- mutate(mutantValues,
                         stepSize = offspring_fitness - progenitor_fitness,
                         jobName = name_batchString(funcBase = ref_batchString,
                                                    func_setID = Set,
                                                    func_jobID = Job,
                                                    func_sepString = "_")) %>%
                left_join(.,
                          all_parmInfo[,c("jobName","distFactor","shape_simModel")],
                          by = "jobName") %>%
                mutate(., colString = paste(distFactor, shape_simModel, sep="_")) %>%
                # Now, for each set, job, rep (ie: unique run) I will find the 
                # ordered stepwise transitions and then returned the numeric
                left_join(.,
                           group_by(., Set, Job, Rep) %>%
                            mutate(., orderedTransition = as.numeric(factor(transitionStep, levels = sort(unique(transitionStep))))) %>%
                            ungroup(.) %>%
                            select(., Set, Job, Rep, orderedTransition),
                          by = c("Set","Job","Rep"))
# This is a mixed effects model for stepSize
model_stepSize <- lmer(stepSize ~ orderedTransition * colString + (1|jobName), data = plot_stepsData)

main_stepPlot <- ggplot(data = plot_stepsData, 
                        aes(x = orderedTransition, y = stepSize,
                            colour = colString, group = colString)) +
                theme_classic(base_size = plotCex["baseMain"]) +
                theme(legend.position = "none",
                      axis.line = element_line(size = plotCex["lwdAxis"],
                                               lineend = "square"),
                      axis.ticks = element_line(size = plotCex["lwdAxis"],
                                                lineend = "butt"),
                      axis.text.y.left = element_text(margin = margin(r = plotCex["baseMain"]/4,
                                                                      l = plotCex["baseMain"]/2)),
                      axis.text.x.bottom = element_text(margin = margin(b = plotCex["baseMain"]/2,
                                                                        t = plotCex["baseMain"]/4)),
                      plot.margin = margin(plotCex["baseMain"]/8,
                                           plotCex["baseMain"]/8,
                                           plotCex["baseMain"]/8,
                                           plotCex["baseMain"]/8)) +
                geom_smooth(mapping = aes(x = orderedTransition, y = stepSize,
                                          colour = colString, group = colString,
                                          fill = colString),
                            method = "glm", 
                            formula = str(model_stepSize@call), 
                            data = plot_stepsData) +
                scale_fill_manual(values = setNames(alpha(plot_namedCols, 0.2),names(plot_namedCols))) +
                scale_colour_manual(values = plot_namedCols) +
                scale_x_continuous(name = "Evolutionary step",
                                   breaks = seq(1,20,by=2)) +
                scale_y_continuous(name = expression("Step size (W"[o]*" - W"[p]*")"))
    


fullPlot <- ggdraw() + 
    draw_plot(main_stepPlot, x = 0 , y = 0,
              width = 1 , height = 1) +
    draw_plot(plotLegend, x = 0.95, y = 0.25,
              width = 0.25, height = 0.25,
              hjust = 1, vjust = 0)

ggsave(paste(outDir,"plot_rSHAPE_validation_stepSize.png",sep=""), 
      plot = fullPlot, dpi = 240, 
      width = 320, height = 180, units = "mm")

################### PLOT - closseness to optima ############################
# Using the prevbiously created set with sequential steps we can calculate the time between steps

# Next plot will be to show how close steps are the a local optima


############################ END BODY ###################################
q(save="no")
