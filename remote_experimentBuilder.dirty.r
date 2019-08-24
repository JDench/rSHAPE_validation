# We load rSHAPE
library(rSHAPE)
# This is my working directory and the template location
expDir <- "R:/Research_Data/SimulationRuns/rSHAPE_validation/"
templateDir <- paste(expDir,"SHAPE_templates/",sep="")
# This is where R is installed
filePath_R <- "\\'C:/Program Files/R/R-3.5.3/bin/R\\'"

# We setup the rSHAPE space
defineSHAPE()
# We build the experiment passing the requirements to run locally
shapeExperiment(func_filepath_toDesign = paste(templateDir,"SHAPE_experimentalDesign.v.1.r",sep=""), 
                func_templateDir = templateDir, 
                func_maxGrouped_perShell = 4, 
                func_filePath_R = "\\'C:/Program Files/R/R-3.5.3/bin/R\\'", 
                func_baseCall = "CMD BATCH --vanilla --slave", 
                func_processingCores = 4, 
                func_suppressOld_summaryFiles = FALSE)





shapeExperiment(func_filepath_toDesign = paste(templateDir,"SHAPE_experimentalDesign.v.1.r",sep=""), 
                func_templateDir = templateDir, 
                func_maxGrouped_perShell = 4, 
                func_filePath_R = "\\'C:/Program Files/R/R-3.5.3/bin/R\\'", 
                func_baseCall = "CMD BATCH --vanilla --slave", 
                func_processingCores = 4, 
                func_suppressOld_summaryFiles = FALSE)
				
				
shapeExperiment(func_filepath_toDesign = "/home/jdenc017/project/jdenc017/shapeValidation/SHAPE_templates/SHAPE_experimentalDesign.v.1.r",
func_templateDir = "/home/jdenc017/project/jdenc017/shapeValidation/SHAPE_templates/",
func_maxGrouped_perShell = 2,
func_filePath_R = "/home/jdenc017/project/jdenc017/R-3.4.1/bin/R",
func_baseCall = "CMD BATCH" ,
func_rArgs = '\"--args shape_thisRep=1 shape_outDir=\'fake_serverPath/fakeDir/\'\"',
func_remoteLocation = "$SLURM_TMPDIR", 
func_submitArgs = c(number_ofCores = "-c 1",
					memory = "--mem=2048", 
					jobName = "--job-name=fakeJob", 
					wallTime = "--time 14-00:00:00",
					fileOut = "--output=fakeOut"),
func_processingCores = 1,
func_suppressOld_summaryFiles = FALSE)
