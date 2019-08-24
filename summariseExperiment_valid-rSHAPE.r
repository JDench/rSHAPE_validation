library(rSHAPE)
defineSHAPE(shape_save_batchSet = 1,
 shape_save_batchJob = 1,
 shape_sepString = "_",
 shape_sepLines = "__and__",
 shape_workDir = "R:/Research_Data/SimulationRuns/rSHAPE_validation/",
 shape_postDir = "R:/Research_Data/SimulationRuns/rSHAPE_validation/postAnal/",
 shape_save_batchBase = "valid_rSHAPE",
 shape_string_lineDescent = "_->_")

# Now I miss wrote my name, I included in the batch the unacceptable character
# of the sepString.... so we change this
allFiles <- list.files(getOption("shape_workDir"), pattern = getOption("shape_save_batchBase"), full.names = TRUE, recursive = TRUE)
# We remove any sqlite files
allFiles <- allFiles[which(!grepl("\\.sqlite",allFiles))]
# we now rename the files
file.rename(from = allFiles, to = gsub(getOption("shape_save_batchBase"),sub(getOption("shape_sepString"),"-",getOption("shape_save_batchBase")),allFiles))

# I now reset the batch base
defineSHAPE(shape_save_batchSet = 1,
            shape_save_batchJob = 1,
            shape_sepString = "_",
            shape_sepLines = "__and__",
            shape_workDir = "R:/Research_Data/SimulationRuns/rSHAPE_validation/",
            shape_postDir = "R:/Research_Data/SimulationRuns/rSHAPE_validation/postAnal/",
            shape_save_batchBase = "valid-rSHAPE",
            shape_string_lineDescent = "_->_")
# Now we can summarise the experiment
summariseExperiment(func_numCores = 4, func_suppressOld = FALSE)


q(save="no")
