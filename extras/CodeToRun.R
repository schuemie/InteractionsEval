library(InteractionsEval)

maxCores <- max(24, parallel::detectCores())

source("extras/DatabaseDetails.R")

for (i in 1:length(databases)) {
  database <- databases[[i]]
  message(sprintf("***** Running on %s *****", database$databaseId))
  execute(
    connectionDetails = database$connectionDetails,
    cdmDatabaseSchema = database$cdmDatabaseSchema,
    cohortDatabaseSchema = database$cohortDatabaseSchema,
    cohortTable = database$cohortTable,
    verifyDependencies = TRUE,
    outputFolder = database$outputFolder,
    databaseId = database$databaseId,
    createCohorts = TRUE,
    runCohortMethod = TRUE,
    maxCores = maxCores
  )
}


for (i in 4:length(databases)) {
  database <- databases[[i]]
  message(sprintf("***** Computing diagnostics on %s *****", database$databaseId))
  computeDiagnostics(database$outputFolder)
}

dataFolders <- sapply(databases, function(x) x$outputFolder)
dataFolders <- dataFolders[dataFolders != "d:/interactionEval/OptumEhr"]
dataSet <- createAnalysisDataSet(dataFolders)
saveRDS(dataSet, file.path(rootFolder, "DataSet.rds"))
dataSet <- readRDS(file.path(rootFolder, "DataSet.rds"))

evaluate2dPadePoisson(dataSet = dataSet, folder = rootFolder)
evaluate2dPadeevaluate2dPadeCox(dataSet = dataSet, folder = rootFolder)