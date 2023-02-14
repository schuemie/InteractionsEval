library(InteractionsEval)

maxCores <- max(24, parallel::detectCores())

source("extras/DatabaseDetails.R")

for (i in 2:length(databases)) {
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
