source("extras/DatabaseDetails.R")

database <- databases[[3]]
# Diclofenac
pathToApp <- "d:/temp/ce_11038"
CohortExplorer::createCohortExplorerApp(
  connectionDetails = database$connectionDetails,
  cdmDatabaseSchema = database$cdmDatabaseSchema,
  cohortDatabaseSchema = database$cohortDatabaseSchema,
  cohortTable = database$cohortTable,
  cohortDefinitionId = 11038,
  databaseId = database$databaseId,
  exportFolder = pathToApp
)
usethis::create_project(pathToApp)

# Celecoxib
pathToApp <- "d:/temp/ce_11037"
CohortExplorer::createCohortExplorerApp(
  connectionDetails = database$connectionDetails,
  cdmDatabaseSchema = database$cdmDatabaseSchema,
  cohortDatabaseSchema = database$cohortDatabaseSchema,
  cohortTable = database$cohortTable,
  cohortDefinitionId = 11037,
  databaseId = database$databaseId,
  exportFolder = pathToApp
)
usethis::create_project(pathToApp)
