# Code for running database diagnostics to select the databases to include in this study

folder <- "d:/temp/dbDiagnostics"
dir.create(folder)

# Create settings --------------------------------------------------------------
ROhdsiWebApi::authorizeWebApi(
  baseUrl = Sys.getenv("baseUrl"),
  authMethod = "windows")

targetCohortDefinition <- ROhdsiWebApi::getCohortDefinition(
  cohortId = 11038, 
  baseUrl = Sys.getenv("baseUrl")
)
targetConceptSet <- targetCohortDefinition$expression$ConceptSets[[1]]
targetConceptSet$name
targetConceptIds <- ROhdsiWebApi::resolveConceptSet(
  conceptSetDefinition = targetConceptSet, 
  baseUrl = Sys.getenv("baseUrl")
)

comparatorCohortDefinition <- ROhdsiWebApi::getCohortDefinition(
  cohortId = 11037, 
  baseUrl = Sys.getenv("baseUrl")
)
comparatorConceptSet <- comparatorCohortDefinition$expression$ConceptSets[[1]]
comparatorConceptSet$name
comparatorConceptIds <- ROhdsiWebApi::resolveConceptSet(
  conceptSetDefinition = comparatorConceptSet, 
  baseUrl = Sys.getenv("baseUrl")
)

outcomeCohortDefinition <-  jsonlite::fromJSON(PhenotypeLibrary::getPlCohortDefinitionSet(77)$json, simplifyVector = FALSE)
outcomeConceptSet <- outcomeCohortDefinition$ConceptSets[[2]]
outcomeConceptSet$name
outcomeConceptIds <- ROhdsiWebApi::resolveConceptSet(
  conceptSetDefinition = outcomeConceptSet, 
  baseUrl = Sys.getenv("baseUrl"))

# indicationCohortDefinition <- ROhdsiWebApi::getCohortDefinition(
#   cohortId = 5907, 
#   baseUrl = Sys.getenv("baseUrl")
# )
# indicationConceptSet <- indicationCohortDefinition$expression$ConceptSets[[1]]
# indicationConceptSet$name
# indicationConceptIds <- ROhdsiWebApi::resolveConceptSet(
#   conceptSetDefinition = indicationConceptSet,
#   baseUrl = Sys.getenv("baseUrl")
# )

analysisSettings1 <- DbDiagnostics::createDataDiagnosticsSettings(
  analysisId = 1,
  analysisName = "A1",
  requiredDurationDays = 365,
  requiredDomains = c("condition","drug"),
  targetName = targetConceptSet$name,
  targetConceptIds = targetConceptIds,
  comparatorName = comparatorConceptSet$name,
  comparatorConceptIds = comparatorConceptIds,
  outcomeName = outcomeConceptSet$name,
  outcomeConceptIds = outcomeConceptIds,
  desiredVisits = c("IP")
)

settingsList <- list(analysisSettings1)
ParallelLogger::saveSettingsToJson(settingsList, file.path(folder, "dbDiagnosticsSettingsList.json"))

# Run database diagnostics -----------------------------------------------------
settingsList <- ParallelLogger::loadSettingsFromJson(file.path(folder, "dbDiagnosticsSettingsList.json"))

dbProfileConnectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  user = keyring::key_get("strategusUser"),
  password = keyring::key_get("strategusPassword"),
  connectionString = keyring::key_get("strategusConnectionString")
)
dbDiagnosticResults <- DbDiagnostics::executeDbDiagnostics(
  connectionDetails = dbProfileConnectionDetails,
  resultsDatabaseSchema = "dp_temp",
  resultsTableName = "dp_achilles_results_augmented",
  outputFolder = folder,
  dataDiagnosticsSettingsList = settingsList
)

library(dplyr)
library(tidyr)
dbDiagnosticSummary <- DbDiagnostics::createDataDiagnosticsSummary(dataDiagnosticResults)
CohortGenerator::writeCsv(dbDiagnosticSummary, file.path(folder,"data_diagnostics_summary.csv"))
