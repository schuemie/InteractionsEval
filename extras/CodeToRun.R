library(InteractionsEval)

options(andromedaTempFolder = "d:/andromedaTemp")
maxCores <- parallel::detectCores()


connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "redshift",
                                                                connectionString = keyring::key_get("redShiftConnectionStringMdcd"),
                                                                user = keyring::key_get("redShiftUserName"),
                                                                password = keyring::key_get("redShiftPassword"))
cdmDatabaseSchema <- "cdm"
cohortDatabaseSchema <- "scratch_mschuemi2"
cohortTable <- "interactions_eval_mdcd"
outputFolder <- "d:/InteractionsEval_MDCD"
