# Code for setting the connection details and schemas for the various databases,
# as well as some folders in the local file system.
options(andromedaTempFolder = "d:/andromedaTemp")
databases <- list()
rootFolder <- "d:/interactionEval"

# IBM_CCAE-20221023 ------------------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "CCAE",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaCcae"),
    user = keyring::key_get("redShiftUserName"),
    password = keyring::key_get("redShiftPassword")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_truven_ccae_v2324"
)

# OPTUM_Extended_DOD-20221111 --------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "OptumDod",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaOptumDod"),
    user = keyring::key_get("redShiftUserName"),
    password = keyring::key_get("redShiftPassword")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_optum_extended_dod_v2323"
)

# Optum_EHR-20221205 -----------------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "OptumEhr",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaOptumEhr"),
    user = keyring::key_get("temp_user"),
    password = keyring::key_get("temp_password")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_optum_ehr_v2247"
)

# PharMetrics-20220831 ---------------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "PharMetrics",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaPharmetrics"),
    user = keyring::key_get("redShiftUserName"),
    password = keyring::key_get("redShiftPassword")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_iqvia_pharmetrics_plus_v2286"
)

# IBM_MDCR-20221021 ------------------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "MDCR",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaMdcr"),
    user = keyring::key_get("redShiftUserName"),
    password = keyring::key_get("redShiftPassword")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_truven_mdcr_v2322"
)

# IBM_MDCD-20220802 ------------------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "MDCD",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaMdcd"),
    user = keyring::key_get("redShiftUserName"),
    password = keyring::key_get("redShiftPassword")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_truven_mdcd_v2321"
)

# JMDC-20221106 ----------------------------------------------------------------
databases[[length(databases) + 1]] <- list(
  databaseId = "JMDC",
  connectionDetails = DatabaseConnector::createConnectionDetails(
    dbms = "redshift",
    connectionString = keyring::key_get("redShiftConnectionStringOhdaJmdc"),
    user = keyring::key_get("redShiftUserName"),
    password = keyring::key_get("redShiftPassword")
  ),
  cohortDatabaseSchema = "scratch_mschuemi",
  cdmDatabaseSchema = "cdm_jmdc_v2325"
)

# Set cohort table and folders -------------------------------------------------
for (i in seq_along(databases)) {
  databases[[i]]$cohortTable <- sprintf(
    "interaction_eval_cohort_%s", 
    databases[[i]]$databaseId
  )
  databases[[i]]$outputFolder <- file.path(
    rootFolder,
    databases[[i]]$databaseId
  )
}

