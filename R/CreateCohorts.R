# Copyright 2023 Observational Health Data Sciences and Informatics
#
# This file is part of InteractionsEval
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# library(dplyr)

createCohorts <- function(connectionDetails,
                          cdmDatabaseSchema,
                          cohortDatabaseSchema,
                          cohortTable,
                          outputFolder) {

  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection))

  cohortTableNames <- CohortGenerator::getCohortTableNames(cohortTable)
  CohortGenerator::createCohortTables(
    connection = connection,
    cohortDatabaseSchema = cohortDatabaseSchema,
    cohortTableNames = cohortTableNames
  )
  
  # Exposure cohorts (from ATLAS) ----------------------------------------------
  cohortDefinitionSet1 <- CohortGenerator::getCohortDefinitionSet(
    settingsFileName = "Cohorts.csv",
    packageName = "InteractionsEval"
  )
  CohortGenerator::generateCohortSet(
    connection = connection,
    cdmDatabaseSchema = cdmDatabaseSchema,
    cohortDatabaseSchema = cohortDatabaseSchema,
    cohortTableNames = cohortTableNames,
    cohortDefinitionSet = cohortDefinitionSet1
  )
  
  # Outcome of interest cohort (from the Phenotype Library) --------------------
  # 77 = GI Bleed
  cohortDefinitionSet2 <- PhenotypeLibrary::getPlCohortDefinitionSet(77)
  CohortGenerator::generateCohortSet(
    connection = connection,
    cdmDatabaseSchema = cdmDatabaseSchema,
    cohortDatabaseSchema = cohortDatabaseSchema,
    cohortTableNames = cohortTableNames,
    cohortDefinitionSet = cohortDefinitionSet2
  )
  
  # Negative controls ----------------------------------------------------------
  negativeControls <- readr::read_csv(
    file = system.file("NegativeControls.csv", package = "InteractionsEval"),
    show_col_types = FALSE
  ) %>%
    select(
      cohortId = "outcomeId",
      cohortName = "outcomeName",
      outcomeConceptId = "outcomeId"
    )
  CohortGenerator::generateNegativeControlOutcomeCohorts(
    connection = connection,
    cdmDatabaseSchema = cdmDatabaseSchema,
    cohortDatabaseSchema = cohortDatabaseSchema,
    cohortTable = cohortTable,
    occurrenceType = "first",
    negativeControlOutcomeCohortSet = negativeControls
  )
  
  # Count cohorts --------------------------------------------------------------
  cohortCounts <- CohortGenerator::getCohortCounts(
    connection = connection,
    cohortDatabaseSchema = cohortDatabaseSchema,
    cohortTable = cohortTable
  )
  cohortCounts <- cohortCounts %>%
    inner_join(
      bind_rows(
        cohortDefinitionSet1 %>%
          select("cohortId", "cohortName"),
        cohortDefinitionSet2 %>%
          select("cohortId", "cohortName"),
        negativeControls %>%
          select("cohortId", "cohortName")
      ),
      by = "cohortId"
    )
  readr::write_csv(cohortCounts, file.path(outputFolder, "CohortCounts.csv"))
}
