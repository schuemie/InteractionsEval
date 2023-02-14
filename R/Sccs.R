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



runSccs <- function(connectionDetails,
                    cdmDatabaseSchema,
                    cohortDatabaseSchema,
                    cohortTable,
                    outputFolder,
                    maxCores) {
  folder <- file.path(outputFolder, "sccs")
  diclofenac <- 1124300
  ppis <- c(911735, 929887, 923645, 904453, 948078, 19039926)

  if (!file.exists(folder)) {
    dir.create(folder)
  }

  sccsData <- SelfControlledCaseSeries::getDbSccsData(
    connectionDetails = connectionDetails,
    cdmDatabaseSchema = cdmDatabaseSchema,
    outcomeDatabaseSchema = cohortDatabaseSchema,
    outcomeTable = cohortTable,
    outcomeIds = 77,
    exposureDatabaseSchema = cdmDatabaseSchema,
    exposureTable = "drug_era",
    exposureIds = c(diclofenac, ppis),
    studyEndDate = "20151231"
  )
  SelfControlledCaseSeries::saveSccsData(sccsData, file.path(folder, "data1.zip"))
  sccsData <- SelfControlledCaseSeries::loadSccsData(file.path(folder, "data1.zip"))
  studyPop <- SelfControlledCaseSeries::createStudyPopulation(
    sccsData = sccsData,
    outcomeId = 77,
    firstOutcomeOnly = FALSE,
    naivePeriod = 180
  )
  covarDiclofenac <- SelfControlledCaseSeries::createEraCovariateSettings(
    label = "Exposure of interest",
    includeEraIds = diclofenac,
    start = 0,
    end = 0,
    endAnchor = "era end"
  )

  covarPreDiclofenac <- SelfControlledCaseSeries::createEraCovariateSettings(
    label = "Pre-exposure",
    includeEraIds = diclofenac,
    start = -60,
    end = -1,
    endAnchor = "era start"
  )

  seasonalityCovariateSettings <- SelfControlledCaseSeries::createSeasonalityCovariateSettings()

  calendarTimeSettings <- SelfControlledCaseSeries::createCalendarTimeCovariateSettings()

  covarPpis <- SelfControlledCaseSeries::createEraCovariateSettings(
    label = "PPIs",
    includeEraIds = ppis,
    stratifyById = FALSE,
    start = 1,
    end = 0,
    endAnchor = "era end"
  )

  sccsIntervalData <- SelfControlledCaseSeries::createSccsIntervalData(
    studyPopulation = studyPop,
    sccsData = sccsData,
    eraCovariateSettings = list(
      covarDiclofenac,
      covarPreDiclofenac,
      covarPpis
    ),
    seasonalityCovariateSettings = seasonalityCovariateSettings,
    calendarTimeCovariateSettings = calendarTimeSettings
  )
  SelfControlledCaseSeries::saveSccsIntervalData(sccsIntervalData, file.path(folder, "sccsIntervalData.zip"))
  sccsIntervalData <- SelfControlledCaseSeries::loadSccsIntervalData(file.path(folder, "sccsIntervalData.zip"))

  model <- SelfControlledCaseSeries::fitSccsModel(sccsIntervalData, control = createControl(
    cvType = "auto",
    selectorType = "byPid",
    startingVariance = 0.1,
    noiseLevel = "quiet",
    threads = 20
  ))
  saveRDS(model, file.path(folder, "ppiModel.rds"))
  
  
  sccsIntervalData2 <- SelfControlledCaseSeries::createSccsIntervalData(
    studyPopulation = studyPop,
    sccsData = sccsData,
    eraCovariateSettings = list(
      covarDiclofenac,
      covarPreDiclofenac,
      covarPpis
    ),
    eventDependentObservation = TRUE
  )
  SelfControlledCaseSeries::saveSccsIntervalData(sccsIntervalData2, file.path(folder, "sccsIntervalData2.zip"))
  sccsIntervalData2 <- SelfControlledCaseSeries::loadSccsIntervalData(file.path(folder, "sccsIntervalData2.zip"))
  
  model <- SelfControlledCaseSeries::fitSccsModel(sccsIntervalData2, control = createControl(
    cvType = "auto",
    selectorType = "byPid",
    startingVariance = 0.1,
    noiseLevel = "quiet",
    threads = 20
  ))
  saveRDS(model, file.path(folder, "ppiModel2.rds"))
  
}
