# Copyright 2022 Observational Health Data Sciences and Informatics
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

runCohortMethod <- function(connectionDetails,
                            cdmDatabaseSchema,
                            cohortDatabaseSchema,
                            cohortTable,
                            outputFolder,
                            maxCores) {
  # Create settings ------------------------------------------------------------
  outcomeOfInterest <- CohortMethod::createOutcome(
    outcomeId = 69,
    outcomeOfInterest = TRUE
  )
  negativeControls <- readr::read_csv(
    file = system.file("NegativeControls.csv", package = "InteractionsEval"),
    show_col_types = FALSE
  )
  negativeControlOutcomes <- lapply(
    negativeControls$outcomeId,
    function(outcomeId) {
      CohortMethod::createOutcome(
        outcomeId = outcomeId,
        outcomeOfInterest = FALSE,
        trueEffectSize = 1
      )
    }
  )

  tcos <- CohortMethod::createTargetComparatorOutcomes(
    targetId = 1,
    comparatorId = 2,
    outcomes = append(
      list(outcomeOfInterest),
      negativeControlOutcomes
    )
  )
  targetComparatorOutcomesList <- list(tcos)

  conceptIdsToExclude <- c(1335471, 1340128, 1341927, 1363749, 1308216, 1310756, 1373225, 1331235, 1334456, 1342439, 1395058, 974166, 978555, 907013)

  covarSettings <- FeatureExtraction::createDefaultCovariateSettings(
    excludedCovariateConceptIds = conceptIdsToExclude,
    addDescendantsToExclude = TRUE
  )

  getDbCmDataArgs <- CohortMethod::createGetDbCohortMethodDataArgs(
    washoutPeriod = 365,
    restrictToCommonPeriod = TRUE,
    firstExposureOnly = TRUE,
    maxCohortSize = 100000,
    covariateSettings = covarSettings
  )

  createStudyPopArgsAcuteTar <- CohortMethod::createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = TRUE,
    removeDuplicateSubjects = "remove all",
    minDaysAtRisk = 1,
    riskWindowStart = 0,
    startAnchor = "cohort start",
    riskWindowEnd = 7,
    endAnchor = "cohort start"
  )

  createPsArgs <- CohortMethod::createCreatePsArgs(
    maxCohortSizeForFitting = 100000,
    control = Cyclops::createControl(
      cvType = "auto",
      startingVariance = 0.01,
      tolerance = 1E-5,
      noiseLevel = "quiet",
      cvRepetitions = 1
    )
  )

  stratifyByPsArgs <- CohortMethod::createStratifyByPsArgs(numberOfStrata = 10)

  computeSharedCovBalArgs <- CohortMethod::createComputeCovariateBalanceArgs()

  fitOutcomeModelArgsLr <- CohortMethod::createFitOutcomeModelArgs(
    modelType = "logistic",
    stratified = TRUE
  )

  cmAnalysis1 <- CohortMethod::createCmAnalysis(
    analysisId = 1,
    description = "Logistic regression, TAR = 0-7",
    getDbCohortMethodDataArgs = getDbCmDataArgs,
    createStudyPopArgs = createStudyPopArgsAcuteTar,
    createPsArgs = createPsArgs,
    stratifyByPsArgs = stratifyByPsArgs,
    computeSharedCovariateBalanceArgs = computeSharedCovBalArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgsLr
  )

  createStudyPopArgsOnTreatment <- CohortMethod::createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = TRUE,
    removeDuplicateSubjects = "remove all",
    minDaysAtRisk = 1,
    riskWindowStart = 0,
    startAnchor = "cohort start",
    riskWindowEnd = 0,
    endAnchor = "cohort end"
  )

  fitOutcomeModelArgsCox <- CohortMethod::createFitOutcomeModelArgs(
    modelType = "cox",
    stratified = TRUE
  )

  cmAnalysis2 <- CohortMethod::createCmAnalysis(
    analysisId = 2,
    description = "Cox regression, TAR = on treatment",
    getDbCohortMethodDataArgs = getDbCmDataArgs,
    createStudyPopArgs = createStudyPopArgsOnTreatment,
    createPsArgs = createPsArgs,
    stratifyByPsArgs = stratifyByPsArgs,
    computeSharedCovariateBalanceArgs = computeSharedCovBalArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgsCox
  )


  black <- 8516004

  fitOutcomeModelArgsLrInteraction <- CohortMethod::createFitOutcomeModelArgs(
    modelType = "logistic",
    stratified = TRUE,
    interactionCovariateIds = black
  )

  fitOutcomeModelArgsCoxInteraction <- CohortMethod::createFitOutcomeModelArgs(
    modelType = "cox",
    stratified = TRUE,
    interactionCovariateIds = black
  )

  cmAnalysis3 <- CohortMethod::createCmAnalysis(
    analysisId = 3,
    description = "Logistic regression, TAR = 0-7, interaction",
    getDbCohortMethodDataArgs = getDbCmDataArgs,
    createStudyPopArgs = createStudyPopArgsAcuteTar,
    createPsArgs = createPsArgs,
    stratifyByPsArgs = stratifyByPsArgs,
    computeSharedCovariateBalanceArgs = computeSharedCovBalArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgsLrInteraction
  )

  cmAnalysis4 <- CohortMethod::createCmAnalysis(
    analysisId = 4,
    description = "Cox regression, TAR = on treatment, interaction",
    getDbCohortMethodDataArgs = getDbCmDataArgs,
    createStudyPopArgs = createStudyPopArgsOnTreatment,
    createPsArgs = createPsArgs,
    stratifyByPsArgs = stratifyByPsArgs,
    computeSharedCovariateBalanceArgs = computeSharedCovBalArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgsCoxInteraction
  )
  cmAnalysisList <- list(cmAnalysis1, cmAnalysis2, cmAnalysis3, cmAnalysis4)

  CohortMethod::saveCmAnalysisList(cmAnalysisList, file.path(outputFolder, "cmAnalysisList.json"))
  CohortMethod::saveTargetComparatorOutcomesList(targetComparatorOutcomesList, file.path(outputFolder, "targetComparatorOutcomesList.json"))

  # Run analyses -----------------------------------------------------------------
  cmAnalysisList <- CohortMethod::loadCmAnalysisList(file.path(outputFolder, "cmAnalysisList.json"))
  targetComparatorOutcomesList <- CohortMethod::loadTargetComparatorOutcomesList(file.path(outputFolder, "targetComparatorOutcomesList.json"))
  multiThreadingSettings <- CohortMethod::createDefaultMultiThreadingSettings(maxCores)

  result <- CohortMethod::runCmAnalyses(
    connectionDetails = connectionDetails,
    cdmDatabaseSchema = cdmDatabaseSchema,
    exposureDatabaseSchema = cohortDatabaseSchema,
    exposureTable = cohortTable,
    outcomeDatabaseSchema = cohortDatabaseSchema,
    outcomeTable = cohortTable,
    outputFolder = outputFolder,
    cmAnalysisList = cmAnalysisList,
    targetComparatorOutcomesList = targetComparatorOutcomesList,
    multiThreadingSettings = multiThreadingSettings
  )

  resultsSummary <- CohortMethod::getResultsSummary(outputFolder)
  readr::write_csv(resultsSummary, file.path(outputFolder, "ResultsSummary.csv"))
  interactionResultsSummary <- CohortMethod::getInteractionResultsSummary(outputFolder)
  readr::write_csv(interactionResultsSummary, file.path(outputFolder, "InteractionResultsSummary.csv"))
}
