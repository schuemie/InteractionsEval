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

runCohortMethod <- function(connectionDetails,
                            cdmDatabaseSchema,
                            cohortDatabaseSchema,
                            cohortTable,
                            outputFolder,
                            maxCores) {
  # Create settings ------------------------------------------------------------
  outcomeOfInterest <- CohortMethod::createOutcome(
    outcomeId = 77,
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
        outcomeOfInterest = TRUE,
        trueEffectSize = 1
      )
    }
  )

  tcos <- CohortMethod::createTargetComparatorOutcomes(
    targetId = 11038,
    comparatorId = 11037,
    outcomes = append(
      list(outcomeOfInterest),
      negativeControlOutcomes
    )
  )
  targetComparatorOutcomesList <- list(tcos)

  conceptIdsToExclude <- 21603933 # NSAIDs

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

  createStudyPopArgsOnTreatment <- CohortMethod::createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = TRUE,
    removeDuplicateSubjects = "keep first",
    minDaysAtRisk = 1,
    riskWindowStart = 0,
    startAnchor = "cohort start",
    riskWindowEnd = 0,
    endAnchor = "cohort end",
    censorAtNewRiskWindow = TRUE
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

  fitOutcomeModelArgsPoisson <- CohortMethod::createFitOutcomeModelArgs(
    modelType = "poisson",
    stratified = TRUE
  )

  cmAnalysis1 <- CohortMethod::createCmAnalysis(
    analysisId = 1,
    description = "Poisson regression",
    getDbCohortMethodDataArgs = getDbCmDataArgs,
    createStudyPopArgs = createStudyPopArgsOnTreatment,
    createPsArgs = createPsArgs,
    stratifyByPsArgs = stratifyByPsArgs,
    computeSharedCovariateBalanceArgs = computeSharedCovBalArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgsPoisson
  )

  fitOutcomeModelArgsCox <- CohortMethod::createFitOutcomeModelArgs(
    modelType = "cox",
    stratified = TRUE
  )

  cmAnalysis2 <- CohortMethod::createCmAnalysis(
    analysisId = 2,
    description = "Cox regression",
    getDbCohortMethodDataArgs = getDbCmDataArgs,
    createStudyPopArgs = createStudyPopArgsOnTreatment,
    createPsArgs = createPsArgs,
    stratifyByPsArgs = stratifyByPsArgs,
    computeSharedCovariateBalanceArgs = computeSharedCovBalArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgsCox
  )


  warfarinDrugEraOverlap <- 1310149413

  fitOutcomeModelArgsPoissonInteraction <- CohortMethod::createFitOutcomeModelArgs(
    modelType = "poisson",
    stratified = TRUE,
    interactionCovariateIds = warfarinDrugEraOverlap
  )

  fitOutcomeModelArgsCoxInteraction <- CohortMethod::createFitOutcomeModelArgs(
    modelType = "cox",
    stratified = TRUE,
    interactionCovariateIds = warfarinDrugEraOverlap
  )

  cmAnalysis3 <- CohortMethod::createCmAnalysis(
    analysisId = 3,
    description = "Poisson regression",
    getDbCohortMethodDataArgs = getDbCmDataArgs,
    createStudyPopArgs = createStudyPopArgsOnTreatment,
    createPsArgs = createPsArgs,
    stratifyByPsArgs = stratifyByPsArgs,
    computeSharedCovariateBalanceArgs = computeSharedCovBalArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgsPoissonInteraction
  )

  cmAnalysis4 <- CohortMethod::createCmAnalysis(
    analysisId = 4,
    description = "Cox regression",
    getDbCohortMethodDataArgs = getDbCmDataArgs,
    createStudyPopArgs = createStudyPopArgsOnTreatment,
    createPsArgs = createPsArgs,
    stratifyByPsArgs = stratifyByPsArgs,
    computeSharedCovariateBalanceArgs = computeSharedCovBalArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgsCoxInteraction
  )
  
  matchOnPsArgs <- CohortMethod::createMatchOnPsArgs(maxRatio = 1)
  fitOutcomeModelArgsPoissonNonStratified <- CohortMethod::createFitOutcomeModelArgs(
    modelType = "poisson",
    stratified = FALSE
  )
  cmAnalysis5 <- CohortMethod::createCmAnalysis(
    analysisId = 5,
    description = "Poisson regression, matching",
    getDbCohortMethodDataArgs = getDbCmDataArgs,
    createStudyPopArgs = createStudyPopArgsOnTreatment,
    createPsArgs = createPsArgs,
    matchOnPsArgs = matchOnPsArgs,
    computeSharedCovariateBalanceArgs = computeSharedCovBalArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgsPoissonNonStratified
  )
  fitOutcomeModelArgsCoxNonStratified <- CohortMethod::createFitOutcomeModelArgs(
    modelType = "cox",
    stratified = FALSE
  )
  cmAnalysis6 <- CohortMethod::createCmAnalysis(
    analysisId = 6,
    description = "Cox regression, matching",
    getDbCohortMethodDataArgs = getDbCmDataArgs,
    createStudyPopArgs = createStudyPopArgsOnTreatment,
    createPsArgs = createPsArgs,
    matchOnPsArgs = matchOnPsArgs,
    computeSharedCovariateBalanceArgs = computeSharedCovBalArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgsCoxNonStratified
  )
  cmAnalysisList <- list(cmAnalysis1, cmAnalysis2, cmAnalysis3, cmAnalysis4, cmAnalysis5, cmAnalysis6)
  # cmAnalysisList <- list(cmAnalysis2, cmAnalysis4)

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
