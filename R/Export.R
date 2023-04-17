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
#' @export
computeDiagnostics <- function(outputFolder) {
  fileReference <- CohortMethod::getFileReference(outputFolder)
  
  # Negative controls ----------------------------------------------------------
  resultsSummary <- CohortMethod::getResultsSummary(outputFolder)
  interactionResultsSummary <- CohortMethod::getInteractionResultsSummary(outputFolder)
  # subset = subsets[[6]]
  plotNegativeControls <- function(subset, interactions = FALSE) {
    fileName <- file.path(outputFolder, sprintf("ncs%s_a%d.png", if(interactions) "_interactions" else "", subset$analysisId[1]))
    ncs <- subset %>%
      filter(.data$trueEffectSize == 1)
    hois <- subset %>%
      filter(is.na(.data$trueEffectSize))
    EmpiricalCalibration::plotCalibrationEffect(
      logRrNegatives = ncs$logRr,
      seLogRrNegatives = ncs$seLogRr,
      logRrPositives = hois$logRr,
      seLogRrPositives = hois$seLogRr,
      fileName = fileName,
      showCis = TRUE,
      showExpectedAbsoluteSystematicError = TRUE
    )
    invisible(NULL)
  }
  subsets <- split(resultsSummary, resultsSummary$analysisId)
  invisible(lapply(subsets, plotNegativeControls))
  subsets <- split(interactionResultsSummary, interactionResultsSummary$analysisId)
  invisible(lapply(subsets, plotNegativeControls, interactions = TRUE))
  
  # Propensity score distribution ----------------------------------------------
  ps <- readRDS(file.path(outputFolder, fileReference$sharedPsFile[1]))
  fileName <- file.path(outputFolder, "Ps.png")
  CohortMethod::plotPs(
    data = ps, 
    targetLabel = "Diclofenac",
    comparatorLabel = "Celexocib",
    showCountsLabel = TRUE,
    showEquiposeLabel = TRUE,
    fileName = fileName
  )
  
  # Covariate balance ----------------------------------------------------------
  balanceFiles <- unique(fileReference$sharedBalanceFile)
  for (balanceFile in balanceFiles) {
    fileName <- file.path(outputFolder, gsub(".rds", ".png", balanceFile))
    balance <- readRDS(file.path(outputFolder, balanceFile))
    CohortMethod::plotCovariateBalanceScatterPlot(balance, fileName = fileName)
  }
  for (balanceFile in balanceFiles) {
    fileName <- file.path(outputFolder, gsub("Balance", "BalTop", gsub(".rds", ".png", balanceFile)))
    balance <- readRDS(file.path(outputFolder, balanceFile))
    CohortMethod::plotCovariateBalanceOfTopVariables(balance, fileName = fileName)
  }
  
  # Concurrent medication ------------------------------------------------------
  # cmData <- CohortMethod::loadCohortMethodData(file.path(outputFolder, fileReference$cohortMethodDataFile[1]))
}

#' Create analysis dataset
#'
#' @param dataFolders A vector of folder names where the cohort method analyses have 
#'                    been executed.
#'
#' @return 
#' A tibble with one row per person, across all outcomes, analyses, and databases
#' 
#' @export
createAnalysisDataSet <- function(dataFolders) {
  data <- list()
  for (i in seq_along(dataFolders)) {
    dataFolder <- dataFolders[i]
    message(sprintf("Collecting data from %s", dataFolder))
    fileReference <- CohortMethod::getFileReference(dataFolder)
    fileReference <- fileReference %>%
      filter(analysisId %in% c(3, 4))
    preFilteredCmData <- CohortMethod::loadCohortMethodData(file.path(dataFolder, fileReference$prefilteredCovariatesFile[1]))
    interactionCov <- preFilteredCmData$covariates %>%
      select("rowId", subgroup = "covariateValue") %>%
      collect()
    for (j in seq_len(nrow(fileReference))) {
      strataPop <- readRDS(file.path(dataFolder, fileReference$strataFile[j]))
      pop <- strataPop %>%
        left_join(interactionCov, by = join_by(rowId)) %>%
        mutate(subgroup = if_else(is.na(.data$subgroup), 0, .data$subgroup)) %>%
        select(y = "outcomeCount", time = "timeAtRisk", "survivalTime", "treatment", "stratumId", "subgroup") %>%
        mutate(
          siteId = i,
          analysisId = fileReference$analysisId[j],
          outcomeId = fileReference$outcomeId[j]
        )
      data[[length(data)+ 1]] <- pop
    }
  }
  data <- bind_rows(data)
  negativeControls <- readr::read_csv(
    file = system.file("NegativeControls.csv", package = "InteractionsEval"),
    show_col_types = FALSE
  )
  data <- data %>%
    mutate(negativeControl = .data$outcomeId %in% negativeControls$outcomeId)
  return(data)
}

