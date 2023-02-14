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
computeDiagnostics <- function(outputFolder) {
  fileReference <- CohortMethod::getFileReference(outputFolder)
  
  # Negative controls ----------------------------------------------------------
  resultsSummary <- CohortMethod::getResultsSummary(outputFolder)
  interactionResultsSummary <- CohortMethod::getInteractionResultsSummary(outputFolder)
  # subset = subsets[[1]]
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
      fileName = fileName
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
  
  # Concurrent medication ------------------------------------------------------
  cmData <- CohortMethod::loadCohortMethodData(file.path(outputFolder, fileReference$cohortMethodDataFile[1]))
}