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
# source("R/load_functions_demo_biProfile.R")

evaluate2dPade <- function(dataSet, folder) {
  alpha <- 0.05
  analysisId <- 3 # Poisson
  outcomeIds <- unique(dataSet$outcomeId)
  # i = 1
  for (i in seq_along(outcomeIds)) {
    fileName <- file.path(folder, sprintf("PadeArtifacts_o%d.rdata", outcomeIds[i]))
    if (!file.exists(fileName)) {
      message(sprintf("Processing outcome %d of %d", i, length(outcomeIds)))
      ppdata <- dataSet %>%
        filter(.data$outcomeId == outcomeIds[i],
               .data$analysisId == !!analysisId) %>%
        mutate(x1 = treatment * (1 - subgroup),
               x2 = treatment * subgroup,
               z1 = 1,
               z2 = subgroup) %>%
        select(period = "time",
               siteID = "siteId",
               stratumID = "stratumId",
               "x1", 
               "x2",
               "z1",
               "z2",
               "y") %>%
        as.data.frame()
      message("- Computing pooled estimate")
      start <- Sys.time()
      ResPool <- estimatePool(ppdata)
      delta <- Sys.time() - start
      message("  computing pooled estimate took ", signif(delta, 3), " ", attr(delta, "units"))
      
      message("- Computing meta-analytic estimate")
      ResMeta <- estimateMeta(ppdata)
      ebar <- ResMeta[[1]]
      
      message("- Computing local derivatives, and combining to global derivative")
      deriv<-GetGlobalDeriv(ebar, ppdata)
      
      message("- Computing local derivatives")
      J<- max(ppdata$siteID)
      LocalDeriv <- lapply(1:J, function(j){
        ipdata <- ppdata[ppdata$siteID==j,]
        return(GetLocalDeriv(ebar,ipdata))
      })
      
      PadeCoef <- EstPadeCoef(deriv)
      PECR <- PadeEstCR(ebar,PadeCoef)
      save(ResPool, ResMeta, deriv, LocalDeriv, PadeCoef, PECR, file = fileName)
      
    } 
    # e <- new.env()
    # load(file = file.path(folder, sprintf("PadeArtifacts_o%d.rdata", outcomeIds[i])), env = e)
    # ls(envir = e)
  }
}