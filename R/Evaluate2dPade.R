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
# source("R/load_functions_demo_coxProfile.R")

evaluate2dPadePoisson <- function(dataSet, folder) {
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

computeRrsFrom2dPade <- function(folder) {
  # alpha <- 0.05
  results <- tibble()
  padeFiles <- list.files(folder, "PadeArtifacts_o")
  for (padeFile in padeFiles) {
    load(file.path(folder, padeFile))
    # Recompute Pade estimate to get CIs:
    ebar <- ResMeta[[1]]
    ResPade <- PadeEstCR(ebar = ebar, PadeCoef = PadeCoef)
    coef <- ResPade$CR
    coef.CI <- ResPade$CI
    plotCI.pade.ellipse <- ellipse.plot.pade(coef.CI, ebar)
    plotCI.pade <- ci.plot.pade(coef.CI, ebar)
    ci <- plotCI.pade
    outcomeId <- as.numeric(gsub("^PadeArtifacts_o", "", gsub(".rdata$", "", padeFile)))
    result <- tibble(
      outcomeId = outcomeId,
      var = paste0("x", c(1, 2)),
      pooled = ResPool$est,
      pade = PECR$Est,
      logRr = PECR$Est,
      logCi95Lb = c(ci$CI1[1], ci$CI2[1]),
      logCi95Ub = c(ci$CI1[2], ci$CI2[2]),
      seLogRr = c(ci$CI1[2] - ci$CI1[1], ci$CI2[2] - ci$CI2[1]) / (2*qnorm(0.975))
    )
    results <- bind_rows(results, result)
  }
  library(ggplot2)
  results$negative <- results$outcomeId != 77
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 8)
  ggplot(results, aes(x = exp(pooled), y = exp(pade), group = var, color = negative)) +
    geom_abline(slope = 1) + 
    geom_point(alpha = 1) +
    scale_x_log10(breaks = breaks, limits = c(0.1, 10)) +
    scale_y_log10(breaks = breaks, limits = c(0.1, 10)) +
    facet_grid(~var)
  ggsave(file.path(folder, "PadeVsPooledPoisson.png"), width = 7, height = 3.5)
  ncs1 <- results %>%
    filter(var == "x1" & negative) 
  ncs2 <- results %>%
    filter(var == "x2" & negative) 
  hoi1 <- results %>%
    filter(var == "x1" & !negative) 
  hoi2 <- results %>%
    filter(var == "x2" & !negative)
  EmpiricalCalibration::plotCalibrationEffect(
    logRrNegatives = ncs1$logRr,
    seLogRrNegatives = ncs1$seLogRr,
    logRrPositives = hoi1$logRr,
    seLogRrPositives = hoi1$seLogRr,
    xLabel = "Incidence Rate Ratio",
    title = "Variable x1",
    showCis = TRUE,
    showExpectedAbsoluteSystematicError = TRUE,
    fileName = file.path(folder, "PadePoissonX1.png")
  )
  EmpiricalCalibration::plotCalibrationEffect(
    logRrNegatives = ncs2$logRr,
    seLogRrNegatives = ncs2$seLogRr,
    logRrPositives = hoi2$logRr,
    seLogRrPositives = hoi2$seLogRr,
    xLabel = "Incidence Rate Ratio",
    title = "Variable x2",
    showCis = TRUE,
    showExpectedAbsoluteSystematicError = TRUE,
    fileName = file.path(folder, "PadePoissonX2.png")
  )
}

evaluate2dPadeCox <- function(dataSet, folder) {
  alpha <- 0.05
  analysisId <- 4 # Cox
  outcomeIds <- unique(dataSet$outcomeId)
  # i = 1
  for (i in seq_along(outcomeIds)) {
    fileName <- file.path(folder, sprintf("PadeArtifactsCox_o%d.rdata", outcomeIds[i]))
    if (!file.exists(fileName)) {
      message(sprintf("Processing outcome %d of %d", i, length(outcomeIds)))
      ppdata <- dataSet %>%
        filter(.data$outcomeId == outcomeIds[i],
               .data$analysisId == !!analysisId) %>%
        mutate(x1 = treatment * (1 - subgroup),
               x2 = treatment * subgroup,
               z1 = subgroup,
               status = as.numeric(y > 0)) %>%
        select("time",
               "status",
               "x1", 
               "x2",
               "z1",
               siteID = "siteId",
               stratumID = "stratumId",
               male = "subgroup"
              ) %>%
        as.data.frame()
      message("- Computing pooled estimate")
      start <- Sys.time()
      ResPool <- estimatePoolCox(ppdata)
      delta <- Sys.time() - start
      message("  computing pooled estimate took ", signif(delta, 3), " ", attr(delta, "units"))
      
      message("- Computing meta-analytic estimate")
      ResMeta <- estimateMetaCox(ppdata)
      ebar <- ResMeta[[1]]
      
      message("- Computing local derivatives, and combining to global derivative")
      deriv<-GetGlobalDerivCox(ebar, ppdata)
      
      message("- Computing local derivatives")
      J<- max(ppdata$siteID)
      LocalDeriv <- lapply(1:J, function(j){
        ipdata <- ppdata[ppdata$siteID==j,]
        return(GetLocalDerivCox(ebar,ipdata))
      })
      
      PadeCoef <- EstPadeCoefCox(deriv)
      PECR <- PadeEstCRCox(ebar,PadeCoef)
      save(ResPool, ResMeta, deriv, LocalDeriv, PadeCoef, PECR, file = fileName)
    } 
    # e <- new.env()
    # load(file = file.path(folder, sprintf("PadeArtifacts_o%d.rdata", outcomeIds[i])), env = e)
    # ls(envir = e)
  }
}

computeRrsFrom2dPadeCox <- function(folder) {
  # alpha <- 0.05
  results <- tibble()
  padeFiles <- list.files(folder, "PadeArtifactsCox_o")
  for (padeFile in padeFiles) {
    load(file.path(folder, padeFile))
    # Recompute Pade estimate to get CIs:
    ebar <- ResMeta[[1]]
    ResPade <- PadeEstCRCox(ebar = ebar, PadeCoef = PadeCoef)
    coef <- ResPade$CR
    coef.CI <- ResPade$CI
    plotCI.pade.ellipse <- ellipse.plot.pade(coef.CI, ebar)
    plotCI.pade <- ci.plot.pade(coef.CI, ebar)
    ci <- plotCI.pade
    outcomeId <- as.numeric(gsub("^PadeArtifactsCox_o", "", gsub(".rdata$", "", padeFile)))
    result <- tibble(
      outcomeId = outcomeId,
      var = paste0("x", c(1, 2)),
      pooled = ResPool$est,
      pade = PECR$Est,
      logRr = PECR$Est
    )
    results <- bind_rows(results, result)
  }
  library(ggplot2)
  results$negative <- results$outcomeId != 77
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 8)
  ggplot(results, aes(x = exp(pooled), y = exp(pade), color = negative)) +
    geom_abline(slope = 1) + 
    geom_point(alpha = 1) +
    scale_x_log10(breaks = breaks, limits = c(0.1, 10)) +
    scale_y_log10(breaks = breaks, limits = c(0.1, 10)) +
    facet_grid(~var)
  ggsave(file.path(folder, "PadeVsPooledCox.png"), width = 7, height = 3.5)
}
