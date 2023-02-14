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

# Format and check code:
styler::style_pkg()
OhdsiRTools::checkUsagePackage("InteractionsEval")
OhdsiRTools::updateCopyrightYearFolder()
devtools::spell_check()

# Create manual:
unlink("extras/InteractionsEval.pdf")
shell("R CMD Rd2pdf ./ --output=extras/InteractionsEval.pdf")

# Insert cohort definitions from ATLAS into package -----------------------
remotes::install_github("ohdsi/ROhdsiWebApi")
ROhdsiWebApi::authorizeWebApi(
  baseUrl = Sys.getenv("baseUrl"),
  authMethod = "windows"
)
ROhdsiWebApi::insertCohortDefinitionSetInPackage(
  fileName = "inst/Cohorts.csv",
  baseUrl = Sys.getenv("baseUrl"),
  insertTableSql = FALSE,
  insertCohortCreationR = FALSE,
  generateStats = FALSE,
  packageName = "InteractionsEval"
)
