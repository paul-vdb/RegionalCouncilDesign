##########################
# Have a crack at sample size for
# regional councils based on "heterogeneity"
# of their LCDB classification.
#
# Paul van Dam-Bates
# Decemeber 28, 2018
##########################
# Logic here is that heterogeneity is not a function of "area"
# but a function of heterogeneity of landcover. 
# We will model birds as a function of different forest types and
# then act as if we will model them as such.

setwd("..")
library(data.table)
library(dplyr)
library(readxl)
library(tidyverse)
library(simr)
library(rgdal)
library(sp)
library(rgeos)
library(sf)

# Start with Tier 1 data:
load("Data/bdi_meta.Rda")		# Object: metadata
load("Data/FiveMinBird.Rda")	# Object: FiveMinBird

metadata
FiveMinBird

# Read in the Wellington bird data
well <- data.table(read_excel("Data/Tier I bird count data summary.xlsx"))
setnames(well, "Site", "PlotID")
allLocs <- fread("Data/Appendix_1_static_master_data_July2016.csv", skip = 10)

well <- merge(well, allLocs, all.x = TRUE, all.y = FALSE)

# Read in lcdb data:
# Keep the data elswhere, too big for github.
lcdb <- readOGR("../../Examples/Data/lris-lcdb-v41-land-cover-database-version-41-mainland-new-zealand-SHP", "lcdb-v41-land-cover-database-version-41-mainland-new-zealand")

