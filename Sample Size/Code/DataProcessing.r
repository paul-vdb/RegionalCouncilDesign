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

# Step 1: Add LCDB class to all sites.

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
library(reshape2)

# Start with Tier 1 data:
load("Data/bdi_meta.Rda")		# Object: metadata
load("Data/FiveMinBird.Rda")	# Object: FiveMinBird

# Read in the Wellington bird data
well <- data.table(read_excel("Data/Tier I bird count data summary.xlsx"))
setnames(well, "Site", "Place")
allLocs <- fread("Data/Appendix_1_static_master_data_July2016.csv", skip = 10)
setnames(allLocs, "PlotID", "Place")

# Read in lcdb use map data.
# Keep the data elswhere, too big for github.
lcdb <- readOGR("../../Examples/Data/mfe-lucas-nz-land-use-map-1990-2008-2012-v018-SHP", "lucas-nz-land-use-map-1990-2008-2012-v018")
rc <- readOGR("../../Examples/Data/RC_Boundaries", "RegionalCouncilBoundaries")
pcl <- readOGR("../../Examples/Data/PCL", "PCL_08July2016")

# Lookup table for LUC_NAME aggregation
lookup <- fread("Data/LandMapUseLookUp.csv")

# nztm <- proj4string(lcdb)
nztm <- "+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +units=m +no_defs"
locs <- SpatialPointsDataFrame(SpatialPoints(cbind(allLocs$NZTM_Easting, allLocs$NZTM_Northing), proj4string = CRS(nztm)), data = data.frame(Place = allLocs$Place, Region = allLocs$Region))

dat.lcdb <- over(locs, lcdb)
locs@data <- cbind(locs@data, dat.lcdb)
save(locs, file = "Data/LocationsLCDB.Rda")
load("Data/LocationsLCDB.Rda")

dat.pcl <- over(locs, pcl)
locs@data <- cbind(locs@data, data.frame(PCL2016 = !is.na(dat.pcl$Type)))

#####################
#Process Tier 1 data as we do in annual report.
#####################
#Hash out the dates
metadata[, "StartDate" := as.Date(Start.Date, format = "%d/%m/%Y")]

#Need just the most recent of the observed plot data
metadata[, "Season" := as.numeric(format(StartDate, '%Y'))]
metadata[as.numeric(format(StartDate, '%m')) <= 06, "Season" := Season - 1]
metadata = metadata[,.SD[Season == max(Season)], by = Place]
metadata <- metadata[!is.na(Wood)]

#Quick averages of 5MBC
avg.birds <- dcast.data.table(data = FiveMinBird, Place + Season + Station ~ CommonName, value.var = "TotalCount", fun = sum, fill = 0)
avg.birds[,"NA" := NULL]

##########################################
# Do the equivallent with Wellington Data
##########################################
well[, "Season" := as.numeric(format(Date, '%Y'))]
well[as.numeric(format(Date, '%m')) <= 06, "Season" := Season - 1]

avg.birds.w <- dcast.data.table(data = well, Place + Season + Point ~ Species_name, value.var = "Number", fun = sum, fill = 0)

avg.birds.w <- merge(avg.birds.w, locs@data, all.x = TRUE, all.y = FALSE, by = "Place")
avg.birds <- merge(avg.birds, locs@data, all.x = TRUE, all.y = FALSE, by = "Place")
avg.birds <- merge(avg.birds, lookup, by = "LUC_NAME")

t1.birds.w <- avg.birds[Region == "Wellington Region"]

wellyFans <- rbind(t1.birds.w[,.(Place, Station, Season, Region, Fantail, Bellbird, GreyWarbler = `Grey Warbler`, LUC_NAME)], 
	avg.birds.w[,.(Place, Station = Point, Season, Region, Fantail, Bellbird, GreyWarbler = `Warbler_Grey`, LUC_NAME)])
wellyFans[, "time" := Season - min(Season)]

wellyFans <- merge(wellyFans, lookup)
wellyFans[, .(fan = sd(Fantail), bell = sd(Bellbird), greyWarb = sd(GreyWarbler)), by = "Level"]
avg.birds[, .(fan = sd(Fantail), bell = sd(Bellbird), greyWarb = sd(`Grey Warbler`))]

m <- glmer.nb(Fantail ~ -1 + Level + (1|Station), data = wellyFans)

# Build scenarios:
# Scenario 1. Can we find the proportion of effort to put into off of Native cover types?
# Does equal knowledge occur if 
#############
results.nat <- data.table()
results.imp <- data.table()
for(i in seq(0, 80, by = 5))
{
	X <- expand.grid(Level = c(rep("Native", 90 - i), rep("Impacted", 10 + i)), Station = c("A", "B", "C", "D", "E"))
	getData(m) <- X
	
	tmpImp <- data.table(getWidths(powerCurve(m, test = ciWidth("LevelImpacted"), along = "Level", nsim = 10)), NI = 10 + i)
	results.imp <- rbind(results.imp, tmpImp)
	tmpNat <- data.table(getWidths(powerCurve(m, test = ciWidth("LevelNative"), along = "Level", nsim = 10)), NI = 10 + i)
	results.nat <- rbind(results.nat, tmpNat)
}

# For stratified sampling need to know the true proportion of Native and Impacted:
lcdb <- merge(lcdb, lookup, by = "LUC_NAME", all.x = TRUE, all.y = FALSE, sort = FALSE)	# Need to merge with sp merge for the sake of order.
locs <- merge(locs, lookup, by = "LUC_NAME", all.x = TRUE, all.y = FALSE, sort = FALSE)	# Need to merge with sp merge for the sake of order.

lcdb.well <- lcdb[lcdb$LUM_REG_NA == "Wellington",]

dt.well <- data.table(lcdb.well@data)
dt.well <- dt.well[, .(totArea = sum(AREA_HA)), by = "Level"]
dt.well[, prop := totArea/sum(totArea)]

results.nat[, "Level" := "Native"]
results.imp[, "Level" := "Impacted"]

results <- rbind(results.nat, results.imp)
ggplot(data = results, aes(x = NI, y = mean, colour = Level)) + geom_point() + geom_smooth()

getData(m) <- X
X <- data.table(X)
X$Y <- doSim(m)
X[, .(v = var(Y),), by = "Level"]