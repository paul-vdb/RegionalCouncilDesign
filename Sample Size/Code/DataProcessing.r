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
library(raster) # intersect function to try.

# library(devtools)
# install_github("ogansell/MSampNZ")
library(MSampNZ)

# Helpful fast clip function.
gClip <- function(shp, bb){
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else b_poly <- as(extent(bb), "SpatialPolygons")
  proj4string(b_poly) <- proj4string(shp)
  gIntersection(shp, b_poly, byid = T)
}

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

# Speed things up by not "rerunning" this bit.
if(!any("LocationsLCDB.Rda" %in% dir("Data/"))){
	dat.lcdb <- over(locs, lcdb)
	locs@data <- cbind(locs@data, dat.lcdb)

	dat.pcl <- over(locs, pcl)
	locs@data <- cbind(locs@data, data.frame(PCL2016 = !is.na(dat.pcl$Type)))

	save(locs, file = "Data/LocationsLCDB.Rda")
}
load("Data/LocationsLCDB.Rda")



# I also want to know what proportion of rc is pcl.
# Super slow so don't do it if you don't have to.
if(!any("AreaPCL.Rda" %in% dir("Data/"))){
	area.pcl <- data.table()
	for( i in unique(rc$NAME) )
	{
		tmp.shp <- rc[rc$NAME == i,]
		tmp.int <- gIntersection(pcl, tmp.shp, byid = TRUE, drop_lower_td = TRUE)
		area.i <- sum(area(tmp.int))
		area.rc <- sum(area(tmp.shp))
		area.pcl <- rbind(area.pcl, data.table(NAME = i, areaPCL = area.i, areaRC = area.rc))
	}
	area.pcl[, "PropPCL" := areaPCL/areaRC]
save(area.pcl, file = "Data/AreaPCL.Rda")
}else{
	load("Data/AreaPCL.Rda")
}
area.pcl[, "LUM_REG_NA" := gsub(" Region|'", "", NAME)]

# Build tables for the report.
# For stratified sampling need to know the true proportion of Native and Impacted:
# We also need a table with area, proportion of each strata and then proportion of that which is PCL.
#####################################################################################################
lcdb <- merge(lcdb, lookup, by = "LUC_NAME", all.x = TRUE, all.y = FALSE, sort = FALSE)	# Need to merge with sp merge for the sake of order.
locs <- merge(locs, lookup, by = "LUC_NAME", all.x = TRUE, all.y = FALSE, sort = FALSE)	# Need to merge with sp merge for the sake of order.

dat.lcdb <- data.table(lcdb@data)
dat.cov <- dat.lcdb[,.(area = sum(AREA_HA)), by = c("LUM_REG_NA", "Level")]
dat.cov[, "prop" := area/sum(area), by = "LUM_REG_NA"]
tab.cover <- dcast.data.table(dat.cov, LUM_REG_NA ~ Level, value.var = "prop")
tab.cover <- merge(tab.cover, area.pcl[,.(LUM_REG_NA, areaRC, PropPCL)], all.x = TRUE, all.y = FALSE)
tab.cover <- tab.cover[,.(Region = LUM_REG_NA, Area = round(areaRC/1000000,1), PCL = round(PropPCL*100, 1), Impacted = round(Impacted*100,1), Native = round(Native*100, 1))]
write.csv(tab.cover, "Data/CoverTypesTable.csv", row.names = FALSE)

dat.samp <- data.table(locs@data)
dat.samp <- dat.samp[, .N, by = c("LUM_REG_NA", "PCL2016", "Level")]
tab.sites <- dcast.data.table(dat.samp, LUM_REG_NA ~ Level + PCL2016, value.var = "N")
tab.sites[!is.na(Impacted_TRUE), Native_TRUE := Native_TRUE + Impacted_TRUE]
tab.sites[, Impacted_TRUE := NULL]
tab.sites[, "Native" := Native_TRUE + Native_FALSE]
tab.sites <- tab.sites[,.(Region = LUM_REG_NA, Native = Native,  NativePCL = Native_TRUE, Impacted = Impacted_FALSE)]
tab.sites[is.na(Impacted), Impacted := 0]
tab.sites[, "Total" := Native + Impacted ]
write.csv(tab.sites, "Data/SampleSizeGridTable.csv", row.names = FALSE)

tab.total <- merge(tab.sites, tab.cover, by = "Region")
tab.total[, "n_Total" := ceiling(Native.x/(2*Native.y/(2*Native.y + 1*Impacted.y)))]
tab.total[,"Native_Ratio" := (2*Native.y/(2*Native.y + 1*Impacted.y))]
tab.total[, "n_Impacted" := n_Total - Native.x]
tab.total <- tab.total[, .(Region, Total = n_Total, Native_PCL = NativePCL, Native_RC = Native.x - NativePCL, Impacted = n_Impacted)]
tab.total[, RC_Total := Impacted + Native_RC]
write.csv(tab.total, "Data/SampleSizeTable.csv", row.names = FALSE)

########################################################
# Now create the map of Southland for Tier 1:
########################################################
rc.south <- rc[rc$NAME == "Southland Region",]
southland.df <- fortify(rc.south)
locs.south <- data.frame(locs[locs$Region == "Southland Region",])
pcl.south <- gClip(pcl, bbox(rc.south))
pcl.df <- fortify(pcl.south)

ggplot(data = southland.df, aes(long, lat, group = group)) + geom_polygon(fill = "white", colour = "black", size = 1) +
	geom_polygon(data = pcl.df, fill = "grey", size = .5, alpha = 0.8) + 
	geom_point(data = locs.south, aes(x = coords.x1, y = coords.x2, shape = PCL2016), group = 1) + 
	scale_shape_manual(name = "", values = c("TRUE" = 4, "FALSE" = 16), labels = c( "Non-PCL", "PCL")) + 
	theme_bw() + coord_fixed() + xlab("Easting (m)") + ylab("Northing (m)")
ggsave("Data/Tier1Southland.png")

locs.y1 <- locs.south[locs.south$PCL2016 == "FALSE",]
locs.y1 <- locs.y1[sample(1:nrow(locs.y1), ceiling(nrow(locs.y1)/5)),]

ggplot(data = southland.df, aes(long, lat, group = group)) + geom_polygon(fill = "white", colour = "black", size = 1) +
	geom_polygon(data = pcl.df, fill = "grey", size = .5, alpha = 0.8) + 
	geom_point(data = rbind(locs.south[locs.south$PCL2016 == TRUE,], locs.y1), aes(x = coords.x1, y = coords.x2, shape = PCL2016), group = 1) + 
	scale_shape_manual(name = "", values = c("TRUE" = 4, "FALSE" = 16), labels = c( "Non-PCL", "PCL")) + 
	theme_bw() + coord_fixed() + xlab("Easting (m)") + ylab("Northing (m)")
ggsave("Data/Tier1SouthlandYear1.png")

ms.pts.native <- masterSample(N = 500, shp = lcdb[lcdb$LUM_REG_NA == "Southland" & lcdb$Level == "Native",], island = "South", J = c(4,3))
ms.pts.impact <- masterSample(N = 500, shp = lcdb[lcdb$LUM_REG_NA == "Southland" & lcdb$Level == "Impacted",], island = "South", J = c(4,3))
ms.pts.impact$Level <- "Impacted"
ms.pts.native$Level <- "Native"
ms.pts.impact <- spTransform(ms.pts.impact, CRS(proj4string(lcdb)))
ms.pts.native <- spTransform(ms.pts.native, CRS(proj4string(lcdb)))

pts1 <- ms.pts.impact[pcl,]
pts2 <- ms.pts.native[pcl,]

ms.pts.impact <- ms.pts.impact[!ms.pts.impact$SiteID %in% pts1$SiteID,]
ms.pts.native <- ms.pts.native[!ms.pts.native$SiteID %in% c(pts2$SiteID, "South2521"),]	# One point is ending up in Milford sound and I think it's a spatial processing issue.

ms.pts <- rbind(ms.pts.impact[1:55,], ms.pts.native[1:112,])

ggplot(data = southland.df, aes(long, lat, group = group)) + geom_polygon(fill = "white", colour = "black", size = 1) +
	geom_polygon(data = pcl.df, fill = "grey", size = .5, alpha = 0.8) + 
	geom_point(data = locs.south[locs.south$PCL2016 == TRUE,], aes(x = coords.x1, y = coords.x2), group = 1, shape = 4) + 
	geom_point(data = data.frame(ms.pts), aes(x = coords.x1, y = coords.x2, shape = Level, colour = Level), group = 1, size = 1.5) +
	scale_shape_manual(name = "Land Use Type", values = c("Impacted" = 16, "Native" = 17), labels = c("Impacted", "Native")) + 
	scale_colour_manual(name = "Land Use Type", values = c("Impacted" = "red", "Native" = "blue"), labels = c("Impacted", "Native")) + 
	theme_bw() + coord_fixed() + xlab("Easting (m)") + ylab("Northing (m)") 
ggsave("Data/StratifiedSouthland.png")

ms.pts <- rbind(ms.pts.impact[1:11,], ms.pts.native[1:23,])
ggplot(data = southland.df, aes(long, lat, group = group)) + geom_polygon(fill = "white", colour = "black", size = 1) +
	geom_polygon(data = pcl.df, fill = "grey", size = .5, alpha = 0.8) + 
	geom_point(data = locs.south[locs.south$PCL2016 == TRUE,], aes(x = coords.x1, y = coords.x2), group = 1, shape = 4) + 
	geom_point(data = data.frame(ms.pts), aes(x = coords.x1, y = coords.x2, shape = Level, colour = Level), group = 1, size = 1.5) +
	scale_shape_manual(name = "Land Use Type", values = c("Impacted" = 16, "Native" = 17), labels = c("Impacted", "Native")) + 
	scale_colour_manual(name = "Land Use Type", values = c("Impacted" = "red", "Native" = "blue"), labels = c("Impacted", "Native")) + 
	theme_bw() + coord_fixed() + xlab("Easting (m)") + ylab("Northing (m)")
ggsave("Data/StratifiedSouthlandYear1.png")



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


