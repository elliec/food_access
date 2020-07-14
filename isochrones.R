###########################
# Calculate isochrones ased on point data using ORS API 
# Eleanor Cervigni
#

options(scipen=999)

library(openrouteservice)
library(rgdal)
library(rgeos)
library(sp)
library(geojsonio)

ors_api_key("your api key")
ors_api_key()

setwd("C:/Users/Ellie/Desktop/Demographics_paper/trying_stuff/Mandurah")

#calculate centroids of polygons
polys <- readOGR("./output/mandmb.gpkg") #read in polygons spatial data
cents <- gCentroid(polys, byid=T) #calculate centroids
cents <- SpatialPointsDataFrame(cents,polys@data) #append original attributes

writeOGR(cents, dsn="./output/centroids.gpkg", layer="centroids", driver="GPKG") #save centroids

#import centroid points data and check projection
cents <- readOGR("./output/centroids.gpkg")
proj4string(cents)

n <- dim(cents)[1]

i <- 1
while(i<n+1){
	print(i)
	
	# API allows 5 points at once - first line calculates isochrones in groups of 5, 
	# second line calls individual points. if running ind. points must also change
	# to append attributes for only one point under "add in ID and coordinate atttributes" below
	# and the incrementing of i at the end

	coordinates <- list(c(cents$coords.x1[i], cents$coords.x2[i]), c(cents$coords.x1[i+1], cents$coords.x2[i+1]), c(cents$coords.x1[i+2], cents$coords.x2[i+2]), c(cents$coords.x1[i+3], cents$coords.x2[i+3]), c(cents$coords.x1[i+4], cents$coords.x2[i+4]))
	#coordinates <- c(cents$coords.x1[i], cents$coords.x2[i])

	x <- ors_isochrones(coordinates, profile = "foot-walking", range = 15*60)

	#set class geojson
	class(x) <- "geo_list"

	#convert geojson to spatialpolygon
	x_sp <- geojson_sp(x)

	#add in ID and coordinates attributes
	
	#for 5 points
	x_sp$mb_id <- as.character(cents$MB_CODE16[i:(i+4)])
	x_sp$long <- as.character(cents$coords.x1[i:(i+4)])
	x_sp$lat <- as.character(cents$coords.x2[i:(i+4)])

	#for 1 point
	#x_sp$mb_id <- as.character(cents$MB_CODE16[i])
	#x_sp$long <- as.character(cents$coords.x1[i])
	#x_sp$lat <- as.character(cents$coords.x2[i])

	#merge with previous output
	if(i<2){
	w_15min <- x_sp
	}else{
	w_15min <- rbind(w_15min, x_sp) 
	}

	# for 5 points
	i=i+4

	#for 1 point
	#i=i+1
}
head(w_15min@data)

writeOGR(w_15min, dsn="./output/w_15min.gpkg", layer="w_15min", driver="GPKG")


