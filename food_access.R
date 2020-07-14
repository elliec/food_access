###############################################
# Food access clusters and profiles
# Eleanor Cervigni
# set copyright
#
###############################################

options(scipen=999)
options(stringsAsFactors = FALSE)
library(rgdal)
library(tmap)
library(openrouteservice)

setwd("C:/Users/Ellie/Desktop/Demographics_paper/3_accessibility")

################################################
# Nearest grocery
# Uses the OpenRouteService API to calculate the nearest grocery from each centroid based on a given food outlet dataset
# Due to API limits the script will stop and save the outputs if limits are reached and must be restarted 24hrs later at the
# point where it stopped
#################################################

API.key <- "5b3ce3597851110001cf624868d86ecb470b4a28b00aea82dc491ade"
cents <- readOGR("./input/perthmbcent.gpkg")
outletdata <- readOGR("./input/open_data.gpkg")
output.csv <- "./output/od_nearest_grocery.csv"
groceryCat <- c("Supermarket","Fruit and veg") # list of categories assigned to grocery stores
catColumn <- "level.1" # name of column containing categories
nCents <- 26296 #number of centroids

#set API key
ors_api_key(API.key)
ors_api_key()

#extract grocery stores
grocery <- outletdata[which(outletdata[[catColumn]] %in% groceryCat), ]
grocery$rowid <- row.names(grocery) #assign row names of original dataframe

#create empty data frame for results
ng <- NULL

#calculate number of origin points per call, (limit is 2500 combinations per API call)
n <- as.integer(2500/nrow(grocery))
n2 <- n-1

#extract coordinates of grocery stores
df <- data.frame(coordinates(grocery))

i=1
while(i<nCents+1){
	print(i)
	points <- cents[i:(i+n2),]

	#extract long and lat
	coords <- data.frame(coordinates(points))

	# join long and lat to top of grocery stores
	df_mb <- rbind(coords, df)

	# query for duration and distance in meters (ORS API call)
	res = ors_matrix(df_mb, sources=c(0:n2), metrics = c("duration", "distance"), units = "m")

	r<- 1
	for(r in 1:n){
		#print(r)
		dur <- res$durations[r,-c(1:n)]
		dist <- res$distances[r,-c(1:n)]
		
		m <- which(dur == min(dur, na.rm = TRUE))
		t <- which(dist == min(dist, na.rm = TRUE))

		results <- data.frame(mb_id=points$MB_CODE16[r], grocery_m=grocery$rowid[m[1]], meters=dist[m[1]], grocery_t=grocery$rowid[t[1]], time_sec=dur[t[1]])
		ng <- rbind(ng, results)
	}

	i=i+n
}

write.csv(ng,"./output/ng.csv")

########################################
# Abundance
# Calculates the abundance of each food outlet type within each isochrone
#######################################
library(rgdal)
library(tmap)

iso <- readOGR("./input/grperthiso.gpkg") #isochrone polygons
outlet.points <- readOGR("./input/open_data.gpkg") #food outlet data points
catColumn <- "level.1" # name of column containing categories

#convert category column to factor for counting outlets
outlet.points[[catColumn]] <- as.factor(outlet.points[[catColumn]])

#create empty dataframe for appending counts
abund <- NULL

#count outlets within each isochrone
for(i in 1:length(iso)){
#for(i in 1:10){

	print(i)

	#select current polygon
	poly <- iso[i,]

	#subset points based on polygon
	outlets <- outlet.points[poly,]

	#Get counts per category and then transpose into dataframe
	df <- data.frame(table(outlets$level.1))
	trans_df = setNames(data.frame(t(df[,-1])), df[,1])

	#add MB id
	trans_df$id <- poly$mb_id
	
	#append results
	abund <- rbind(abund,trans_df)
}


#append total category counts
df <- as.data.frame(table(outlet.points$level.1))
trans_df = setNames(data.frame(t(df[,-1])), df[,1])

trans_df$id <- "99999999999"
abund <- rbind(abund,trans_df)

write.csv(abund, "./output/od_abundance.csv")

####################################################
# Accessibility
# Calculates different measures of accessibility based on abundance data
####################################################

abund <- read.csv("./output/od_abundance.csv")
abund <- abund[,-1]
mb <- readOGR("./input/perthmbcent.gpkg")
healthy <- c(2,3,6) #indexes of columns containing "healthy" outlets. Required if calculating MRFEI
unhealthy <- c(1,4,5) #indexes of columns containing "unhealthy" outlets. Required if calculating MRFEI
totalOutlets <- c(1:6) #indexes of columns containing abundance data - for calculating total abundance(number of outlets) around a MB

groceryCat <- c(3,7) #indexes of columns containing fruit and veg (eg supermarkets, fruit and veg store)
fastfoodCat <- 2 #index of columns containing fast food categories

## Total, healthy, unhealthy

abund$total <- rowSums(abund[,totalOutlets])
abund$healthy <- rowSums(abund[,healthy])
abund$unhealthy <- rowSums(abund[,unhealthy])

## load and merge population data - required for calculating population weighted density
mbpop <- read.csv("./input/2016censusmeshblockcounts.csv")
abund <- merge(abund,mbpop, by.x="id", by.y="MB_CODE_2016")
dim(abund) #important to check merging worked - can be issues if data has been imported from a csv

# calculate fast food and grocery outlets per person
abund$dens_fv <- rowSums(abund[,groceryCat])/abund$Person
abund$dens_ff <- rowSums(abund[,fastfoodCat])/abund$Person

## mrfei

abund$mrfei <- (abund$unhealthy/abund$total)*100

write.csv(abund, "./output/access.csv")

###############################################################################
#  Diversity measures
##############################################################################

abund <- read.csv("./output/access.csv")
names(abund)
abund <- abund[,-1]
names(abund)

spe.pres <- apply(abund[,c(2:7)] >0, 2, sum)
spe.pres <- sort(spe.pres)
spe.relf <- 100*spe.pres/nrow(abund)
spe.relf <- round(sort(spe.relf), 1)

spec <- cbind(spe.pres = spe.pres, spe.relf = spe.relf)
write.csv(spec, "./output/odmb_presence.csv")

N0 <- rowSums(abund[,c(2:7)] > 0) 		#Species richness
H <- diversity(abund[,c(2:7)])		#Shannon Entropy
N1 <- exp(H)			#Shannon diversity number
N2 <- diversity(abund[,c(2:7)], "inv")	#Simpson diversity number
J <- H/log(N0)			#Pielou evenness
E1 <- N1/N0				#Shannon evenness (Hill's ratio)
E2 <- N2/N0				#Simpson's evenness (HIll's ratio)

div <- data.frame(N0, H, N1, N2, E1, E2, J)
dim(div)
dim(abund)
abund <- cbind(abund, div)
dim(abund)

write.csv(abund, "./output/access_div.csv")

#######
# add in nearest grocery

## nearest grocery

#Mesh block
ng <- read.csv("./output/ng_driving.csv")

abund <- merge(abund, ng, by.x="id", by.y="mb_id")

save(abund, file="./output/access_div.RData")

###########################################################
## Clustering MB SA1
############################################################

library(vegan)

(load("./output/access_div.RData"))

#Clustering SA1 
names(abund)
spe <- abund[, c(1:6)]
spe$dummy <- 1
names(spe)

#Calculate Bray-curtis dissimilarities on log-transformed data. Data is log transformed to
#minimise dominance of restaurants and cafes.
bc <- vegdist(log1p(spe))

#ward's clustering
wardbc <- hclust(bc, method="ward.D2")

#plot dendrogram, kills r sometimes, save before running
windows()
plot(wardbc)
rect.hclust(wardbc, k=8, border="red")

#plot heights
windows()
plot(wardbc$height, nrow(spe):2, type="S", main="Fusion levels - BC - Ward", ylab="k (number of clusters)",
					xlab="h (node height)", 
					ylim=c(0, 30), xlim=c(2, 20), 
					col="grey")
text(wardbc$height, nrow(spe):2, nrow(spe):2, col="red", cex=0.8)

#Cut dendrogram at defined number of groups and assign group number to each MB
k <- 7
wardg7 <- cutree(wardbc, k)
table(wardg7)
abund<- cbind(abund,data.frame(wardg7))
head(abund)

#plot dendrogram with clusters identified - hcoplot available here https://github.com/JoeyBernhardt/NumericalEcology
source("scripts/hcoplot.R")
hcoplot(wardbc, bc, k=7)

save(abund, file="./output/od_cluster.RData")
(load(file="./output/od_cluster.RData"))

#reclassify in order of total abundance
boxplot(sa1$total+1 ~ sa1$wardg5,
		log="y",
		main = "Food Outlet Type Richness",
		xlab = "Cluster",
		ylab = "# of Food Outlet Types",)

sa1$cluster <- NA

sa1$cluster[which(sa1$wardg5 == 1)] <- 5
sa1$cluster[which(sa1$wardg5 == 2)] <- 4
sa1$cluster[which(sa1$wardg5 == 3)] <- 2
sa1$cluster[which(sa1$wardg5 == 4)] <- 3
sa1$cluster[which(sa1$wardg5 == 5)] <- 1

save(abund,sa1, file="./output/od_cluster.RData")

###########################################################
## Clustering SA1 urban
############################################################

library(vegan)

sa1 <- readOGR("C:/Users/Ellie/Desktop/Demographics_paper/4_demo_analysis/output/sa1_alldat.gpkg")
sa1 <- sa1[which(sa1$urb ==1),]
dim(sa1)

#Clustering SA1 
names(sa1)
spe <- sa1[, c(3:8)]
spe$dummy <- 1
names(spe)

spe <- spe@data

#Calculate Bray-curtis dissimilarities on log-transformed data. Data is log transformed to
#minimise dominance of restaurants and cafes.
bc <- vegdist(log1p(spe))

#ward's clustering
wardbc <- hclust(bc, method="ward.D2")

#plot dendrogram, kills r sometimes
windows()
plot(wardbc)
rect.hclust(wardbc, k=6, border="red")

#plot heights
windows()
plot(wardbc$height, nrow(spe):2, type="S", main="Fusion levels - BC - Ward", ylab="k (number of clusters)",
					xlab="h (node height)", 
					ylim=c(0, 30), xlim=c(2, 20), 
					col="grey")
text(wardbc$height, nrow(spe):2, nrow(spe):2, col="red", cex=0.8)

#Cut dendrogram at defined number of groups and assign group number to each SA1
k <- 3
uwardg3 <- cutree(wardbc, k)
table(uwardg3)
sa1$uwardg3 <- uwardg3
View(head(sa1))
k <- 6
uwardg6 <- cutree(wardbc, k)
table(uwardg6)
sa1$uwardg6 <- uwardg6
head(sa1$uwardg6)

#plot dendrogram with clusters identified
source("scripts/hcoplot.R")
hcoplot(wardbc, bc, k=3)
hcoplot(wardbc, bc, k=6)


save(abund,sa1, file="./output/od_cluster.RData")
(load(file="./output/od_cluster.RData"))

#reclassify in order of total abundance
boxplot(sa1$total+1 ~ sa1$uwardg3,
		log="y",
		main = "Food Outlet Type Richness",
		xlab = "Cluster",
		ylab = "# of Food Outlet Types",)

#reclassify in order of total abundance
boxplot(sa1$total+1 ~ sa1$uwardg6,
		log="y",
		main = "Food Outlet Type Richness",
		xlab = "Cluster",
		ylab = "# of Food Outlet Types",)

sa1$ucluster3 <- NA

sa1$ucluster3[which(sa1$uwardg3 == 1)] <- 3
sa1$ucluster3[which(sa1$uwardg3 == 2)] <- 2
sa1$ucluster3[which(sa1$uwardg3 == 3)] <- 1

sa1$ucluster6 <- NA
sa1$ucluster6[which(sa1$uwardg6 == 1)] <- 6
sa1$ucluster6[which(sa1$uwardg6 == 2)] <- 3
sa1$ucluster6[which(sa1$uwardg6 == 3)] <- 5
sa1$ucluster6[which(sa1$uwardg6 == 4)] <- 4
sa1$ucluster6[which(sa1$uwardg6 == 5)] <- 2
sa1$ucluster6[which(sa1$uwardg6 == 6)] <- 1

save(sa1, file="./output/sa1_urban_cluster.RData")


#############################
# nMDS
##############################
library(vegan)
bcnmds <- metaMDS(bc, k=2, parallel=2, maxit=1000, trymax=150, sratmax=0.9999999)
bcnmds <- metaMDS(bc, k=2, parallel=2, maxit=1000, trymax=150, sratmax=0.9999999, sfgrmin = 1e-7, previous.best=bcnmds )
env.spec <- envfit(bcnmds, log1p(spe@data))
bcnmds$stress

save(bcnmds, file="./output/urban/unmds.RData")

(load(file="./output/nmds.RData"))

#plot nmds with 3 and 6 groups
cols <- c("black","orange","blue")
groups <- as.factor(sa1$ucluster3)
#ordiplot(samp.nmds,col="white")
#orditorp(samp.nmds,display="species",col="red",air=0.01)
par(mai=rep(0.1,4))
plot(bcnmds$points, col = cols[groups],pch=as.numeric(groups),xlab="",ylab="",xaxt="n",yaxt="n")
legend('topleft', col=cols, legend=levels(groups), pch = 1:7, cex = 1,inset=0.02,title="Group")
#ordihull(samp.nmds,groups=groups,draw="polygon",col="grey90",label=F)
ordiellipse(bcnmds,groups=groups, col = cols ,lty=1, lwd=2)
plot(env.spec)

cols <- c("black","orange","blue","green3","red", "purple")
groups <- as.factor(sa1$ucluster6)
#ordiplot(samp.nmds,col="white")
#orditorp(samp.nmds,display="species",col="red",air=0.01)
windows()
par(mai=rep(0.1,4))
plot(bcnmds$points, col = cols[groups],pch=as.numeric(groups),xlab="",ylab="",xaxt="n",yaxt="n")
legend('topleft', col=cols, legend=levels(groups), pch = 1:7, cex = 1,inset=0.02,title="Group")
#ordihull(samp.nmds,groups=groups,draw="polygon",col="grey90",label=F)
ordiellipse(bcnmds,groups=groups, col = cols ,lty=1, lwd=2)
plot(env.spec)

#Calculate mean abundances for each species in each category for chosen method

names(sa1)
spe <- sa1[, c(3:8)]
names(spe)

groups <-as.factor(sa1$ucluster6)
spe.means <- matrix(0,ncol(spe), length(levels(groups)))
row.names(spe.means) <- colnames(spe@data)

for(i in 1:ncol(spe)) {
	spe.means[i,] <- tapply(spe@data[,i], sa1$ucluster6, mean)
}

View(spe.means)

#order the groups and list the dominant species in each group

group1 <- round(sort(spe.means[,1], decreasing=TRUE), 2)
group2 <- round(sort(spe.means[,2], decreasing=TRUE), 2)
group3 <- round(sort(spe.means[,3], decreasing=TRUE), 2)
group4 <- round(sort(spe.means[,4], decreasing=TRUE), 2)
group5 <- round(sort(spe.means[,5], decreasing=TRUE), 2)
group6 <- round(sort(spe.means[,6], decreasing=TRUE), 2)

write.csv(group1, "./output/urban/group1.csv")
write.csv(group2, "./output/urban/group2.csv")
write.csv(group3, "./output/urban/group3.csv")
write.csv(group4, "./output/urban/group4.csv")
write.csv(group5, "./output/urban/group5.csv")
write.csv(group6, "./output/urban/group6.csv")

#calculate presence
spe.pres <- apply(spe@data >0, 2, sum)
spe.pres <- sort(spe.pres)
spe.relf <- 100*spe.pres/nrow(spe@data)
spe.relf <- round(sort(spe.relf), 1)

spec <- cbind(spe.pres = spe.pres, spe.relf = spe.relf)
write.csv(spec, "./output/urban/odsa1_presence.csv")

#run adonis and save results to R file
names(spe)
spe$dummy <- 1
adonis.res <- adonis(spe@data ~ as.factor(ucluster6), data = sa1@data)
adonis.res # note down R2 and P values
save(adonis.res, file= "./output/urban/od_adonis.RData")

######################################
# Mapping and figures
#####################################

options(scipen = 999)
library(sp)
library(rgeos)
library(rgdal)
library(sf)
library(ggspatial)
library(rnaturalearth)

sa4 <- "./input/ABS_SA4/SA4_2016_AUST.shp"
aust = "./input/coastline.gpkg"

sa4 <- readOGR(sa4)
sa <- sa4[which(sa4$GCC_NAME16 == "Greater Perth"),]

#project data
sa <- spTransform(sa, "+init=epsg:28350")
ss_centroids <- gCentroid(sa, byid=T)
ss_centroids <- spTransform(ss_centroids, CRS("+init=epsg:4326"))
ss_centroids <- as.data.frame(ss_centroids@coords)
row.names(ss_centroids) <- c(1:6)
ss_centroids$labels <- sa@data$SA4_NAME16
coordinates(ss_centroids) <- ~x+y
sa_unp <- spTransform(sa, CRS("+init=epsg:4326"))
proj4string(ss_centroids) <- proj4string(sa_unp)

sa_sf <- st_as_sf(sa_unp)

#move labels
ss_centroids@coords[1,1] <- 115.677338
ss_centroids@coords[1,2] <- -31.964003
ss_centroids@coords[4,1] <- 116.138043
ss_centroids@coords[4,2] <- -32.184545

ss_cent_sf <- st_as_sf(ss_centroids)

aoi <- sa
aust <- ne_countries(country = "australia")
b <- bbox(aust)
aoi_sf <- st_as_sf(aoi)

#map of clusters
sa1poly <- readOGR("./input/1270055001_sa1_2016_aust_shape/SA1_2016_AUST.shp")
perth <- readOGR("./input/perth.gpkg")
oddata <- sa1

dim(sa1poly)
dim(oddata)
alldatasp <- merge(sa1poly, oddata, by.x="SA1_MAIN16", by.y="id", all.x=F)
dim(alldatasp)
alldatasp$cluster <- as.factor(alldatasp$cluster)
alldata_sf <- st_as_sf(alldatasp)

alldatasp <- sa1
alldatasp$ucluster3 <- as.factor(alldatasp$ucluster3)
alldata_sf <- st_as_sf(alldatasp)

str(alldatasp@data)

ss_cent_df <- as.data.frame(ss_centroids)

library(viridis)

od <- ggplot(ss_cent_sf) +
	layer_spatial(alldata_sf, aes(fill=ucluster3), color=NA)+ 
	scale_fill_viridis(discrete=T,alpha=0.8, option="magma", name="Cluster")+
	layer_spatial(sa_sf, fill=NA, col="grey50", lty=3, size=1) +
  	geom_sf(col="red")+ 
	#geom_sf_label(aes(label=labels))+
	annotate("text", label = "Open Data", x = 115.5, y = -31.47, vjust=-1.5, size = 6, colour="black", fontface="bold")+
	annotate("text", label = "Rockingham", x = 115.63, y = -32.28, size = 3, colour="black", fontface="bold")+
	annotate("text", label = "Armadale", x = 116.1, y = -32.15, size = 3, colour="black", fontface="bold")+
	annotate("text", label = "Midland", x = 116.01, y = -31.86, size = 3, colour="black", fontface="bold")+
	annotate("text", label = "Fremantle", x = 115.67, y = -32.05, size = 3, colour="black", fontface="bold")+
	annotate("text", label = "Joondalup", x = 115.65, y = -31.72, size = 3, colour="black", fontface="bold")+
	annotate("text", label = "Perth", x = 115.84, y = -31.945, size = 3, colour="white", fontface="bold")+
	annotation_scale(location = "br", width_hint = 0.3, style="ticks") +
  	annotation_north_arrow(location = "bl", which_north = "true")+
	theme(axis.title = element_blank(),
		panel.background = element_rect(fill = "white"),
		#panel.grid.major = element_line(color = "grey86", linetype = 3, size = 0.3),
		panel.border = element_rect(colour = "black", fill=NA),
		legend.position = c(.1, .8)
	) 
	#annotate("text",x=115.75,y=-32.5,vjust=3.5,label="ASGS MB & SA4 Boundaries Source: ABS 2016 (CC BY 4.0)",
		#size=2.5)

plot(od)

### Boxplots

#round those less than 50m to 50m

sa1$mean_ng[which(sa1$mean_ng < 50)] <- 50

cols <- c("grey23","orange","blue","green3","red")

par(mfrow=c(1,3), cex=1.1, cex.main=1.2, pch=20)
boxplot(sa1$total+1 ~ sa1$cluster,
		col = cols,
		log = "y",
		main = "Food Outlet Abundance",
		#xlab = "Cluster",
		ylab = "# of Food Outlets",
		range = 0)
boxplot(sa1$N0 ~ sa1$cluster,
		col = cols,
		main = "Food Outlet Type Richness",
		#xlab = "Cluster",
		ylab = "# of Food Outlet Types",
		range = 0)
boxplot(sa1$N2 ~ sa1$cluster,
		col = cols,
		main = "Simpson Diversity",
		#xlab = "Cluster",
		ylab = "Simpson Diversity",
		range = 0)

par(mfrow=c(2,2), cex=1.1, cex.main=0.8, pch=20, cex.axis=0.8)
boxplot(sa1$densff ~ sa1$cluster,
		col = cols,
		#log = "y",
		main = "Pop. Weighted Density Fast-Food",
		#xlab = "Cluster",
		#yaxt = "n",		
		ylab = "# per Person",
		range = 0)
#axis(side=2,at=c(0.001,0.01,0.1,1), labels=c("0.001","0.01","0.1","1"))
boxplot(sa1$densfv ~ sa1$cluster,
		col = cols,
		#log = "y",
		main = "Pop. Weighted Density Grocery Stores",
		#xlab = "Cluster",
		#yaxt = "n",
		ylab = "# per Person",
		range = 0)
#axis(side=2,at=c(0.001,0.01,0.1,1), labels=c("0.001","0.01","0.1","1"))
boxplot(sa1$mrfei ~ sa1$cluster,
		col = cols,
		main = "MRFEI",
		#ylim = c(30,100),
		#xlab = "Cluster",
		ylab = "% less healthy",
		range = 0)
boxplot(sa1$mean_ng ~ sa1$cluster,
		col = cols,
		#log = "y",
		ylim = c(50,6000),
		main = "Distance to Nearest Grocery (m)",
		#xlab = "Cluster",
		ylab = "Meters",
		range = 0)

cols <- c("grey23","orange","blue","green3","red", "purple")

par(mfrow=c(1,3), cex=1.1, cex.main=1.2, pch=20)
boxplot(sa1$total+1 ~ sa1$ucluster6,
		col = cols,
		log = "y",
		main = "Food Outlet Abundance",
		#xlab = "Cluster",
		ylab = "# of Food Outlets",
		range = 0)
boxplot(sa1$N0 ~ sa1$ucluster6,
		col = cols,
		main = "Food Outlet Type Richness",
		#xlab = "Cluster",
		ylab = "# of Food Outlet Types",
		range = 0)
boxplot(sa1$N2 ~ sa1$ucluster6,
		col = cols,
		main = "Simpson Diversity",
		#xlab = "Cluster",
		ylab = "Simpson Diversity",
		range = 0)

par(mfrow=c(2,2), cex=1.1, cex.main=0.8, pch=20, cex.axis=0.8)
boxplot(sa1$pden_ff ~ sa1$ucluster6,
		col = cols,
		#log = "y",
		main = "Pop. Weighted Density Fast-Food",
		ylim = c(0,15),
		#xlab = "Cluster",
		#yaxt = "n",		
		ylab = "Density",
		range = 0)
#axis(side=2,at=c(0.001,0.01,0.1,1), labels=c("0.001","0.01","0.1","1"))
boxplot(sa1$pden_fv ~ sa1$ucluster6,
		col = cols,
		#log = "y",
		main = "Pop. Weighted Density Grocery Stores",
		ylim = c(0,15),
		#xlab = "Cluster",
		#yaxt = "n",
		ylab = "Density",
		range = 0)
#axis(side=2,at=c(0.001,0.01,0.1,1), labels=c("0.001","0.01","0.1","1"))
boxplot(sa1$mrfei ~ sa1$ucluster6,
		col = cols,
		main = "MRFEI",
		#ylim = c(30,100),
		#xlab = "Cluster",
		ylab = "% less healthy",
		range = 0)
boxplot(sa1$mean_ng ~ sa1$ucluster6,
		col = cols,
		#log = "y",
		ylim = c(50,6000),
		main = "Distance to Nearest Grocery (m)",
		#xlab = "Cluster",
		ylab = "Meters",
		range = 0)


