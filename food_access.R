###############################################
# Food access clusters and metrics for profiling
# Eleanor Cervigni
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

API.key <- "Your API Key"
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

save(abund, file="./output/cluster.RData")
(load(file="./output/cluster.RData"))

#check boxplot and reclassify in order of total abundance
boxplot(abund$total+1 ~ abund$wardg7,
		log="y",
		main = "Food Outlet Abundance",
		xlab = "Cluster",
		ylab = "# of Food Outlets",)

abund$cluster <- NA

abund$cluster[which(abund$wardg7 == 1)] <- 5
abund$cluster[which(abund$wardg7 == 2)] <- 4
abund$cluster[which(abund$wardg7 == 3)] <- 2
abund$cluster[which(abund$wardg7 == 4)] <- 3
abund$cluster[which(abund$wardg7 == 5)] <- 1
abund$cluster[which(abund$wardg7 == 6)] <- 6
abund$cluster[which(abund$wardg7 == 7)] <- 7

save(abund,file="./output/cluster.RData")


#############################
# nMDS
##############################
library(vegan)

#run nmds
bcnmds <- metaMDS(bc, k=2, parallel=2, maxit=1000, trymax=150, sratmax=0.9999999)

#fit abundance data
env.spec <- envfit(bcnmds, log1p(spe@data))

save(bcnmds, file="./output/nmds.RData")

(load(file="./output/nmds.RData"))

#plot nmds with points coloured according to assigned clusters
cols <- c("black","orange","blue", "purple", "green", "red", "yellow")
groups <- as.factor(abund$ucluster)
par(mai=rep(0.1,4))
plot(bcnmds$points, col = cols[groups],pch=as.numeric(groups),xlab="",ylab="",xaxt="n",yaxt="n")
legend('topleft', col=cols, legend=levels(groups), pch = 1:7, cex = 1,inset=0.02,title="Group")
ordiellipse(bcnmds,groups=groups, col = cols ,lty=1, lwd=2)
#add abundance vectors
plot(env.spec)

#Calculate mean abundances for each food outlet type in each cluster

names(abund)
spe <- abund[, c(3:8)]
names(spe)

groups <-as.factor(abund$ucluster7)
spe.means <- matrix(0,ncol(spe), length(levels(groups)))
row.names(spe.means) <- colnames(spe@data)

for(i in 1:ncol(spe)) {
	spe.means[i,] <- tapply(spe@data[,i], abund$cluster7, mean)
}

View(spe.means)

#order the groups and list the dominant food outlet type in each group

group1 <- round(sort(spe.means[,1], decreasing=TRUE), 2)
group2 <- round(sort(spe.means[,2], decreasing=TRUE), 2)
group3 <- round(sort(spe.means[,3], decreasing=TRUE), 2)
group4 <- round(sort(spe.means[,4], decreasing=TRUE), 2)
group5 <- round(sort(spe.means[,5], decreasing=TRUE), 2)
group6 <- round(sort(spe.means[,6], decreasing=TRUE), 2)

write.csv(group1, "./output/group1.csv")
write.csv(group2, "./output/group2.csv")
write.csv(group3, "./output/group3.csv")
write.csv(group4, "./output/group4.csv")
write.csv(group5, "./output/group5.csv")
write.csv(group6, "./output/group6.csv")
write.csv(group7, "./output/group7.csv")

#calculate presence
spe.pres <- apply(spe@data >0, 2, sum)
spe.pres <- sort(spe.pres)
spe.relf <- 100*spe.pres/nrow(spe@data)
spe.relf <- round(sort(spe.relf), 1)

spec <- cbind(spe.pres = spe.pres, spe.relf = spe.relf)
write.csv(spec, "./output/presence.csv")

#run adonis and save results to R file
names(spe)
spe$dummy <- 1
adonis.res <- adonis(spe ~ as.factor(cluster), data = abund)
adonis.res # note down R2 and P values
save(adonis.res, file= "./output/adonis_clust.RData")

