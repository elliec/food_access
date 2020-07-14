# Food Access Clustering #

These workflow in these scripts identifies clusters of distinct food access based on a given good outlet dataset.

_Inputs required_
* food outlets point data
* points or polygons covering geographic area of interest at desired resolution. These could be eg census areas or a simple grid of points

_Inputs Optional_
* population data if calculating #stores/person

Steps

1. Calculate isochrones around points (isochrones.R)
	1. Calculate centroids if geographic area is polygons
	1. Calculate isochrones using [OpenRouteService API](https://openrouteservice.org/)

1. Identify clusters and profile them based food outlets, accessibility and diversity metrics (food_access.R)
	1. Calculate distance to the nearest grocery store from each centroid
	1. Calculate MRFEI, density of fast food and grocery stores/person for each isochrone
	1. Calculate Simpson's diversity, total abundance, and richness of food outlets for each isochrone
	1. Identify clusters based on food outlet type and abundance
	1. Verify clusters with non-metric multidimensional scaling and permanova (adonis)

This is a workflow - check and save outputs as you go

For more information on clustering and multivariate methods see [Numerical Ecology with R](https://www.springer.com/gp/book/9783319714035)
