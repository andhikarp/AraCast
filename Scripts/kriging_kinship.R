#!/usr/bin/env Rscript
#Kriging genomic similarity parallelized

args = commandArgs(trailingOnly = TRUE) #arguments are number of cores available and nth job.

library(RcppCNPy)
library(raster)
library(rgdal)
library(tictoc)
library(tidyverse)
library(gstat)
library(stringr)



library(sp)
library(dplyr) 
library(ggplot2)
library(scales) 
library(magrittr)
library(tmap)
library(maptools)
library(spatstat)
library(automap)

setwd("~/Clim_GWAS/Clim_GWAS_2")


gsm=npyLoad("1001genomes_k01wg.npy") #2029 * 5623 GSM , whole genome for SNPs trimmed at 0.1 pairwise correlation threshold
data=read.csv("k2029_accessions_coord_tidied.csv") #2171 * 5 , information on coordinates for each accession


#Matching GSM entries to accession coordinates

idx=match(unique(data$ecotype_id),data$ecotype_id) 
	data=data[idx,]
filter=is.na(data$longitude) #logical for data without coordinates

lats=data$latitude
lons=data$longitude
coords<-cbind(lons,lats)
coords=coords[!filter,] #subset of data with coordinates

crs=CRS("+proj=robin") # coordinate reference system


setwd("~/Clim_GWAS/E-OBSv19.0HOM")
map=raster("tn_ens_mean_0.1deg_reg_v19.0eHOM.nc")

#Europe outline

# a <- map[[1]] > -Inf
# pp <- rasterToPolygons(a, dissolve=TRUE) #this step takes a while....
# shapefile(pp,"europe_shape.shp") #saves the polygon

pp <- readOGR("europe_shape.shp")

land  <-readOGR('~/Clim_GWAS/Clim_GWAS_2/ne_110m_land.shp')
proj4string(land)<-proj4string(a) #define extent of projections
land@bbox <- pp@bbox 




indices=seq(1:ncol(gsm))
chunk <- function(x, n) split(x, sort(rank(x) %% n))
ncores<-as.numeric(args[1])
chunks=as.data.frame(chunk(indices,ncores)[as.numeric(args[2])])
print(paste("Chunk no.:",args[2]))


grd <- as.data.frame(spsample(pp, cellsize=c(1),type="regular",fn=mean,pretty=TRUE,offset=c(0.5,0.5))) #samples in center of cell

# Kriging each K 
tic()
fin_pred=matrix(NA,nrow=3015)
fin_var=matrix(NA,nrow=3015)
fin_stdev=matrix(NA,nrow=3015)

for (i in rownames(chunks)){
	tic()
	i=as.numeric(chunks[i,1])
	print(paste("Kriging GSM entry no.",i))

	datas<-gsm[,i]
	datas<-datas[!filter]
	datas=as.data.frame(datas)
	spdf<-SpatialPointsDataFrame(coords=coords,data=datas,proj4string=crs)

	
	coords_binded=paste0(spdf@coords[,1],"_",spdf@coords[,2])
	ag=aggregate(spdf$datas~coords_binded,FUN=function(x)c(mean=mean(x),count=length(x))) # aggregate for the unique coordinates
	ag_vals=matrix(NA,ncol=3)
	

	#split the aggregated coordinates

	for (j in (1:nrow(ag))){
		split_coords=str_split_fixed(ag[j,1],pattern="_",n=2)
		coords1=as.numeric(split_coords[1])
		coords2=as.numeric(split_coords[2])
		mean_val=as.numeric(ag$spdf[j])
		combined=cbind(coords1,coords2,mean_val)
		ag_vals=rbind(ag_vals,combined)
	}

	ag_vals=ag_vals[-1,]
	ag_coords=cbind(ag_vals[,1],ag_vals[,2])
	ag_spdf<-SpatialPointsDataFrame(coords=data.frame(ag_coords),data=data.frame(ag_vals[,3]),proj4string=crs) #create aggregated spdf for kriging

	v=grd
	names(v)       <- c("X", "Y")
	coordinates(v) <- c("X", "Y")
	gridded(v)     <- TRUE  # Create SpatialPixel object
	fullgrid(v)    <- TRUE  # Create SpatialGrid object
	proj4string(v) <- proj4string(pp) # grid is 3015
	ag_spdf@bbox<-v@bbox
	crs(v)<-crs
	gsm.kriged<-autoKrige(ag_vals...3.~1,input_data=ag_spdf,new_data=v)
	prediction_spdf = gsm.kriged$krige_output
	pred_vals=prediction_spdf$var1.pred
	var_vals=prediction_spdf$var1.var
	stdev_vals=prediction_spdf$var1.stdev

	fin_pred=cbind(fin_pred,pred_vals)
	fin_var=cbind(fin_var,var_vals)
	fin_stdev=cbind(fin_stdev,stdev_vals)
	# sample_variogram = gsm.kriged$exp_var
	# variogram_model = gsm.kriged$var_model
	toc()
	
}

fin_pred=fin_pred[,-1]
fin_var=fin_var[,-1]
fin_stdev=fin_stdev[,-1]
toc()


#Partial results

setwd("~/Clim_GWAS/Clim_GWAS_2/Temp_Files")
write.csv(fin_pred,paste0('krigval',args[2],'.csv'))
write.csv(fin_var,paste0('krigvar',args[2],'.csv'))
write.csv(fin_var,paste0('krigstdev',args[2],'.csv'))

#Then load and concatenate them all together into a single dataframe