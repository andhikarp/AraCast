#Resampling Exposito-Alonso (2020)'s tif files

#!/usr/bin/env Rscript
library(raster)
library(viridis)
b=raster("germinationpc2raw.tiff")
tmin=raster("~/Documents/multitab_aracast/tmin_reduced.nc")
crs(b)<-crs(tmin)
b_crop=crop(b,tmin)
extent(b)<-extent(tmin[[1]])
boints=rasterToPoints(b)
crs<-CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") #obtained by running proj4string(map)
spdf<-SpatialPointsDataFrame(coords=as.data.frame(cbind(boints[,1],boints[,2])),data=as.data.frame(boints[,3]),proj4string=crs)
grd <- as.data.frame(spsample(spdf, cellsize=c(1.035821,1.008889),type="regular",fn=mean,pretty=TRUE,offset=c(0.5,0.5))) #samples in center of cell
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object
proj4string(grd) <- proj4string(tmin) # grid is 3015
v=grd 
grd=raster(grd)
resampled_raster=resample(b_crop,grd,"ngb")
library(automap)

boints2=rasterToPoints(resampled_raster)

crs2=CRS("+proj=robin")
gridded(v)     <- TRUE  # Create SpatialPixel object
fullgrid(v)    <- TRUE  # Create SpatialGrid object
proj4string(v) <- proj4string(tmin) # grid is 3015
crs(v)<-crs2

ag_spdf<-SpatialPointsDataFrame(coords=data.frame(cbind(boints2[,1],boints2[,2])),data=data.frame(boints2[,3]),proj4string=crs2)
gsm.kriged<-autoKrige(boints2...3.~1,input_data=ag_spdf,new_data=v)

dtharvestraw=tmin[[1]]
values(dtharvestraw)<-gsm.kriged$krige_output$var1.pred
dtharvestraw=mask(dtharvestraw,tmin[[1]])


germinationharvestraw=tmin[[1]]
values(germinationharvestraw)<-gsm.kriged$krige_output$var1.pred
germinationharvestraw=mask(germinationharvestraw,tmin[[1]])


germinationharvestraw2=tmin[[1]]
values(germinationharvestraw2)<-gsm.kriged$krige_output$var1.pred
germinationharvestraw2=mask(germinationharvestraw2,tmin[[1]])


#Clustering
library(caret)
vals=cbind(values(dtharvestraw),values(germinationharvestraw),values(germinationharvestraw2))
vals=scale(vals)
clus=kmeans(na.omit(vals),4,10)

clusmap=tmin[[1]]
clusmap[!is.na(clusmap)]<-clus$cluster
plot(clusmap,col=viridis(4))


clus1=clusmap
clus1[clus1!=1]<-NA
clus2=clusmap
clus2[clus2!=2]<-NA 
clus3=clusmap
clus3[clus3!=3]<-NA
clus4=clusmap
clus4[clus4!=4]<-NA

clusmaps=stack(clus1,clus2,clus3,clus4)
writeRaster(clusmaps,'lifehistory_clusters.nc')

#replacing outlier life history cells
lh=clusmaps
lhm=merge(lh)

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

for (i in seq(1,3015)){
	x=values(lhm)[i]
	if (is.na(x))
	{
	}else{
		nbs=values(lhm)[adjacent(lhm,i,8,pairs=FALSE)] #values of neighboring cells
		bb=sum(nbs==x,na.rm=TRUE)
		if (bb>0)
		{} else{
			z=getmode(nbs)
			values(lhm)[i]<-z

		}
	}
}

library(rgdal)
world=readOGR("~/Clim_GWAS/Clim_GWAS_2/ne_110m_admin_0_countries.shp")
europe=world[world$CONTINENT=='Europe',]
lhm=mask(lhm,europe)

lhm1=lhm
lhm1[lhm1!=1]<-NA
lhm2=lhm
lhm2[lhm2!=2]<-NA 
lhm3=lhm
lhm3[lhm3!=3]<-NA
lhm4=lhm
lhm4[lhm4!=4]<-NA
lhm_fin=stack(lhm1,lhm2,lhm3,lhm4)
writeRaster(lhm_fin,"~/Clim_GWAS/Clim_GWAS_2/lifehistory_clusters2.nc",overwrite=TRUE)