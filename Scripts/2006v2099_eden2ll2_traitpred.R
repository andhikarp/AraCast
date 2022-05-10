# Comparing DTB and SN for Eden-2 and Ll-2


edenid=933
ll2id=1124

library(RcppCNPy)
library(raster)
library(viridis)
library(rgdal)

setwd("~/Clim_GWAS/Clim_GWAS_2")

set.seed(12)



tmin26=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP2.6/tasmin_rcp26_20410101-21001231.nc")
tmax26=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP2.6/tasmax_rcp26_20410101-21001231.nc")
# tmin45=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP4.5/tasmin_rcp45_20410101-21001231.nc")
# tmax45=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP4.5/tasmax_rcp45_20410101-21001231.nc")
# tmin60=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP6.0/tasmin_rcp60_20410101-21001231.nc")
# tmax60=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP6.0/tasmax_rcp60_20410101-21001231.nc")
tmin85=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP8.5/tasmin_rcp85_20410101-21001231.nc")
tmax85=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP8.5/tasmax_rcp85_20410101-21001231.nc")

germinationdates=raster('pred_germdate.nc')
germdates=as.data.frame(round(values(germinationdates)))
coords=as.data.frame(coordinates(germinationdates))

# Shapes and present day temperature data
	setwd("~/Clim_GWAS/E-OBSv19.0HOM")
	tmin=stack("tmin_reduced.nc") #2000-2010 data maximum temperature
	tmax=stack("tmax_reduced.nc") #2000-2010 data minimum temperature
	world=readOGR("~/Clim_GWAS/Clim_GWAS_2/ne_110m_admin_0_countries.shp")
	europe=world[world$CONTINENT=='Europe',]

#
accessiondata=read.csv("~/Clim_GWAS/Clim_GWAS_2/k2029_accessioninfo.csv")
edencoords=accessiondata[edenid,c('longitude','latitude')]
llcoords=accessiondata[ll2id,c('longitude','latitude')]

#ADMIXTURE assignments 

setwd("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/ADMIXTURE/admixture_linux-1.3.0")
admixture=read.csv("k2029_trimmed.4.Q",header=F,sep=' ')
admixture=cbind(rep(1,2029),admixture)
admixture_eden=admixture[edenid,]
admixture_ll2=admixture[ll2id,]

admix=read.csv("~/Clim_GWAS/Clim_GWAS_2/krigclus4.csv") #kriged admixture proportions
admix=admix[,-1]

K_values=npyLoad('kriged_valuesnoNA.npy') # Kinship matrix applied to the landscape

meanIgnoringZeroes <- function(x) {
  mean(x[x!=0],na.rm=T)
}


dtb_intercept=-0.578511150335065
sn_intercept=29.54953719927532

#intercepts were obtained from the fitted G x E model


# Load LMMLASSO model weights for seed number
	setwd("~/Clim_GWAS/LMMLASSO/Weights_Corrected")
	snw_ridge=npyLoad('SNBASE2Weightsridge_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
	snw_ridge=as.matrix(snw_ridge,ncol=1)
	snbetas=npyLoad('SNBASE2Weights_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
	snbetas=as.matrix(snbetas,ncol=1)




setwd("~/Clim_GWAS/E-OBSv19.0HOM")

	tmin=stack("tmin_reduced.nc")
	tmax=stack("tmax_reduced.nc")
	world=readOGR("~/Clim_GWAS/Clim_GWAS_2/ne_110m_admin_0_countries.shp")
	europe=world[world$CONTINENT=='Europe',]

setwd("~/Clim_GWAS/Clim_GWAS_2/")
	k2029=npyLoad('1001genomes_k01wg.npy')
	k1=k2029[edenid,]
	k2=k2029[ll2id,]



k1=rep(k1,times=3015)
k1=matrix(k1,nrow=3015,byrow=TRUE)
k2=rep(k2,times=3015)
k2=matrix(k2,nrow=3015,byrow=TRUE)
# k3=rep(k3,times=3015)
# k3=matrix(k3,nrow=3015,byrow=TRUE)


#### Eden-2 SN predictions

#### Eden-2 seed number in 2006

## Get environmental predictors for 2006

	tmin2006=c()
		for (i in (1:nrow(germdates))){
			print(i)
			germdate=germdates[i,]
			if(is.na(germdate)==TRUE){
					b=rep(NA,203)
					tmin2006=rbind(tmin2006,b)
				}
				else{
			coordinate=coords[i,]
			a=tmin[[(germdate+(365*6)):((germdate+365*6)+202)]]
			b=extract(a,coordinate)
			tmin2006=rbind(tmin2006,b)
		}
		}

		tmax2006=c()
		for (i in (1:nrow(germdates))){
			print(i)
			germdate=germdates[i,]
			if(is.na(germdate)==TRUE){
					b=rep(NA,203)
					tmax2006=rbind(tmax2006,b)
				}
				else{
			coordinate=coords[i,]
			a=tmax[[(germdate+(365*6)):((germdate+365*6)+202)]]
			b=extract(a,coordinate)
			tmax2006=rbind(tmax2006,b)
		}
		}

	predictors_2006_eden=cbind(tmin2006,tmax2006)
	predictors_2006_eden=cbind((predictors_2006_eden*1),(predictors_2006_eden*admixture_eden[,1]),(predictors_2006_eden*admixture_eden[,2]),(predictors_2006_eden*admixture_eden[,3]),(predictors_2006_eden*admixture_eden[,4]))

	sn_2006_eden=(predictors_2006_eden%*%snbetas) + (k1 %*% snw_ridge) + sn_intercept
	sn_2006_eden[sn_2006_eden<0]=0

################################################ RCP 2.6 #################################
##########################################################################################
##########################################################################################


		############# 2099 predictions RCP2.6 ############################################

			# Get the environmental predictors for 2099 RCP2.6

			tmin2099rcp26=c()
			for (i in (1:nrow(germdates))){
				print(i)
				germdate=germdates[i,]
				if(is.na(germdate)==TRUE){
						b=rep(NA,203)
						tmin2099rcp26=rbind(tmin2099rcp26,b)
					}
					else{
				coordinate=coords[i,]
				a=tmin26[[(germdate+(21184)):((germdate+21184)+202)]]
				b=extract(a,coordinate)
				b=b-273.15 #Kelvin to Celsius
				tmin2099rcp26=rbind(tmin2099rcp26,b)
			}
			}

			tmax2099rcp26=c()
			for (i in (1:nrow(germdates))){
				print(i)
				germdate=germdates[i,]
				if(is.na(germdate)==TRUE){
						b=rep(NA,203)
						tmax2099rcp26=rbind(tmax2099rcp26,b)
					}
					else{
				coordinate=coords[i,]
				a=tmax26[[(germdate+(21184)):((germdate+21184)+202)]]
				b=extract(a,coordinate)
				b=b-273.15 #Kelvin to Celsius
				tmax2099rcp26=rbind(tmax2099rcp26,b)
			}
			}

			# get predicted seed number for the baseline, eden-2, and ll-2

			# BASELINE
			
			predictors_2099rcp26=cbind(tmin2099rcp26,tmax2099rcp26)
			baseline_predictors=cbind((predictors_2099rcp26*1),(predictors_2099rcp26*admix[,1]),(predictors_2099rcp26*admix[,2]),(predictors_2099rcp26*admix[,3]),(predictors_2099rcp26*admix[,4]))
			sn_2099rcp26=(baseline_predictors%*%snbetas) + (K_values %*% snw_ridge) + sn_intercept
			sn_2099rcp26[sn_2099rcp26<0]=0

			# EDEN-2 

			eden_predictors=cbind((predictors_2099rcp26*1),(predictors_2099rcp26*admixture_eden[,2]),(predictors_2099rcp26*admixture_eden[,3]),(predictors_2099rcp26*admixture_eden[,4]),(predictors_2099rcp26*admixture_eden[,5]))
			sn_2099rcp26_eden=(eden_predictors%*%snbetas) + (k1 %*% snw_ridge) + sn_intercept
			sn_2099rcp26_eden[sn_2099rcp26_eden<0]=0

			# LL-2

			ll2_predictors=cbind((predictors_2099rcp26*1),(predictors_2099rcp26*admixture_ll2[,2]),(predictors_2099rcp26*admixture_ll2[,3]),(predictors_2099rcp26*admixture_ll2[,4]),(predictors_2099rcp26*admixture_ll2[,5]))
			sn_2099rcp26_ll2=(ll2_predictors%*%snbetas) + (k2 %*% snw_ridge) + sn_intercept
			sn_2099rcp26_ll2[sn_2099rcp26_ll2<0]=0



			# compare differences

			# EDEN-2 VS BASELINE

			eden_baseline = sn_2099rcp26_eden - sn_2099rcp26
			eden_baseline_raster = tmin[[1]]
			values(eden_baseline_raster) = eden_baseline
			eden_baseline_raster=disaggregate(eden_baseline_raster,fact=25)
			eden_baseline_raster=mask(eden_baseline_raster,europe)

			# LL-2 VS BASELINE

			ll2_baseline = sn_2099rcp26_ll2 - sn_2099rcp26
			ll2_baseline_raster= tmin[[1]]
			values(ll2_baseline_raster) = ll2_baseline
			ll2_baseline_raster = disaggregate(ll2_baseline_raster,fact=25)
			ll2_baseline_raster = mask(ll2_baseline_raster,europe)

			# EDEN-2 - LL-2
			ll2_eden = sn_2099rcp26_eden - sn_2099rcp26_ll2
			ll2_eden_raster = tmin[[1]]
			values(ll2_eden_raster) = ll2_eden
			ll2_eden_raster = disaggregate(ll2_eden_raster,fact=25)
			ll2_eden_raster = mask(ll2_eden_raster,europe)

			# making comparison plots

				# color palette

				pos_colors=rev(viridis(100))
				neg_colors=(magma(50))
				sn_negbreaks=seq(-20000,0,length.out=33)
				sn_posbreaks=seq(0,40000,length.out=67)
				sn_bps=c(sn_negbreaks,sn_posbreaks[-1])
				sn_col_list=c(neg_colors[1:32],'#FFFFFF',pos_colors[1:66])

				# life history outlines
				shapes=readOGR("ashape3.shp") #life history borders




			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_eden_baseline_rcp26.svg")
			plot(eden_baseline_raster,col=sn_col_list,breaks=sn_bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="EDEN-2 - BASELINE; 2099; RCP 2.6")
			plot(eden_baseline_raster,legend.only=TRUE,col=sn_col_list,breaks=sn_bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-40000,-20000,0,20000,40000),labels=c("-40","-20",0,"+20","+40")),legend.args=list(text=expression(paste(Delta,"SN x 1000")),line=1.3))
			plot(shapes,border=c("#3DDC97","#A67DB8","#0D1821","#3F612D"),col="transparent",lwd=4,add=T)
			dev.off()


			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_ll2_baseline_rcp26.svg")
			plot(ll2_baseline_raster,col=sn_col_list,breaks=sn_bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="LL-2 - BASELINE; 2099; RCP 2.6")
			plot(ll2_baseline_raster,legend.only=TRUE,col=sn_col_list,breaks=sn_bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-20000,0,20000,40000),labels=c("-20",0,"+20","+40")),legend.args=list(text=expression(paste(Delta,"SN x 1000")),line=1.3))
						plot(shapes,border=c("#3DDC97","#A67DB8","#0D1821","#3F612D"),col="transparent",lwd=4,add=T)
			dev.off()


			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_ll2_eden_rcp26.svg")
			plot(ll2_eden_raster,col=sn_col_list,breaks=sn_bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="LL-2 - EDEN-2; 2099; RCP 2.6")
			plot(ll2_eden_raster,legend.only=TRUE,col=sn_col_list,breaks=sn_bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-20000,0,20000,40000),labels=c("-20",0,"+20","+40")),legend.args=list(text=expression(paste(Delta,"SN x 1000")),line=1.3))
						plot(shapes,border=c("#3DDC97","#A67DB8","#0D1821","#3F612D"),col="transparent",lwd=4,add=T)
			dev.off()


################################################ RCP 8.5 #################################
##########################################################################################
##########################################################################################


		############# 2099 predictions RCP8.5 ############################################

			# Get the environmental predictors for 2099 RCP8.5

			tmin2099rcp85=c()
			for (i in (1:nrow(germdates))){
				print(i)
				germdate=germdates[i,]
				if(is.na(germdate)==TRUE){
						b=rep(NA,203)
						tmin2099rcp85=rbind(tmin2099rcp85,b)
					}
					else{
				coordinate=coords[i,]
				a=tmin85[[(germdate+(21184)):((germdate+21184)+202)]]
				b=extract(a,coordinate)
				b=b-273.15 #Kelvin to Celsius
				tmin2099rcp85=rbind(tmin2099rcp85,b)
			}
			}

			tmax2099rcp85=c()
			for (i in (1:nrow(germdates))){
				print(i)
				germdate=germdates[i,]
				if(is.na(germdate)==TRUE){
						b=rep(NA,203)
						tmax2099rcp85=rbind(tmax2099rcp85,b)
					}
					else{
				coordinate=coords[i,]
				a=tmax85[[(germdate+(21184)):((germdate+21184)+202)]]
				b=extract(a,coordinate)
				b=b-273.15 #Kelvin to Celsius
				tmax2099rcp85=rbind(tmax2099rcp85,b)
			}
			}

			# get predicted seed number for the baseline, eden-2, and ll-2

			# BASELINE
			
			predictors_2099rcp85=cbind(tmin2099rcp85,tmax2099rcp85)
			baseline_predictors=cbind((predictors_2099rcp85*1),(predictors_2099rcp85*admix[,1]),(predictors_2099rcp85*admix[,2]),(predictors_2099rcp85*admix[,3]),(predictors_2099rcp85*admix[,4]))
			sn_2099rcp85=(baseline_predictors%*%snbetas) + (K_values %*% snw_ridge) + sn_intercept
			sn_2099rcp85[sn_2099rcp85<0]=0

			# EDEN-2 

			eden_predictors=cbind((predictors_2099rcp85*1),(predictors_2099rcp85*admixture_eden[,2]),(predictors_2099rcp85*admixture_eden[,3]),(predictors_2099rcp85*admixture_eden[,4]),(predictors_2099rcp85*admixture_eden[,5]))
			sn_2099rcp85_eden=(eden_predictors%*%snbetas) + (k1 %*% snw_ridge) + sn_intercept
			sn_2099rcp85_eden[sn_2099rcp85_eden<0]=0

			# LL-2

			ll2_predictors=cbind((predictors_2099rcp85*1),(predictors_2099rcp85*admixture_ll2[,2]),(predictors_2099rcp85*admixture_ll2[,3]),(predictors_2099rcp85*admixture_ll2[,4]),(predictors_2099rcp85*admixture_ll2[,5]))
			sn_2099rcp85_ll2=(ll2_predictors%*%snbetas) + (k2 %*% snw_ridge) + sn_intercept
			sn_2099rcp85_ll2[sn_2099rcp85_ll2<0]=0



			# compare differences

			# EDEN-2 VS BASELINE

			eden_baseline = sn_2099rcp85_eden - sn_2099rcp85
			eden_baseline_raster = tmin[[1]]
			values(eden_baseline_raster) = eden_baseline
			eden_baseline_raster=disaggregate(eden_baseline_raster,fact=25)
			eden_baseline_raster=mask(eden_baseline_raster,europe)

			# LL-2 VS BASELINE

			ll2_baseline = sn_2099rcp85_ll2 - sn_2099rcp85
			ll2_baseline_raster= tmin[[1]]
			values(ll2_baseline_raster) = ll2_baseline
			ll2_baseline_raster = disaggregate(ll2_baseline_raster,fact=25)
			ll2_baseline_raster = mask(ll2_baseline_raster,europe)

			# EDEN-2 - LL-2
			ll2_eden = sn_2099rcp85_eden - sn_2099rcp85_ll2
			ll2_eden_raster = tmin[[1]]
			values(ll2_eden_raster) = ll2_eden
			ll2_eden_raster = disaggregate(ll2_eden_raster,fact=25)
			ll2_eden_raster = mask(ll2_eden_raster,europe)

			# making comparison plots

				# color palette

				pos_colors=rev(viridis(100))
				neg_colors=(magma(100))
				sn_negbreaks=seq(-40000,0,length.out=67)
				sn_posbreaks=seq(0,40000,length.out=67)
				sn_bps=c(sn_negbreaks,sn_posbreaks[-1])
				sn_col_list=c(neg_colors[1:66],'#FFFFFF',pos_colors[1:66])


			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_eden_baseline_rcp85.svg")
			plot(eden_baseline_raster,col=sn_col_list,breaks=sn_bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="EDEN-2 - BASELINE; 2099; RCP 8.5")
			plot(eden_baseline_raster,legend.only=TRUE,col=sn_col_list,breaks=sn_bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-40000,-20000,0,20000,40000),labels=c("-40","-20",0,"+20","+40")),legend.args=list(text=expression(paste(Delta,"SP x 1000")),line=1.3))
						plot(shapes,border=c("#3DDC97","#A67DB8","#0D1821","#3F612D"),col="transparent",lwd=4,add=T)
						points(edencoords,pch=19,cex=1.5,col="blue")
			dev.off()


			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_ll2_baseline_rcp85.svg")
			plot(ll2_baseline_raster,col=sn_col_list,breaks=sn_bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="LL-2 - BASELINE; 2099; RCP 8.5")
			plot(ll2_baseline_raster,legend.only=TRUE,col=sn_col_list,breaks=sn_bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-40000,-20000,0,20000,40000),labels=c("-40","-20",0,"+20","+40")),legend.args=list(text=expression(paste(Delta,"SP x 1000")),line=1.3))
									plot(shapes,border=c("#3DDC97","#A67DB8","#0D1821","#3F612D"),col="transparent",lwd=4,add=T)
									points(llcoords,pch=19,cex=1.5,col="blue")
			dev.off()


			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_ll2_eden_rcp85.svg")
			plot(ll2_eden_raster,col=sn_col_list,breaks=sn_bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="EDEN-2 - LL-2; 2099; RCP 8.5")
			plot(ll2_eden_raster,legend.only=TRUE,col=sn_col_list,breaks=sn_bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-40000,-20000,0,20000,40000),labels=c("-40","-20",0,"+20","+40")),legend.args=list(text=expression(paste(Delta,"SP x 1000")),line=1.3))
									plot(shapes,border=c("#3DDC97","#A67DB8","#0D1821","#3F612D"),col="transparent",lwd=4,add=T)
			dev.off()


##########################################


# Just the life strategies
			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_lifehistory_clusters.svg")
			plot(ll2_eden_raster,col="grey",colNA='grey',legend=FALSE,axes=FALSE,box=FALSE)
			plot(shapes,col=c("#3DDC97","#A67DB8","#0D1821","#3F612D"),border=c("#3DDC97","#A67DB8","#0D1821","#3F612D"),add=T)
			dev.off()

##########################################
##########################################
##########################################
##########################################


# Doing the same for DTB
	
	# Load LMMLASSO model weights for seed number
		setwd("~/Clim_GWAS/LMMLASSO/Weights_Corrected")
		dtbw_ridge=npyLoad('Weightsridge_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
		dtbw_ridge=as.matrix(dtbw_ridge,ncol=1)
		dtbbetas=npyLoad('Weights_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
		dtbbetas=as.matrix(dtbbetas,ncol=1)


# get predicted DTB for the baseline, eden-2, and ll-2

			# BASELINE
			
			predictors_2099rcp85=cbind(tmin2099rcp85,tmax2099rcp85)
			baseline_predictors=cbind((predictors_2099rcp85*1),(predictors_2099rcp85*admix[,1]),(predictors_2099rcp85*admix[,2]),(predictors_2099rcp85*admix[,3]),(predictors_2099rcp85*admix[,4]))
			dtb_2099rcp85=(baseline_predictors%*%dtbbetas) + (K_values %*% dtbw_ridge) + dtb_intercept
			dtb_2099rcp85[dtb_2099rcp85<0]=0

			# EDEN-2 

			eden_predictors=cbind((predictors_2099rcp85*1),(predictors_2099rcp85*admixture_eden[,2]),(predictors_2099rcp85*admixture_eden[,3]),(predictors_2099rcp85*admixture_eden[,4]),(predictors_2099rcp85*admixture_eden[,5]))
			dtb_2099rcp85_eden=(eden_predictors%*%dtbbetas) + (k1 %*% dtbw_ridge) + dtb_intercept
			dtb_2099rcp85_eden[dtb_2099rcp85_eden<0]=0

			# LL-2

			ll2_predictors=cbind((predictors_2099rcp85*1),(predictors_2099rcp85*admixture_ll2[,2]),(predictors_2099rcp85*admixture_ll2[,3]),(predictors_2099rcp85*admixture_ll2[,4]),(predictors_2099rcp85*admixture_ll2[,5]))
			dtb_2099rcp85_ll2=(ll2_predictors%*%dtbbetas) + (k2 %*% dtbw_ridge) + dtb_intercept
			dtb_2099rcp85_ll2[dtb_2099rcp85_ll2<0]=0



			# compare differences

			# EDEN-2 VS BASELINE

			eden_baseline = dtb_2099rcp85_eden - dtb_2099rcp85
			eden_baseline_raster = tmin[[1]]
			values(eden_baseline_raster) = eden_baseline
			eden_baseline_raster=disaggregate(eden_baseline_raster,fact=25)
			eden_baseline_raster=mask(eden_baseline_raster,europe)

			# LL-2 VS BASELINE

			ll2_baseline = dtb_2099rcp85_ll2 - dtb_2099rcp85
			ll2_baseline_raster= tmin[[1]]
			values(ll2_baseline_raster) = ll2_baseline
			ll2_baseline_raster = disaggregate(ll2_baseline_raster,fact=25)
			ll2_baseline_raster = mask(ll2_baseline_raster,europe)

			# EDEN-2 - LL-2
			ll2_eden = dtb_2099rcp85_eden - dtb_2099rcp85_ll2
			ll2_eden_raster = tmin[[1]]
			values(ll2_eden_raster) = ll2_eden
			ll2_eden_raster = disaggregate(ll2_eden_raster,fact=25)
			ll2_eden_raster = mask(ll2_eden_raster,europe)

			# making comparison plots

				# color palette

				pos_colors=rev(viridis(100))
				neg_colors=(magma(100))
				dtb_negbreaks=seq(-250,0,length.out=67)
				dtb_posbreaks=seq(0,250,length.out=67)
				dtb_bps=c(dtb_negbreaks,dtb_posbreaks[-1])
				dtb_col_list=c(neg_colors[1:66],'#FFFFFF',pos_colors[1:66])


			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_dtbeden_baseline_rcp85.svg")
			plot(eden_baseline_raster,col=dtb_col_list,breaks=dtb_bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="EDEN-2 - BASELINE; 2099; RCP 8.5")
			plot(eden_baseline_raster,legend.only=TRUE,col=dtb_col_list,breaks=dtb_bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-200,-100,0,100,200),labels=c("-200","-100",0,"+100","+200")),legend.args=list(text=expression(paste(Delta,"DTB")),line=1.3))
						plot(shapes,border=c("#3DDC97","#A67DB8","#0D1821","#3F612D"),col="transparent",lwd=4,add=T)
						points(edencoords,pch=19,cex=1.5,col="blue")
			dev.off()


			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_dtbll2_baseline_rcp85.svg")
			plot(ll2_baseline_raster,col=dtb_col_list,breaks=dtb_bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="LL-2 - BASELINE; 2099; RCP 8.5")
			plot(ll2_baseline_raster,legend.only=TRUE,col=dtb_col_list,breaks=dtb_bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-200,-100,0,100,200),labels=c("-200","-100",0,"+100","+200")),legend.args=list(text=expression(paste(Delta,"DTB")),line=1.3))
									plot(shapes,border=c("#3DDC97","#A67DB8","#0D1821","#3F612D"),col="transparent",lwd=4,add=T)
									points(llcoords,pch=19,cex=1.5,col="blue")
			dev.off()


			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_dtbll2_eden_rcp85.svg")
			plot(ll2_eden_raster,col=dtb_col_list,breaks=dtb_bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="EDEN-2 - LL-2; 2099; RCP 8.5")
			plot(ll2_eden_raster,legend.only=TRUE,col=dtb_col_list,breaks=dtb_bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-200,-100,0,100,200),labels=c("-200","-100",0,"+100","+200")),legend.args=list(text=expression(paste(Delta,"DTB")),line=1.3))
									plot(shapes,border=c("#3DDC97","#A67DB8","#0D1821","#3F612D"),col="transparent",lwd=4,add=T)
			dev.off()
























