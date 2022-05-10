# Predicting DTB across the Landscape in 2099

#intercepts were obtained from the fitted G x E model

dtb_intercept=-0.578511150335065
sn_intercept=29.54953719927532

library(RcppCNPy)
library(raster)
library(rgdal)
library(viridis)
set.seed(12)


setwd("~/Clim_GWAS/Clim_GWAS_2")


tmin26=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP2.6/tasmin_rcp26_20410101-21001231.nc")
tmax26=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP2.6/tasmax_rcp26_20410101-21001231.nc")

tmin85=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP8.5/tasmin_rcp85_20410101-21001231.nc")
tmax85=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP8.5/tasmax_rcp85_20410101-21001231.nc")

K_values=npyLoad('kriged_valuesnoNA.npy') # Kinship matrix applied to the landscape
germinationdates=raster('pred_germdate.nc')

# Load LMMLASSO model weights for both DTB and seed number
	setwd("~/Clim_GWAS/LMMLASSO/Weights_Corrected")
	w_ridge=npyLoad('Weightsridge_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
	w_ridge=as.matrix(w_ridge,ncol=1)
	betas=npyLoad('Weights_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
	betas=as.matrix(betas,ncol=1)
	snw_ridge=npyLoad('SNBASE2Weightsridge_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
	snw_ridge=as.matrix(snw_ridge,ncol=1)
	snbetas=npyLoad('SNBASE2Weights_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
	snbetas=as.matrix(snbetas,ncol=1)

# Shapes and temperature data
	setwd("~/Clim_GWAS/E-OBSv19.0HOM")
	tmin=stack("tmin_reduced.nc") #2000-2010 data maximum temperature
	tmax=stack("tmax_reduced.nc") #2000-2010 data minimum temperature
	world=readOGR("~/Clim_GWAS/Clim_GWAS_2/ne_110m_admin_0_countries.shp")
	europe=world[world$CONTINENT=='Europe',]

admix=read.csv("~/Clim_GWAS/Clim_GWAS_2/krigclus4.csv") #kriged admixture proportions
admix=admix[,-1]

germdates=as.data.frame(round(values(germinationdates)))
coords=as.data.frame(coordinates(germinationdates))

############# 2006 predictions ############################################

	# Get the environmental predictors

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

	predictors_2006=cbind(tmin2006,tmax2006)
	predictors_2006=cbind((predictors_2006*1),(predictors_2006*admix[,1]),(predictors_2006*admix[,2]),(predictors_2006*admix[,3]),(predictors_2006*admix[,4]))

	dtb_2006=(predictors_2006 %*% betas) + (K_values %*% w_ridge) + dtb_intercept
	dtb_2006[dtb_2006<16]=NA 
	dtb_2006[dtb_2006>246]=NA 
	sn_2006=(predictors_2006%*%snbetas) + (K_values %*% snw_ridge) + sn_intercept
	sn_2006[sn_2006<0]=0



################################################ RCP 2.6 #################################
##########################################################################################
##########################################################################################


		############# 2099 predictions RCP2.6 ############################################

			# Get the environmental predictors

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

			predictors_2099rcp26=cbind(tmin2099rcp26,tmax2099rcp26)
			predictors_2099rcp26=cbind((predictors_2099rcp26*1),(predictors_2099rcp26*admix[,1]),(predictors_2099rcp26*admix[,2]),(predictors_2099rcp26*admix[,3]),(predictors_2099rcp26*admix[,4]))

			dtb_2099rcp26=(predictors_2099rcp26 %*% betas) + (K_values %*% w_ridge) + dtb_intercept
			dtb_2099rcp26[dtb_2099rcp26<16]=NA 
			dtb_2099rcp26[dtb_2099rcp26>246]=NA 
			sn_2099rcp26=(predictors_2099rcp26%*%snbetas) + (K_values %*% snw_ridge) + sn_intercept
			sn_2099rcp26[sn_2099rcp26<0]=0


		############# Difference between 2006 and 2099 RCP2.6

		dtb_diff= dtb_2099rcp26 - dtb_2006
		dtb_diff_raster = tmin[[1]]
		values(dtb_diff_raster) = dtb_diff

		sn_diff = sn_2099rcp26 - sn_2006
		sn_diff_raster = tmin[[1]]
		values(sn_diff_raster) = sn_diff

		# Making the figure

			pos_colors=rev(viridis(100))
			neg_colors=(magma(50))
			negbreaks=seq(-200,0,length.out=50)
			posbreaks=seq(0,200,length.out=50)
			bps=c(negbreaks,posbreaks[-1])
			col_list=c(neg_colors[1:49],'#FFFFFF',pos_colors[1:49])

			sn_negbreaks=seq(-20000,0,length.out=33)
			sn_posbreaks=seq(0,40000,length.out=67)
			sn_bps=c(sn_negbreaks,sn_posbreaks[-1])
			sn_col_list=c(neg_colors[1:32],'#FFFFFF',pos_colors[1:66])


			# resample the raster to make smoother boundaries
			dtb_diff_raster=disaggregate(dtb_diff_raster,fact=25)
			dtb_diff_raster=mask(dtb_diff_raster,europe)

			sn_diff_raster=disaggregate(sn_diff_raster,fact=25)
			sn_diff_raster=mask(sn_diff_raster,europe)

		
			#plot
			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_DTB2099vs2006_rcp26.svg")
			plot(dtb_diff_raster,col=col_list,breaks=bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="2099; RCP 2.6")
			plot(dtb_diff_raster,legend.only=TRUE,col=col_list,breaks=bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-200,-100,0,100,200),labels=c(-200,-100,0,"+100","+200")),legend.args=list(text=expression(paste(Delta,"DTB")),line=1.3))
			dev.off()

			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_SN2099vs2006_rcp26.svg")
			plot(sn_diff_raster,col=sn_col_list,breaks=sn_bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="2099; RCP 2.6")
			plot(sn_diff_raster,legend.only=TRUE,col=sn_col_list,breaks=sn_bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-20000,0,20000,40000),labels=c("-20",0,"+20","+40")),legend.args=list(text=expression(paste(Delta,"SN x 1000")),line=1.3))
			dev.off()



################################################ RCP 8.5 #################################
##########################################################################################
##########################################################################################

		############# 2099 predictions RCP8.5 ############################################

			# Get the environmental predictors

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

			predictors_2099rcp85=cbind(tmin2099rcp85,tmax2099rcp85)
			predictors_2099rcp85=cbind((predictors_2099rcp85*1),(predictors_2099rcp85*admix[,1]),(predictors_2099rcp85*admix[,2]),(predictors_2099rcp85*admix[,3]),(predictors_2099rcp85*admix[,4]))

			dtb_2099rcp85=(predictors_2099rcp85 %*% betas) + (K_values %*% w_ridge) + dtb_intercept
			dtb_2099rcp85[dtb_2099rcp85<16]=NA 
			dtb_2099rcp85[dtb_2099rcp85>246]=NA 
			sn_2099rcp85=(predictors_2099rcp85%*%snbetas) + (K_values %*% snw_ridge) + sn_intercept
			sn_2099rcp85[sn_2099rcp85<0]=0


		############# Difference between 2006 and 2099 RCP8.5

		dtb_diff= dtb_2099rcp85 - dtb_2006
		dtb_diff_raster = tmin[[1]]
		values(dtb_diff_raster) = dtb_diff

		sn_diff = sn_2099rcp85 - sn_2006
		sn_diff_raster = tmin[[1]]
		values(sn_diff_raster) = sn_diff

		# Making the figure

			pos_colors=rev(viridis(100))
			neg_colors=(magma(100))
			negbreaks=seq(-200,0,length.out=50)
			posbreaks=seq(0,200,length.out=50)
			bps=c(negbreaks,posbreaks[-1])
			col_list=c(neg_colors[1:49],'#FFFFFF',pos_colors[1:49])

			sn_negbreaks=seq(-40000,0,length.out=67)
			sn_posbreaks=seq(0,40000,length.out=67)
			sn_bps=c(sn_negbreaks,sn_posbreaks[-1])
			sn_col_list=c(neg_colors[1:66],'#FFFFFF',pos_colors[1:66])


			# resample the raster to make smoother boundaries
			dtb_diff_raster=disaggregate(dtb_diff_raster,fact=25)
			dtb_diff_raster=mask(dtb_diff_raster,europe)

			sn_diff_raster=disaggregate(sn_diff_raster,fact=25)
			sn_diff_raster=mask(sn_diff_raster,europe)

			# life history outlines
				shapes=readOGR("ashape3.shp") #life history borders
				
			#plot
			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_DTB2099vs2006_rcp85.svg")
			plot(dtb_diff_raster,col=col_list,breaks=bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="2099; RCP 8.5")
			plot(dtb_diff_raster,legend.only=TRUE,col=col_list,breaks=bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-200,-100,0,100,200),labels=c(-200,-100,0,"+100","+200")),legend.args=list(text=expression(paste(Delta,"DTB")),line=1.3))
			dev.off()

			svg("~/Clim_GWAS/Clim_GWAS_2/Figures/MER_SP2099vs2006_rcp85.svg")
			plot(sn_diff_raster,col=sn_col_list,breaks=sn_bps,colNA='grey',legend=FALSE,axes=FALSE,box=FALSE,main="2099; RCP 8.5")
			plot(sn_diff_raster,legend.only=TRUE,col=sn_col_list,breaks=sn_bps,legend.width=1,legend.shrink=0.6,cex=1.5,axis.arg=list(at=c(-40000,-20000,0,20000,40000),labels=c("-40","-20",0,"+20","+40")),legend.args=list(text=expression(paste(Delta,"SP x 1000")),line=1.3))
			dev.off()