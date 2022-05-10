#Seeing how seed number changes from year to year
#27 February- Spring planting
#25 May     - Summer planting
#3  Oct     - Fall planting
dtb_intercept=-0.578511150335065
sn_intercept=29.54953719927532

library(RcppCNPy)
library(raster)
set.seed(12)
setwd("~/Clim_GWAS/Clim_GWAS_2")
rasterOptions(tmpdir='~/Clim_GWAS/r_tmp')
values=npyLoad('kriged_valuesnoNA.npy')
clusters=stack('lifehistory_clusters.nc')
values=as.matrix(values)
tmin26=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP2.6/tasmin_rcp26_20410101-21001231.nc")
tmax26=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP2.6/tasmax_rcp26_20410101-21001231.nc")
tmin45=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP4.5/tasmin_rcp45_20410101-21001231.nc")
tmax45=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP4.5/tasmax_rcp45_20410101-21001231.nc")
tmin60=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP6.0/tasmin_rcp60_20410101-21001231.nc")
tmax60=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP6.0/tasmax_rcp60_20410101-21001231.nc")
tmin85=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP8.5/tasmin_rcp85_20410101-21001231.nc")
tmax85=stack("~/Clim_GWAS/Clim_GWAS_2/CMIP5/RCP8.5/tasmax_rcp85_20410101-21001231.nc")
admix=read.csv("~/Clim_GWAS/Clim_GWAS_2/krigclus4.csv")
admix=admix[,-1]


setwd("~/Clim_GWAS/LMMLASSO/Weights_Corrected")
w_ridge=npyLoad('Weightsridge_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
w_ridge=as.matrix(w_ridge,ncol=1)
betas=npyLoad('Weights_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
betas=as.matrix(betas,ncol=1)
snw_ridge=npyLoad('SNBASE2Weightsridge_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
snw_ridge=as.matrix(snw_ridge,ncol=1)
snbetas=npyLoad('SNBASE2Weights_k01wg_minmax_4clusters_ADMIXTUREfitint.npy')
snbetas=as.matrix(snbetas,ncol=1)
setwd("~/Clim_GWAS/E-OBSv19.0HOM")
library(raster)
tmin=stack("tmin_reduced.nc")
tmax=stack("tmax_reduced.nc")
library(rgdal)
world=readOGR("~/Clim_GWAS/Clim_GWAS_2/ne_110m_admin_0_countries.shp")
europe=world[world$CONTINENT=='Europe',]


years=seq(2041,2099)


nacount26=matrix(NA,nrow=1,ncol=1)
CEURspring26=matrix(NA,nrow=3015,ncol=1)
dtbCEURspring26=matrix(NA,nrow=3015,ncol=1)
zeron26=matrix(NA,nrow=1,ncol=1)	
# nacount45=matrix(NA,nrow=1,ncol=1)
# CEURspring45=matrix(NA,nrow=3015,ncol=1)
# dtbCEURspring45=matrix(NA,nrow=3015,ncol=1)
# zeron45=matrix(NA,nrow=1,ncol=1)
# nacount60=matrix(NA,nrow=1,ncol=1)
# CEURspring60=matrix(NA,nrow=3015,ncol=1)
# dtbCEURspring60=matrix(NA,nrow=3015,ncol=1)
# zeron60=matrix(NA,nrow=1,ncol=1)
nacount85=matrix(NA,nrow=1,ncol=1)
CEURspring85=matrix(NA,nrow=3015,ncol=1)
dtbCEURspring85=matrix(NA,nrow=3015,ncol=1)
zeron85=matrix(NA,nrow=1,ncol=1)
for (i in c(26,85)){
	startdate= 58 #sowing date
	for (j in (1:length(years))){
	ns=seq(startdate,(startdate+202))
	a=eval(as.name(paste0('tmin',i)))[[ns]]
	a=mask(a,clusters[[2]])
	a=mask(a,europe)
	b=eval(as.name(paste0('tmax',i)))[[ns]]
	b=mask(b,clusters[[2]])
	b=mask(b,europe)
	currentyear=years[j]
	vmin=values(a)-273.15
	vmax=values(b)-273.15
	vvb=cbind(vmin,vmax)
	vvb=cbind((vvb*1),(vvb*admix[,1]),(vvb*admix[,2]),(vvb*admix[,3]),(vvb*admix[,4]))
	dtb_preds=(vvb %*% betas) + (values %*% w_ridge) + dtb_intercept #intercept
	#dtb_preds[dtb_preds<16]=NA 
	#dtb_preds[dtb_preds>246]=NA 
	sn_preds=(vvb%*%snbetas) + (values %*% snw_ridge) + sn_intercept
	sn_preds[sn_preds<0]=0
	sn_preds=matrix(sn_preds,nrow=3015,ncol=1)
	dtb_preds=matrix(dtb_preds,nrow=3015,ncol=1)
	count_na=sum(is.na(sn_preds))
	zerocount=sum(na.omit(sn_preds)==0)
	assign(paste0('CEURspring',i),cbind(eval(as.name(paste0('CEURspring',i))),sn_preds))
	assign(paste0('dtbCEURspring',i),cbind(eval(as.name(paste0('dtbCEURspring',i))),dtb_preds))
	assign(paste0('zeron',i),cbind(eval(as.name(paste0('zeron',i))),zerocount))
	assign(paste0('nacount',i),cbind(eval(as.name(paste0('nacount',i))),count_na))
	remainder=currentyear%%4
	if (remainder==0 && i==45|i == 85){
		startdate=startdate+366
	} else{
	startdate=startdate+365
}}}

#2006 predictions

tmins=tmin[[(58+(365*6)):((58+365*6)+202)]]
tmins=mask(tmins,clusters[[2]])
tmins=mask(tmins,europe)
tmaxs=tmax[[(58+(365*6)):((58+365*6)+202)]]
tmaxs=mask(tmaxs,clusters[[2]])
tmaxs=mask(tmaxs,europe)
vmins=values(tmins)
vmaxs=values(tmaxs)
vvbs=cbind(vmins,vmaxs)
vvbs=cbind((vvbs*1),(vvbs*admix[,1]),(vvbs*admix[,2]),(vvbs*admix[,3]),(vvbs*admix[,4]))
dtb_preds=(vvbs %*% betas) + (values %*% w_ridge) + dtb_intercept
#dtb_preds[dtb_preds<16]=NA 
#dtb_preds[dtb_preds>246]=NA 
CEURspring2006=(vvbs%*%snbetas) + (values %*% snw_ridge) + sn_intercept
CEURspring2006[CEURspring2006<0]=0
CEURspring2006=matrix(CEURspring2006,nrow=3015,ncol=1)
dtbCEURspring2006=matrix(dtb_preds,nrow=3015,ncol=1)
nacount=sum(is.na(sn_preds))
zeron=sum(na.omit(sn_preds)==0)


CEURspring26=CEURspring26[,-1]
# CEURspring45=CEURspring45[,-1]
# CEURspring60=CEURspring60[,-1]
CEURspring85=CEURspring85[,-1]
dtbCEURspring26=dtbCEURspring26[,-1]
# dtbCEURspring45=dtbCEURspring45[,-1]
# dtbCEURspring60=dtbCEURspring60[,-1]
dtbCEURspring85=dtbCEURspring85[,-1]

save(dtbCEURspring2006,dtbCEURspring26,dtbCEURspring85,CEURspring2006,CEURspring26,CEURspring85,file="~/Clim_GWAS/Clim_GWAS_2/CEUR_springleapadmixfitint.RData")




