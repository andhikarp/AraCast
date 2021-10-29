#Environmental blocking figure (DTB)
library(RcppCNPy)
setwd("~/Clim_GWAS/Clim_GWAS_2")
plantings=c("HalleFall2006","NorwichFall2006","NorwichSpring2007","NorwichSummer2006","NorwichSummer2007","OuluFall2007","ValenciaFall2006")
titles=c("Halle, FAL 2006","Norwich, FAL 2006","Norwich, SPR 2007","Norwich, SUM 2006","Norwich, SUM 2007","Oulu, FAL 2007","Valencia, FAL 2006")


plantings=a$Planting
levels(plantings)[levels(plantings)=='HalleFall2006']<-'#8f0000'
levels(plantings)[levels(plantings)=='NorwichSummer2006']<-'#c54259'
levels(plantings)[levels(plantings)=='NorwichSummer2007']<-'#e881ab'
levels(plantings)[levels(plantings)=='NorwichSpring2007']<-'#ffc2f3'
levels(plantings)[levels(plantings)=='NorwichFall2006']<-'#d890e7'
levels(plantings)[levels(plantings)=='OuluFall2007']<-'#9f66e1'
levels(plantings)[levels(plantings)=='ValenciaFall2006']<-'#3a47de'

plantingmains=c("HF","NF","NSP",expression("NSU"^"06"),expression("NSU"^"07"),"OF","VF")




svg("~/Clim_GWAS/Clim_GWAS_2/Figures/envblock_dtb.svg",width=7,height=7)
par(mfrow=c(3,3),pty='s',oma=c(3,3,0,0),mai=c(0.2,0.2,0.3,0.1)) #oma = outer margin
	i=1
	plantingname=plantings[i]
	base=read.csv(paste0("LOO_",plantingname,"_fitint.csv")) #no gxe
	admix=read.csv(paste0("LOO_",plantingname,"_admix2.csv")) #admixture-based gxe

	setwd("~/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")

	r2_base=npyLoad(paste0("LOOresults_",plantingname,"_fitint.npy"))[4]
	r2_admix=npyLoad(paste0("LOOresults_",plantingname,"_admix2.npy"))[4]

	plot(Y_hat~Y_test,data=base,main=as.character(plantingmains[i]),font.main=2,xlim=c(0,250),xaxt='n',yaxt='n',ylim=c(0,400),pch=19,col='#000000',xlab="",ylab="")
	axis(1, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
	axis(2, at = seq(0, 400, by = 100), las=1,labels=c("0",NA,"200",NA,"400"))
	grid(lty=2,col='grey')
	abline(a=0,b=1)
	points(Y_hat~Y_test,data=admix,pch=17,col='#4974a5')
	legend("topright",	bty="n",cex=1.0,text.col=c('#000000','#4974a5'),legend = c(as.expression(bquote(bold(italic(r)^2 == .(format(r2_base,digits=3,nsmall=3))))),as.expression(bquote(bold(italic(r)^2 == .(format(r2_admix,digits=3,nsmall=3)))))))
	mtext(side=1,expression(bold("Observed Days-to-bolting")),outer=TRUE,line=1.25)
	mtext(side=2,expression(bold("Predicted Days-to-bolting")),outer=TRUE,line=1.25)


for (i in (2:3)){
		setwd("~/Clim_GWAS/Clim_GWAS_2")
		plantingname=plantings[i]
		base=read.csv(paste0("LOO_",plantingname,"_fitint.csv")) #no gxe
		admix=read.csv(paste0("LOO_",plantingname,"_admix2.csv")) #admixture-based gxe
		plot(Y_hat~Y_test,data=base,main=as.character(plantingmains[i]),font.main=2,xlim=c(0,250),xaxt='n',yaxt='n',ylim=c(0,250),pch=19,col='#000000',xlab="",ylab="")
		axis(1, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
		axis(2, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
		grid(lty=2,col='grey')

		setwd("~/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")
		r2_base=npyLoad(paste0("LOOresults_",plantingname,"_fitint.npy"))[4]
		r2_admix=npyLoad(paste0("LOOresults_",plantingname,"_admix2.npy"))[4]

		abline(a=0,b=1)
		points(Y_hat~Y_test,data=admix,pch=17,col='#4974a5')
		legend("topright",	bty="n",cex=1.0,text.col=c('#000000','#4974a5'),legend = c(as.expression(bquote(bold(italic(r)^2 == .(format(round(r2_base,3),digits=3,nsmall=3))))),as.expression(bquote(bold(italic(r)^2 == .(format(round(r2_admix,3),digits=3,nsmall=3)))))))

	}

	i=4
	setwd("~/Clim_GWAS/Clim_GWAS_2")
	plantingname=plantings[i]
	base=read.csv(paste0("LOO_",plantingname,"_fitint.csv")) #no gxe
	admix=read.csv(paste0("LOO_",plantingname,"_admix2.csv")) #admixture-based gxe

	setwd("~/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")

	r2_base=npyLoad(paste0("LOOresults_",plantingname,"_fitint.npy"))[4]
	r2_admix=npyLoad(paste0("LOOresults_",plantingname,"_admix2.npy"))[4]

	plot(Y_hat~Y_test,data=base,main=as.expression(bquote(bold("NSU"^"06"))),xlim=c(0,250),xaxt='n',yaxt='n',ylim=c(0,250),pch=19,col='#000000',xlab="",ylab="")
	axis(1, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
	axis(2, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
	grid(lty=2,col='grey')
	abline(a=0,b=1)
	points(Y_hat~Y_test,data=admix,pch=17,col='#4974a5')
	legend("topright",	bty="n",cex=1.0,text.col=c('#000000','#4974a5'),legend = c(as.expression(bquote(bold(italic(r)^2 == .(format(r2_base,digits=3,nsmall=3))))),as.expression(bquote(bold(italic(r)^2 == .(format(r2_admix,digits=3,nsmall=3)))))))
	mtext(side=1,expression(bold("Observed Days-to-bolting")),outer=TRUE,line=1.25)
	mtext(side=2,expression(bold("Predicted Days-to-bolting")),outer=TRUE,line=1.25)

	i=5
	setwd("~/Clim_GWAS/Clim_GWAS_2")
	plantingname=plantings[i]
	base=read.csv(paste0("LOO_",plantingname,"_fitint.csv")) #no gxe
	admix=read.csv(paste0("LOO_",plantingname,"_admix2.csv")) #admixture-based gxe

	setwd("~/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")

	r2_base=npyLoad(paste0("LOOresults_",plantingname,"_fitint.npy"))[4]
	r2_admix=npyLoad(paste0("LOOresults_",plantingname,"_admix2.npy"))[4]

	plot(Y_hat~Y_test,data=base,main=as.expression(bquote(bold("NSU"^"07"))),xlim=c(0,250),xaxt='n',yaxt='n',ylim=c(0,250),pch=19,col='#000000',xlab="",ylab="")
	axis(1, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
	axis(2, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
	grid(lty=2,col='grey')
	abline(a=0,b=1)
	points(Y_hat~Y_test,data=admix,pch=17,col='#4974a5')
	legend("topright",	bty="n",cex=1.0,text.col=c('#000000','#4974a5'),legend = c(as.expression(bquote(bold(italic(r)^2 == .(format(r2_base,digits=3,nsmall=3))))),as.expression(bquote(bold(italic(r)^2 == .(format(r2_admix,digits=3,nsmall=3)))))))
	mtext(side=1,expression(bold("Observed Days-to-bolting")),outer=TRUE,line=1.25)
	mtext(side=2,expression(bold("Predicted Days-to-bolting")),outer=TRUE,line=1.25)

for (i in (6:7)){
		setwd("~/Clim_GWAS/Clim_GWAS_2")
		plantingname=plantings[i]
		base=read.csv(paste0("LOO_",plantingname,"_fitint.csv")) #no gxe
		admix=read.csv(paste0("LOO_",plantingname,"_admix2.csv")) #admixture-based gxe
		plot(Y_hat~Y_test,data=base,main=as.character(plantingmains[i]),font.main=2,xlim=c(0,250),xaxt='n',yaxt='n',ylim=c(0,250),pch=19,col='#000000',xlab="",ylab="")
		axis(1, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
		axis(2, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
		grid(lty=2,col='grey')

		setwd("~/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")
		r2_base=npyLoad(paste0("LOOresults_",plantingname,"_fitint.npy"))[4]
		r2_admix=npyLoad(paste0("LOOresults_",plantingname,"_admix2.npy"))[4]

		abline(a=0,b=1)
		points(Y_hat~Y_test,data=admix,pch=17,col='#4974a5')
		legend("topright",	bty="n",cex=1.0,text.col=c('#000000','#4974a5'),legend = c(as.expression(bquote(bold(italic(r)^2 == .(format(round(r2_base,3),digits=3,nsmall=3))))),as.expression(bquote(bold(italic(r)^2 == .(format(round(r2_admix,3),digits=3,nsmall=3)))))))

	}



	setwd("~/Clim_GWAS/Clim_GWAS_2")

	korves=read.csv("Korves_pvo2.csv")
	korves_admix=read.csv("Korves_admix_pvo2.csv")

	r2_base=cor(korves$Y_test,korves$Y_hat)**2
	r2_admix=cor(korves_admix$Y_test,korves_admix$Y_hat)**2

	plot(Y_hat~Y_test,data=korves,main="RS",xaxt='n',yaxt='n',xlim=c(0,250),ylim=c(0,250),pch=19,col='#000000',xlab="",ylab="")
	grid(lty=2,col='grey')
	points(Y_hat~Y_test,data=korves_admix,pch=17,col='#4974a5')
	axis(1, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
	axis(2, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
	abline(a=0,b=1)
	legend("topright",	bty="n",cex=1.0,text.col=c('#000000','#4974a5'),legend = c(as.expression(bquote(bold(italic(r)^2 == .(format(round(r2_base,3),digits=3,nsmall=3))))),as.expression(bquote(bold(italic(r)^2 == .(format(round(r2_admix,3),digits=3,nsmall=3)))))))

	plot(NA)
	legend('topleft',legend=c(expression(italic("G + E")),expression(italic("G x E"))),bg="transparent",pch=c(19,17),col=c('#000000','#4974a5'))
dev.off()



