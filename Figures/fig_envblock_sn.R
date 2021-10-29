#Environmental blocking figure (SN)

library(RcppCNPy)
setwd("~/Clim_GWAS/Clim_GWAS_2")
plantings=c("HalleFall2006","NorwichFall2006","NorwichSpring2007","NorwichSummer2006","NorwichSummer2007","OuluFall2007","ValenciaFall2006")
titles=c("Halle, FAL 2006","Norwich, FAL 2006","Norwich, SPR 2007","Norwich, SUM 2006","Norwich, SUM 2007","Oulu, FAL 2007","Valencia, FAL 2006")
plantingmains=c("HF","NF","NSP",expression("NSU"^"06"),expression("NSU"^"07"),"OF","VF")




svg("~/Clim_GWAS/Clim_GWAS_2/Figures/envblock_sn.svg",width=7,height=7)
par(mfrow=c(3,3),pty='s',oma=c(3,3,0,0),mai=c(0.2,0.2,0.3,0.1)) #oma = outer margin
setwd("~/Clim_GWAS/Clim_GWAS_2")
	i=1
	plantingname=plantings[i]
	base=read.csv(paste0("SNLOO_",plantingname,".csv"))/1000 #no gxe
	admix=read.csv(paste0("SNLOO_",plantingname,"_admix2.csv"))/1000 #admixture-based gxe

	setwd("~/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")

	r2_base=read.csv(paste0("SNLOOresults_",plantingname,".csv"))[5]
	r2_admix=read.csv(paste0("SNLOOresults_",plantingname,"_admix2.csv"))[5]

	plot(Y_hat~Y_test,data=base,main=as.character(plantingmains[i]),font.main=2,xlim=c(0,250),xaxt='n',yaxt='n',ylim=c(-10,50),pch=19,col='#000000',xlab="",ylab="")
	axis(1, at = seq(0, 250, by = 50), las=1,labels=c("0",NA,"100",NA,"200",NA))
	axis(2, at = seq(-10, 50, by = 10), las=1,labels=c(NA,"0",NA,"20",NA,"40",NA))
	grid(lty=2,col='grey')
	abline(a=0,b=1)
	points(Y_hat~Y_test,data=admix,pch=17,col='#4974a5')
	legend("topright",	bty="n",cex=1.0,text.col=c('#000000','#4974a5'),legend = c(as.expression(bquote(bold(italic(r)^2 == .(format(as.numeric(round(r2_base,3)),nsmall=3))))),as.expression(bquote(bold(italic(r)^2 == .(format(as.numeric(round(r2_admix,3)),digits=3,nsmall=3)))))))
	mtext(side=1,expression(bold("Observed Seed Proxy / 1000")),outer=TRUE,line=1.25)
	mtext(side=2,expression(bold("Predicted Seed Proxy / 1000")),outer=TRUE,line=1.25)


for (i in (2:3)){
		setwd("~/Clim_GWAS/Clim_GWAS_2")
		plantingname=plantings[i]
		base=read.csv(paste0("SNLOO_",plantingname,".csv"))/1000
		admix=read.csv(paste0("SNLOO_",plantingname,"_admix2.csv"))/1000
		plot(Y_hat~Y_test,data=base,main=as.character(plantingmains[i]),font.main=2,xlim=c(0,50),xaxt='n',yaxt='n',ylim=c(-10,50),pch=19,col='#000000',xlab="",ylab="")
		grid(lty=2,col='grey')
		axis(1, at = seq(0, 50, by = 10), las=1,labels=c("0",NA,"20",NA,"40",NA))
		axis(2, at = seq(-10, 50, by = 10), las=1,labels=c(NA,"0",NA,"20",NA,"40",NA))


		setwd("~/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")
		r2_base=read.csv(paste0("SNLOOresults_",plantingname,".csv"))[5]
		r2_admix=read.csv(paste0("SNLOOresults_",plantingname,"_admix2.csv"))[5]

		abline(a=0,b=1)
		points(Y_hat~Y_test,data=admix,pch=17,col='#4974a5')
		legend("topright",	bty="n",cex=1.0,text.col=c('#000000','#4974a5'),legend = c(as.expression(bquote(bold(italic(r)^2 == .(format(as.numeric(round(r2_base,3)),nsmall=3))))),as.expression(bquote(bold(italic(r)^2 == .(format(as.numeric(round(r2_admix,3)),digits=3,nsmall=3)))))))

	}

i=4
	setwd("~/Clim_GWAS/Clim_GWAS_2")

	plantingname=plantings[i]
	base=read.csv(paste0("SNLOO_",plantingname,".csv"))/1000 #no gxe
	admix=read.csv(paste0("SNLOO_",plantingname,"_admix2.csv"))/1000 #admixture-based gxe

	setwd("~/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")

	r2_base=read.csv(paste0("SNLOOresults_",plantingname,".csv"))[5]
	r2_admix=read.csv(paste0("SNLOOresults_",plantingname,"_admix2.csv"))[5]

	plot(Y_hat~Y_test,data=base,main=as.expression(bquote(bold("NSU"^"06"))),font.main=2,xlim=c(0,50),xaxt='n',yaxt='n',ylim=c(-10,50),pch=19,col='#000000',xlab="",ylab="")
	axis(1, at = seq(0, 50, by = 10), las=1,labels=c("0",NA,"20",NA,"40",NA))
	axis(2, at = seq(-10, 50, by = 10), las=1,labels=c(NA,"0",NA,"20",NA,"40",NA))
	grid(lty=2,col='grey')
	abline(a=0,b=1)
	points(Y_hat~Y_test,data=admix,pch=17,col='#4974a5')
	legend("topright",	bty="n",cex=1.0,text.col=c('#000000','#4974a5'),legend = c(as.expression(bquote(bold(italic(r)^2 == .(format(as.numeric(round(r2_base,3)),nsmall=3))))),as.expression(bquote(bold(italic(r)^2 == .(format(as.numeric(round(r2_admix,3)),digits=3,nsmall=3)))))))

i=5
	setwd("~/Clim_GWAS/Clim_GWAS_2")
	plantingname=plantings[i]
	base=read.csv(paste0("SNLOO_",plantingname,".csv"))/1000 #no gxe
	admix=read.csv(paste0("SNLOO_",plantingname,"_admix2.csv"))/1000 #admixture-based gxe

	setwd("~/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")

	r2_base=read.csv(paste0("SNLOOresults_",plantingname,".csv"))[5]
	r2_admix=read.csv(paste0("SNLOOresults_",plantingname,"_admix2.csv"))[5]

	plot(Y_hat~Y_test,data=base,main=as.expression(bquote(bold("NSU"^"07"))),font.main=2,xlim=c(0,50),xaxt='n',yaxt='n',ylim=c(-10,50),pch=19,col='#000000',xlab="",ylab="")
	axis(1, at = seq(0, 50, by = 10), las=1,labels=c("0",NA,"20",NA,"40",NA))
	axis(2, at = seq(-10, 50, by = 10), las=1,labels=c(NA,"0",NA,"20",NA,"40",NA))
	grid(lty=2,col='grey')
	abline(a=0,b=1)
	points(Y_hat~Y_test,data=admix,pch=17,col='#4974a5')
	legend("topright",	bty="n",cex=1.0,text.col=c('#000000','#4974a5'),legend = c(as.expression(bquote(bold(italic(r)^2 == .(format(as.numeric(round(r2_base,3)),nsmall=3))))),as.expression(bquote(bold(italic(r)^2 == .(format(as.numeric(round(r2_admix,3)),digits=3,nsmall=3)))))))


for (i in (6:7)){
		setwd("~/Clim_GWAS/Clim_GWAS_2")
		plantingname=plantings[i]
		base=read.csv(paste0("SNLOO_",plantingname,".csv"))/1000
		admix=read.csv(paste0("SNLOO_",plantingname,"_admix2.csv"))/1000
		plot(Y_hat~Y_test,data=base,main=as.character(plantingmains[i]),font.main=2,xlim=c(0,50),xaxt='n',yaxt='n',ylim=c(-10,50),pch=19,col='#000000',xlab="",ylab="")
		grid(lty=2,col='grey')
		axis(1, at = seq(0, 50, by = 10), las=1,labels=c("0",NA,"20",NA,"40",NA))
		axis(2, at = seq(-10, 50, by = 10), las=1,labels=c(NA,"0",NA,"20",NA,"40",NA))


		setwd("~/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")
		r2_base=read.csv(paste0("SNLOOresults_",plantingname,".csv"))[5]
		r2_admix=read.csv(paste0("SNLOOresults_",plantingname,"_admix2.csv"))[5]

		abline(a=0,b=1)
		points(Y_hat~Y_test,data=admix,pch=17,col='#4974a5')
		legend("topright",	bty="n",cex=1.0,text.col=c('#000000','#4974a5'),legend = c(as.expression(bquote(bold(italic(r)^2 == .(format(as.numeric(round(r2_base,3)),nsmall=3))))),as.expression(bquote(bold(italic(r)^2 == .(format(as.numeric(round(r2_admix,3)),digits=3,nsmall=3)))))))

	}

	plot(NA)
	legend('topleft',legend=c(expression(italic("G + E")),expression(italic("G x E"))),bg="transparent",pch=c(19,17),col=c('#000000','#4974a5'))
dev.off()



