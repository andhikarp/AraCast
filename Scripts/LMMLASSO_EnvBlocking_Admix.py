#Leave One Planting Out Cross-Validation (environmental blocking) with ADMIXTURE proportions 
import sys
import limix.vardec.vardec as va 
import scipy.linalg as LA
from sklearn.linear_model import Lasso
import numpy as np 
from sklearn.metrics import mean_squared_error
import math
from sklearn.model_selection import KFold 
import gc
import scipy as sp
import os
import csv
import time
import pandas as pd
import matplotlib
ncluster=4
exclplanting=["ValenciaFall2006","NorwichSummer2006","NorwichSummer2007","NorwichSpring2007","NorwichFall2006","OuluFall2007","HalleFall2006"]

os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
plantings=np.load("Plantings_NEW.npy")
ids=np.load('Ecotype_ids_NEW_ORDERED.npy')
ids_2029=pd.read_csv("k2029_accessioninfo.csv")
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/ADMIXTURE/admixture_linux-1.3.0") 
admixture=pd.read_csv("k2029_trimmed."+str(ncluster)+".Q",header=None,sep=' ')
design_matrix=np.empty((0,ncluster))
for i in range(ids.shape[0]):
	id_no=ids[i]
	idx_2029=np.where(ids_2029['ecotype_id']==id_no)
	d=admixture.loc[idx_2029[0],:]
	dd=np.reshape(d,(-1,ncluster))
	design_matrix=np.append(design_matrix,dd,axis=0)


ones=np.ones(5623)
ones=np.reshape(ones,(-1,1))
design_matrix=np.append(ones,design_matrix,axis=1)
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")
K=np.load('k01wg.npy') #Load the preferred K-matrix
#K=np.load("K_MATRIX_PLINKIBS.npy")
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")
env=pd.read_csv('Microclimate_minmaxDaily_threshold_0_0.csv',sep=',')
logical=env.columns.str.startswith(('PAR','TT','PTT','daylength') )
env=env.iloc[:,~logical] #get columns that don't start with PAR
environment=np.array(env)
X=np.kron(design_matrix,environment) #not enough memory?
X=np.empty((5623,0))
for i in range(design_matrix.shape[1]):
	dd=design_matrix[:,i]
	dd=np.reshape(dd,(5623,-1))
	preds=environment*dd
	X=np.append(X,preds,axis=1)

os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
y=np.load("DTB_NEW.npy")


#Leaving one planting out
for j in exclplanting:
	exc=plantings!=j
	K_idx=K[exc,:]
	K_idx=K_idx[:,exc]
	X_idx=X[exc,:]
	y_idx=y[exc]
	K_in=K[~exc,:]
	K_in=K_in[:,exc]
	X_in=X[~exc,:]
	y_in=y[~exc]
	#Running LMMLASSO
	alphas = 2.**(sp.linspace(-10,10,10)) #list of alphas to test
	n_splits=10
	N = X_idx.shape[0]
	kf = KFold(n_splits,shuffle=True,random_state=12)
	n_alphas = len(alphas)
	MSE_train = np.zeros((n_splits,n_alphas))
	MSE_test  = np.zeros((n_splits,n_alphas))
	W_nonzero = np.zeros((n_splits,n_alphas))
	rsquared  = np.zeros((n_splits,n_alphas))
	os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
	import lmmlasso
	lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=True,tol=0.5) #note the tolerance value
	MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,X_idx,y_idx,alphas,n_splits=10,K=K_idx,verbose=True)
	MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)#Interpolating the values for alphas within -2 to 10...
	MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)
	alphas_inter = 2.**(sp.linspace(-10,10,100))
	idx_train = sp.argmin(MSE_train.mean(axis=0))
	idx_test = sp.argmin(MSE_test.mean(axis=0))
	alpha_cv = (float(alphas[idx_train])+float(alphas[idx_test]))/2
	#Model fitting
	os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
	lasso.set_params(alpha=alpha_cv)
	model_fit=lasso.fit(X_idx,y_idx,K=K_idx) #Fit the model to the training set
	Y_hat = lasso.predict(X_in, K_in)
	pd_Yhat=pd.DataFrame(Y_hat)
	pd_Yhat.reset_index(drop=True,inplace=True)
	pd_ytest=pd.DataFrame(y_in)
	pd_plantings=pd.DataFrame(plantings[~exc])
	for_plotting=pd.concat([pd_plantings,pd_ytest,pd_Yhat],axis=1)
	for_plotting.columns=("Planting","Y_test","Y_hat")
	for_plotting.to_csv("LOO_"+j+"_admix2.csv",header=True,sep=',')
	weights = lasso.coef_
	alpha_cv_np=np.array(alpha_cv)
	c=np.array(sum(weights!=0))
	W_nonzero=np.array(c)
	rmse=np.sqrt(mean_squared_error(y_in,Y_hat))
	r2_model=sp.stats.pearsonr(Y_hat,y_in)[0]**2
	kendall_model=sp.stats.kendalltau(y_in,Y_hat)[0]
	os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")
	results=np.array((alpha_cv_np,W_nonzero,rmse,r2_model,kendall_model))
	#np.save("LOOresults_"+j+"_admix2.npy",results)
	weights = lasso.coef_
	wridge= lasso.w_ridge
	os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Weights_Corrected")
	np.save("Weights_LOO"+str(j)+"_admix2.npy",weights)
	np.save("Weightsridge_LOO"+str(j)+"_admix2.npy",wridge)



