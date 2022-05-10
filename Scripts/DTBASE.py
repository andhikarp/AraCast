# Base G + E model for testing different K-matrices

import sys
#import limix.vardec.vardec as va 
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
import scipy.stats

os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
plantings=np.load("Plantings_NEW.npy")

import matplotlib
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
K=np.load("k01wg.npy")

os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/New_Microclim")
env=pd.read_csv('Microclimate_minmaxDaily_threshold_0_0.csv',sep=',')
logical=env.columns.str.startswith(('PAR','TT','PTT','daylength') )
env=env.iloc[:,~logical] 
environment=np.array(env)
X=environment
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
y=np.load("DTB_NEW.npy")
# X=np.concatenate((dtb,environment),axis=1)

#Running LMMLASSO
alphas = 2.**(np.linspace(-10,10,10)) #list of alphas to test
n_splits=10
N = X.shape[0]
n_alphas = len(alphas)
MSE_train = np.zeros((n_splits,n_alphas))
MSE_test  = np.zeros((n_splits,n_alphas))
W_nonzero = np.zeros((n_splits,n_alphas))
rsquared  = np.zeros((n_splits,n_alphas))
os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2/Temp_Files")

import lmmlasso
lasso = lmmlasso.LmmLasso(warm_start=True,fit_intercept=True,tol=0.5,random_state=12)
MSE_train,MSE_test,W_nonzero, rsquared = lmmlasso.runCrossValidation(lasso,X,y,alphas,n_splits=10,K=K,verbose=True)
MSE_train_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_train.mean(axis=0))).derivative(n=2)
MSE_test_inter=sp.interpolate.UnivariateSpline(x=alphas, y=(MSE_test.mean(axis=0))).derivative(n=2)
alphas_inter = 2.**(np.linspace(-10,10,100))
idx_train = np.argmin(MSE_train.mean(axis=0))
idx_test = np.argmin(MSE_test.mean(axis=0))

alpha_cv = (float(alphas[idx_train])+float(alphas[idx_test]))/2


os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/Clim_GWAS_2")
kf12 = KFold(n_splits,shuffle=True,random_state=12)
kf12.get_n_splits(X)

#Model fitting step
N = X.shape[0]
kf = KFold(n_splits,shuffle=True,random_state=12)
MSE_train_final = np.zeros((n_splits,))
MSE_test_final  = np.zeros((n_splits,))
W_nonzero_final = np.zeros((n_splits,))
rsquared_final  = np.zeros((n_splits,)) 
kendall_final   = np.zeros((n_splits,))

kf.get_n_splits(X)
lasso.set_params(alpha=alpha_cv)

ifold = 0 
for train_index,test_index in kf.split(X):
    print(('running fold %d'%ifold)) 
    X_train, X_test, y_train, y_test= X[train_index], X[test_index], y[train_index], y[test_index]
    K_train = K[train_index][:,train_index]
    K_test  = K[test_index][:,train_index]
    model_fit=lasso.fit(X_train,y_train,K=K_train) 
    ytrain_star=model_fit.predict(X_train,K_train)
    ytest_star=model_fit.predict(X_test,K_test)
    MSE_train_final[ifold]=mean_squared_error(ytrain_star,y_train)
    MSE_test_final[ifold]=mean_squared_error(ytest_star,y_test)
    W_nonzero_final[ifold]=sp.sum(model_fit.coef_!=0)
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(y_test, ytest_star)
    rsquared_final[ifold]=r_value**2
    kendall_final[ifold]=sp.stats.kendalltau(y_test,ytest_star)[0]
    ifold +=1

for train_index,test_index in kf12.split(X):
	X_train, X_test, y_train, y_test, plantings_train, plantings_test= X[train_index], X[test_index], y[train_index], y[test_index], plantings[train_index], plantings[test_index] #Split into training and testing data sets
	K_train = K[train_index][:,train_index]
	K_test  = K[test_index][:,train_index]

lasso.set_params(alpha=alpha_cv)
lasso = lasso.fit(X_train,y_train,K=K_train)
weights = lasso.coef_
Y_hat = lasso.predict(X_test, K_test)
pd_Yhat=pd.DataFrame(Y_hat)
pd_Yhat.reset_index(drop=True,inplace=True)
pd_ytest=pd.DataFrame(y_test)
pd_plantings=pd.DataFrame(plantings_test)
for_plotting=pd.concat([pd_plantings,pd_ytest,pd_Yhat],axis=1)
for_plotting.columns=("Planting","Y_test","Y_hat")

HalleFall2006=for_plotting.loc[for_plotting['Planting']=="HalleFall2006"]
NorwichSummer2006=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2006"]
NorwichSummer2007=for_plotting.loc[for_plotting['Planting']=="NorwichSummer2007"]
NorwichSpring2007=for_plotting.loc[for_plotting['Planting']=="NorwichSpring2007"]
NorwichFall2006=for_plotting.loc[for_plotting['Planting']=="NorwichFall2006"]
OuluFall2007=for_plotting.loc[for_plotting['Planting']=="OuluFall2007"]
ValenciaFall2006=for_plotting.loc[for_plotting['Planting']=="ValenciaFall2006"]
for_plotting.to_csv('DTBASE_pvo.csv',header=True)


#final fit

lasso.set_params(alpha=alpha_cv)
lasso = lasso.fit(X,y,K=K)
wride = lasso.w_ridge
intercept = lasso.intercept_
weights = lasso.coef_

os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Weights_Corrected")
np.save("Weights_DTBASE.npy",weights)
np.save("Weightsridge_DTBASE.npy",wride)

alpha_cv_np=np.array(alpha_cv)
c=np.array(sum(weights!=0))
W_nonzero=np.array(c)
rmse_test=math.sqrt(np.mean(MSE_test_final))
r2_model=np.mean(rsquared_final)
kendall_model=np.mean(kendall_final)

os.chdir("/home/student.unimelb.edu.au/andhikap/Clim_GWAS/LMMLASSO/Summary_Results_Corrected")

results=np.array((intercept,alpha_cv_np,W_nonzero,rmse_test,r2_model,kendall_model))
resultspd=pd.DataFrame(results)
resultspd=resultspd.T
resultspd.columns=("intercept","alpha","W_nonzero","rmse_test","R2","kendalltau")
resultspd.to_csv("Results_DTBASE.csv")

