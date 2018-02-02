'''
Created on 2018/01/15

@author: Peter Grossmann
'''
import numpy as np
from HSICLasso import *
#from kernel_Gaussian import *
from pylab import *
import scipy.io as spio
import csv
from pandas import read_csv
import pandas as pd

#Read  file
data = read_csv('../../../data/measurements/Summary human sample data.csv', sep = '\t', header = 0, decimal = ',')
data = data[data['Day'] == 0]
data = data[data['CAP / FP'] != '-']
Yin  = np.asfarray(data['Survival'] == 'S') + 1
Yin.shape = [1,Yin.shape[0]]
Xin  = data.loc[:, 'C0':]
Xin  = Xin.dropna(axis = 1)
Xin  = Xin[Xin.columns.values[np.where(Xin.std(0) > 1e-15)]]
Xfields = Xin.columns.values
Xin  = Xin.T.values

#Generate label data
path, beta, A, lam = hsiclasso(Xin, Yin, numFeat=20, ykernel='Delta')

#Print selected column names
ind = np.asarray(A) + 5
dc = Xfields[ind]
print dc

nonzero_ind = beta.nonzero()[0]
t = path.sum(0)
figure()
hold(True)
for ind in range(0,len(A)):
    plot(t,path[A[ind],:])

show()


