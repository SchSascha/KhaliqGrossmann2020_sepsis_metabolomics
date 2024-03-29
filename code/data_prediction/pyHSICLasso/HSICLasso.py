__author__ = 'makotoy'

import numpy as np
from kernel_tools import *

#We used the Nonnegative LARS solver
#http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/5523/zip/imm5523.zip
def hsiclasso(Xin,Yin,numFeat=10,ykernel='Gauss'):
    
    d,n = Xin.shape

    #Centering matrix
    H = np.eye(n) - 1.0/float(n)*np.ones(n)

    #Normalization
    Xin = Xin - Xin.mean(1)[:,None] #Added by PG
    XX = Xin/(Xin.std(1)[:,None])*np.sqrt(float(n-1)/(n))
    
    if ykernel == 'Gauss':
        YY = Yin/(Yin.std(1)[:,None]+10e-20)*np.sqrt(float(n-1)/(n))
        L = kernel_Gaussian(YY,YY,1.0)
    else:
        L = kernel_Delta_norm(Yin,Yin)

    L = np.dot(H,np.dot(L,H))

    #Preparing design matrix for HSIC Lars
    X = np.zeros((n*n,d))
    Xty = np.zeros((d,1))
    for ii in range(0,d):
        Kx = kernel_Gaussian(XX[ii,None],XX[ii,None],1.0)
        tmp = np.dot(np.dot(H,Kx),H)
        X[:,ii] = tmp.flatten()
        Xty[ii] = (tmp*L).sum()

    #Nonnegative Lars
    #Inactive set
    I = range(0,d)
    A = []
    beta = np.zeros((d,1))
    
    XtXbeta = np.dot(X.transpose(),np.dot(X,beta))
    c = Xty - XtXbeta
    j = c.argmax()
    C = c[j]
    A.append(I[j])
    I.remove(I[j])

    inc_path = True
    k = 0
    if inc_path:
        path = np.zeros((d,4*d))
        lam = np.zeros((1,4*d))
        
        path[:,k] = beta.transpose()
    
    
    while float(sum(c[A]))/float(len(A)) >= 1e-9 and len(A) < numFeat+1:
        tmp = float(sum(c[A]))/float(len(A))
        print tmp
        
        s = np.ones((len(A),1))
        #w = np.linalg.solve(np.dot(X[:,A].transpose(),X[:,A])+1e-10*np.eye(len(A)),s)
        w = np.linalg.solve(np.dot(X[:,A].transpose(),X[:,A]),s)
        XtXw = np.dot(X.transpose(),np.dot(X[:,A],w))
        
        gamma1 = (C - c[I])/(XtXw[A[0]] - XtXw[I])
        gamma2 = -beta[A]/(w)
        gamma3 = np.zeros((1,1))
        gamma3[0] = c[A[0]]/(XtXw[A[0]])
        gamma = np.concatenate((np.concatenate((gamma1,gamma2)),gamma3))
        
        gamma[gamma <= 1e-9] = np.inf
        t = gamma.argmin()
        mu = min(gamma)
        
        beta[A] = beta[A] + mu*w
        #print beta[A]
        #print t, len(gamma1), (len(gamma1) + len(gamma2) + 1)
        
        if t >= len(gamma1) and t < (len(gamma1) + len(gamma2)):
            lassocond = 1
            j = t - len(gamma1)
            I.append(A[j])
            A.remove(A[j])
            
        else:
            lassocond = 0
        
        XtXbeta = np.dot(X.transpose(),np.dot(X,beta))
        c = Xty - XtXbeta
        j = np.argmax(c[I])
        C = max(c[I])
        
        k += 1
        if inc_path:
            path[:,k] = beta.transpose()
            
            if len(C) == 0:
                lam[k] = 0
            else:
                lam[0,k] = C[0]
        #print mu,t,len(I)
        if lassocond is 0:
            A.append(I[j])
            I.remove(I[j])
            
        #print tmp
    
    if inc_path:
        path_final = path[:,0:(k+1)]
        lam_final = lam[0:(k+1)]
    
    return path_final,beta,A,lam_final
