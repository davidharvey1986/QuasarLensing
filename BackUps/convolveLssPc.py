'''
I want to convolve the Pl from TurboGL with the 
Pc from probablitymagnification
'''
import probabilityMagnification as pm
import numpy as np
from matplotlib import pyplot as plt
import ipdb as pdb
import deltaTable as dt
import os as os
import pickle as pkl
def main(z=1.0, alpha=0.83,nMu=1000):

    pklFile = 'pickles/PL_%0.1f_%0.2f_%i.pkl' % (z,alpha,nMu)
    if os.path.isfile( pklFile ):
        magConvolve, pdfConvolve = \
          pkl.load( open( pklFile, 'rb'))
        print len(magConvolve)
    else:
        magConvolve, pdfConvolve = totalPl(z=z, nMu=nMu, alpha=alpha)
        pkl.dump( [ magConvolve, pdfConvolve], open( pklFile, 'wb'))


    meanMag = pm.getMeanMag( z)
    dMag = magConvolve[1] - magConvolve[0]
    #pdfConvolve /= np.sum(pdfConvolve)*dMag

    print np.sum(dMag*pdfConvolve)
    print np.sum(dMag*magConvolve*pdfConvolve), pm.getMeanMag( z )

    plt.plot(np.append(0.,magConvolve), np.append(1e-3,pdfConvolve))
    
    plt.yscale('log')
    plt.ylim(0.05,40)
   

    
def matchPBHpdf(z):

    magLSS, pdfLSS = getTurboGL()
    dMag = magLSS[1] - magLSS[0]
    magPBH, pdfPBH = getPBHpdf(z, dMag=dMag)

    return magPBH, pdfPBH, \
      magLSS, pdfLSS

def getTurboGL():


    TurboGLfile = '/Users/DavidHarvey/Documents/Work/QuasarLensing/turboGL3.0/results/hpdf_3.dat'

    return  np.loadtxt(TurboGLfile, unpack=True)
    


def getPBHpdf( z, dMag ):
    mag = np.arange(0.,1.,dMag)[1:]

    
    delta = dt.getDelta(z)
    
    mag, pdf =  pm.magPDF( delta, mag=mag )
    mag = np.append(0.,mag)
    pdf = np.append(0.,pdf)

    mag -=  pm.getMeanMag( z )
    return mag, pdf



def convolutionPl( mu, alpha, z ):
    
    P_L = 1e-5
    endInt = mu / (1.-alpha)
    endInt = np.min([endInt,100.])

    #dMuPrime = endInt/nInt
    #dMuPrime should be the P_LSS one as this defines min
    dMuPrime = get_dMuPrime()
    meanMag = pm.getMeanMag(z)

    for iMu, MuPrime in enumerate(np.arange(0., endInt, dMuPrime)):
        print("PLSS")
        P_LSS = getPLSS( MuPrime-meanMag, z)
        print("PC")
        P_C = getPC( mu - MuPrime*(1.-alpha), \
                         alpha*MuPrime)
        
        dP_L = dMuPrime*P_LSS*P_C
        if (dP_L == 0) & (P_L > 1e-5):
            return P_L
    
        P_L += dP_L
        
    return P_L


def getPLSS( mu, z):

    mag, pdf = getTurboGL()

    
    correspondingMag = mag[ np.argmin( np.abs(mu - mag)) ]
   
   
    #if the required mu is outside the pdf of TurboGL return 0
    if np.abs(correspondingMag - mu ) > 2.*(mag[1]-mag[0]):
        return 0.
    else:
        return pdf[ np.argmin( np.abs(mu - mag)) ]

def getPC( mu, muBar):

    delta = dt.getDelta( meanMag=muBar )

    mag, pdf =  pm.magPDF( delta )
    dMag = mag[1] - mag[0]
    #print muBar, np.sum(pdf*mag*dMag), pdf[  np.argmin( np.abs(mu - mag)) ]
    
    return pdf[  np.argmin( np.abs(mu - mag)) ] 



def totalPl(z=1.0, nMu=100, alpha=0.):
    meanMu = pm.getMeanMag( z )
    MuList = np.linspace(0., 1.0, nMu)[1:]
    P_l = np.zeros(nMu-1)
   
    for i, iMu in enumerate(MuList):
        print("%i/%i" %(i,nMu-1))
        P_l[i] = convolutionPl(  iMu, alpha, z )

    return MuList, P_l


def get_dMuPrime():
    mag, pdf = getTurboGL()

    return mag[1]-mag[0]
