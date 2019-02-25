'''
I want to convolve the Pl from TurboGL with the 
Pc from probablitymagnification
'''
import probabilityMagnification as pm
import numpy as np
from matplotlib import pyplot as plt
import deltaTable as dt
import os as os
import pickle as pkl
import lensingProbabilityDistribution as lpd

def main(z=1.0, alpha=0.83,nMu=1000):

    lensingPDF = \
      lpd.lensingProbabilityDistribution( redshift=z, alpha=alpha, nMagnitudeBins=nMu)
    

    lensingPDF.plotTotalProbabilityDistribution()
    print lensingPDF.convolvedPbhPdfWithLssPdf['y']
def getTurboGL():


    TurboGLfile = 'TurboGL.dat'

    return  np.loadtxt(TurboGLfile, unpack=True)
    


def getPBHpdf( z, dMag ):
    mag = np.arange(0.,1.,dMag)[1:]

    
    delta = dt.getDelta(z)
    
    mag, pdf =  pm.magPDF( delta, mag=mag )
    mag = np.append(0.,mag)
    pdf = np.append(0.,pdf)

    mag -=  pm.getMeanMag( z )


    return mag, pdf



def convolutionPl( mu, alpha, z, pdfMags, MeanMagnitudes, PDFarray ):
    
   
    endInt = mu / (1.-alpha)
    endInt = np.min([endInt,10.])

    #dMuPrime = endInt/nInt
    #dMuPrime should be the P_LSS one as this defines min
    #Does this need to be changed?
    dMuPrime = get_dMuPrime()
    meanMag = pm.getMeanMag(z)

    MuPrimeList = np.arange(0., endInt, dMuPrime) 


    inputPdfPBH = \
      {    'pdfMagnitudes':pdfMags, \
        'meanMagnitudes':MeanMagnitudes, \
      'pdfDistributionArray':PDFarray}
      
    P_C =  dt.getPdfPBH(  mu - MuPrimeList*(1.-alpha), \
                                  alpha*MuPrimeList,\
                                  inputPdfPBH)



    P_LSS = getPLSS_New(MuPrimeList-meanMag, z)
        
    dP_L = dMuPrime*P_LSS*P_C

    
    P_L = np.sum(dP_L)

    return P_L

def getPLSS_New( mu, z):
    print mu
    mag, pdf = getTurboGL()


    pdfMagsMatrix = np.matrix(mag).T*np.matrix(np.ones(len(mu)))

    pdfMagsIndex =\
      np.argmin(np.abs((pdfMagsMatrix - np.matrix(mu))),axis=0)
    


    
    correspondingMag = mag[pdfMagsIndex]

    returnPDF = pdf[pdfMagsIndex]
    returnPDF[ correspondingMag-mu > 2.*(mag[1]-mag[0])] = 0.

    return returnPDF
    
  
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



def totalPl(z=1.0, nMu=100, alpha=0.83):

    #pklFile = 'pickles/PL_%0.1f_%0.2f_%i.pkl' % (z,alpha,nMu)
    #if os.path.isfile( pklFile ):
    #    return pkl.load( open( pklFile, 'rb'))

    meanMu = pm.getMeanMag( z )
    MuList = np.linspace(0., 1.0, nMu)[1:]
    P_l = np.zeros(nMu-1)

    #Get the array of pdfArrays for all possible meanMags (speed up)
    #PBH PDF array showing for P_C(mu|MeanMu)
    pdfMags, MeanMagnitudes, PDFarray = \
      dt.generatePDFs()
    
    #The 
    for i, iMu in enumerate(MuList):
        print("%i/%i" %(i,nMu-1))
        P_l[i] = convolutionPl(  iMu, alpha, z, \
                    pdfMags, MeanMagnitudes, PDFarray )

    #pkl.dump( [ MuList, P_l], open( pklFile, 'wb'))
    return MuList, P_l


def get_dMuPrime():
    mag, pdf = getTurboGL()

    return (mag[1]-mag[0])



if __name__ == '__main__':
    main()
