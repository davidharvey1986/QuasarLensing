'''
When creating the PDFs for the PBH it takes a while to get the right deltas

so i need a table outlining expected mean mag for a given delta.

'''

import probabilityMagnification as pm
import pickle as pkl
import numpy as np
import sys
import ipdb as pdb
from predictDeltaAndNormalisation import *

def deltaTable():


    deltaList = 10**np.linspace(-2.,2.,1000)


    expecMag = []
    normalistation = []
    
    for i, delta in enumerate(deltaList):
        sys.stdout.write("PROGRESS: %i/%i\r" \
                             %(i+1, len(deltaList)))
        sys.stdout.flush()
        iPDF = pbhPDF( delta )
    
        expecMag.append(iPDF.mean)
        normalistation.append(iPDF.norm)
        
        
    pkl.dump([deltaList, np.array(expecMag), np.array(normalistation)],\
                 open('deltaTable.pkl','wb'))


def getDelta( z=1.0, meanMag=None):
    
    '''
    Derive the delta for the PBH 
    '''

    if meanMag is None:
        meanMag = pm.getMeanMag( z )

    deltaList, magnitudes, norm = pkl.load(open('deltaTable.pkl','rb'))

    delta = deltaList[ np.argmin( np.abs( magnitudes - meanMag))]

    return delta
    
    

def generatePDFs():

    '''
    Generate PBH PDFs for each given delta (i.e. mean Mag)

    '''

    
    deltaList, MeanMagnitudes, norm = \
      pkl.load(open('deltaTable.pkl','rb'))

    #mag, pdf =  pm.magPDF( deltaList[0] )
    #Genereate an array that is magnitude across and meanMag
    #PDFarray = np.zeros((len(mag),len(MeanMagnitudes)), float)


   # for i, iDelta in enumerate(deltaList):
   #     print("%i/%i" % (i,len(deltaList)))
   #     mag, pdf =  pm.magPDF( iDelta )

#        PDFarray[:,i] = pdf
#    print np.mean(PDFarray[:,0])
    mag, PDFarray =  pm.magPDFarray( deltaList )
    
    return mag, MeanMagnitudes, PDFarray

        

        
def getPdfPBH( Magnitudes, MeanMags):

    
    #Get the deltas that correspond to the meanMags
    chosenNorms, chosenDeltas = \
      getNormAndDeltaForMagnitude(MeanMags)
      
    probablityGivenMagnitudes = \
      probPbhGivenMag( Magnitudes, chosenDeltas, chosenNorms )

    probablityGivenMagnitudes[ chosenNorms == 0] = 0

    return probablityGivenMagnitudes



class pbhPDF:

    def __init__( self, delta ):

        '''
        Using equation A3 from Z&S
        or eqyuation 8 from Rauch 1990
        
        '''
        self.delta = delta
        self.mag = np.linspace(0.0, 100.0, 10000000)[1:]
        
        self.dMag = self.mag[1]-self.mag[0]
        
       
        P = rauchFunct( self.mag, self.delta)
        
        self.norm = np.sum(P)*self.dMag

        self.pdf = P/self.norm
    
        self.mean = np.sum(self.pdf*self.dMag*self.mag)


def rauchFunct( mag, delta):
    '''
    Unormalisaed PDF from the rauch function
    for lensing of pbg
    '''

    topFrac = 1. - np.exp( -mag/delta)
    botFrac = (mag+1)**2 - 1.
        
    P = ( topFrac/botFrac)**(3./2.)
    return P
        
def probPbhGivenMag( mag, delta, norm):
    '''
    Return the probability of magnitude, given the delta
    and the normalisation
    '''
    
    return rauchFunct( mag, delta)/norm
if __name__ == '__main__':
    deltaTable()
    
