'''
When creating the PDFs for the PBH it takes a while to get the right deltas

so i need a table outlining expected mean mag for a given delta.

'''

import probabilityMagnification as pm
import pickle as pkl
import numpy as np

def deltaTable():

    deltaList = np.linspace(0.,10.,100000)
    expecMag = []
    for delta in deltaList:
        
        mag, iPDF = pm.magPDF( delta )
        dMag = mag[1]-mag[-0]
        expecMag.append(np.sum(iPDF*mag*dMag))
        print np.sum(iPDF*mag*dMag)
    pkl.dump([deltaList, np.array(expecMag)], open('deltaTable.pkl','wb'))


def getDelta( z=1.0, meanMag=None):
    
    '''
    Derive the delta for the PBH 
    '''

    if meanMag is None:
        meanMag = pm.getMeanMag( z )

    deltaList, magnitudes = pkl.load(open('deltaTable.pkl','rb'))

    delta = deltaList[ np.argmin( np.abs( magnitudes - meanMag))]

    return delta
    
    

def generatePDFs():

    '''
    Generate PBH PDFs for each given delta (i.e. mean Mag)

    '''

    
    deltaList, MeanMagnitudes = pkl.load(open('deltaTable.pkl','rb'))

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

        

        
def getPdfPBH( Magnitudes, MeanMags, pdfMags, MeanMagnitudes, PDFarray):

    #This generates a 2x2 array of mea


    pdfMagsMatrix = np.matrix(pdfMags).T*np.matrix(np.ones(len(Magnitudes)))

    pdfMagsIndex =\
      np.argmin(np.abs((pdfMagsMatrix - np.matrix(Magnitudes))),axis=0)


    MeanMagsMatrix = np.matrix(MeanMagnitudes).T*np.matrix(np.ones(len(MeanMags)))

    MeanMagsIndex = np.argmin(np.abs((MeanMagnitudes - np.matrix(MeanMags).T)),axis=1).T

    
    
    return  PDFarray[pdfMagsIndex,MeanMagsIndex]
     
