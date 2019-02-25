
import probabilityMagnification as pm
import deltaTable as dt
import sys
from matplotlib import pyplot as plt
import os
import pickle as pkl
import numpy as np

class lensingProbabilityDistribution():
    '''
    This class contains the functions that will take an output 
    from TurboGL which is the probability density function for
    weak lensing by large scale structure, and convovle it 
    with the probability distribution funcion of compact 
    objects. 

    Currently it reproduces the distributions in 
    https://www.sciencedirect.com/science/article/pii/S2212686418300293
    
    However ultimately it needs to be adapted to Equvalent widths.

    For now it outputs a probablity as a function of magnitude
    '''

    
    def __init__( self, redshift=1.0, alpha=0., nMagnitudeBins=1000):
        '''
        This distribution is a function of two parameters
        redshift : the redshift at which i want to calcuate the probability
                   distribution
        alpha : the fraction of Primordial Black Holes that make up 
                dark matter. Throughout the code it is known as pbh
        nMagnitudeBins : the number of magnitude bins that the code
                 will calcualte te convolution at and the final PDF
                 is sampled at. Default = 1000
        '''
        self.delta = dt.getDelta(redshift)
        self.redshift = redshift
        self.alpha = alpha
        self.nMagnitudeBins = nMagnitudeBins
        #to be changed to include redshift dependence:
        self.getProbabilityLensingByLss()
        self.pickklFileName = \
          'pickles/PL_%0.1f_%0.2f_%i.pkl' % (redshift,alpha,nMagnitudeBins)

        self.getMeanMagnitude()
        self.getDmuPrime()
        self.checkPickleFileExists()

        if not self.boolPickleFileExists:
            self.convolvePbhPdfWithLssPdf()


        
    def checkPickleFileExists( self ):
        if os.path.isfile( self.pickklFileName ):
            
            totalMagnitude, totalProbability = \
            pkl.load( open( self.pickklFileName, 'rb'))

            self.convolvedPbhPdfWithLssPdf = \
            {'x':totalMagnitude, 'y':totalProbability}

        self.boolPickleFileExists =\
          os.path.isfile( self.pickklFileName )

    def getMeanMagnitude( self ):
        self.meanMagnification = pm.getMeanMag( self.redshift )

    def makeArrayOfPbhProbabilityDistributions( self ):

        pdfMags, meanMagnitudes, PDFarray = \
          dt.generatePDFs()

        self.pbhProbabilityDistributionArray = \
          {'pdfMagnitudes':pdfMags ,\
               'meanMagnitudes':meanMagnitudes, \

               'pdfDistributionArray':meanMagnitudes}
        

    def convolvePbhPdfWithLssPdf( self ):
        '''
        This is where the meat of the work is done
        '''
        self.makeArrayOfPbhProbabilityDistributions()
        
        self.convolvedPbhPdfWithLssPdf = \
           np.zeros(self.nMagnitudeBins-1)
           
        self.totalLensingPDFmagnitudes = \
          np.linspace(0., 1.0, self.nMagnitudeBins)[1:]
          
        for howFarThroughList, givenMagnitude \
          in enumerate(self.totalLensingPDFmagnitudes):
            self.reportProgress(howFarThroughList)
            self.convolveLssWithPbhForGivenMagnitude( givenMagnitude )
            
        self.totalProbabilityDistribution = \
          {'x':self.totalLensingPDFmagnitudes, \
               'y':self.convolvedPbhPdfWithLssPdf}

    def reportProgress( self, progress):
        sys.stdout.write("%i/%i" %(progress, self.nMagnitudeBins-1))
        sys.stdout.flush()
        
    def convolveLssWithPbhForGivenMagnitude( self, magnitude ):

        #Get the end of the integral, but to stop
        #it going on forvere i curtail it here.
        endInt = np.min([magnitude / (1.-self.alpha), 10.])
        
        #dMuPrime = endInt/nInt
        #dMuPrime should be the P_LSS one as this defines min
        #Does this need to be changed?
        MuPrimeList = np.arange(0., endInt,  self.dMuPrime) 
     
        self.probabilityLensedByCompactObject =  \
          dt.getPdfPBH(  magnitude - MuPrimeList*(1.-self.alpha), \
                             self.alpha*MuPrimeList,\
                             self.pbhProbabilityDistributionArray )


                                  
        self.getProbabilityLensingByLssForGivenMagnitude(MuPrimeList)
                                
        
        self.dP_L = self.dMuPrime*\
          self.probabilityLensingByLssForGivenMagnitude*\
          self.probabilityLensedByCompactObject
    
    
        self.totalProbabilityForGivenMagnitude = np.sum(self.dP_L)

        index = self.totalLensingPDFmagnitudes == magnitude
        
        self.convolvedProbabilityDistribution[ index ] = \
          totalProbabilityForGivenMagnitude

    def getDmuPrime(self):
        self.getProbabilityLensingByLss()
        self.dMuPrime = self.probabilityLensingByLss['x'][1]-\
                    self.probabilityLensingByLss['x'][0]



    def getProbabilityLensingByLss(self):


        TurboGLfile = 'TurboGL.dat'
        magnitudes, PDF = np.loadtxt(TurboGLfile, unpack=True)

        self.probabilityLensingByLss = \
          {'x':magnitudes, 'y':PDF}

    
    def getProbabilityLensingByLssForGivenMagnitude( self, mu ):
        

        pdfMagsMatrix = np.matrix(self.probabilityLensingLss['x']).T*\
          np.matrix(np.ones(len(mu)))
        
        pdfMagsIndex = \
          np.argmin(np.abs((self.probabilityLensingLss['y'] - \
                                np.matrix(mu))),\
                        axis=0)
    
        correspondingMag = self.probabilityLensingLss['x'][pdfMagsIndex]
        
        returnPDF = pdf[pdfMagsIndex]
        returnPDF[ correspondingMag-mu > 2.*(mag[1]-mag[0])] = 0.
    
        self.probabilityLensingByLssForGivenMagnitude = returnPDF
        
    


    def plotTotalProbabilityDistribution(self):
        #TOTAL PDF
        plotMagnitudesAppendZero = \
          np.append(0., self.convolvedPbhPdfWithLssPdf['x']) - \
          self.meanMagnification
        plotPdfAppendZero = \
          np.append(1e-3,self.convolvedPbhPdfWithLssPdf['y'])
          
        plt.plot(plotMagnitudesAppendZero, plotPdfAppendZero)

        #LSS PDF
        plt.plot(self.probabilityLensingByLss['x'],\
                     self.probabilityLensingByLss['y'], '--')
        
        #PBH PDF
        magPBH, pdfPBH =  pm.magPDF( self.delta )
        plt.plot(np.append(0.,magPBH)-self.meanMagnification, \
                     np.append(1e-5,pdfPBH), ':')
        

        plt.yscale('log')
        plt.ylim(0.05,35)
        plt.xlim(0.,0.6)
        plt.show()
