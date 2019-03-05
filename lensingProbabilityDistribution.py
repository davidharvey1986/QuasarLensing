
import probabilityMagnification as pm
import deltaTable as dt
import sys
from matplotlib import pyplot as plt
import os
import pickle as pkl
import numpy as np
import ipdb as pdb


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
        self.getMeanMagnitude()

        #to be changed to include redshift dependence:
        self.getProbabilityLensingByLss()
        self.pickleFileName = \
          'pickles/PL_%0.1f_%0.2f_%i.pkl' % (redshift,alpha,nMagnitudeBins)

        self.getDmuPrime()
        self.checkPickleFileExists()

        if not self.boolPickleFileExists:
            self.convolvePbhPdfWithLssPdf()

            self.writeToPickleFile()

        self.getConvolvedPbhPdfWithLssPdf()
            

    def getConvolvedPbhPdfWithLssPdf( self ):
        return self.convolvedPbhPdfWithLssPdf['x'], \
          self.convolvedPbhPdfWithLssPdf['y']
          
    def checkPickleFileExists( self ):
        '''
        See if the data file for this PDF exists.
        Load up the one that exists, otherwise 
        let the code know with a bool that we 
        need to do the pdf or not
        '''
        
        if os.path.isfile( self.pickleFileName ):
            
            totalMagnitude, totalProbability = \
            pkl.load( open( self.pickleFileName, 'rb'))

            self.convolvedPbhPdfWithLssPdf = \
            {'x':totalMagnitude, 'y':totalProbability}

        self.boolPickleFileExists = \
          os.path.isfile( self.pickleFileName )

    def getMeanMagnitude( self ):
        '''
        Determine the mean magnification 
        for the given redsihft (the FLRW magnification given 
        the existence of a "empty beam")
        '''
        
        self.meanMagnification = pm.getMeanMag( self.redshift )

    def makeArrayOfPbhProbabilityDistributions( self ):
        '''
        Generate a 2D array of probability distributions where
        one axis is PDF as afuncion of magnitude
        Then this is done for a variety of mean magnitudes
        
        '''

        pdfMags, meanMagnitudes, PDFarray = \
          dt.generatePDFs()

        self.pbhProbabilityDistributionArray = \
          {'pdfMagnitudes':pdfMags ,\
               'meanMagnitudes':meanMagnitudes, \
               'pdfDistributionArray':PDFarray}
        

    def convolvePbhPdfWithLssPdf( self ):
        '''
        This is where the meat of the work is done
        Here I will take a list of magnitudes at which 
        i will calcualte the probabilty of the lensing by
        LSS and compact objects. 
        
        For a given magnitude what is the lensing by PBH and LSS
        
        '''
        self.makeArrayOfPbhProbabilityDistributions()
        
        self.totalConvolvedPbhPdfWithLssPdf = \
           np.zeros(self.nMagnitudeBins-1)
           
        self.totalLensingPDFmagnitudes = \
          np.linspace(0., 1.0, self.nMagnitudeBins)[1:]
          
        if self.alpha == 0:
            self.convolvedPbhPdfWithLssPdf =  \
              {'x':self.probabilityLensingByLss['x'],  \
               'y':self.probabilityLensingByLss['y']}
        else:
            for howFarThroughList, givenMagnitude \
              in enumerate(self.totalLensingPDFmagnitudes):
                self.reportProgress(howFarThroughList)
                self.convolveLssWithPbhForGivenMagnitude( givenMagnitude )

        
            self.convolvedPbhPdfWithLssPdf =  \
            {'x':self.totalLensingPDFmagnitudes,  \
                'y':self.totalConvolvedPbhPdfWithLssPdf}

    def reportProgress( self, progress):
        '''
        Just some stdout stuff to help the user
        '''
        
        sys.stdout.write("PROGRESS: %i/%i\r" \
                             %(progress+1, self.nMagnitudeBins-1))
        sys.stdout.flush()
        
    def convolveLssWithPbhForGivenMagnitude( self, magnitude ):
        '''
        For a given input magnitude what is the 
        probability that you will be lensed by compact objects
        and LSS
        '''

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


        #
        self.getProbabilityLensingByLssForGivenMagnitude(MuPrimeList)
                                
        
        self.dP_L = self.dMuPrime*\
          self.probabilityLensingByLssForGivenMagnitude*\
          self.probabilityLensedByCompactObject
    
          
        self.totalProbabilityForGivenMagnitude = np.sum(self.dP_L)

        index = self.totalLensingPDFmagnitudes == magnitude
        
        self.totalConvolvedPbhPdfWithLssPdf[ index ] = \
          self.totalProbabilityForGivenMagnitude

      

    def getDmuPrime(self):
        '''
        For the convoutions to work the dMagnitudePrime needs
        to be the largest common demonator, i.e. 
        that given by the LSS (turboGL)
        '''
        
        self.getProbabilityLensingByLss()
        self.dMuPrime = (self.probabilityLensingByLss['x'][1]-\
                    self.probabilityLensingByLss['x'][0])/4.



    def getProbabilityLensingByLss(self):
        '''
        Read the asciss file output by TurboGL
        to get the lensing by large scale structure 
        on lines of sight therough the universe

        TurboGL produces a PDF that is relative to 
        zero. I need to add on the mean magnituification
        to match the PBH PDF range.
        '''

        TurboGLfile = 'TurboGL.dat'
        magnitudes, PDF = np.loadtxt(TurboGLfile, unpack=True)
        
        self.probabilityLensingByLss = \
          {'x':magnitudes+self.meanMagnification, 'y':PDF}

    
    def getProbabilityLensingByLssForGivenMagnitude( self, magnitudeList ):
        '''
        For a list of magnitudes, mu find where for each
        the closest magnitude lies in the Lss PDF
        and find the corresponding probablity of each magnitude in
        the list of being lensed.

        magnitudeList : an array of magnitudes I want to find the 
        correspodnig probabilty to in TUrboGL PDF>
        '''
        pdfMagsMatrix = np.matrix(self.probabilityLensingByLss['x']).T*\
          np.matrix(np.ones(len(magnitudeList)))
        
        pdfMagsIndex = \
          np.argmin(np.abs((pdfMagsMatrix - np.matrix(magnitudeList))),\
                        axis=0)
    
        correspondingMag = self.probabilityLensingByLss['x'][pdfMagsIndex]
        
        returnPDF = self.probabilityLensingByLss['y'][pdfMagsIndex]
        returnPDF[ correspondingMag-magnitudeList > 2.*self.dMuPrime] = 0.
    
        self.probabilityLensingByLssForGivenMagnitude = returnPDF
        
    


    def plotTotalProbabilityDistribution(self, show=False):
        '''
        Plot the total convovled PDF
        Plot the LSS PDF
        Plot the PDF by compact objects
        '''
        
        #TOTAL PDF
        plt.plot(self.convolvedPbhPdfWithLssPdf['x'], self.convolvedPbhPdfWithLssPdf['y'])

        #LSS PDF
        plt.plot(self.probabilityLensingByLss['x'],\
                     self.probabilityLensingByLss['y'], '--')
        
        #PBH PDF
        magPBH, pdfPBH =  pm.magPDF( self.delta )
        plt.plot(np.append(0.,magPBH), \
                     np.append(1e-5,pdfPBH), ':')
        

        plt.yscale('log')
        plt.ylim(0.05,35)
        plt.xlim(0.,0.6)
        if show:
            plt.show()


    def writeToPickleFile( self ):
        '''
        Write the determined probability distribution
        to a pickle file so save time in the future

        '''
        
    
        pkl.dump([self.convolvedPbhPdfWithLssPdf['x'],\
                      self.convolvedPbhPdfWithLssPdf['y']],\
                     open( self.pickleFileName, 'wb'))


