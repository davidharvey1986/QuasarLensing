
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
    objects to produce the distributions of equivalent widths

    It is based in the code
    https://www.sciencedirect.com/science/article/pii/S2212686418300293
    
    but adapted for equivalent widths in linear

    Will need to adapted for log

    '''

    
    def __init__( self, redshift=1.0, alpha=0., nEquivalentWidthBins=1000):
        '''
        This distribution is a function of two parameters
        redshift : the redshift at which i want to calcuate the probability
                   distribution
        alpha : the fraction of Primordial Black Holes that make up 
                dark matter. Throughout the code it is known as pbh
        nEquivalentWidthBins : the number of equvalent width bins 
                  that the code will calcualte te convolution 
                  at and the final PDF is sampled at. Default = 1000
        '''
        self.delta = dt.getDelta(redshift)
        self.redshift = redshift
        self.alpha = alpha
        self.nEquivalentWidthBins = nEquivalentWidthBins
        self.getMeanMagnitude()

        #to be changed to include redshift dependence:
        self.getProbabilityLensingByLss()
        self.pickleFileName = \
          'pickles/PL_EW_%0.1f_%0.2f_%i.pkl' % \
          (redshift,alpha,nEquivalentWidthBins)

        self.getDmuPrime()
        self.checkPickleFileExists()

        if not self.boolPickleFileExists:
            self.convolvePbhPdfWithLssPdf()

            self.writeToPickleFile()

        #self.getConvolvedPbhPdfWithLssPdf()
            
        self.normalisePDF()
        self.getPDFmean()
        
    def normalisePDF(self):
        '''
        Normalise the pdf so the integral is 1
        '''

        normalisation = np.sum( self.convolvedPbhPdfWithLssPdf['y'])\
          *self.dEquivalentWidth
           
        self.convolvedPbhPdfWithLssPdf['y'] /= normalisation

    def getPDFmean( self ):

        self.pdfMean = np.sum(self.convolvedPbhPdfWithLssPdf['y']*\
                            self.convolvedPbhPdfWithLssPdf['x'])*\
                            self.dEquivalentWidth
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
            
            totalEquivalentWidth, totalProbability = \
            pkl.load( open( self.pickleFileName, 'rb'))

            self.convolvedPbhPdfWithLssPdf = \
            {'x':totalEquivalentWidth, 'y':totalProbability}
            self.dEquivalentWidth = totalEquivalentWidth[1]-totalEquivalentWidth[0]

        self.boolPickleFileExists = \
          os.path.isfile( self.pickleFileName )

    def getMeanMagnitude( self ):
        '''
        Determine the mean magnification 
        for the given redsihft (the FLRW magnification given 
        the existence of a "empty beam")
        '''
        
        self.meanMagnification = pm.getMeanMag( self.redshift )


    def convolvePbhPdfWithLssPdf( self ):
        '''
        This is where the meat of the work is done
        Here I will take a list of magnitudes at which 
        i will calcualte the probabilty of the lensing by
        LSS and compact objects. 
        
        For a given magnitude what is the lensing by PBH and LSS
        
        '''
        
        self.totalConvolvedPbhPdfWithLssPdf = \
           np.zeros(self.nEquivalentWidthBins-2)
           
        self.totalLensingPDFequivalentWidths = \
          np.linspace(-1., 1.0, self.nEquivalentWidthBins)[1:-1]
          
        self.dEquivalentWidth = \
          self.totalLensingPDFequivalentWidths[1]-\
          self.totalLensingPDFequivalentWidths[0]

        for howFarThroughList, givenEquivalentWidth \
          in enumerate(self.totalLensingPDFequivalentWidths):
            self.reportProgress(howFarThroughList)
            self.convolveLssWithPbhForGivenEquivalentWidth( givenEquivalentWidth )

        
        self.convolvedPbhPdfWithLssPdf =  \
            {'x':self.totalLensingPDFequivalentWidths,  \
                'y':self.totalConvolvedPbhPdfWithLssPdf}

    def reportProgress( self, progress):
        '''
        Just some stdout stuff to help the user
        '''
        
        sys.stdout.write("PROGRESS: %i/%i\r" \
            %(progress+1, self.nEquivalentWidthBins-1))
        sys.stdout.flush()
        
    def convolveLssWithPbhForGivenEquivalentWidth( self, equivalentWidth ):
        '''
        For a given input magnitude what is the 
        probability that you will be lensed by compact objects
        and LSS
        '''
        #Get the start of the integral
        startInt = np.max([0,equivalentWidth/self.alpha])
        
        #End of the integral is formally infinity
        #so will have to curtail some point
        endInt = 10 #Will have to do a convergence test
        

        MinNumEwPrime = np.max([10., (endInt-startInt) / self.dMuPrime])
        MinNEwPrime = 1000
        EwPrimeList = np.linspace(startInt,endInt,MinNEwPrime)[1:]
        
        self.probabilityLensedByCompactObject =  \
          dt.getPdfPBH( self.alpha*EwPrimeList-equivalentWidth,\
                            self.alpha*EwPrimeList)

        
        #
        self.getProbabilityLensingByLssForGivenMagnitude(EwPrimeList)
                                
        dEWPrime = EwPrimeList[1]-EwPrimeList[0]
        self.dP_L = dEWPrime*\
          self.probabilityLensingByLssForGivenMagnitude*\
          self.probabilityLensedByCompactObject
    
          
        self.totalProbabilityForGivenEquivalentWidth = np.sum(self.dP_L)

        index = self.totalLensingPDFequivalentWidths == equivalentWidth
        
        self.totalConvolvedPbhPdfWithLssPdf[ index ] = \
          self.totalProbabilityForGivenEquivalentWidth
    
      

    def getDmuPrime(self):
        '''
        For the convoutions to work the dMagnitudePrime needs
        to be the largest common demonator, i.e. 
        that given by the LSS (turboGL)
        '''
        
        self.getProbabilityLensingByLss()
        self.dMuPrime = (self.probabilityLensingByLss['x'][1]-\
                    self.probabilityLensingByLss['x'][0]) / 8.



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
        self.probabilityLensingByLssForGivenMagnitude = \
          np.interp( magnitudeList, self.probabilityLensingByLss['x'],\
                         self.probabilityLensingByLss['y'])
        
        self.probabilityLensingByLssForGivenMagnitude[ magnitudeList < \
                                        self.probabilityLensingByLss['x'][0]] = 0.

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


