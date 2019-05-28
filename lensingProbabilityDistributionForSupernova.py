
import probabilityMagnification as pm
import deltaTable as dt
import sys
from matplotlib import pyplot as plt
import os
import pickle as pkl
import numpy as np
import ipdb as pdb
'''
This class is used to create the lensing probability distrinutio
assuming lened supernova
i.e. the probability distirubtoi in the papers cited below
'''

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

    Notes:
    - This is now being changed to (1+kappa) instead of mu
    '''

    
    def __init__( self, redshift=1.0, alpha=0., nKappaBins=1000, modelType='Log'):
        '''
        This distribution is a function of two parameters
        redshift : the redshift at which i want to calcuate the probability
                   distribution
        alpha : the fraction of Primordial Black Holes that make up 
                dark matter. Throughout the code it is known as pbh
        nKappaBins : the number of kappa bins that the code
                 will calcualte te convolution at and the final PDF
                 is sampled at. Default = 1000
        modelType : Whether the PDFs have been worked through in log or not.
        '''
        self.delta = dt.getDelta(redshift)
        self.redshift = redshift
        self.alpha = alpha
        self.nKappaBins = nKappaBins
        self.getMeanKappa()
        self.modelType = modelType
        
        #to be changed to include redshift dependence:
        self.getProbabilityLensingByLss()
        
        self.pickleFileName = \
          'pickles/PL_Sn_kappa_%s_%0.1f_%0.2f_%i.pkl' % \
          (self.modelType, redshift, alpha, nKappaBins)

        self.getDkappaPrime()
        self.checkPickleFileExists()
        
        if not self.boolPickleFileExists:
            self.convolvePbhPdfWithLssPdf()

        self.normalisePDF()
        self.getPDFmean()
        
        if not self.boolPickleFileExists:
            self.writeToPickleFile()
            


    def convolvePbhPdfWithLssPdf( self ):
        '''
        This is where the meat of the work is done
        Here I will take a list of magnitudes at which 
        i will calcualte the probabilty of the lensing by
        LSS and compact objects. 
        
        For a given magnitude what is the lensing by PBH and LSS
        
        '''
        #self.makeArrayOfPbhProbabilityDistributions()
        
        self.totalConvolvedPbhPdfWithLssPdf = \
           np.zeros(self.nKappaBins-1)
           
        self.totalLensingPDFamplifications = \
          np.linspace(0., 1.0, self.nKappaBins)[1:]
          
        self.dAmplification = self.totalLensingPDFamplifications[1] - \
          self.totalLensingPDFamplifications[0]
        
        if self.alpha == 0:
            self.convolvedPbhPdfWithLssPdf =  \
              {'x':self.probabilityLensingByLss['x'],  \
               'y':self.probabilityLensingByLss['y']}
        else:
            for howFarThroughList, givenAmplification \
              in enumerate(self.totalLensingPDFamplifications):
                self.reportProgress(howFarThroughList)
                if self.modelType == 'Log':
                    self.convolveLssWithPbhForGivenAmplificationInLog( givenAmplification )
                else:
                    self.convolveLssWithPbhForGivenAmplification( givenAmplification )
                    
        
            self.convolvedPbhPdfWithLssPdf =  \
            {'x':self.totalLensingPDFamplifications,  \
                'y':self.totalConvolvedPbhPdfWithLssPdf}

 
        
    def convolveLssWithPbhForGivenAmplification( self, amplification ):
        '''
        For a given input magnitude what is the 
        probability that you will be lensed by compact objects
        and LSS
        '''

        kappa = amplification  / 2.
        #Get the end of the integral, but to stop
        #it going on forvere i curtail it here.
        endInt = np.min([kappa / (1.-self.alpha), 10.])
        
        #dMuPrime = endInt/nInt
        #dMuPrime should be the P_LSS one as this defines min
        #Does this need to be changed?

        MinNkappaPrime = np.max([10., endInt / self.dKappaPrime])
        KappaPrimeList = np.linspace(0.,endInt,MinNkappaPrime)[1:-1]
        
        self.probabilityLensedByCompactObject =  \
          dt.getPdfPBH(  (kappa - KappaPrimeList*(1.-self.alpha))*2., \
                             self.alpha*KappaPrimeList*2. )

        
        #
        self.getProbabilityLensingByLssForGivenKappa(KappaPrimeList)
                                
        
        self.dP_L = self.dKappaPrime*\
          self.probabilityLensingByLssForGivenKappa*\
          self.probabilityLensedByCompactObject
    
          
        self.totalProbabilityForGivenAmplification = np.sum(2.*self.dP_L)

        index = self.totalLensingPDFamplifications == amplification
        
        self.totalConvolvedPbhPdfWithLssPdf[ index ] = \
          self.totalProbabilityForGivenAmplification
        
        
    def convolveLssWithPbhForGivenAmplificationInLog( self, amplification ):
        '''
        For a given input kappa what is the 
        probability that you will be lensed by compact objects
        and LSS in log form
        '''

        kappa = amplification / 2.
        
        #Get the end of the integral, but to stop
        #it going on forvere i curtail it here.
        startInt = 0.
        endInt = np.min([kappa / (1.-self.alpha), 10.])

        #dMuPrime = endInt/nInt
        #dMuPrime should be the P_LSS one as this defines min
        #Does this need to be changed?

        MinNkappaPrime = np.max([10., endInt / self.dKappaPrime])
        KappaPrimeList = np.linspace(0., endInt, MinNkappaPrime)[1:-1]
        
        #self.probabilityLensedByCompactObject =  \
        #  dt.getPdfPBH(  amplification*(1 - 2.*KappaPrimeList*(1.-self.alpha)), \
        #                     self.alpha*KappaPrimeList*2. )

        self.probabilityLensedByCompactObject =  \
          dt.getPdfPBH(  (kappa - KappaPrimeList*(1.-self.alpha))*2., \
                             self.alpha*KappaPrimeList*2. )
        #
        self.getProbabilityLensingByLssForGivenKappa(KappaPrimeList)
                                
        
        self.dP_L = self.dKappaPrime*\
          self.probabilityLensingByLssForGivenKappa*\
          self.probabilityLensedByCompactObject

        #this extra term comes from the fact that we now use
        #logs
        
        extraTerm = (1.-(1.-self.alpha)*2.*KappaPrimeList)
        self.totalProbabilityForGivenAmplification = np.sum(self.dP_L*extraTerm)

        index = self.totalLensingPDFamplifications == amplification
        
        self.totalConvolvedPbhPdfWithLssPdf[ index ] = \
          self.totalProbabilityForGivenAmplification
        

          
    def getDkappaPrime(self):
        '''
        For the convoutions to work the dMagnitudePrime needs
        to be the largest common demonator, i.e. 
        that given by the LSS (turboGL)
        '''
        
        self.getProbabilityLensingByLss()
        self.dKappaPrime = (self.probabilityLensingByLss['x'][1]-\
                    self.probabilityLensingByLss['x'][0]) / 32.



    def getProbabilityLensingByLss(self):
        '''
        Read the asciss file output by TurboGL
        to get the lensing by large scale structure 
        on lines of sight therough the universe

        TurboGL produces a PDF that is relative to 
        zero. I need to add on the mean magnituification
        to match the PBH PDF range.
        '''

        TurboGLfile = 'TurboGL_Kappa.dat'
        kappa, PDF = np.loadtxt(TurboGLfile, unpack=True)
        
        self.probabilityLensingByLss = \
          {'x':kappa+self.meanKappa, 'y':PDF}

    
    def getProbabilityLensingByLssForGivenKappa( self, kappaList ):
        '''
        For a list of kappas interpolated and find
        the approprtiat
        kappaList : an array of kappas I want to find the 
        correspodnig probabilty to in TUrboGL PDF>
        '''
        self.probabilityLensingByLssForGivenKappa = \
          np.interp( kappaList, self.probabilityLensingByLss['x'],\
                         self.probabilityLensingByLss['y'])
        self.probabilityLensingByLssForGivenKappa[ kappaList < \
                                self.probabilityLensingByLss['x'][0]] = 0.
        



    def writeToPickleFile( self ):
        '''
        Write the determined probability distribution
        to a pickle file so save time in the future

        '''
        self.__model__ = "lensingProbabilityDistributionForSupernova"
    
        pkl.dump( self.__dict__, open( self.pickleFileName, 'wb'))


    def normalisePDF(self):
        '''
        Normalise the pdf so the integral is 1
        '''
        
        normalisation = np.sum( self.convolvedPbhPdfWithLssPdf['y']) *self.dAmplification
           
        self.convolvedPbhPdfWithLssPdf['y'] /= normalisation
        
    def getPDFmean( self ):

        self.pdfMean = np.sum(self.convolvedPbhPdfWithLssPdf['y']*\
                            self.convolvedPbhPdfWithLssPdf['x'])*\
                            self.dAmplification
                            
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
            temp =  pkl.load( open( self.pickleFileName, 'rb'))
            self.__dict__.update(temp)
       
        self.boolPickleFileExists =  os.path.isfile( self.pickleFileName )

    def getMeanKappa( self ):
        '''
        Determine the mean magnification 
        for the given redsihft (the FLRW magnification given 
        the existence of a "empty beam")
        '''
        
        self.meanKappa = pm.getMeanMag( self.redshift )/2.

    def reportProgress( self, progress):
        '''
        Just some stdout stuff to help the user
        '''
        
        sys.stdout.write("PROGRESS: %i/%i\r" \
                             %(progress+1, self.nKappaBins-1))
        sys.stdout.flush()
