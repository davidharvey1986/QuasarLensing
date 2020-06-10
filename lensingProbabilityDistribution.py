
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

    
    def __init__( self, **inputParams):
        '''
        This distribution is a function of two parameters
        redshift : the redshift at which i want to calcuate the probability
                   distribution
        alpha : the fraction of Primordial Black Holes that make up 
                dark matter. Throughout the code it is known as pbh
        nEquivalentWidthBins : the number of equvalent width bins 
                  that the code will calcualte te convolution 
                  at and the final PDF is sampled at. Default = 1000

        It has been adapted to rerturn log(1+mu_ew)
        
        updated: added cosmology and inputParams now contains all
        parameters
        '''
        self.inputParams = {'omega_m':0.3, 'sigma8':0.8288, \
                'redshift':1.0, 'alpha':0.83, \
                'nEquivalentWidthBins':1000, \
                'modelType':'Linear'}
        
        for iParam in inputParams.keys():
            #if iParam not in self.inputParams.keys():
                #print("WARNING: '%s' parameter no recognised" \
                #                     % iParam)

            self.inputParams[iParam] = inputParams[iParam]
        
        self.delta = \
          dt.getDelta(self.inputParams['redshift'],\
                    omega_m=self.inputParams['omega_m'], \
                    sigma_8=self.inputParams['sigma8'])

        self.getMeanMagnitude()

        if  self.inputParams['alpha'] < 0.08:
            raise ValueError("Alpha is too small, valid only >= 0.08")
        
        #to be changed to include redshift dependence:
        self.getProbabilityLensingByLss()
        self.pickleFileName = \
          'pickles/PL_EW_%s_%0.1f_%0.2f_%0.4f_%0.2f_%i.pkl' % \
          (self.inputParams['modelType'],self.inputParams['redshift'],\
               self.inputParams['omega_m'],self.inputParams['sigma8'],\
               self.inputParams['alpha'],\
               self.inputParams['nEquivalentWidthBins'])

        self.getDmuPrime()
        self.checkPickleFileExists()
            
        if not self.boolPickleFileExists:
            self.convolvePbhPdfWithLssPdf()

            self.writeToPickleFile()

            
        self.normalisePDF()
        self.getPDFmean()
    


    def convolvePbhPdfWithLssPdf( self, verbose=False ):
        '''
        This is where the meat of the work is done
        Here I will take a list of magnitudes at which 
        i will calcualte the probabilty of the lensing by
        LSS and compact objects. 
        
        For a given magnitude what is the lensing by PBH and LSS
        
        '''
        
        self.totalConvolvedPbhPdfWithLssPdf = \
           np.zeros(self.inputParams['nEquivalentWidthBins']-2)

        #I want to return in log(1+equivalentwidhts)
        self.totalLensingPDFequivalentWidths = \
          np.linspace(-2., 2., \
                self.inputParams['nEquivalentWidthBins'])[1:-1]
          
        self.dEquivalentWidth = \
          self.totalLensingPDFequivalentWidths[1]-\
          self.totalLensingPDFequivalentWidths[0]

        for howFarThroughList, givenEquivalentWidth \
          in enumerate(self.totalLensingPDFequivalentWidths):
            if verbose:
                self.reportProgress(howFarThroughList)
            self.index = howFarThroughList
            if self.inputParams['modelType'] == 'Log':
                self.convolveLssWithPbhForGivenEquivalentWidthInLog( 10**givenEquivalentWidth -1. )
            else:
                #All the convolution was written in linear space
                #So the input mu needs to be in linear mu not log
                self.convolveLssWithPbhForGivenEquivalentWidth( 10**givenEquivalentWidth - 1. )
       
        #Convert the pdf to log(pdf) with the jacobian
        self.totalConvolvedPbhPdfWithLssPdf *= \
          (10**self.totalLensingPDFequivalentWidths)*np.log(10.)

        self.convolvedPbhPdfWithLssPdf =  \
            {'x':self.totalLensingPDFequivalentWidths,  \
                'y':self.totalConvolvedPbhPdfWithLssPdf}
            


    def convolveLssWithPbhForGivenEquivalentWidth( self, equivalentWidth ):
        '''
        For a given input magnitude what is the 
        probability that you will be lensed by compact objects
        and LSS
        '''
        #Get the start of the integral
        startInt = np.max([0,equivalentWidth/self.inputParams['alpha']])
        
        #End of the integral is formally infinity
        #so will have to curtail some point
        endInt = np.max([10., startInt*10.]) #Will have to do a convergence test
        MinNEwPrime = 10000
        EwPrimeList = np.linspace(startInt,endInt,MinNEwPrime)[1:]
        
        self.probabilityLensedByCompactObject =  \
          dt.getPdfPBH( (self.inputParams['alpha']*EwPrimeList-\
                             equivalentWidth),\
                            self.inputParams['alpha']*EwPrimeList)

        
        #
        self.getProbabilityLensingByLssForGivenMagnitude(EwPrimeList)
                                
        dEWPrime = EwPrimeList[1]-EwPrimeList[0]
        self.dP_L = dEWPrime*\
          self.probabilityLensingByLssForGivenMagnitude*\
          self.probabilityLensedByCompactObject
    
          
        self.totalProbabilityForGivenEquivalentWidth = np.sum(self.dP_L)

        
        self.totalConvolvedPbhPdfWithLssPdf[ self.index ] = \
          self.totalProbabilityForGivenEquivalentWidth
    
    def convolveLssWithPbhForGivenEquivalentWidthInLog( self, equivalentWidth ):
        '''
        For a given input magnitude what is the 
        probability that you will be lensed by compact objects
        and LSS 

        In the maths where we assume log
        '''
        #Get the start of the integral

        startInt = np.max([0,equivalentWidth/self.alpha])
        
        #End of the integral is formally infinity
        #so will have to curtail some point
        endInt = np.max([10., startInt*10.]) #Will have to do a convergence test
        MinNEwPrime = 10000
        EwPrimeList = np.linspace(startInt,endInt,MinNEwPrime)[1:]
        
        self.probabilityLensedByCompactObject =  \
          dt.getPdfPBH( (self.alpha*EwPrimeList-equivalentWidth),\
                            self.alpha*EwPrimeList)

        
        #
        self.getProbabilityLensingByLssForGivenMagnitude(EwPrimeList)
                                
        dEWPrime = EwPrimeList[1]-EwPrimeList[0]

        extraTerm = (1.+2.*EwPrimeList)/ \
          (1.+(1.-self.alpha)*2.*EwPrimeList)
          
        self.dP_L = dEWPrime*\
          self.probabilityLensingByLssForGivenMagnitude*\
          self.probabilityLensedByCompactObject*\
          extraTerm
          
    
          
        self.totalProbabilityForGivenEquivalentWidth = \
          np.sum(self.dP_L)/(1.+2.*equivalentWidth)**2

        
        self.totalConvolvedPbhPdfWithLssPdf[ self.index ] = \
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

        TurboGLfile = \
          '../turboGL3.0/results/hpdf_z_%0.1f_OM_%0.1f_S8_%0.1f.dat' %\
          (self.inputParams['redshift'],\
            self.inputParams['omega_m'], self.inputParams['sigma8'])

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
        
        self.__module__ = "lensingProbabilityDistribution"
        pkl.dump(self.__dict__, open( self.pickleFileName, 'wb'))


    
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
            
            temp =  pkl.load( open( self.pickleFileName, 'rb'), encoding='latin1')
            self.__dict__.update(temp)
            
        self.boolPickleFileExists = \
          os.path.isfile( self.pickleFileName )

    def getMeanMagnitude( self ):
        '''
        Determine the mean magnification 
        for the given redsihft (the FLRW magnification given 
        the existence of a "empty beam")
        '''
        
        self.meanMagnification = \
          pm.getMeanMag( self.inputParams['redshift'] )
        
    def reportProgress( self, progress):
        '''
        Just some stdout stuff to help the user
        '''
        
        sys.stdout.write("PROGRESS: %i/%i\r" \
            %(progress+1, self.inputParams['nEquivalentWidthBins']-1))
        sys.stdout.flush()
        
