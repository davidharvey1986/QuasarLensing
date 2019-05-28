import numpy as np
import pyfits as fits
import convolveLssPc as clp
from scipy.stats import norm
from scipy.stats import rayleigh
import cosmolopy.distance as dist
import lensingProbabilityDistribution as lpd
import ipdb as pdb
from matplotlib import pyplot as plt
class emissionLine:
    '''
    log the all the data from the shen catalopgues
    in to this class which will just make things look a bit cleaner

    also do all the convolution with the lensing pdfs!

    '''
    def __init__( self, emissionLine,\
                      minRedshift=0.1, \
                      maxRedshift=1.0, \
                      nRedshiftBins=5,\
                      catalogue='dr9_weq.fits'):
        '''
        Init all the pars, 
        emissionLine: name of the emission line
        dataArray: the inputs fits from the shen cat
        '''
        dataArray = fits.open(catalogue)[1].data
        #init the keywords
        self.minRedshift=minRedshift
        self.maxRedshift=maxRedshift
        self.nRedshiftBins=nRedshiftBins
        
        #first check that the name exists
        self.emissionLine = emissionLine
        self.completeDataListNames =\
          dataArray.columns.names

          
        self.checkEmissionLine()
        #now set in the class
        #init the redshifts of all the quasars
        
       
        self.redshift = dataArray['Z_PIPE']
        self.getLuminosity( dataArray )

        #init the rest frame equvilent width
        
        self.restFrameEquivilentWidth = \
          dataArray['REW_'+self.emissionLine]
        
        self.cleanEquivilentWidths()
        self.setRedshiftBins( )

    def getLuminosity( self, dataArray ):
        '''
        Get the luminsoity for each line 
        '''
        cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.7}
        cosmo = dist.set_omega_k_0(cosmo)  
        luminosityDistance = \
          dist.luminosity_distance(self.redshift, **cosmo)
        self.logLuminosities = {}
        
        for iName in self.completeDataListNames:
            if ('LOGF' in iName) & ('ERR' not in iName):

                    self.logLuminosities[iName] = dataArray[iName] + \
                      np.log(4.*np.pi*luminosityDistance**2)
  
    

    def checkEmissionLine( self ):
        '''
        Check that the input name exists in the dictionary
        of column names
        Otherwise reaise a flag
        '''

        exactFlag = [ 'REW_'+self.emissionLine == i for i in \
                          self.completeDataListNames]
        flag = [ 'REW_'+self.emissionLine in i for i in  \
                     self.completeDataListNames]
        
        if (not np.any(np.array(flag))) & \
          (not np.any(np.array(exactFlag))):
            raise ValueError("Cant find emission line name")
        elif  (np.any(np.array(flag))) & \
          (not np.any(np.array(exactFlag))):
          proposedName = np.array(self.completeDataListNames)[flag][0]
          raise ValueError("Can't find exact emission line name,"+\
                               " did you mean %s?" % proposedName)
                             

    def applyLuminosityCut( self, luminosityCut ):
        '''
        Apply a luminosity cut at a particularly emission line
        '''
        index =  self.restFrameEquivilentWidth > luminosityCut

        self.restFrameEquivilentWidth = \
            self.restFrameEquivilentWidth[index]
        self.redshift = self.redshift[index]
        print("After Luminosity Cut there are %i quasars" % len(self.redshift))
    def cleanEquivilentWidths( self ):
        '''
        I need to clean out all the rubbish equivilent widths
        as some are a load of rubbish
        '''
        index = (self.restFrameEquivilentWidth>0) & \
          (self.restFrameEquivilentWidth<1e3) & \
          (self.restFrameEquivilentWidth>0) & \
          (self.restFrameEquivilentWidth<1e3)
          
        #init the redshifts of all the quasars
        self.redshift = self.redshift[index]
        #init the rest frame equvilent width
        
        self.restFrameEquivilentWidth = \
          self.restFrameEquivilentWidth[index]

        
    def setRedshiftBins( self ):
        '''
        Set the redshift bins that will be used for all
        the analysis

        '''
        bins = np.linspace(self.minRedshift, \
                               self.maxRedshift, \
                               self.nRedshiftBins+1)
                            
        self.redshiftBins = bins
        self.redshiftCentres = (bins[1:] + bins[:-1])/2.

        
    def getEquvilentWidthMeansAsFunctionOfRedshift(self):
        '''
        Take an in out of x and y and bin it with nBins, from 0 to maxRedshift.
        Returning a vector contaiing hte centers of the bins 
        and an 2xnBins array containing the means and the std deviation 
        
        Inputs : x and y are N length equal vectors to be binned
        Outputs : centres : a vector of nBins long
        MeanWeq : a 2xnBins array of the means and the std deviation

        Note: Should be clean vectors with no issues in thenm
        '''
        MeanWeq = np.zeros((2, self.nRedshiftBins))

        self.binIndices = np.zeros(len(self.redshift)) - 1
        for iBin in xrange(self.nRedshiftBins):
            
            inBin = (self.redshift >self.redshiftBins[iBin]) & \
              (self.redshift <self.redshiftBins[iBin+1])

            
            MeanWeq[0, iBin] = \
              np.nanmean(self.restFrameEquivilentWidth[inBin])
            MeanWeq[1, iBin] = \
              np.nanstd(self.restFrameEquivilentWidth[inBin])/\
              np.sqrt(len(self.restFrameEquivilentWidth[inBin]))
            
            self.binIndices[ inBin ] = iBin
            
        self.meanEquivilentWidth = MeanWeq[0, :]
        self.meanEquivilentError = MeanWeq[1, :]



    def histogramEquivalentWidth(self, nEquivilentWidthBins=10):
        '''
        Hisogram each redshift bin and turn them into
        a dictionary set within a dictionary set 
        i.e. [0.1].x

        '''

        self.equivilentWidthHistograms = {}

        for iBin in xrange(self.nRedshiftBins):
            selectRestFrameEquivilentWidth = \
              self.restFrameEquivilentWidth[self.binIndices==iBin]

            
            yHistogram, xHistogram = \
              np.histogram(selectRestFrameEquivilentWidth, \
                               bins=nEquivilentWidthBins, \
                               density=True)
            redshiftBinName = "%0.2f" % self.redshiftCentres[iBin]
            histogramCentres = (xHistogram[1:] + xHistogram[:-1])/2.
            iRedshiftDict = {'x':histogramCentres, \
                                 'y':yHistogram}
            self.equivilentWidthHistograms[redshiftBinName] = iRedshiftDict
            
    def cumulativeDistributionEquivalentWidth(self):
        '''
        Hisogram each redshift bin and turn them into
        a dictionary set within a dictionary set 
        i.e. [0.1].x

        '''

        self.equivilentWidthCumSum = {}

        for iBin in xrange(self.nRedshiftBins):
            selectRestFrameEquivilentWidth = \
              self.restFrameEquivilentWidth[self.binIndices==iBin]

            probabilites = np.ones(len(selectRestFrameEquivilentWidth))\
              /len(selectRestFrameEquivilentWidth)
            cumSumEquivalentWidths = np.cumsum(probabilites)
            
         
            redshiftBinName = "%0.2f" % self.redshiftCentres[iBin]

            iRedshiftDict = {'x':np.sort(selectRestFrameEquivilentWidth), \
                                 'y':cumSumEquivalentWidths}
            self.equivilentWidthCumSum[redshiftBinName] = iRedshiftDict  

        
            
        
    def getIntrinsicDistribution( self, redshiftCut=0.3, intrinsicDistribution='Gaussian'):
        '''
        In the lensing paper they assume that the PDF of the
        lensed supernova is a Gaussian with a sigma of 0.15
        
        In quasars we do not have the same luxury. 
        So to get the 'intrinsic distribution' i will use the
        low redshift quasars ( z<redshiftCut )

        

        Note that to do this the dEquivalentWidth (increment in EW)
        must be the same as the lensing probability distribution
        that we will use to convovlve it with

        Therefore the bins should be preset

        mimicSuperNova : Instead of using the raw data of quasar distribtion
           mimic the distributions from the supenova papers, and 
           therefor a Gauassian with a variance=0.15
        
        '''

        

        #This is purely for test purposes
        if intrinsicDistribution == 'Gaussian':
            x = np.arange(-3.0,3.0,  self.dEquivalentWidth)
            y = norm.pdf(x, 0., 0.15)
            #these need to be symmetric otherwise, the position of the final
            #convoution is weird since the fft switches things around,
            #I could probably figure it out but for now this will do

            #Also if things go NaN then try making this symmetric in even or
            #odd numbers i.e. (-1,1) or (-2,2)
        if intrinsicDistribution == 'Rayleigh':
            x = np.arange(-5.0, 5.0,  self.dEquivalentWidth)
            y = rayleigh.pdf(x, 0., 0.15)

        if intrinsicDistribution == 'Semi-Gaussian':
            x = np.arange(-31.0,31.0,  self.dEquivalentWidth)
            y = norm.pdf(x, 0., 5.0)
            y[x < 0] = 0.

                
        
        if intrinsicDistribution == 'data':
            #Get the data
            self.intrinsicEquivalentWidths =  \
              self.restFrameEquivilentWidth[ (self.redshift < redshiftCut) &\
                                (np.log10(self.restFrameEquivilentWidth) >-2) &\
                                (np.log10(self.restFrameEquivilentWidth) <2) ]
                                            
            
            params = norm.fit( np.log10(self.intrinsicEquivalentWidths))
            x = np.arange(-3., 3., self.dEquivalentWidth)

            y = norm.pdf( x, *params)
            
            
        self.intrinsicEquivalentWidthDistribution = \
          {'y':y, 'x':x}


    def getLensingProbability( self, z=1.0, alpha=0.83 ):
        '''
        I need to the probability distrubtion that a quasar
        of redshift z is magnified or demagnified by large scale
        structure and primordial black holes. Where the fraction
        of dark matter in primordial black holes is alpha

        z: Quasar Redshift
        alpha: fraction of dark matter is primordial black holes

        NOTES: for now i will use a given redshift of 1. and alpha=0.27
        this will need to be changed such that i can iterate and fit

        '''
        #TO DO This needs to be checked this funcion.
        lensingPDF = \
          lpd.lensingProbabilityDistribution( redshift=z, \
                                              alpha=alpha, \
                                              nEquivalentWidthBins=1000,\
                                              modelType='Linear')
       
        #force some kind of normalisation,
        #although this should alrady be normalised
        #lensingPDF /= np.sum(lensingPDF)*(magnitude[1]-magnitude[0])
        self.dEquivalentWidth = lensingPDF.dEquivalentWidth*2.

        #this probability distribution is the effect of kappe_ew
        #the amplictude will be 1+2kappa
        self.lensingEquivalentWidth = \
          2.*lensingPDF.convolvedPbhPdfWithLssPdf['x'] 
        #times the dEW / dKappa
        self.lensingPDF = lensingPDF.convolvedPbhPdfWithLssPdf['y']*2.
        self.alpha=alpha
        
        self.lensingEquivalentWidth = \
          self.lensingEquivalentWidth[ self.lensingPDF > 0]
          
        self.lensingPDF = self.lensingPDF[ self.lensingPDF > 0 ]
        
    def convolveIntrinsicEquivalentWidthWithLensingProbability(self):
        '''
        In order to get an estimate of the expected distribution of 
        equivalent widths i need to take the intrinsic distrubotion
        and convole it with the lensing PDF

        It would be faster to do this in numpy convolve but i dont trust it
        
        Using tge fft directly is the best way
        '''
        nTotalXbins = len(self.lensingEquivalentWidth)+\
          len(self.intrinsicEquivalentWidthDistribution['x']) - 1

        startValue = self.intrinsicEquivalentWidthDistribution['x'][0] - \
          np.ceil(len(self.lensingEquivalentWidth)/2.)*self.dEquivalentWidth 

        #Create a common axis for the two distributions
        convolvedX = np.arange(0.,nTotalXbins)*self.dEquivalentWidth + \
          startValue
          
        convolvedYintrinsic = np.zeros(len(convolvedX))
        convolvedYlensing = np.zeros(len(convolvedX))

        for i, iX in enumerate(convolvedX):
            matchVectors = np.abs(self.intrinsicEquivalentWidthDistribution['x'] - iX) < 1e-10
            matchVectorsLensing = \
              np.abs(self.lensingEquivalentWidth - iX) < 1e-10

            if len(self.intrinsicEquivalentWidthDistribution['y'][matchVectors]) > 0:
                convolvedYintrinsic[i] = \
                  self.intrinsicEquivalentWidthDistribution['y'][matchVectors ]
            if len(self.lensingPDF[matchVectorsLensing]) > 0:
                convolvedYlensing[i] = self.lensingPDF[matchVectorsLensing]

        #Convolve in fourier space
        convolvedYnumpy = np.fft.ifft( np.fft.fft(convolvedYintrinsic)*np.fft.fft(convolvedYlensing))
        convolvedYnumpy /= np.sum(convolvedYnumpy)*self.dEquivalentWidth
        
        #so the fft switches the negative and positive bits around
        #so switch it back through the middle, however if the intrinsic distribution is not symmetric
        #around 0 this gets complicated
        middleIndex = np.int(np.floor(len(convolvedYnumpy)/2.))
        convolvedYnumpy = np.append( convolvedYnumpy[ middleIndex:], \
                                        convolvedYnumpy[:middleIndex])


        self.predictedLensedEquivalentWidthDistribution = \
          {'x':convolvedX,  'y': np.real(convolvedYnumpy)}
                             
        
