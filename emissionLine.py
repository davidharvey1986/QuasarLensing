import numpy as np
from astropy.io import fits
import convolveLssPc as clp
from scipy.stats import norm
from scipy.stats import rayleigh
import cosmolopy.distance as dist
import lensingProbabilityDistribution as lpd
import ipdb as pdb
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import fitPBHdistribution as fitPBH

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
        #But given i am doing everything in log(1.+mu) do this
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
        #Bolometric luminsoity of the continuum
        A = dataArray['CONTI_FIT'][:,0]
        B = dataArray['CONTI_FIT'][:,1]
        
        bolometricLuminosity =  A/(B+1.)*(1. - (1./3.)**(B+1.))
        self.continuumLuminosity =np.log10( bolometricLuminosity * 4.*np.pi*luminosityDistance**2)


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
                             

    def applyLuminosityCut( self, luminosityCut=10. ):
        '''
        Apply a luminosity cut at a particularly emission line
        '''
        index =  (self.continuumLuminosity > luminosityCut)

        self.restFrameEquivilentWidth = \
            self.restFrameEquivilentWidth[index[self.cleanEquivalentWidths ]]
            
        self.redshift = self.redshift[index[self.cleanEquivalentWidths ]]
        self.continuumLuminosity = self.continuumLuminosity[index & self.cleanEquivalentWidths ]
        print("After Luminosity Cut there are %i quasars" % len(self.redshift))

        
    def cleanEquivilentWidths( self ):
        '''
        I need to clean out all the rubbish equivilent widths
        as some are a load of rubbish
        '''
        self.cleanEquivalentWidths = (self.restFrameEquivilentWidth>0) & \
          (self.restFrameEquivilentWidth<1e3) & \
          (self.restFrameEquivilentWidth>-1) & \
          (self.restFrameEquivilentWidth<1e3) & \
          (self.redshift > 0)
          
          
        #init the redshifts of all the quasars
        self.redshift = self.redshift[self.cleanEquivalentWidths ]
        #init the rest frame equvilent width
        
        self.restFrameEquivilentWidth = \
          self.restFrameEquivilentWidth[self.cleanEquivalentWidths ]

        
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



    def histogramEquivalentWidth(self, nEquivilentWidthBins=None):
        '''
        Hisogram each redshift bin and turn them into
        a dictionary set within a dictionary set 
        i.e. [0.1].x

        '''
        if nEquivilentWidthBins is None:
            nEquivilentWidthBins = 10**np.linspace(-3, 3, 30)

        self.equivilentWidthHistograms = {}

        for iBin in xrange(self.nRedshiftBins):
            selectRestFrameEquivilentWidth = \
              self.restFrameEquivilentWidth[self.binIndices==iBin]
            
            
            yHistogram, xHistogram = \
              np.histogram(selectRestFrameEquivilentWidth, \
                               bins=nEquivilentWidthBins, \
                               density=False)
            yHistogram = yHistogram.astype(float)/np.sum(yHistogram)
            redshiftBinName = "%0.2f" % self.redshiftCentres[iBin]
            histogramCentres = (xHistogram[1:] + xHistogram[:-1])/2.
            iRedshiftDict = {'x':histogramCentres, \
                                 'y':yHistogram, \
                                'samples':selectRestFrameEquivilentWidth}
            self.equivilentWidthHistograms[redshiftBinName] = iRedshiftDict

    def logGaussianFitEquivalentWidths( self ):

        '''
        Fit a gaussian in log(equivalent widths )
        '''
        logNormalParameters = {}
        
        for iRedshiftBin in self.equivilentWidthHistograms.keys():

            params = gaussFit( self.equivilentWidthHistograms[iRedshiftBin] )
            logNormalParameters[ iRedshiftBin ] = params
            print("Fitted lognormal gaussian parameters  "+\
                  "amplitude:%0.2f, mean:%0.2f, scale:%0.24f" %\
                  tuple(params))
        self.logNormalParameters =logNormalParameters

    
    def plotLogGaussianFitEquivalentWidths( self, iRedshiftBin, **plottingParams ):
        logNormalParameters = {}
        
        for iRedshiftBin in self.equivilentWidthHistograms.keys():
            params = norm.fit( np.log10( self.equivilentWidthHistograms[iRedshiftBin]['samples'] ))
            x = 10**np.linspace(np.log10(self.equivilentWidthHistograms[iRedshiftBin]['x'][0]),\
                                np.log10(self.equivilentWidthHistograms[iRedshiftBin]['x'][-1]),\
                                100)
            y = gauss( np.log10(x), *self.logNormalParameters[ iRedshiftBin ])

            ax = plt.gca()
            ax.plot(x,y, **plottingParams)
            
        
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

        
            
        
    def setIntrinsicDistribution( self, **intrinsicDistParams):
        '''
        In the lensing paper they assume that the PDF of the
        lensed supernova is a Gaussian with a sigma of 0.15
        
        In quasars we do not have the same luxury. 
        So to get the 'intrinsic distribution' i will use the
        low redshift quasars (i.e. the first redshift bin)

        

        Note that to do this the dEquivalentWidth (increment in EW)
        must be the same as the lensing probability distribution
        that we will use to convovlve it with

        Therefore the bins should be preset

        mimicSuperNova : Instead of using the raw data of quasar distribtion
           mimic the distributions from the supenova papers, and 
           therefor a Gauassian with a variance=0.15
        
        '''
        #These have been derived from data
        self.intrinsicDistParams = \
          {'type':'LogNormal', 'mean':0.84, 'scale':0.3}
          
        for iParam in intrinsicDistParams.keys():
            self.intrinsicDistParams[iParam] = \
              intrinsicDistParams[iParam]
        

        #This is purely for test purposes
        if self.intrinsicDistParams['type'] == 'LogNormal':
            x = np.arange(-3.0,3.0,  self.dEquivalentWidth)
            y = norm.pdf(x,  self.intrinsicDistParams['mean'],  \
                            self.intrinsicDistParams['scale'])
            #these need to be symmetric otherwise, the position of the final
            #convoution is weird since the fft switches things around,
            #I could probably figure it out but for now this will do

            #Also if things go NaN then try making this symmetric in even or
            #odd numbers i.e. (-1,1) or (-2,2)
        if self.intrinsicDistParams['type'] == 'Rayleigh':
            x = np.arange(-5.0, 5.0,  self.dEquivalentWidth)
            y = rayleigh.pdf(x, 0., self.intrinsicDistParams['scale'])

        if self.intrinsicDistParams['type'] == 'Semi-Gaussian':
            x = np.arange(-31.0,31.0,  self.dEquivalentWidth)
            y = norm.pdf(x, 0., self.intrinsicDistParams['scale'])
            y[x < 0] = 0.
                
 

            dX = x[1]-x[0]
            y /= np.sum(y*dX)
            
        if self.intrinsicDistParams['type'] == 'data':
            #Get the data and fit to it
            self.applyLuminosityCut( )
            self.getEquvilentWidthMeansAsFunctionOfRedshift()
            self.histogramEquivalentWidth()
            self.logGaussianFitEquivalentWidths()
            
            x = np.arange(-3, 3, self.dEquivalentWidth)
            redshift = np.sort(np.array(self.logNormalParameters.keys()))[0]
            print("Intrincsic distrinbution of quasars are "+\
                      "in the redshift range <%s" %redshift)
            y = gauss( x, * self.logNormalParameters[redshift])
            
            dX = x[1]-x[0]
            y /= np.sum(y*dX)

        self.intrinsicEquivalentWidthDistribution = \
          {'y':y, 'x':x}


    def setLensingProbability( self, **inputParams ):
        '''
        I need to the probability distrubtion that a quasar
        of redshift z is magnified or demagnified by large scale
        structure and primordial black holes. Where the fraction
        of dark matter in primordial black holes is alpha

        z: Quasar Redshift
        alpha: fraction of dark matter is primordial black holes

        NOTES: for now i will use a given redshift of 1. and alpha=0.27
        this will need to be changed such that i can iterate and fit

        NOTES: This returns a log(p(1+mu))

        '''
        #TO DO This needs to be checked this funcion.
        lensingPDF = \
          lpd.lensingProbabilityDistribution( **inputParams )
       

        self.dEquivalentWidth = lensingPDF.dEquivalentWidth

        #this probability distribution is the effect of kappe_ew
        #the amplictude will be 1+2kappa
        self.lensingEquivalentWidth = \
          lensingPDF.convolvedPbhPdfWithLssPdf['x'] 
        #times the dEW / dKappa
        self.lensingPDF = lensingPDF.convolvedPbhPdfWithLssPdf['y']
        self.alpha = lensingPDF.inputParams['alpha']

        self.lensingEquivalentWidth = \
          self.lensingEquivalentWidth[ self.lensingPDF > 0] 
          
        self.lensingPDF = self.lensingPDF[ self.lensingPDF > 0 ]

        self.lensingProbability = \
          { 'x':self.lensingEquivalentWidth, 'y':self.lensingPDF }


    def interpolatePDF( self, x, pdf ):

        returnInterp =  np.interp( x, pdf['x'], pdf['y'] )
        returnInterp[ (x < pdf['x'][0]) |  (x > pdf['x'][-1]) ] = 0
        return returnInterp
    
    def convolveIntrinsicEquivalentWidthWithLensingProbability(self):
        '''
        In order to get an estimate of the expected distribution of 
        equivalent widths i need to take the intrinsic distrubotion
        and convole it with the lensing PDF

        It would be faster to do this in numpy convolve but i dont trust it
        
        Using tge fft directly is the best way
        '''
        nTotalXbins = (len(self.lensingEquivalentWidth)+\
          len(self.intrinsicEquivalentWidthDistribution['x']) - 1)

        startValue = self.intrinsicEquivalentWidthDistribution['x'][0] - \
          np.ceil(len(self.lensingEquivalentWidth)/2.)*self.dEquivalentWidth

        #Create a common axis for the two distributions
        observedEquivalentWidth = \
          np.arange(0.,nTotalXbins)*self.dEquivalentWidth + \
          startValue
        dEW = observedEquivalentWidth[1] - observedEquivalentWidth[0]
        totalConvolved = np.zeros(nTotalXbins)
        
        for iBinInConv, iObservedEquivWidth in enumerate(observedEquivalentWidth):

            probabilityIntrinsicEquivalentWidth = \
              self.interpolatePDF( observedEquivalentWidth, \
                        self.intrinsicEquivalentWidthDistribution)

            probabilityLensingKernel = \
              self.interpolatePDF( iObservedEquivWidth - \
                                    observedEquivalentWidth, \
                                      self.lensingProbability)

            
            totalConvolved[iBinInConv] = \
              np.sum( probabilityIntrinsicEquivalentWidth * \
                          probabilityLensingKernel * dEW)

        
        self.predictedLensedEquivalentWidthDistribution = \
          {'x':observedEquivalentWidth,  'y': totalConvolved}


    def getInterplationFromIntrinsic( self ):
        
        fittedDist  = fittedPBHdistribution( logNormalParams=self.logNormalParameters[0])
        fittedDist.loadTrueDistributions()

        fittedDist.interpolateDistribution()
        self.fittedDist = fittedDist

    def predictProbability( self, equivalentWidth, alpha):
        return self.fittedDist.predictDistribution(equivalentWidth, alpha)
        
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def gaussFit(  PDF ):
    '''
    gaussfit the pdf where PDF is a dictionary conatinaing ['x'], ['y']
    '''
    p0 = [1., 0., 1.]

    coeff, var_matrix = \
      curve_fit(gauss, np.log10(PDF['x']),PDF['y'], p0=p0)

    return coeff
