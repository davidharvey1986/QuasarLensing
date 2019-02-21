import numpy as np
import pyfits as fits
import convolveLssPc as clp


class emissionLine:
    '''
    log the all the data from the shen catalopgues
    in to this class which will just make things look a bit cleaner
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
        completeDataListNames =\
          dataArray.columns.names

          
        self.checkEmissionLine(completeDataListNames)
        #now set in the class
        #init the redshifts of all the quasars
        self.redshift = dataArray['Z_PIPE']
        #init the rest frame equvilent width
        
        self.restFrameEquivilentWidth = \
          dataArray['REW_'+self.emissionLine]
        
        self.cleanEquivilentWidths()
        self.setRedshiftBins( )

    def checkEmissionLine( self, dataNames ):
        '''
        Check that the input name exists in the dictionary
        of column names
        Otherwise reaise a flag
        '''

        exactFlag = [ 'REW_'+self.emissionLine == i for i in  dataNames]
        flag = [ 'REW_'+self.emissionLine in i for i in  dataNames]
        
        if (not np.any(np.array(flag))) & \
          (not np.any(np.array(exactFlag))):
            raise ValueError("Cant find emission line name")
        elif  (np.any(np.array(flag))) & \
          (not np.any(np.array(exactFlag))):
          proposedName = np.array(dataNames)[flag][0]
          raise ValueError("Can't find exact emission line name, did you mean %s?" % proposedName)
                             
      
          
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

        
            
        
    def getIntrinsicDistribution( self, redshiftCut=0.3 ):
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
        
        '''

        self.dEquivalentWidth = clp.get_dMuPrime()
        self.intrinsicEquivalentWidths = \
          self.restFrameEquivilentWidth[ self.redshift < redshiftCut ]
        #this might not work as the data is not well sampled enough
        #i might need to bin it coarse and simply sub-sample it into smaller bins
        #TO DO: I NEED TO DO THIS!
        equivalentWidthBins = np.arange(0., np.max(self.intrinsicEquivalentWidths),\
                                       self.dEquivalentWidth)     


        #need to confirm what i am doing here is correct.
        #does it all need to be the same dEW
        
        equivalentWidthBins = \
          np.append(equivalentWidthBins, equivalentWidthBins[-1]+self.dEquivalentWidth)
        


        y, x = np.histogram( self.intrinsicEquivalentWidths, \
                                 bins=equivalentWidthBins, density=True)
                                 
        xc = (x[:-1]+x[1:])/2.
        self.intrinsicEquivalentWidthDistribution = \
          {'y':y, 'x':xc}


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
        
        magnitude, lensingPDF  = \
          clp.totalPl(z=z, alpha=alpha,nMu=1000)

        self.lensingMagnitude = magnitude
        self.lensingPDF = lensingPDF

    def convolveIntrinsicEquivalentWidthWithLensingProbability(self):
        '''
        In order to get an estimate of the expected distribution of 
        equivalent widths i need to take the intrinsic distrubotion
        and convole it with the lensing PDF

        NOTE: Currently the lensing pronbability is in terms of magnitudes
        this will change once we have it in terms of Equivalent widhts
        
        '''

        convolvePDF =\
          np.convolve( self.intrinsicEquivalentWidthDistribution['y'], \
                                self.lensingPDF)

        middleIndex = np.int(np.ceil(len(self.intrinsicEquivalentWidths)/2.))
        magnitudeMiddleValue = self.lensingMagnitude[middleIndex]

                  
        
        self.predictedLensedEquivalentWidthDistribution = \
          {'x':np.arange(0.,self.dEquivalentWidth*len(convolvePDF),\
                             self.dEquivalentWidth) + magnitudeMiddleValue, \
            'y': convolvePDF}
                             
        
