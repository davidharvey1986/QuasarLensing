import numpy as np
import pyfits as fits
import convolveLssPc as clp
from scipy.stats import norm
import lensingProbabilityDistribution as lpd
import ipdb as pdb
from matplotlib import pyplot as plt
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

        
            
        
    def getIntrinsicDistribution( self, redshiftCut=0.3, mimicSuperNova=True):
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

        
        self.intrinsicEquivalentWidths = \
          self.restFrameEquivilentWidth[ self.redshift < redshiftCut ]
          
    
        
        if mimicSuperNova:
            #This is purely for test purposes
            x = np.arange(-1.0,1.0,  self.dEquivalentWidth)
            y = norm.pdf(x, 0., 0.15)
            self.intrinsicEquivalentWidthDistribution = \
              {'y':y, 'x':x}
            return
        else:
            equivalentWidthBins = np.arange(0.,np.max(self.intrinsicEquivalentWidths), \
                                                 self.dEquivalentWidth)
            
        #this might not work as the data is not well sampled enough
        #i might need to bin it coarse and simply sub-sample it into smaller bins
        
        self.nEquivalentWidthBins = len(equivalentWidthBins)
        print("Intrinsic Equivalent Width sampled at %0.3f\n" %  self.dEquivalentWidth)
        
        #need to confirm what i am doing here is correct.
        #does it all need to be the same dEW
        
        equivalentWidthBins = \
          np.append(equivalentWidthBins, \
                        equivalentWidthBins[-1]+self.dEquivalentWidth)
        
        
        y, x = np.histogram( self.intrinsicEquivalentWidths, \
                                 bins=equivalentWidthBins, \
                                 density=True)
        xCenters = (x[1:] + x[:-1])/2.
        self.intrinsicEquivalentWidthDistribution = \
          {'y':y, 'x':xCenters}


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
                                              nEquivalentWidthBins=1000)
        #force some kind of normalisation, although this should alrady be normalised
        #lensingPDF /= np.sum(lensingPDF)*(magnitude[1]-magnitude[0])
        self.dEquivalentWidth = lensingPDF.dEquivalentWidth
        self.lensingEquivalentWidth = lensingPDF.convolvedPbhPdfWithLssPdf['x'] + \
          lensingPDF.dEquivalentWidth/2.
        self.lensingPDF = lensingPDF.convolvedPbhPdfWithLssPdf['y']

        
        self.lensingEquivalentWidth[ self.lensingPDF > 0]
        self.lensingPDF[ self.lensingPDF > 0 ]
       
    def convolveIntrinsicEquivalentWidthWithLensingProbability(self):
        '''
        In order to get an estimate of the expected distribution of 
        equivalent widths i need to take the intrinsic distrubotion
        and convole it with the lensing PDF

        NOTE: Currently the lensing pronbability is in terms of magnitudes
        this will change once we have it in terms of Equivalent widhts
        
        '''
        nTotalXbins = len(self.lensingEquivalentWidth)+\
          len(self.intrinsicEquivalentWidthDistribution['x']) - 1

        startValue = self.intrinsicEquivalentWidthDistribution['x'][0] - \
          len(self.lensingEquivalentWidth)*self.dEquivalentWidth
          
        convolvedX = np.arange(0.,nTotalXbins)*self.dEquivalentWidth + startValue
        convolvedYlss = np.zeros(len(convolvedX))
       
        convolvedYtotal = np.zeros(len(convolvedX))

        for i, iX in enumerate(convolvedX):
            matchVectors = np.abs(self.intrinsicEquivalentWidthDistribution['x'] - iX) < 1e-10
            if len(self.intrinsicEquivalentWidthDistribution['y'][matchVectors]) == 0:
                convolvedYlss[i] = 0
            else:

                convolvedYlss[i] = \
                  self.intrinsicEquivalentWidthDistribution['y'][matchVectors ]

        for i, xPrime in enumerate(convolvedX):
            convolvedYlensing = np.zeros(len(convolvedX))
            for j, jX in enumerate(convolvedX):
                matchVectors = np.abs(self.lensingEquivalentWidth + xPrime - iX) < 1e-10
            
                if len(self.lensingPDF[matchVectors]) == 0:
                    convolvedYlensing[j] = 0
                else:
                    convolvedYlensing[j] = \
                    self.lensingPDF[ matchVectors ]
            
            print i, np.sum( convolvedYlensing*convolvedYlss*self.dEquivalentWidth)
            
            convolvedYtotal[i] = np.sum( convolvedYlensing*convolvedYlss*self.dEquivalentWidth)
            if     convolvedYtotal[i] > 0:
                plt.plot(convolvedX, convolvedYlensing )
                plt.plot(convolvedX, convolvedYlensing )
                


                
        plt.plot(convolvedX,convolvedYlss,'-')

        plt.show()
        pdb.set_trace()


    

        #Need to get the two pdfs on a common frame
        


        self.predictedLensedEquivalentWidthDistribution = \
          {'x':newX,  'y': convolvePDF}
                             
        
