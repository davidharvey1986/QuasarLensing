'''

get a bunch of convoled log normal distributions and
try and fit a analytical prescription so i can use it in an MCMC

'''
from analyseWeq import *

from scipy import interpolate


class fittedPBHdistribution:

    def __init__( self, logNormalParams=[0.8,0.3]):

        self.logNormParams = logNormalParams
        self.alphaList = [0.1, 0.3, 0.5, 0.83]
        
        
    
    def loadTrueDistributions( self ):

        emissionLineCl = emissionLine('NARROW_HB', nRedshiftBins=1)
        emissionLineCl.getLensingProbability( z=1.0, alpha=0.83 )

        x = np.arange(-1.0,4.0,  emissionLineCl.dEquivalentWidth)[1:]
         
        y = norm.pdf(x, self.logNormParams[0], self.logNormParams[1])

         
        emissionLineCl.intrinsicEquivalentWidthDistribution = \
          {'x':x, 'y':y}
        axis = plt.gca()

        #Mapp all the distributiosn to the same x
        self.trueDistributionsEquivWidths = x
        
        trueDistributions = { }
        
        for iAlpha, iAlphaValue in enumerate(self.alphaList):

            emissionLineCl.getLensingProbability( z=1.0, alpha=iAlphaValue )

            emissionLineCl.convolveIntrinsicEquivalentWidthWithLensingProbability()

            mappedTrueDistribution =  \
              emissionLineCl.interpolatePDF( x, emissionLineCl.predictedLensedEquivalentWidthDistribution)
              
            trueDistributions[ iAlphaValue ] = mappedTrueDistribution
            
        self.trueDistributions = trueDistributions



    def getInterpolateFunction( self ):

        interpolationData = np.zeros((len(self.trueDistributionsEquivWidths), len(self.trueDistributions.keys())))
        print self.trueDistributions.keys()
        for iCount, iKey in enumerate(self.trueDistributions.keys()):
            interpolationData[:, iCount] = self.trueDistributions[iKey]


        self.interpolationFunction = \
          interpolate.interp2d(  np.array(self.trueDistributions.keys()), self.trueDistributionsEquivWidths,\
                                    interpolationData, kind='cubic')

    def predictDistribution( self, equivalentWidths, alpha ):

        return self.interpolationFunction( alpha, equivalentWidths)





        


                                     
        
        

if __name__ == '__main__':
    main()
