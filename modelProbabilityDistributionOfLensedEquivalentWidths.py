from emissionLine import emissionLine
from matplotlib import pyplot as plt
import numpy as np
import ipdb as pdb
def main():
    '''
    Run some tests on this class
    '''

    newModel = modelProbabilityDistributionOfLensedEquivalentWidths()
    ax = plt.gca()
    for iAlpha in np.linspace(0.1,0.9,9):
        newDistribution = newModel.getNewDistribution(lensingParams={'alpha':iAlpha})
        plt.plot(newDistribution['x'],newDistribution['y'],'-')
    
    ax.set_xlim(-2,3.5)
    ax.set_ylim(0.0001,2.)
    ax.set_yscale('log')
    plt.show()
    
class modelProbabilityDistributionOfLensedEquivalentWidths:
    def __init__( self ):
        '''
        This class is set up to learn the distriubtion
        of magnitudes for a given distribution and return
        the predicted probabilityh distribution
        for given parameters
            '''

        self.setDefaultDistriubtion()

    def setDefaultDistriubtion( self):
        self.emissionLineCl = emissionLine('NARROW_HB', nRedshiftBins=2 )
        
        self.emissionLineCl.setLensingProbability( )
        
        self.emissionLineCl.setIntrinsicDistribution( )

        self.emissionLineCl.convolveIntrinsicEquivalentWidthWithLensingProbability()

    def getNewDistribution( self, lensingParams={}, intrinsicDistParams={}):

        self.emissionLineCl.setLensingProbability( **lensingParams )
        
        self.emissionLineCl.setIntrinsicDistribution(**intrinsicDistParams)
        
        self.emissionLineCl.convolveIntrinsicEquivalentWidthWithLensingProbability()

        return self.emissionLineCl.predictedLensedEquivalentWidthDistribution

if __name__ == '__main__':
    main()
