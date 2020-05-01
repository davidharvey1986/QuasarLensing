from emissionLine import emissionLine
from matplotlib import pyplot as plt
import numpy as np
import ipdb as pdb
from itertools import product
from scipy.interpolate import interpn
from scipy.interpolate import NearestNDInterpolator
import pickle as pkl
from matplotlib import gridspec
import progressbar
def main():
    '''
    Run some tests on this class, and show how the predictive works
    '''

    newModel = modelProbabilityDistributionOfLensedEquivalentWidths()
    newModel.setInterpolatorFunctionParamGrid()
    newModel.fillParamGrid()
    newModel.fitInterpolator(loadPklFile=False)
    
    predictAlphas = np.linspace(0.2,0.4, 8)
    fig=plt.figure()
    gs = gridspec.GridSpec(5,1)
    ax1 = plt.subplot(gs[:3,0])
    ax2 = plt.subplot(gs[3:,0])
    fig.subplots_adjust(hspace=0)
    for iAlpha in predictAlphas:
        predictThese = {'alpha':0.2, 'scale':iAlpha}
        
        trueDist = newModel.getNewDistribution(predictThese)

        predictedPDF = \
          newModel.predictPDFforParameterCombination(predictThese, \
                                    xVector=trueDist['x'] )

        ax1.plot(trueDist['x'], predictedPDF['y'])
        
        ax2.plot(trueDist['x'], predictedPDF['y']/(1+trueDist['y']))
        
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim(1e-4,1)
    ax1.set_ylim(1e-4,3)
        
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

    def getNewDistribution( self, inputParams={}):

        self.emissionLineCl.setLensingProbability( **inputParams )
        
        self.emissionLineCl.setIntrinsicDistribution(**inputParams)
        
        self.emissionLineCl.convolveIntrinsicEquivalentWidthWithLensingProbability()

        return self.emissionLineCl.predictedLensedEquivalentWidthDistribution


    def setInterpolatorFunctionParamGrid( self ):

        self.interpolateParams = { \
            'alpha': np.linspace(0.1,0.83, 100),
            'scale' : np.linspace(0.2,0.5, 100)}
                
        self.paramKeys = self.interpolateParams.keys()

    def fitInterpolator( self, loadPklFile=False ):

        pklFile = 'pickles/interpolator.func'
        
        if loadPklFile:
            self.interpolatorFunction = \
              pkl.load(open(pklFile,'rb'))
        else:

            self.interpolatorFunction = \
              NearestNDInterpolator( self.points, self.values)
            pkl.dump(self.interpolatorFunction, \
                        open(pklFile, 'wb'))

    def predictPDFforParameterCombination( self, inputParams, xVector=None):
        '''
        inputParams needs to be a dict so i can ensure i put them in 
        the correct order
        '''
        reOrderedParams = []
        
        for iKey in self.paramKeys:
            reOrderedParams.append(inputParams[iKey])

        if xVector is None:
            xVector = \
              self.emissionLineCl.predictedLensedEquivalentWidthDistribution['x']
            
        predict = [ np.append(i,np.array(reOrderedParams)) \
                        for i in xVector ]


        predictedPDF = \
          self.interpolatorFunction( predict )

        pdf = {'x':xVector, 'y': predictedPDF}
        
        return pdf
    

        
        
        
    def fillParamGrid( self ):
        '''
        Using the scipy nearest neihbour which requires 
        points : the value of the parameter combineation (equivWidth, alpha, scale) for ecample
        values: the y value at eac of these parameter combinations
        '''
      

        nParamCombos = np.prod([ len(self.interpolateParams[i]) \
                    for i in self.paramKeys])
        
        interpolateParamCube = np.array([])

        listsOfParams = \
          [ list(self.interpolateParams[i]) for \
                i in self.paramKeys]

        self.points = []
        self.values = np.array([])

        bar = progressbar.ProgressBar(maxval=nParamCombos, \
            widgets=[progressbar.Bar('=', '[', ']'), ' ', \
                         progressbar.Percentage()])
        bar.start()

        for progress, iParamCombination in enumerate(product(*listsOfParams)):

            bar.update(progress+1)
            
            
            inputParams = {}
            for iParam, iKey in enumerate(self.paramKeys):
                inputParams[iKey] = iParamCombination[iParam]

            iPDF = self.getNewDistribution(inputParams)
            self.points = self.points + \
              [ np.append(i,np.array(iParamCombination)) for i in iPDF['x'] ]
            self.values = np.append(self.values, iPDF['y'])

        bar.finish()
if __name__ == '__main__':
    main()
