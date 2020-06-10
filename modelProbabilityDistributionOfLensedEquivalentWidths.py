from emissionLine import emissionLine
from matplotlib import pyplot as plt
import numpy as np
import ipdb as pdb
from itertools import product
from scipy.interpolate import interpn
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import LinearNDInterpolator
import pickle as pkl
from matplotlib import gridspec
import progressbar
from sklearn.linear_model import LinearRegression


import os

def main():
    '''
    Run some tests on this class, and show how the predictive works
    '''
    print("Initiatilising")
    newModel = modelProbabilityDistributionOfLensedEquivalentWidths()
    print("Setting Grid")
    newModel.setInterpolatorFunctionParamGrid()
    if os.path.isfile(newModel.paramGridFile):
        newModel.loadParamGrid()
    else:
        print("Filling Grid")
        newModel.fillParamGrid()
        print("Saving param grid")
        newModel.saveParamGrid()
    print("Fitting interpolator")
    newModel.fitInterpolator()
    
    predictAlphas = np.linspace(0.2,0.4, 8)
    
    fig=plt.figure()
    gs = gridspec.GridSpec(5,1)
    ax1 = plt.subplot(gs[:3,0])
    ax2 = plt.subplot(gs[3:,0])
    fig.subplots_adjust(hspace=0)
    for iAlpha in predictAlphas:
        predictThese = {'alpha':iAlpha, 'scale':0.3, 'redshift':1.0}
        
        trueDist = newModel.getNewDistribution(predictThese)

        predictedPDF = \
          newModel.predictPDFforParameterCombination(predictThese, \
                                    xVector=trueDist['x'] )

        ax1.plot(trueDist['x'], predictedPDF['y'])
        
        ax2.plot(trueDist['x'], \
            np.cumsum(predictedPDF['y'])/np.sum(predictedPDF['y']) - \
            np.cumsum(trueDist['y'])/np.sum(trueDist['y']))
        
    #ax1.set_yscale('log')
    #ax2.set_yscale('log')
    ax2.set_ylim(-2e-3,2e-3)
    ax1.set_ylim(1e-4,2.2)
    ax1.set_xlim(-1,2.2)
    ax2.set_xlim(-1,2.2)

    plt.show()
    
class modelProbabilityDistributionOfLensedEquivalentWidths:
    def __init__( self, paramGridFile='pickles/paramGrid.pkl' ):
        '''
        This class is set up to learn the distriubtion
        of magnitudes for a given distribution and return
        the predicted probabilityh distribution
        for given parameters
            '''
        self.paramGridFile = paramGridFile
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
            'alpha': np.linspace(0.1,0.83, 10),
            'scale' : np.linspace(0.2,0.5, 10),\
            'redshift':np.linspace(0.5,5, 10) }
                
        self.paramKeys = self.interpolateParams.keys()

    def fitInterpolator( self, loadPklFile=False ):

        pklFile = 'pickles/interpolatorWithRedshift.func'
        
        if loadPklFile:
            self.interpolatorFunction = \
              pkl.load(open(pklFile,'rb'), encoding='latin1') 
        else:

            self.interpolatorFunction = \
              LinearNDInterpolator( self.points, self.values)
              
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
        
    def saveParamGrid( self):
        paramGrid = [self.points, self.values]
        pkl.dump(paramGrid, open(self.paramGridFile,'wb'))

    def loadParamGrid( self):
        self.points, self.values = pkl.load(open(self.paramGridFile,'rb'), encoding='latin1') 
            
    def generateDistributionOfEquivalentWidths( self, inputParameters, nObservations):
        '''
        Generate a selection of equivalent widths, as predicted by this
        model predictor for the given inputParams and the number of nObservations

        inputParams : a dict of the input params
        nObservations: an int of the nuymber of observations
        '''

        #gfenerate a vector of many possible EW for which i will get their
        #probabilities (many so non get double selected) reality = inifite number of chocies
        xVector = np.linspace(-3, 3, nObservations*100)
        theoreticalPDF = \
          self.predictPDFforParameterCombination(inputParameters, xVector=xVector)
          
        probabilities = theoreticalPDF['y'] / np.sum(theoreticalPDF['y'])
        selectedEquivalentWidths = \
          np.random.choice( xVector,  size=nObservations,p=probabilities )
          
        probabiltiyDistribution, x = \
          np.histogram( selectedEquivalentWidths, \
                    bins=np.linspace(-3,3,nObservations/100), density=True)
                    
        binCenters = (x[1:] + x[:-1])/2.
                
        pdf = {'x': binCenters, 'y':probabiltiyDistribution}
            
        return pdf

            
            
            
if __name__ == '__main__':
    main()
