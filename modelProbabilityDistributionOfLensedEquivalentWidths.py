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
from sklearn.gaussian_process.kernels import  Matern
from sklearn.decomposition import PCA
import os
from sklearn.gaussian_process import GaussianProcessRegressor

def main():
    '''
    Run some tests on this class, and show how the predictive works
    '''
    newModel = modelProbabilityDistributionOfLensedEquivalentWidths()
    newModel.setInterpolatorFunctionParamGrid()
    if os.path.isfile(newModel.paramGridFile):
        newModel.loadParamGrid()
    else:
        newModel.fillParamGrid()
        newModel.saveParamGrid()

    newModel.fitInterpolator(loadPklFile=True)
    
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
    
    def __init__( self, paramGridFile='pickles/paramGrid.pkl', \
                    nPrincipalComponents=10 ,\
                    regressorNoiseLevel=1e-3):
        '''
        This class is set up to learn the distriubtion
        of magnitudes for a given distribution and return
        the predicted probabilityh distribution
        for given parameters

        To do this it
        1. Generates a parameter grid and determines the PDF for PBH and
           each point in that parameter grid
        2. Then compresses these PDFs using a PCA
        3. The learns the PCAs with a Gaussian Processor
        
        '''
        print("Initiatilising Class")

        self.paramGridFile = paramGridFile
        self.setDefaultDistribution()
        self.nPrincipalComponents = nPrincipalComponents
        self.regressorNoiseLevel = regressorNoiseLevel
    def setDefaultDistribution( self):
        self.emissionLineCl = emissionLine('NARROW_HB', nRedshiftBins=2 )
        
        self.emissionLineCl.setLensingProbability( )
        
        self.emissionLineCl.setIntrinsicDistribution( )

        self.emissionLineCl.convolveIntrinsicEquivalentWidthWithLensingProbability()
        self.interpolateToTheseEqWidth = \
          self.emissionLineCl.predictedLensedEquivalentWidthDistribution['x']

    def getNewDistribution( self, inputParams={}):

        self.emissionLineCl.setLensingProbability( **inputParams )
        
        self.emissionLineCl.setIntrinsicDistribution(**inputParams)
        
        self.emissionLineCl.convolveIntrinsicEquivalentWidthWithLensingProbability()

        return self.emissionLineCl.predictedLensedEquivalentWidthDistribution


    def setInterpolatorFunctionParamGrid( self ):
        print("Setting Parameter Grid")

        self.interpolateParams = { \
            'alpha': np.linspace(0.1,0.83, 10),
            'scale' : np.linspace(0.2,0.5, 10),\
            'redshift':np.linspace(0.5,5, 10) }
                
        self.paramKeys = self.interpolateParams.keys()

    def fitInterpolator( self, loadPklFile=False ):

     
        pklFile = 'pickles/interpolatorWithRedshift.func'
        
        self.compressPDFwithPCA()

        if loadPklFile:
            print("Loading Gaussian Processor")
            self.predictor = \
              pkl.load(open(pklFile,'rb'), encoding='latin1') 
        else:
            self.learnPrincipalComponents()             
            pkl.dump(self.predictor, open(pklFile, 'wb'))

    def compressPDFwithPCA( self ):
        '''
        Compress all the PDFs use PCA analysis

        '''
        print("Compressing PDFs with a PCA")
        self.pca = PCA(n_components=self.nPrincipalComponents)
        self.pca.fit( self.targetPDF )
        self.principalComponents = \
          self.pca.transform( self.targetPDF )
          
    def learnPrincipalComponents( self, length_scale=1., nu=3./2.):
        '''
        Using a mixture of Gaussian Processes 
        predict the distributions of compoentns
        It will return a list of nPrincipalComponent models that try to
        learn each component to the data
        '''
        print("Learning PCAs with Gaussian Processor")
        self.predictor = []

        kernel =  Matern(length_scale=length_scale, nu=nu)
       
        for i in range(self.nPrincipalComponents):
            print("Fitting Component %i/%i" % (i, self.nPrincipalComponents))
            gaussProcess = \
              GaussianProcessRegressor( alpha=self.regressorNoiseLevel, kernel=kernel,\
                                            n_restarts_optimizer=10)

            gaussProcess.fit(self.features, self.principalComponents[:,i])

                
            self.predictor.append(gaussProcess)
        
    def predictPDFforParameterCombination( self, inputParams, xVector=None):
        '''
        inputParams needs to be a dict so i can ensure i put them in 
        the correct order
        '''
        reOrderedParams = []

        for iKey in self.paramKeys:
            reOrderedParams.append(inputParams[iKey])
        reOrderedParams = np.array(reOrderedParams)
        
        if xVector is None:
            xVector = \
              self.emissionLineCl.predictedLensedEquivalentWidthDistribution['x']
              

        predictedComponents = \
          np.zeros(self.nPrincipalComponents)

        for iComponent in range(self.nPrincipalComponents):
            #there are now many preditors for the subsamples
            #So the predictor for this subsample is
            predictor = self.predictor[iComponent]


            predictedComponents[iComponent] = \
              predictor.predict(reOrderedParams.reshape(1,-1))

            
        predictedTransformCDF =  \
          self.pca.inverse_transform( predictedComponents)

        if xVector is not None:
            pdfDict = \
              {'x': self.interpolateToTheseEqWidth, \
               'y':predictedTransformCDF}
            finalPDF =  \
              self.interpolatePDF(pdfDict, xVector=xVector)
        else:
            xVector = self.interpolateToTheseEqWidth
            finalPDF = predictedTransformCDF


        finalPDF[ finalPDF < 0 ] = 0
        return {'x': xVector, 'y':finalPDF}
        
        
        
    def fillParamGrid( self ):
        '''
        Using the scipy nearest neihbour which requires 
        points : the value of the parameter combineation (equivWidth, alpha, scale) for ecample
        values: the y value at eac of these parameter combinations
        '''
        print("Filling Parameter Grid")


        nParamCombos = np.prod([ len(self.interpolateParams[i]) \
                    for i in self.paramKeys])
        
        interpolateParamCube = np.array([])

        listsOfParams = \
          [ list(self.interpolateParams[i]) for \
                i in self.paramKeys]

        self.targetPDF = None 
        self.features = None 
          
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

            interpolatedPDF = self.interpolatePDF( iPDF )

            if self.targetPDF is None:
                self.features = iParamCombination
                self.targetPDF = interpolatedPDF
            else:
                self.features = np.vstack((self.features, iParamCombination))
                self.targetPDF = np.vstack((self.targetPDF, interpolatedPDF))

                
        bar.finish()
    def interpolatePDF(self, pdfDict, xVector=None):
        '''
        interpolte the input pdf to the x values
        '''
        if xVector is None:
            xVector = self.interpolateToTheseEqWidth
        return np.interp( xVector,pdfDict['x'], pdfDict['y'])
        
        
    def saveParamGrid( self):
        paramGrid = [self.features, self.targetPDF]
        pkl.dump(paramGrid, open(self.paramGridFile,'wb'))

    def loadParamGrid( self):
        self.features, self.targetPDF =\
          pkl.load(open(self.paramGridFile,'rb'), encoding='latin1') 
            
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
