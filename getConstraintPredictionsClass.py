from fitDistributionToData import fitEquivalentWidthDistribution
import modelProbabilityDistributionOfLensedEquivalentWidths as getModel
import os
import numpy as np
import ipdb as pdb
import pickle as pkl
import corner
from matplotlib import pyplot as plt
class getConstraintPredictions:
    '''
    This will genetate the predictions we will get from a samle of
    observed quasars
    '''
    def __init__( self, inputParams, nObservations, pklRootName, \
                    nBootStraps=5):

        #The number of simualted observed quasars
        self.nObservations = nObservations
        #The number of bootstraps that i will carr out
        #to get the covariance
        self.nBootStraps = nBootStraps
        #The ground truth
        self.inputParameters = inputParams
        self.nParameters = len(self.inputParameters.keys())
        #pklRootName for the saved samples
        self.pklRootName = '%s_%i' % (pklRootName,nObservations)
        
    def getGroundTruthDistribution( self ):
        self.groundTruth = \
          getModel.modelProbabilityDistributionOfLensedEquivalentWidths()
        self.groundTruth.setInterpolatorFunctionParamGrid()
        self.groundTruth.loadParamGrid()

        self.groundTruth.fitInterpolator(loadPklFile=True)


        #Get an input model from the trainer (in future change this to a true
        #model)
        self.groundTruthDist = \
          self.groundTruth.predictPDFforParameterCombination(self.inputParameters)

        
    def bootStrapSamplingForCovariance(self ):
        '''
        Loop over getSamplesForNumberOfObservations, nBootStrap times
        and get the distribution of estimated parameters
        '''

        self.listOfBootstrappedParams = []
        for iBootStrap in range(self.nBootStraps):
            samplesSaveFile = "%s_%i.pkl" % (self.pklRootName, iBootStrap+1)
            iBootStrapParams = \
              self.getSamplesForNumberOfObservations(samplesSaveFile=samplesSaveFile)
            self.listOfBootstrappedParams.append( iBootStrapParams )


    def getSamplesForNumberOfObservations( self, samplesSaveFile=None, verifyDistribution=False ):
        '''
        For a given number of observations, create a mock pdf and 
        sample from this to see what the constraints could be
        '''
        if os.path.isfile(samplesSaveFile):
            samples = pkl.load(open(samplesSaveFile, 'rb'))
            params = self.getParamsFromSamples(samples, statType='median')
            return params
       
        #so the posteior was a delta function on the truth so
        #add some noise by putting in a sampled for a given
        #number of observations
        observedSampledDist = \
          self.groundTruth.generateDistributionOfEquivalentWidths(self.inputParameters, \
                                            np.int(self.nObservations))
      
        if verifyDistribution:
        #Plot to see if the two distributions are the same

            plt.plot(self.groundTruthDist['x'], \
                    self.groundTruthDist['y']/np.max(self.groundTruthDist['y']),'r',label='Truth')
            plt.plot(observedSampledDist['x'], observedSampledDist['y']/\
                    np.max(observedSampledDist['y']),'b', label='Observed')
            plt.legend()
            plt.show()
        
        fitModelClass = \
          fitEquivalentWidthDistribution( observedSampledDist, self.groundTruth)
        fitModelClass.fitProbabilityDistribution()
        fitModelClass.saveSamples(samplesSaveFile)

        return fitModelClass.params

    def getParamsFromSamples( self, samples, statType='median'):
        #And work out some simple statistics (need to add maxLike)
        errorLower, median, errorUpper = \
          np.percentile(samples, [16, 50, 84], axis=0)

        error = np.mean([median - errorLower, errorUpper - median], axis=0)

        return {'params':median, 'error':error}


    def reformBootStrappedParamsForCorner( self ):

        self.bestFitParams = np.zeros((self.nBootStraps,self.nParameters))
        for thisParamSet, iParamSet in enumerate(self.listOfBootstrappedParams):
            
            self.bestFitParams[thisParamSet, :] = iParamSet['params']

    def plotBootStrappedResults( self, figcorner=None, color='red' ):


        if figcorner is None:
            figcorner, axarr = \
              plt.subplots(self.nParameters,self.nParameters,figsize=(7,7))

        truths = [self.inputParameters[i] for i in self.inputParameters.keys()]
        labels = [i for i in self.inputParameters.keys()]
        
        nsamples = self.nBootStraps
        
        corner.corner( self.bestFitParams, truths=truths, labels=labels, \
                fig=figcorner, color=color, \
                smooth=True, levels=(0.20,), plot_datapoints=False,\
                weights=np.ones(nsamples)/nsamples,\
                hist_kwargs={'linewidth':2.},\
                contour_kwargs={'linewidths':2.},\
                bins=20)
                           
if __name__ == '__main__':
    main()
