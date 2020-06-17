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
    def __init__( self, inputParams, nObservations, pklRootName ):

        #The number of simualted observed quasars
        self.nObservations = nObservations
        #The ground truth
        self.inputParameters = inputParams
        self.nParameters = \
          len([i for i in self.inputParameters.keys() if not 'red' in i]) 
        #pklRootName for the saved samples
        self.pklRootName = '%s_%i' % (pklRootName,nObservations)
        self.samplesSaveFile = '%s.pkl' % self.pklRootName
        
    def getGroundTruthDistribution( self ):
        self.groundTruth = \
          getModel.modelProbabilityDistributionOfLensedEquivalentWidths()
        self.groundTruth.setInterpolatorFunctionParamGrid()
        self.groundTruth.loadParamGrid()

        self.groundTruth.fitInterpolator(loadPklFile=True)


        #Get an input model from the trainer (in future change this to a true
        #model)
        self.groundTruthDist = {}
        for iRedshiftBin in self.inputParameters['redshiftBins']:
            self.inputParameters['redshift'] = iRedshiftBin
            self.groundTruthDist[iRedshiftBin] = \
              self.groundTruth.predictPDFforParameterCombination(self.inputParameters)

        
    def getSamplesForNumberOfObservations( self, verifyDistribution=False ):
        '''
        For a given number of observations, create a mock pdf and 
        sample from this to see what the constraints could be
        '''
        
        if os.path.isfile(self.samplesSaveFile):
            samples = pkl.load(open(self.samplesSaveFile, 'rb'))
            return samples
       
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
        fitModelClass.saveSamples(self.samplesSaveFile)

    def getParamsFromSamples( self, samples, statType='median'):
        #And work out some simple statistics (need to add maxLike)
        if statType == 'median':
            errorLower, median, errorUpper = \
              np.percentile(samples, [16, 50, 84], axis=0)

            error = np.mean([median - errorLower, errorUpper - median], axis=0)

            return {'params':median, 'error':error}
        elif statType == 'maxLike':
            errorLower, median, errorUpper = \
              np.percentile(samples, [16, 50, 84], axis=0)
            error = np.mean([median - errorLower, errorUpper - median], axis=0)

            maxLike = []
            for iPar in range(samples.shape[1]):
                y, x = np.histogram( samples[:, iPar], bins=50)
                xc = (x[:-1] + x[1:])/2.
                maxLike.append(xc[ np.argmax( y) ])

            return {'params':np.array(maxLike), 'error':error}    


    def plotBootStrappedResults( self,  figcorner=None, color='red', **kwargs ):
        
        samples = pkl.load(open(self.samplesSaveFile, 'rb'))

        if figcorner is None:
            figcorner, axarr = \
              plt.subplots(self.nParameters,self.nParameters,figsize=(7,7))

        truths = [self.inputParameters[i] for i in self.inputParameters.keys()]
        labels = [i for i in self.inputParameters.keys()]
        
        nsamples = samples.shape[0]

        corner.corner( samples, truths=truths, labels=labels, \
                fig=figcorner, color=color, \
                smooth=True, plot_datapoints=False,\
                weights=np.ones(nsamples)/nsamples,\
                hist_kwargs={'linewidth':2.},\
                contour_kwargs={'linewidths':2.},\
                bins=20, **kwargs)
                           
if __name__ == '__main__':
    main()
