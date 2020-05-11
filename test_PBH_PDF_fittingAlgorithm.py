
from fitDistributionToData import fitEquivalentWidthDistribution
import modelProbabilityDistributionOfLensedEquivalentWidths as getModel
import numpy as np
import corner
from matplotlib import pyplot as plt
import ipdb as pdb
import os
import pickle as pkl

def main(inputParameters = {'alpha':0.4, 'scale':0.4}):
    '''
    Test the mcmc fitting algorithm by simulating a pdf and fitting to it
    '''
    ndim = 2
    figcorner, axarr = plt.subplots(ndim,ndim,figsize=(7,7))

    plottingColours = ['red','b','g','c','y','k','grey']
    
    #Do it over a range of mock samples sizes
    listOfMockSampleSizes = 10**np.linspace(3,6,4)
    for iColor, nObservations in enumerate(listOfMockSampleSizes):

        #Save the samples in a pickle file, check to see if exists
        pklFile = 'pickles/testDataFitterPickles_logPost_%i.pkl' \
                                  %np.log10(nObservations)
        print("Saving in %s" % pklFile)
        if os.path.isfile(pklFile):
            samples = pkl.load(open(pklFile, 'rb'))
        else:
            samples = \
              getSamplesForNumberOfObservations( nObservations, pklFile, \
                                    inputParameters=inputParameters)

                                    
        truths = [inputParameters[i] for i in inputParameters.keys()]
        labels = [i for i in inputParameters.keys()]
        nsamples = len(samples)
        corner.corner( samples, truths=truths, labels=labels, \
            fig=figcorner, color=plottingColours[iColor], \
            smooth=True, levels=(0.20,), plot_datapoints=False,\
            weights=np.ones(nsamples)/nsamples,\
                           hist_kwargs={'linewidth':2.},\
                        contour_kwargs={'linewidths':2.},\
                           bins=20)


        axarr[1,1].plot(0,0, label='log(Nobs)=%i' % \
                          np.log10(nObservations), \
                        color=plottingColours[iColor])
    axarr[1,1].legend(loc=2)
    plt.savefig('../plots/testDataFitter.pdf')
    plt.show()



def getSamplesForNumberOfObservations( nObservations, pklFile, \
                                           inputParameters=None):
    '''
    For a given number of observations, create a mock pdf and 
    sample from this to see what the constraints could be
    '''
    
    #What combination of parameters I want to test?
    if inputParameters is None:
        inputParameters =  {'alpha':0.4, 'scale':0.4}
      
    pbdPDFmodel = \
      getModel.modelProbabilityDistributionOfLensedEquivalentWidths()
    pbdPDFmodel.setInterpolatorFunctionParamGrid()
    pbdPDFmodel.fitInterpolator(loadPklFile=True)


    #Get an input model from the trainer (in future change this to a true
    #model)
    trueDist = \
      pbdPDFmodel.predictPDFforParameterCombination(inputParameters)
    #so the posteior was a delta function on the truth so add some noise by putting in a sampled for a given
    #number of observations
    observedSampledDist = \
      pbdPDFmodel.generateDistributionOfEquivalentWidths(inputParameters, \
                                            np.int(nObservations))
      
    '''
    #Plot to see if the two distributions are the same

    plt.plot(trueDist['x'], trueDist['y']/np.max(trueDist['y']),'r')
    plt.plot(observedSampledDist['x'], \
                 observedSampledDist['y']/\
                 np.max(observedSampledDist['y']),'r')
    plt.show()
    '''
    fitModelClass = \
      fitEquivalentWidthDistribution( observedSampledDist, pbdPDFmodel)
    fitModelClass.fitProbabilityDistriubtion()
    fitModelClass.saveSamples(pklFile)

    return fitModelClass.samples
    
if __name__ == '__main__':
    main()
