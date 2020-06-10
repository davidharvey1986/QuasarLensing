
import numpy as np
import corner
from matplotlib import pyplot as plt
import ipdb as pdb
import os
import pickle as pkl
from getConstraintPredictionsClass import getConstraintPredictions

def main(inputParameters = {'alpha':0.4, 'scale':0.4, 'redshift':1.0}):
    '''
    Test the mcmc fitting algorithm by simulating a pdf and fitting to it
    '''
    ndim = len(inputParameters.keys())
    figcorner, axarr = plt.subplots(ndim,ndim,figsize=(7,7))

    plottingColours = ['red','b','g','c','y','k','grey']
    
    #Do it over a range of mock samples sizes
    listOfMockSampleSizes = 10**np.linspace(3,6,4)
    for iColor, nObservations in enumerate(listOfMockSampleSizes):
        pklRootName = 'pickles/predictConstraints'
        predictionClass = \
          getConstraintPredictions( inputParameters, nObservations, pklRootName, \
                                        nBootStraps=10)
        
        predictionClass.getGroundTruthDistribution()
        predictionClass.bootStrapSamplingForCovariance()
        predictionClass.reformBootStrappedParamsForCorner()
        predictionClass.plotBootStrappedResults(figcorner=figcorner, \
                                color=plottingColours[iColor] )
                                
    plt.show()


if __name__ == '__main__':
    main()
