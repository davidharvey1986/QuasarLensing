
import numpy as np
import corner
from matplotlib import pyplot as plt
import ipdb as pdb
import os
import pickle as pkl
from getConstraintPredictionsClass import getConstraintPredictions

def main(inputParameters = {'alpha':0.4, 'scale':0.4, 'redshiftBins':[1.0,2.0]}):
    '''
    Test the mcmc fitting algorithm by simulating a pdf and fitting to it
    '''
    ndim = len(inputParameters.keys()) - 1
    figcorner, axarr = plt.subplots(ndim,ndim,figsize=(2,2))

    plottingColours = ['red','b','g','c','y','k','grey']
    
    #Do it over a range of mock samples sizes
    listOfMockSampleSizes = 10**np.linspace(3,5,3)
    for iColor, nObservations in enumerate(listOfMockSampleSizes):
        pklRootName = 'pickles/predictConstraints'
        predictionClass = \
          getConstraintPredictions( inputParameters, nObservations, pklRootName )
    
        predictionClass.getGroundTruthDistribution()
        predictionClass.getSamplesForNumberOfObservations()
        axarr[0,1].plot(0,0,color=plottingColours[iColor],\
                    label = r'log(nObs)=%i' % np.log10(nObservations))
        predictionClass.plotBootStrappedResults(figcorner=figcorner, \
                                color=plottingColours[iColor] )

        plt.show()
    axarr[0,1].legend()
    plt.savefig('../plots/test_PDB.pdf')
    plt.show()


if __name__ == '__main__':
    main()
