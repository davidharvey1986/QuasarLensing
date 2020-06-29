
import numpy as np
import corner
from matplotlib import pyplot as plt
import ipdb as pdb
import os
import pickle as pkl
from getConstraintPredictionsClass import getConstraintPredictions
from astro_tools import colourFromRange
from matplotlib import gridspec

def main(inputParameters = {'alpha':0.4, 'scale':0.4}, nRedshiftCombos=10, nSampleSizes=3):
    '''
    Test the mcmc fitting algorithm by simulating a pdf and fitting to it
    '''

    fittingParams =  {'chain_len':1000, 'burn_len':100}
    
    #ndim = len(inputParameters.keys()) 
    #figcorner, axarr = plt.subplots(ndim,ndim,figsize=(10,10))
    fig = plt.figure()
    fig.subplots_adjust(hspace=0)
    gs = gridspec.GridSpec( 2, 1)
    axAlpha = plt.subplot(gs[0,0])
    axScale = plt.subplot(gs[1,0])
    
    plottingColours = ['red','b','g','c','y','k','grey']
    cRange = colourFromRange([0, nSampleSizes])
    
    
    #Do it over a range of mock samples sizes
    listOfMockSampleSizes = 10**np.linspace(3,5,nSampleSizes)

    #loop through the different number of redshift bins to see how this helps
    redshiftBins = np.linspace(0.5, 5.5, nRedshiftCombos+1)
    
    for iColor, nObservations in enumerate(listOfMockSampleSizes):
        error = None
        #Number of redshift bin combinations
        for iRedshiftBinCombination in range(nRedshiftCombos):

            #add the redshift combinations to the input params
            inputParameters['redshiftBins'] = [redshiftBins[iRedshiftBinCombination]]
        
            #the save file
            pklRootName = 'pickles/predictConstraintsConserveNgalaxies'
            for i in inputParameters['redshiftBins']:
                pklRootName +=  "_z%0.1f" % (i) 

            #initialise the prediction class, conserve total number of quasars
            nObsPerRedshiftBins = nObservations/(iRedshiftBinCombination+1)
            predictionClass = \
              getConstraintPredictions( inputParameters, nObsPerRedshiftBins, pklRootName )

            #check for the save file with the mcmc samples in
            if not predictionClass.checkForSampleFile():
                predictionClass.getGroundTruthDistribution()
                predictionClass.getSamplesForNumberOfObservations(fittingParams=fittingParams)
            predictionClass.plotBootStrappedResults()
            
            plt.show()
            #get a dict with the median params in
            params = predictionClass.getBestFitParams()

            if error is None:
                error = params['error']/params['params']
            else:
                error = np.vstack((error, params['error']))
                

        axAlpha.plot( redshiftBins[:-1], error[:,0], color=plottingColours[iColor], \
                          label='log(nObs)=%i' % np.log10(nObservations))
        axScale.plot( redshiftBins[:-1], error[:,1], color=plottingColours[iColor])

    axAlpha.set_yscale('log')
    axScale.set_yscale('log')
    axAlpha.legend(loc=3)
    axAlpha.set_ylabel(r'$\sigma_\alpha/\alpha$')
    axScale.set_ylabel(r'$\sigma_{\sigma{\rm int}}/\sigma_{\rm int}$')

    axScale.set_xlabel(r'Redshift')

        
    plt.savefig('../plots/sensitivity_PBH.pdf')
    plt.show()


if __name__ == '__main__':
    main()
