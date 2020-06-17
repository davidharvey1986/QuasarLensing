

import ipdb as pdb
import emcee
import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt 
import pickle as pkl
def lnProbOfDataGivenModel( theta, pdfDict, modelPredictor ):
    '''
    yError in this case is teh error in the cumulative sum
    '''
    #theta is a just a list, need to turn it in to a dict
   
    posterior = 0
    for iRedshiftBin in pdfDict.keys():
        thetaDict = {'alpha': theta[0], 'scale':theta[1], 'redshift':iRedshiftBin}

        xTrue = pdfDict[iRedshiftBin]['x']
        yTrue = pdfDict[iRedshiftBin]['y']
        yError = pdfDict[iRedshiftBin]['error']
        
        theoreticalPrediction = modelPredictor.predictPDFforParameterCombination( thetaDict, xVector=xTrue)

        priorOnParameters = getPriorOnParameters( thetaDict, modelPredictor )

        cumsumTheory = np.cumsum(theoreticalPrediction['y'])/np.sum(theoreticalPrediction['y'])
        cumsumTrue = np.cumsum(yTrue)/np.sum(yTrue)

        #for the uneven error bars
        upperIndex = (cumsumTrue >= cumsumTheory) & (yError[1] != 0)
        likelihoodUpper = norm.logpdf( cumsumTheory[upperIndex], \
                    cumsumTrue[upperIndex], scale=yError[1][upperIndex])
                    
        lowerIndex = (cumsumTrue < cumsumTheory) & (yError[0] != 0)
        likelihoodLower = norm.logpdf( cumsumTheory[lowerIndex], \
                        cumsumTrue[lowerIndex], scale=yError[0][lowerIndex])

        likelihood = np.nansum(likelihoodUpper) + np.nansum(likelihoodLower)

        posterior += priorOnParameters*likelihood


    if np.isfinite(posterior) == False:
        posterior = -np.inf
    return posterior


def getPriorOnParameters( thetaDict, modelPredictor):
    '''
    For now the prior is just a flat one within the bounds
    of the trained interpolator
    '''

    prior = 1
    
    for iKey in thetaDict.keys():
        

        minParam = \
          np.min(modelPredictor.interpolateParams[iKey])
        maxParam = \
          np.max(modelPredictor.interpolateParams[iKey])
          
        if thetaDict[iKey] < minParam:
            return -np.inf
        if thetaDict[iKey] > maxParam:
            return -np.inf
        
    return prior
          
    
    
    return 1.


class fitEquivalentWidthDistribution:

    def __init__( self, inputProbabliltyDistribution, modelClass ):

          
        '''
        This class is going to fit the predicted lensing
        distribution to the data. 

        For now we are using simply 2 parameters:
        alpha: the fraction of PbH in the Universe
        intrinsicWidth: the width of the (assumed) 
        log normal distribution of intrinsic equivalent widths
        
        input: inputProbabilityDistribution: a dict of 'x' and 'y'
        where 'x': the list of equivalent widths
              'y': the correspodning probabilities
              'error': error in the y values

        modelClass: a class of the input model from 
           "modelProbabilityDistributionfOfLensedEquivalentWidths.py"

        '''
                   
        self.pdf = inputProbabliltyDistribution
        self.modelClass = modelClass
      
    def fitProbabilityDistribution( self, nthreads=4, **kwargs):

        ### options for the sampling
        nwalkers = 20
        ndim = len(self.modelClass.interpolateParams.keys())-1
        burn_len=1000
        chain_len=10000
        #####
        #initial set up of samplers
        pos0 = np.zeros((nwalkers,ndim))

        for iPos, iDim in \
          enumerate(self.modelClass.interpolateParams.keys()):
          if iDim == 'redshift':
            continue
          pos0[:,iPos] = \
              np.random.uniform(self.modelClass.interpolateParams[iDim][0],\
                                self.modelClass.interpolateParams[iDim][-1], nwalkers)
       
        args = (self.pdf, self.modelClass )

        dmsampler = \
          emcee.EnsembleSampler(nwalkers, ndim, \
                            lnProbOfDataGivenModel, \
                            args=args, \
                            threads=nthreads)

        #DO a burn in
        print("Burning In")
        pos, prob, state  = dmsampler.run_mcmc(pos0, burn_len, \
                                                   progress=True)


        #Now burnt in, do the proper sampling
        print("Sampling")
        pos, prob, state  = dmsampler.run_mcmc(pos, chain_len,\
                                        progress=True)

        #Get the samples of the chain
        self.samples = dmsampler.flatchain

        #And work out some simple statistics
        errorLower, median, errorUpper = \
          np.percentile(self.samples, [16, 50, 84], axis=0)

        error = np.mean([median - errorLower, errorUpper - median], axis=0)

        self.params = {'params':median, 'error':error}


    def getPredictedProbabilities( self, xInput=None ):
        if xInput is None:
            xInput = self.xNoZeros

        return 10**self.fitFunc( xInput, *self.params['params'])
    
   

    def saveSamples( self, pklFile):
        pkl.dump( self.samples, open(pklFile, 'wb'))
