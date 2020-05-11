

import ipdb as pdb
import emcee
import numpy as np
from scipy.stats import norm

import pickle as pkl
def lnProbOfDataGivenModel( theta, xTrue, yTrue,  modelPredictor ):

    #theta is a just a list, need to turn it in to a dict
    thetaDict = {'alpha': theta[0], 'scale':theta[1]}
    theoreticalPrediction = modelPredictor.predictPDFforParameterCombination( thetaDict, xVector=xTrue)

    priorOnParameters = \
      getPriorOnParameters( theta, modelPredictor )

    cumsumTheory = np.cumsum(theoreticalPrediction['y'])
    cumsumTrue = np.cumsum(yTrue)
    
    likelihood = 1./np.sum( (cumsumTheory - cumsumTrue)**2)*len(cumsumTrue)
    
    posterior = priorOnParameters*likelihood
    #pdb.set_trace()
    if np.isfinite(posterior) == False:
        posterior = -np.inf
    return posterior


def getPriorOnParameters( theta, modelPredictor):
    '''
    For now the prior is just a flat one within the bounds
    of the trained interpolator
    '''
    thetaDict = {'alpha': theta[0], 'scale':theta[1]}

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

    def fitProbabilityDistriubtion( self, nthreads=4, **kwargs):

        ### options for the sampling
        nwalkers = 20
        ndim = len(self.modelClass.interpolateParams.keys())
        burn_len=100
        chain_len=2000
        #####
        #initial set up of samplers
        pos0 = np.zeros((nwalkers,ndim))
    
        for iPos, iDim in \
          enumerate(self.modelClass.interpolateParams.keys()):
            pos0[:,iPos] = \
              np.random.uniform(self.modelClass.interpolateParams[iDim][0],\
                                self.modelClass.interpolateParams[iDim][-1], nwalkers)
       
        args = (self.pdf['x'], self.pdf['y'], \
             self.modelClass )

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
