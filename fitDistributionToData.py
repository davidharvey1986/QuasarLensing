


import emcee

import numpy as np

def lnProbOfDataGivenModel( theta, xTrue, yTrue, error ):

    theoreticalPrediction = \
      getModelOfProbabilityDistribution( xTrue, **theta)
      
    priorOnParameters = getPrior( theta )

    likelihood = \
      np.sum(norm.logpdf( theoreticalPrediction, yTrue, scale=error))
    
    posterior = priorOnParameters*likelihood

       
    return posterior


def priorOnParameters():
    return 1.


class fitEquivalentWidthDistribution:

    def __init__( self, inputProbabliltyDistirubtoin ):

          
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

        '''
                   
        self.pdf = pdf
       
        self.fitProbabilityDistriubtion()

    def fitProbabilityDistriubtion( self, nthreads=4, **kwargs):

        ### options for the sampling
        nwalkers = 20
        ndim = 3
        burn_len=500
        chain_len=1000
        #####
        
        pos0 = np.random.rand(nwalkers,ndim)
      
        args = (self.pdf['x'], self.pdf['y'], \
            self.pdf['error'], self.hubbleInterpolator )

        dmsampler = \
          emcee.EnsembleSampler(nwalkers, ndim, \
                            lnprob, args=args, \
                                    threads=nthreads)

        #DO a burn in 
        pos, prob, state  = dmsampler.run_mcmc(pos0, burn_len)


        #Now burnt in, do the proper sampling
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
