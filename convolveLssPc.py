'''
I want to convolve the Pl from TurboGL with the 
Pc from probablitymagnification
'''
import probabilityMagnification as pm
import numpy as np
from matplotlib import pyplot as plt
import deltaTable as dt
import os as os
import pickle as pkl
import lensingProbabilityDistributionLog as Loglpd
import lensingProbabilityDistribution as lpd
import lensingProbabilityDistributionForSupernova as supernova
from matplotlib import gridspec 
from matplotlib import colors as colors
import matplotlib.cm as cmx
from matplotlib import rc

def main(z=1.0, nMu=1000, nOmegaM=10):
    '''
    Under the linear taylor expansion assumption
    this works out the equivalent widths for different
    alphas for a given redshift of 1.0
    '''
    #seems to be valid only above 0.08
    alphaList = [0.10, 0.30, 0.50, 0.83]
    color = ['r','b','g','c']

    colorMaps = ['Reds','Greens','Blues','Purples']
                                                                          


    sigmaEightList = [0.7, 0.8288, 0.9]
    omegaMatterList = np.linspace(0.1,1.,nOmegaM)
    omegaMatterList = [0.19, 0.3, 0.43]

    for j, alpha in enumerate(alphaList):
           #For aesthetics                                                         
        jet = cm = plt.get_cmap(colorMaps[j])
        cNorm  = colors.Normalize(vmin=-1, vmax=3+1)
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
        ##### 
        for i, omega_m in enumerate(omegaMatterList):
       
                        

      
            lensingPDF = \
              lpd.lensingProbabilityDistribution( redshift=z, \
                                              alpha=alpha, \
                                              omega_m=omega_m, \
                                                nEquivalentWidthBins=nMu,\
                                                  modelType='Linear')
                                                
                                                
            plt.plot( lensingPDF.convolvedPbhPdfWithLssPdf['x'], \
                     lensingPDF.convolvedPbhPdfWithLssPdf['y'],\
                     '-', color=scalarMap.to_rgba(i), alpha=1., \
                      label=r'$\Omega_m=%0.2f$, $\alpha=%0.2f$' % \
                          (omega_m, alpha))
                     
        print lensingPDF.pdfMean
        plt.yscale('log')
        plt.ylim(0.05,130)
        plt.xlim(-0.3,0.3)

    plt.xlabel(r'log$(1+\mu_{\rm EW})$')
    plt.ylabel(r'$P($log$(1.+\mu_{\rm EW}))$')
    plt.legend()
    plt.savefig('../plots/quasarConvolutionOmegaM.pdf')
    plt.show()

    
def diffLogLinear(z=1.0, nMu=1000):
    '''
    This script will look at the difference between the pdf's
    if we assume log or linear convestaion

    They are identical, so therefore identical mathematically
    so probably should talk to Jon about this
    '''


    alphaList = [ 0.1, 0.3,0.5,0.83]
    gs = gridspec.GridSpec( 4,1 )
    mainAx = plt.subplot(gs[0:3,0])
    diffAx = plt.subplot(gs[3:,0])
    gs.update(hspace=0.0)
    for alpha in alphaList:

        
        lensingPDF = \
          lpd.lensingProbabilityDistribution( redshift=z, \
                                              alpha=alpha, \
                                                nEquivalentWidthBins=nMu)
                                                
        logLensingPDF = \
          Loglpd.lensingProbabilityDistribution( redshift=z, \
                                              alpha=alpha, \
                                                nEquivalentWidthBins=nMu)


                                                
        mainAx.plot(lensingPDF.convolvedPbhPdfWithLssPdf['x'], \
                     lensingPDF.convolvedPbhPdfWithLssPdf['y'],\
                     label='%0.2f' % alpha)

        diffAx.plot(lensingPDF.convolvedPbhPdfWithLssPdf['x'], \
                     lensingPDF.convolvedPbhPdfWithLssPdf['y']/
                     logLensingPDF.convolvedPbhPdfWithLssPdf['y']-1,\
                     label='%0.2f' % alpha)
        
        print lensingPDF.pdfMean - logLensingPDF.pdfMean
    mainAx.set_yscale('log')
    mainAx.set_xscale('log')

    mainAx.set_ylim(0.05,130)
    mainAx.set_xlim(-0.6,0.6)
    diffAx.set_xlim(-0.6,0.6)
    diffAx.set_ylim(-0.01,0.01)
    diff = np.sum(lensingPDF.convolvedPbhPdfWithLssPdf['y']/\
                     logLensingPDF.convolvedPbhPdfWithLssPdf['y']-1)
    print("Diff is %0.3f"% diff)
    mainAx.legend()
    mainAx.set_xticklabels([])
    plt.show()


def testSupernovaDistribution(z=1.0,nMu=1000 ):
    alphaList = [ 0.1, 0.3,0.5,0.83]
    gs = gridspec.GridSpec( 4,1 )
    mainAx = plt.subplot(gs[0:,0])
    
    color = ['r','b','g','c']
    for i, alpha in enumerate(alphaList):

        
        lensingPDFlog = \
          supernova.lensingProbabilityDistribution( redshift=z, \
                                              alpha=alpha, \
                                                nKappaBins=nMu,\
                                                modelType='Log')
   
                                                
        mainAx.plot(lensingPDFlog.convolvedPbhPdfWithLssPdf['x'], \
                     lensingPDFlog.convolvedPbhPdfWithLssPdf['y'],\
                     label='%0.2f' % alpha,  color=color[i])
        
         
        lensingPDF = \
          supernova.lensingProbabilityDistribution( redshift=z, \
                                              alpha=alpha, \
                                                nKappaBins=nMu,\
                                                modelType='Linear')
   
                                                
        mainAx.plot(lensingPDF.convolvedPbhPdfWithLssPdf['x'], \
                     lensingPDF.convolvedPbhPdfWithLssPdf['y'],\
                     label='%0.2f' % alpha, ls='--',  color=color[i])
        
    mainAx.set_yscale('log')
    mainAx.set_xscale('log')

    mainAx.set_ylim(0.05,130)
    mainAx.set_xlim(0.,0.6)
  
    mainAx.legend()

    plt.show()

    
if __name__ == '__main__':
    main()
