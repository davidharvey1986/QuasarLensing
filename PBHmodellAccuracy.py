from modelProbabilityDistributionOfLensedEquivalentWidths import modelProbabilityDistributionOfLensedEquivalentWidths
import numpy as np
from matplotlib import pyplot as plt
import os
from matplotlib import gridspec
from astro_tools.getColourFromRange import colourFromRange
import ipdb as pdb
def main():
    '''
    What is the accuracy of the GPR to predict the PBH PDF
    '''
    newModel = modelProbabilityDistributionOfLensedEquivalentWidths()
    newModel.setInterpolatorFunctionParamGrid()
    if os.path.isfile(newModel.paramGridFile):
        newModel.loadParamGrid()
    else:
        newModel.fillParamGrid()
        newModel.saveParamGrid()
    newModel.targetPDF[newModel.targetPDF==0] = 1e-80
    newModel.targetPDF = np.log10(newModel.targetPDF)
    
    newModel.fitInterpolator(loadPklFile=False, \
                pklFile='pickles/logInterplator.pkl')
    
    predictAlphas = np.linspace(0.2,0.4, 8)
    
    fig=plt.figure()
    gs = gridspec.GridSpec(5,1)
    ax1 = plt.subplot(gs[:3,0])
    ax2 = plt.subplot(gs[3:,0])
    fig.subplots_adjust(hspace=0)
    for iAlpha in predictAlphas:
        predictThese = {'alpha':iAlpha, 'scale':0.3, 'redshift':1.0}
        
        trueDist = newModel.getNewDistribution(predictThese)

        predictedPDF = \
          newModel.predictPDFforParameterCombination(predictThese, \
                                    xVector=trueDist['x'] )

        ax1.plot(trueDist['x'], predictedPDF['y'])
        
        ax2.plot(trueDist['x'], \
            np.cumsum(predictedPDF['y'])/np.sum(predictedPDF['y']) - \
            np.cumsum(trueDist['y'])/np.sum(trueDist['y']))
        
    #ax1.set_yscale('log')
    #ax2.set_yscale('log')
    ax2.set_ylim(-2e-3,2e-3)
    ax1.set_ylim(1e-4,2.2)
    ax1.set_xlim(-1,2.2)
    ax2.set_xlim(-1,2.2)

    plt.show()
        
 
def extrapolationToLowAlpha(nAlphas=10):
    '''
    If i extrapolate the PDF does it behave okay?

    '''
    newModel = modelProbabilityDistributionOfLensedEquivalentWidths()
    newModel.setInterpolatorFunctionParamGrid()
    if os.path.isfile(newModel.paramGridFile):
        newModel.loadParamGrid()
    else:
        newModel.fillParamGrid()
        newModel.saveParamGrid()

    newModel.fitInterpolator(loadPklFile=True)

    cRange = colourFromRange( [0, nAlphas], 'rainbow')
    predictAlphas = np.linspace(0.01,0.83, nAlphas)

    fig=plt.figure()
    gs = gridspec.GridSpec(5,1)
    ax1 = plt.subplot(gs[:3,0])
    ax2 = plt.subplot(gs[3:,0])
    fig.subplots_adjust(hspace=0)
    for iColour, iAlpha in enumerate(predictAlphas):
        predictThese = {'alpha':iAlpha, 'scale':0.3, 'redshift':1.0}
        
        trueDist = newModel.getNewDistribution(predictThese)

        predictedPDF = \
          newModel.predictPDFforParameterCombination(predictThese, \
                                    xVector=trueDist['x'] )

                                
        label=r'$\alpha=%0.2f$' % iAlpha
        ax1.plot(trueDist['x'], trueDist['y'], '-', \
                     color=cRange.getColour(iColour), label=label)
        ax1.plot(trueDist['x'], predictedPDF['y'], '--',\
                     color=cRange.getColour(iColour), label=label)

        ax2.plot(trueDist['x'], \
            np.cumsum(predictedPDF['y'])/np.sum(predictedPDF['y']) - \
            np.cumsum(trueDist['y'])/np.sum(trueDist['y']))
        
    ax1.set_yscale('log')
    ax1.legend()
    #ax2.set_yscale('log')
    ax2.set_ylim(-2e-3,2e-3)
    ax1.set_ylim(1e-4,2.2)
    ax1.set_xlim(-1,2.2)
    ax2.set_xlim(-1,2.2)

    plt.show()
                   
if __name__ == '__main__':
    main()
