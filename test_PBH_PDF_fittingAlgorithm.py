
from fitDistributionToData import fitEquivalentWidthDistribution
import modelProbabilityDistributionOfLensedEquivalentWidths as getModel
import numpy as np
import corner
from matplotlib import pyplot as plt
def main():
    '''
    Test the mcmc fitting algorithm by simulating a pdf and fitting to it
    '''

    #What combination of parameters I want to test?
    inputParameters = \
      {'alpha':0.2, 'scale':0.4}
      
    pbdPDFmodel = \
      getModel.modelProbabilityDistributionOfLensedEquivalentWidths()
    pbdPDFmodel.setInterpolatorFunctionParamGrid()
    #pbdPDFmodel.fillParamGrid()
    pbdPDFmodel.fitInterpolator(loadPklFile=True)


    #Get an input model from the trainer (in future change this to a true
    #model)
    trueDist = \
      pbdPDFmodel.predictPDFforParameterCombination(inputParameters)

    trueDist['error'] = np.ones(len(trueDist['x']))/\
      len(trueDist['x'])
    fitModelClass = fitEquivalentWidthDistribution( trueDist, pbdPDFmodel)
    fitModelClass.fitProbabilityDistriubtion()

    truths = [inputParameters[i] for i in inputParameters.keys()]
    labels = [i for i in inputParameters.keys()]
    corner.corner( fitModelClass.samples, truths=truths, labels=labels)
    fitModelClass.saveSamples('pickles/testDataFitterPickles.pkl')
    plt.savefig('../plots/testDataFitter.pdf')
    plt.show()

if __name__ == '__main__':
    main()
