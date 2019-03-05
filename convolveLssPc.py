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
import lensingProbabilityDistribution as lpd

def main(z=1.0, nMu=1000):

    alphaList = [0.,0.3,0.5,0.83]
    for alpha in alphaList:
        lensingPDF = \
        lpd.lensingProbabilityDistribution( redshift=z, \
                                              alpha=alpha, \
                                                nMagnitudeBins=nMu)
        dMu = (lensingPDF.convolvedPbhPdfWithLssPdf['x'][1] - lensingPDF.convolvedPbhPdfWithLssPdf['x'][0])


        lensingPDF.convolvedPbhPdfWithLssPdf['y'] /= np.sum(lensingPDF.convolvedPbhPdfWithLssPdf['y']*dMu)
        
        print np.sum(lensingPDF.convolvedPbhPdfWithLssPdf['y']*dMu)
        print np.sum(lensingPDF.convolvedPbhPdfWithLssPdf['y']*dMu*lensingPDF.convolvedPbhPdfWithLssPdf['x'])

        print lensingPDF.meanMagnification
                #TOTAL PDF
        plt.plot(lensingPDF.convolvedPbhPdfWithLssPdf['x'], \
                     lensingPDF.convolvedPbhPdfWithLssPdf['y'],\
                     label='%0.2f' % alpha)

      
        
                     
        plt.yscale('log')
        plt.ylim(0.05,35)
        plt.xlim(0.,0.6)
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()
