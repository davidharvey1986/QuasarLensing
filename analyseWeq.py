from matplotlib import pyplot as plt
from matplotlib import gridspec
import numpy as np
import pyfits as fits
import lensing as lensing
from emissionLine import *
import sys

def exampleEquivalentWidthConvolutedWithLensingAndCompared():
    '''
    This example will take an emission line
    HBeta for example

    Determine the intrinsic EW distribution from low redshift
    quasars. Assuming now evolution, it will convolve with the
    lensing probability distribution of redshift one quasars 
    aand then compare to the redshift one quasars.

    Lots of to dos noted in the emissionLine class
    '''

    #get the lensing probability from convolvelssPc.
    #currently only for z=1., and alpha=0.83
    #future models will sample these two paraemeters to find
    #the best fitting for the observed tomorgraphic sampels
    EquivilentWidthBins =  10**np.linspace(-3, 3, 30)

    alphaList = [0.1,0.3,0.5,0.83]
    gs = gridspec.GridSpec( 5, 1)
    axis1 = plt.subplot(gs[:,0])
    
    
    for index, iAlpha in enumerate(alphaList):
    #initiate the class
        emissionLineCl = emissionLine('NARROW_HB', nRedshiftBins=2macs )

        emissionLineCl.getLensingProbability( z=1.0, alpha=iAlpha )
        emissionLineCl.applyLuminosityCut( luminosityCut=10.)
        emissionLineCl.getEquvilentWidthMeansAsFunctionOfRedshift()

        
        emissionLineCl.histogramEquivalentWidth(nEquivilentWidthBins=EquivilentWidthBins)
        emissionLineCl.logGaussianFitEquivalentWidths()
        
    #get what i am calling the intrinsic distribution of HB narrow EQ
    #this can be discussed.the redshift cut determines thoise quasars
    #that havent been lensed, but still need enough quasars to get a good
    #distribution
        emissionLineCl.getIntrinsicDistribution( intrinsicDistribution= 'data' )
    

    #now convolve the the two together to get the expected
    #distribtuion of EW at a redshift of 1.0 given an intrinsic
    #distribution of quasars
        emissionLineCl.convolveIntrinsicEquivalentWidthWithLensingProbability()


    #plot the resulting convolution
        axis1.plot(emissionLineCl.predictedLensedEquivalentWidthDistribution['x'],\
                    emissionLineCl.predictedLensedEquivalentWidthDistribution['y'], \
                    label=r'$\alpha=%0.2f$' % iAlpha)
        #Since the pdfs have different x values if i want to compare, i need
        #to interpolate between them
     

        if index == len(alphaList)-1:
            axis1.plot(emissionLineCl.intrinsicEquivalentWidthDistribution['x'],\
                         emissionLineCl.intrinsicEquivalentWidthDistribution['y'], label='Intrinsic Distribution')
    plt.legend(loc=2, prop={'size': 6})
#    plt.ylim(0,4.3)

    axis1.set_xlim(-2,2.5)
    axis1.set_ylim(0.0001,2.)
    axis1.set_yscale('log')
    plt.xlabel(r'log(Equivalent Width)')
    plt.ylabel(r'p(log(Equivalent Width))')
    plt.savefig('../plots/effectOfLensingOnEquivalentWidths.pdf')
    plt.show()
    
    


def examplePlotEquivalentWidthHistograms(  ):
    '''
    A script to analyse and bin the OIII, HBeta and the MgIII lines 

    Note that i should be careful about the equivalent widths
    on the web page it states that the only narrow line is in 
    fact the HBeta
    
    So now that seems like the one we trust.

    These three lines therefore show interesting things.
    1. That the MGII has some strange evolution with redshift
    2. The OIII does show any
    3. THe Narrow HB shows some shift. Which could be interesting.

    '''


    emissionLineList = ['MGII','OIII_5007','NARROW_HB']
    emissionLineList = ['NARROW_HB']

    fig, axarr = plt.subplots(len(emissionLineList),1)
    if len(emissionLineList) == 1:
        axarr = [axarr]
    nEquivilentWidthBins = 30
    maxEquivilentWidth = np.array([1000.,100., 50., 20.])
    luminosityCuts = np.linspace( 8,10,5)
    color=['r','b','g','c','orange']
    for iColor, iLuminosityCut in enumerate(luminosityCuts):
        for iAxis, iEmissionLine in enumerate(emissionLineList):
        
            emissionLineCl = emissionLine(iEmissionLine, nRedshiftBins=1)
            emissionLineCl.applyLuminosityCut( luminosityCut=iLuminosityCut)

            emissionLineCl.getEquvilentWidthMeansAsFunctionOfRedshift()
            EquivilentWidthBins = \
              10**np.linspace(-3, np.log10(maxEquivilentWidth[iAxis]), nEquivilentWidthBins)

        
            emissionLineCl.histogramEquivalentWidth(nEquivilentWidthBins=EquivilentWidthBins)
            emissionLineCl.logGaussianFitEquivalentWidths()

            
            redshiftList = emissionLineCl.equivilentWidthHistograms.keys()
            redshiftListSorted = \
              list((np.sort(np.array(redshiftList).astype(float))).astype(str))
            
            for iRedshift in redshiftListSorted:
            
                axarr[iAxis].plot(emissionLineCl.equivilentWidthHistograms[iRedshift]['x'], \
                    emissionLineCl.equivilentWidthHistograms[iRedshift]['y'], \
                    label=iLuminosityCut, color=color[iColor])
                
                emissionLineCl.plotLogGaussianFitEquivalentWidths(iRedshift,\
                                                    **{'color':color[iColor],'ls':'--'})

                      
            axarr[iAxis].text(0.1,0.5,iEmissionLine,  transform=axarr[iAxis].transAxes)
            axarr[iAxis].legend()
            axarr[iAxis].set_xscale('log')

            axarr[iAxis].set_xlim(0,maxEquivilentWidth[iAxis])
    axarr[iAxis].set_xlabel('Equivalent Width')
    axarr[iAxis].set_ylabel('P(Equivalent Width)')
    
    plt.savefig('../plots/emissionLineHistogram.pdf')
    plt.show()

def main( catalogue='dr9_weq.fits' ):
    '''
    Analyse the dr9 boss quasar catalogue from 
    http://quasar.astro.illinois.edu/BH_mass/dr9.
    '''

    data = fits.open(catalogue)[1].data

    dataKeys = data.columns.names
    
    ax = plt.gca()
    for iKey in dataKeys:
        
        if (not 'REW' in iKey) | ('ERR' in iKey):
            continue
        print iKey


        
        zBin, weqBin = \
          redshiftBinWeq( data['Z_PIPE'], data[iKey])
        lineName = ' '.join(iKey.split('_')[1:3])
        print data[iKey]
        angdist = lensing.ang_distance(zBin)
        ax.errorbar(zBin, weqBin[0,:], \
                        yerr=weqBin[1,:], label=lineName)
        
    ax.legend()
    ax.set_xlim(0.,1.)
    ax.set_xlabel(r'z')
    ax.set_ylabel(r'\Delta\mu')


        
    

    
def compareCats():


    cat1 = 'DR12Q.fits'
    cat2 = 'dr9_weq.fits'
    
    dataA = fits.open(cat1)[1].data
    dataB = fits.open(cat2)[1].data

    dataAKeys = dataA.columns.names
    dataBKeys = dataB.columns.names


    REW_keys = [ i for i in dataAKeys if ('REW' in i) & (not 'ERR' in i)]

    fig, axarr = plt.subplots( 5, 1 )
    plt.subplots_adjust(wspace=0, hspace=0)
    iPlot = 0
    for i, iKey in enumerate(REW_keys):

        

        

        lineName = iKey.split('_')[1]
        print lineName
        if (lineName == 'BROAD') | (lineName == 'NARROW'):
            lineName = iKey.split('_')[2]
        keyB = [ j for j in dataBKeys if (lineName in j) & ('REW' in j) & (not 'ERR' in j) & (not 'FPG' in j)]
        print keyB, iKey
        if len(keyB) == 1:
          
          ax = axarr[iPlot]

          if not 'REWE' in iKey:
            label = lineName+' absorption'
          else:
            label = lineName
          iPlot += 1
          zBinA, weqBinA = \
            redshiftBinWeq( dataA['Z_PIPE'], dataA[iKey])

        
          zBinB, weqBinB = \
            redshiftBinWeq( dataB['Z_PIPE'], dataB[keyB[0]])

          ratio = (weqBinA[0,:]/weqBinB[0,:])
          error = np.sqrt((weqBinA[1,:]/weqBinA[0,:])**2 +\
                            (weqBinB[1,:]/weqBinB[0,:]))*ratio
          ax.errorbar(zBinA, ratio, yerr=error, label=label)
        
          ax.set_ylim(0.1,3)
          
          ax.plot([0.,6.],[1.,1.],'--')
          ax.legend()
          if iPlot == 3:
              ax.set_ylabel(r'DR12 PARIS ET AL / DR9 SHEN ET AL')

    ax.set_xlabel('Redshift')


    plt.savefig('CompareDR12toDR9.pdf')


if __name__ == '__main__':
    if len(sys.argv) == 1:
        exampleEquivalentWidthConvolutedWithLensingAndCompared()
    else:
        if sys.argv[1] == 'hist':
            examplePlotEquivalentWidthHistograms()
