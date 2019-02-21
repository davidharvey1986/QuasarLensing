from matplotlib import pyplot as plt
import numpy as np
import pyfits as fits
import lensing as lensing
from emissionLine import *

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
    #initiate the class
    emissionLineCl = emissionLine('NARROW_HB', nRedshiftBins=2)

    #get the lensing probability from convolvelssPc.
    #currently only for z=1., and alpha=0.83
    #future models will sample these two paraemeters to find
    #the best fitting for the observed tomorgraphic sampels
    emissionLineCl.getLensingProbability( z=1.0, alpha=0.83 )

    #get what i am calling the intrinsic distribution of HB narrow EQ
    #this can be discussed.the redshift cut determines thoise quasars
    #that havent been lensed, but still need enough quasars to get a good
    #distribution
    emissionLineCl.getIntrinsicDistribution( redshiftCut=0.3 )


    #now convolve the the two together to get the expected
    #distribtuion of EW at a redshift of 1.0 given an intrinsic
    #distribution of quasars
    emissionLineCl.convolveIntrinsicEquivalentWidthWithLensingProbability()

    #plot the initial distrubtion
    plt.plot(emissionLineCl.intrinsicEquivalentWidthDistribution['x'],\
                 emissionLineCl.intrinsicEquivalentWidthDistribution['y'], 'b', label='Intrinsic Distribution')
                 
    #plot the resulting convolution
    plt.plot(emissionLineCl.predictedLensedEquivalentWidthDistribution['x'],\
                 emissionLineCl.predictedLensedEquivalentWidthDistribution['y'],'r',label=r'Predicted Distribution ($\alpha=0.83$)')

    plt.legend(loc=2, prop={'size': 6})

    plt.xlim(-0.7,0.7)

    plt.xlabel(r'$\Delta$m')
    plt.savefig('Figure2Bellido.pdf')
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
    fig, axarr = plt.subplots(len(emissionLineList),1)

    nEquivilentWidthBins = 30
    maxEquivilentWidth = np.array([100., 50., 20.])
    for iAxis, iEmissionLine in enumerate(emissionLineList):
        emissionLineCl = emissionLine(iEmissionLine, nRedshiftBins=2)
        emissionLineCl.getEquvilentWidthMeansAsFunctionOfRedshift()
        EquivilentWidthBins = np.linspace(0., maxEquivilentWidth[iAxis], nEquivilentWidthBins)
        emissionLineCl.cumulativeDistributionEquivalentWidth()
        
        redshiftList = emissionLineCl.equivilentWidthCumSum.keys()
        redshiftListSorted = list((np.sort(np.array(redshiftList).astype(float))).astype(str))
        for i in redshiftListSorted:

            axarr[iAxis].plot(emissionLineCl.equivilentWidthCumSum[i]['x'], \
                emissionLineCl.equivilentWidthCumSum[i]['y'], \
                      label=i)
                      

        axarr[iAxis].legend()
        axarr[iAxis].set_xlim(0,maxEquivilentWidth[iAxis])
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
