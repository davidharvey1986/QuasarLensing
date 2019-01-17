from matplotlib import pyplot as plt
import numpy as np
import pyfits as fits
import lensing as lensing



def OIIIwithMGII( catalogue='dr9_weq.fits' ):
    shen = fits.open(catalogue)[1].data
    ax = plt.gca()
    index = (shen['REW_MGII']>0) & \
      (shen['REW_MGII']<1e3) & \
      (shen['REW_OIII_5007']>0) & \
      (shen['REW_OIII_5007']<1e3)


    zBin, weqBin = \
      redshiftBinWeq( shen['Z_PIPE'][index], \
                          shen['REW_MGII'][index])
    ax.errorbar(zBin, weqBin[0,:], \
                        yerr=weqBin[1,:], label='MGII')
    zBin, weqBin = \
      redshiftBinWeq( shen['Z_PIPE'][index], \
                          shen['REW_OIII_5007'][index])

    ax.errorbar(zBin, weqBin[0,:], \
                    yerr=weqBin[1,:], label='OIII')
    ax.legend()
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

def redshiftBinWeq( x, y, nBins=5 ):
    '''
    bin up the equviliant width
    in to redshift
    '''

    #First remove all the zeros

    xCut = x[(y>0) & (y<1e2)]
    yCut = y[(y>0) & (y<1e2)]

    bins = np.linspace(0., 1., nBins+1)

    MeanWeq = np.zeros((2, nBins))
    for iBin in xrange(nBins):
        inBin = (xCut>bins[iBin]) & (xCut<bins[iBin+1])

        
        MeanWeq[0, iBin] = np.mean(yCut[inBin])

        MeanWeq[1, iBin] = np.std(yCut[inBin])/np.sqrt(len(yCut[inBin]))

    centres = (bins[1:] + bins[:-1])/2.

    MeanWeq /= MeanWeq[0,0]
    
    return centres, MeanWeq

        
        
    

    
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
