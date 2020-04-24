'''
plot the mean magnitification and mean kappa as a function
of the cosmological parameter

'''
from matplotlib import pyplot as plt
import numpy as np
import probabilityMagnification as pm
import ipdb as pdb
def main(nOmegaM=20, nSigma8=20, z=1.0):

    omegaMlist = np.linspace(0.1,1.0, nOmegaM)
    sigma8list = np.linspace(0.6,1.5, nSigma8)
    
    meanKappa = np.zeros(nOmegaM)

    norm = pm.getMeanMag( z )/2
    for iOmegaM in xrange(nOmegaM):
           
        meanKappa[iOmegaM] = \
              pm.getMeanMag( z, omega_m =omegaMlist[iOmegaM])/2 - norm

    omegaMLow = omegaMlist[ np.argmin(np.abs(meanKappa+1e-2)) ]
    omegaMHi = omegaMlist[ np.argmin(np.abs(meanKappa-1e-2)) ]

    plt.plot([omegaMLow, omegaMLow], [-0.1,-1e-2],'k--')
    plt.plot([omegaMHi, omegaMHi], [-0.1,1e-2],'k--')
    
    plt.plot([0., omegaMLow], [-1e-2,-1e-2],'k--')
    plt.plot([0.1, omegaMHi], [1e-2,1e-2],'k--')
    print('Omega M hi=%0.2f and Omega M low = %0.2f' % (omegaMHi,omegaMLow))
    plt.xlim(omegaMlist[0],omegaMlist[-1])
    plt.ylim(-0.02,0.04)
    pdb.set_trace()
    plt.plot(omegaMlist, meanKappa)
    plt.xlabel(r'$\Omega_M$')
    plt.ylabel(r'$\langle \kappa \rangle$ -$\langle \kappa(\Omega_M=0.3) \rangle$ ')
    plt.savefig('../plots/meanKappaFunctionCosmo.pdf')
    plt.show()

if __name__ == '__main__':
    main()
