'''
Calculate the probaility of a point source being magnified

'''
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import colors as colors
import matplotlib.cm as cmx
from matplotlib import rc
from cosmolopy import distance as dist
import deltaTable as dt    

def getMeanMagOther(z, omegaM=0.31):

    '''
    From http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1993ApJ...404..436P&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf
    

    mean mag = equation 27, or 28
    '''
    beta = np.sqrt(1+24.*omegaM)
    
    tau = 3./beta*omegaM*( ( (1.+z)**(beta/2.)+1)/( (1.+z)**(beta/2.)-1) * np.log(1+z) - 4./beta)


    eqn26 =np.exp(2.*tau)

    eqn28 = beta**(-2)*(np.sinh(0.25*beta*np.log(1.+z)))**2*\
      (np.sinh(0.25*np.log(1.+z)))**(-2)
      
    return eqn28 - 1

def getInverseAnalPDF( z):

    delta = dt.getDelta(z)
    mag, pdf = magInvPDF( delta )

    return mag, pdf

def getInversePDF( z ):

    mag, pdf = getPDF(z)

    pdf /= np.sum(pdf[1:])*(mag[1]-mag[0])
    
    InvMag = 1./(1.+mag)
    inverse = pdf/InvMag**2

    dMag =   np.roll(InvMag,1) - InvMag 

    return InvMag[1:], inverse[1:], dMag[1:]

def getDelta(z, meanMag=None, omegaM=0.30):
    '''
    get the delta of the pdf in the normal version
    '''
    deltaList = np.linspace(0.,2,1000)
    if meanMag is None:
        meanMag = getMeanMag( z, omegaM=omegaM )

    expecMag = []
    
    for delta in deltaList:
        
        mag, iPDF = magPDF( delta )
        dMag = mag[1]-mag[-0]
        expecMag.append(np.sum(iPDF*mag*dMag))



    indexClosest = np.argmin( np.abs(np.array(expecMag) - meanMag))

    closestDelta = deltaList[indexClosest]
    closestExpec = np.array(expecMag)[indexClosest]
    diff = ((closestExpec -  meanMag)/(1.+meanMag))*100.
    if diff > 10.:
        print("WARNING DIFFERENCE IN EXPEC MAG AND MEAN MAG TOO HIGH at %0.2f" % diff)
        #pdb.set_trace()
    return closestDelta
        
def getPDF( z ):
    delta = dt.getDelta(z)
    mag, pdf =  magPDF( delta )
    return mag, pdf



def magPDF( delta, mag=None ):
    '''
    Using equation A3 from Z&S
    or eqyuation 8 from Rauch 1990

    '''
    if mag is None:
        mag = np.linspace(0.0, 10000.0, 1000000)[1:]


    dMag = mag[1]-mag[0]
        
    topFrac = 1. - np.exp( -mag/delta)
    botFrac = (mag+1)**2 - 1.
        
    P = ( topFrac/botFrac)**(3./2.)
        
    P /= np.sum(P)*dMag
    
    return mag, P

def magPDFarray( delta, mag=None ):
    '''
    Using equation A3 from Z&S
    or eqyuation 8 from Rauch 1990

    This generates an array for the PDF of PBH that for a given 
    delta (which then corresponds to a given mean magnitude)

    It is done in an array such that it generates all the PDFs
    for all possible deltas (delta is an array) at the same time

    so it outputs an NXM dimension array where 
    N is the number of elements in a signle PDF and M is the 
    number of deltas (or meanMags)

    '''
    magVec = np.linspace(0.0, 100.0, 10000)[1:]
    
            
    mag = np.matrix(magVec).T*np.ones(len(delta))


    dMag = mag[1,0]-mag[0,0]
        
    topFrac = np.array(1. - np.exp( -mag/delta))
    
    botFrac = (np.array(mag)+1)**2 - 1.
        
    P = ( topFrac/botFrac)**(3./2.)
        
    P /= np.sum(P,axis=0)*dMag

    return magVec, np.array(P)

    
    
def multiRedshiftPDF():

    redshiftList = np.linspace(0.01,2.,10.)
    jet = cm = plt.get_cmap('rainbow')
    cNorm  = colors.Normalize(vmin=0, vmax=2.)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    fig, [ax1, ax2] = plt.subplots( 2, 1 )
    plt.subplots_adjust(hspace=0.25)
    meanMagTot = []
    InvFLRWmeanMag = []
    for iRedshift in redshiftList:

        mag, PDF, dMag = \
          getInversePDF( iRedshift )

        iInvFLRWmeanMag = 1./(1.+getMeanMag(iRedshift))
        InvFLRWmeanMag.append( iInvFLRWmeanMag)
                               
        meanMag =  np.sum(dMag*mag*PDF)
        meanMagTot.append(meanMag)
        color = scalarMap.to_rgba(iRedshift)

        diff = mag - iInvFLRWmeanMag
        diff = np.append(1.-iInvFLRWmeanMag, diff)
        PDF = np.append(1e-10,PDF)
        ax2.plot( diff,PDF,color=color, label='z=%0.2f' % iRedshift)
        ax2.plot([meanMag- iInvFLRWmeanMag,meanMag- iInvFLRWmeanMag],[1e-5,100],'--', \
                     color=color)

   
    
    ax1.plot( redshiftList, np.array(meanMagTot) - np.array(InvFLRWmeanMag))


    ax1.set_xlabel(r'z')
    ax1.set_ylabel(r'$\langle (1+\mu)^{-1}\rangle$ - $(1+\bar{\mu})^{-1}$')
    
    ax2.set_yscale('log')
    ax2.set_xscale('linear')
    ax2.set_ylim(1e-5,100)
    ax2.set_xlabel(r'$\Delta(1.+\mu)^{-1}$')
    ax2.set_ylabel(r'$P_C$')
   
    plt.savefig('PDFinvesePBH.pdf')

    
def plotPBHpdf():
    redshiftList = np.linspace(0.01,2.,10.)
    jet = cm = plt.get_cmap('rainbow')
    cNorm  = colors.Normalize(vmin=0, vmax=2.)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    redshiftList = np.linspace(0.1,2.1,5)
    ax = plt.gca()

    for z in redshiftList:
        mag, pdf = getPDF(z)
        pdf = np.append(0., pdf)
        mag = np.append(0., mag)

        color = scalarMap.to_rgba(z)
        pdf/=np.max(pdf)
        ax.plot(mag-getMeanMag(z), pdf, color=color,label="z=%0.2f" % z, lw=2)

    ax.legend()
    ax.set_ylim(0.,1.1)
    ax.set_ylabel(r'Normalised Probability P($\Delta\mu_{\rmPBH}$)', fontsize=15)
    ax.set_xlabel(r'$\Delta\mu_{\rmPBH}$', fontsize=15)
    ax.axes.get_yaxis().set_ticks([])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.savefig('probabilityPBH.pdf')
    plt.show()



def getMeanMag( z, dz=1e-4, omega_m=0.30, sigma_8=0.9):
    #return 0.1311
    cosmo = {'omega_M_0' : omega_m, \
            'omega_lambda_0' : 0.6914, \
            'h' : 0.6777, \
            'sigma_8':sigma_8}
    cosmo = dist.set_omega_k_0(cosmo)
    distanceEB = 0
    distanceFB = 0
    for i in np.arange(0.,z,dz):
        dist_hz = dist.hubble_distance_z(i, **cosmo)
        
        distanceEB += dz*dist_hz/(1+i)**2

        distanceFB += dz*dist_hz

    distanceFB /= 1.+z

    return  (distanceEB/distanceFB)**2-1.
        
