'''
I did this already but apparently got delted

I need to create a function that predicsts the normalisation and
delta in the compact object probabilty funciton

To do this i will load up the delta table and see if i can make some
predictions

so given a mean magnitdue what is the delta and normlsiation

'''

from scipy.optimize import curve_fit
from deltaTable import *
from matplotlib import pyplot as plt
from scipy.optimize import leastsq

def singleRegion( x, p1, p2, p3, p4, p5, p6, p7, p8, p9):
    return \
      p1 + p2*x + p3*x**2 +p4*x**3 + \
      p5*x**4 +p6*x**5 +p7*x**6 +p8*x**7 + \
      p9*x**8# +p10*x**9 +p11*x**10 +p12*x**11


def magnitude2param( x, y, split, returnParams=False):

    if (split > 10.) | (split < 1.0):
        return np.zeros(len(y))+10000

    xLow = x[ x < split]
    yLow = y[ x < split]
    
    xHigh = x[ x>=split]
    yHigh = y[ x>=split]
    
    paramsLow, cov = curve_fit( singleRegion, \
                    np.log10(xLow), np.log10(yLow), maxfev=100000)
                    
    paramsHi, cov = curve_fit( singleRegion, \
                    np.log10(xHigh), np.log10(yHigh), maxfev=100000)
    
    yModel = np.zeros(len(x))
    yModel[  x < split] = 10**singleRegion( np.log10(xLow), *paramsLow)
    yModel[  x >= split]= 10**singleRegion( np.log10(xHigh), *paramsHi)

    if returnParams:
        return paramsLow, paramsHi
    else:
        return yModel

def lossFunction( split, x, delta, norm ):
    '''
    My ansatz of how to go from mean magnitude to delta or the norm
    '''
    print split
    yDelta = magnitude2param( x, delta, split[0])
    yNorm = magnitude2param( x, norm, split[1])

    return (yDelta - delta)  + (yNorm-norm)
    


#These are the results from the predictNormAndDelta function
def getNormAndDeltaFromMeanMag( meanMag ):
    deltaCoefficients = np.array([0.4966286 , 1.95295421])
    normCoefficients = np.array([ 0.78158726, -1.0664505 ])

    delta = magnitude2param( meanMag, *deltaCoefficients)
    norm = magnitude2param( meanMag, *normCoefficients)

    return norm, delta



def getDeltaTable():
    '''
    Read in the pickle file for the delta table
    '''

    deltaList, \
      magnitudes, \
      norm = pkl.load(open('deltaTable.pkl','rb'))

    return pkl.load(open('deltaTable.pkl','rb'))

    
def predictNormAndDelta( ):
    '''
    Fit some function to the mean magnitudes
    and the delta
    '''

    
    deltaList, \
      magnitudes, \
      normList = getDeltaTable()
    print np.max(deltaList)
    '''
    Used for linear deltatable
    index = np.floor(10**np.linspace(0,np.log10(len(magnitudesAll)-1),5000)).astype(int)

    normList = normListAll[index]
    deltaList=deltaListAll[index]
    magnitudes = magnitudesAll[index]
    '''
    initGuess = [1.0,1.0]
    deltaArgs =  (magnitudes, deltaList, normList)
    #split, cov, infodict, mesg, ier = \
    #  leastsq( lossFunction, initGuess, args = deltaArgs, \
    #               full_output=True, epsfcn=0.001, xtol=0.000001)
    #split = [8.0, 8.0]
    #found if i run the fit
    split = [3.15491079, 2.54681736]
 
    fig, axarr = plt.subplots(2)
    axarr[0].plot( magnitudes, \
                magnitude2param(magnitudes, deltaList, split[0]),'--')
    axarr[0].plot( magnitudes, deltaList,'-')

    axarr[1].plot( magnitudes, \
                magnitude2param(magnitudes, normList, split[1]),'--')
    axarr[1].plot( magnitudes, normList,'-')
    axarr[1].plot( magnitudes, normList,'-')

    
    axarr[0].set_xscale('log')
    axarr[1].set_xscale('log')

    print("The splits is %0.2f for Delta and %0.2f for Norm" % tuple(split))

    paramsLoDelta, paramsHiDelta = \
      magnitude2param( magnitudes, deltaList, split[0], returnParams=True)
    print("The fitted params for delta low are", paramsLoDelta)
    print("The fitted params for delta high are", paramsLoDelta)
      
    paramsLoNorm, paramsHiNorm = \
      magnitude2param( magnitudes, normList, split[1], returnParams=True)
      
    print("The fitted params for norm low are", paramsLoNorm)
    print("The fitted params for norm high are", paramsHiNorm)

    pkl.dump( [split,paramsLoDelta, paramsHiDelta, \
                   paramsLoNorm, paramsHiNorm],\
                  open('getNormDeltaParams.pkl','wb'))
    plt.show()


def getNormAndDeltaForMagnitude( magnitudes ):

    split, paramsDeltalo, paramsDeltaHi, paramsNormLo, paramsNormHi = \
      pkl.load(open('getNormDeltaParams.pkl','rb'))

    norm = getParams( magnitudes, split[1], paramsNormLo, paramsNormHi)
    delta = getParams( magnitudes, split[0], paramsDeltalo, paramsDeltaHi)

    return norm, delta

def getParams( magnitude, split, paramsLo, paramsHi):
    
    indexLo = magnitude < split
    indexHi = magnitude >= split
    
    returnY = np.zeros(len(magnitude))
    
    returnY[ indexLo ] = \
      10**singleRegion( np.log10(magnitude[ indexLo ]), *paramsLo)
    returnY[ indexHi ] = \
      10**singleRegion( np.log10(magnitude[ indexHi ]), *paramsHi)
      

    return returnY


def testGetNormAndDeltaForMagnitude():

    deltaList, \
      magnitudes, \
      normList = getDeltaTable()

    normPredict, deltaPredict = getNormAndDeltaForMagnitude(magnitudes)

    fig, axarr= plt.subplots(2)

    axarr[0].plot( magnitudes, deltaList / deltaPredict)
    axarr[1].plot( magnitudes, normList / normPredict)
    axarr[0].set_xscale('log')
    axarr[1].set_xscale('log')
    plt.show()
if __name__ == '__main__':
     predictNormAndDelta( )
    

