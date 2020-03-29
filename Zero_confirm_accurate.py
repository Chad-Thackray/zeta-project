import sympy as sym
import numpy as np
import pynverse
import time
import concurrent.futures
import math

###############################################################
##############                    #############################
##############     Error Terms    #############################
##############                    #############################
###############################################################

p = sym.symbols("p")
t = sym.symbols("t")
N = sym.symbols("N")

### Error term constants
psi0 = sym.cos(2*sym.pi *(p**2 -p -0.0625))/sym.cos(2*sym.pi*p)
psi0 = sym.series(psi0, p, x0 = 0.5, n = 20).removeO()




    
c1 = -1/(2**5 * 3 * sym.pi**2)* sym.diff( psi0, p, 3)
c2 = 1/(2**11 * 3**2 * sym.pi**4) * sym.diff(psi0,p,6) + 1/(2**6 * sym.pi**2 )* sym.diff(psi0, p, 2)

cut_off = 12

psi0 = sym.series(psi0, p, x0 = 0.5, n = cut_off).removeO()
c1 = sym.series(c1, p, x0 = 0.5, n = cut_off).removeO()
c2 = sym.series(c2, p, x0 = 0.5, n = cut_off).removeO()
    
##### Estimation of remainder term   
remfunc = ( (-1)**(N-1)*(t/(2*sym.pi))**(-1/4) ) * (psi0+ c1*(t/(2*sym.pi))**(-1/2) + c2*(t/(2*sym.pi))**(-2/2) ) 


###############################################################
##############                       ##########################
##############   Basic methods       ##########################
##############                       ##########################
###############################################################

def giveRange(n):
    """
    Given a value n, returns a range of t values such that floor(sqrt(t/2pi)) = n
    """
    return [  (n)**2 *2*np.pi ,  (n+1)**2 *2*np.pi ]


def tau(x):
    return np.ceil(np.sqrt(x/(2*np.pi)))


def theta(t):
    """
    Calculates the value of the Riemann-Siegel theta
    function at t by approximation using Stirling's
    series
    accurate for large values of t
    """
    
    return 0.5*t* np.log(t/(2*np.pi)) - (0.5*t) - (np.pi*0.125) + 1/(48*t) + 7/(5760*(t**3))


def gram(n): 
    """
    returns g_n, then nth gram point
    """
    invtheta = pynverse.inversefunc(theta, domain=10, open_domain=True, accuracy=2)

    n = np.array(n)
    return invtheta(n*np.pi)


def Zlist(x):
    """
    Given x, a list of t values, returns the z values. 
    """

    x = np.array(x)
    zlist = np.array([])

    for n in range(int(tau(x[0]))-1, int(tau(x[-1]))):
        ranges = giveRange(n)
        part =  x[ (x > ranges[0]) & (x < ranges[1]) ]

        zlist = np.append(zlist, ZBulk(part))

    return zlist
        

def ZBulk(x):
    """
    Uses first x terms of the Riemann Siegel formulas
    to calculate Z(t). Ensure all terms have the same N value
    (the number of terms to be added in the sum)
    """

    global remfunc

    if x.size == 0:
        return []
   
    integer = int(np.floor(np.sqrt(x[0]/(2*np.pi))))

    fractional = np.sqrt(x/(2*np.pi)) - integer 

    total = 0
    for n in range(1, integer+1):
        total += n**(-0.5) * np.cos( theta(x) - x*np.log(n))

    remainder = sym.lambdify([N,p,t], remfunc, "numpy")

    return 2*total + remainder(integer, fractional, x)


def law(x, extra = 0):
    """
    given a list of gram points, returns the gram points that fail
    Gram's law (-1)^n Z(g_n)>0
    """

    gramlist = np.array(gram(x))
    
    zlist = Zlist(gramlist + extra)

    signs = [ (-1)**n for n in x]
    errors = 0.77*(gramlist+ extra)**(-7/4)

    toSort = np.multiply(zlist, signs)
    toSort2 = toSort[ toSort < errors]

    difficulties = []
    for z in toSort2:
        difficulties.append(np.where(toSort == z)[0][0])
    
    difficulties = np.take(x, difficulties)

    return difficulties

def findh(z, step = 0.1):
    """
    Given a list of n's. Find corresponding h_n's for g_n. 
    """
    h = 0
    final = {}
    current = np.array(z)
    while current.size != 0:

        failed = law(current, h)
        
        passed = np.setdiff1d(current,failed)
        for x in passed:
            final[x] = h

        if failed.size == 0:
            break

        current = failed

        failed = law(current, -h)
        
        passed = np.setdiff1d(current,failed)
        for x in passed:
            final[x] = -h

        current = failed

        h += step

    final2 = []
    for x in z:
        final2.append(final[x])
        
    return final2
        
###############################################################
##############                    #############################
##############  Turing's Method   #############################
##############                    #############################
###############################################################

def findk(kmin, n, step = 0.1):
    """
    Finds a value of k such that h_{n+k} and h_{n-k} = 0
    """
    
    k=kmin
    while findh([n+k], step)[0] != 0 or findh([n-k], step)[0] != 0:
        k += 1
        if k > 30:
            return 0
    return k


def turingtest(k,n, step = 0.1):
    """
    Attempts to verify the turing bound at g_n
    that is, there are n+1 zeros of zeta below g_n. 
    """

    upper = np.sum(findh( [j for j in range(n+1, n+k)], step))
    lower = np.sum(findh( [n-j for j in range(1,k)], step) )
    
    if ( (1.91 + 0.114*np.log(gram(n+k)/(2*np.pi)) + upper)/(gram(n+k)-gram(n)) >= 1    ):
        return False

    if ( (1.91 + 0.114*np.log(gram(n)/(2*np.pi)) + lower)/(gram(n)-gram(n-k)) >= 1    ):
        return False

    return True
    
    

def turing_bound(n, step = 0.1):
    """
    Given a value of n, seeks to find whether
    Gram's law applies at the value g_n,
    i.e. N(g_n) = n-1
    """

    h_n = findh([n], step)[0]
    
    if h_n != 0:
        return "h_n != 0"

    k = findk(3,n, step)
    if k==0:
        return "no k"

    while turingtest(k,n, step) != True:
        k = findk(k+1,n, step)
        if k ==0:
            return "no k"
    
    return turingtest(k,n, step)

###############################################################
##############                    #############################
##############  Verification      #############################
##############                    #############################
###############################################################

def Bulklehmans(z, step = 0.1):
    """
    given a list of gram points that do no follow grams law, determine if there can be found suitable h_n's such that
    (-1)^n Z(g_n + h_n) > err and g_n-1 + h_n-1 < g_n + h_n < g_n+1 + h_n+1
    """
    
    gn = gram(z)
    gnhigh = gram(z+1)
    gnlow = gram(z-1)

    for x in range(10):
        print("go")
        step = step/2
        hn_guess = findh(z,step)
        
        if (  False not in (gn + hn_guess > gnlow + findh(z-1, step)) and False not in (gn + hn_guess < gnhigh + findh(z+1, step))   ):
            return True

    return False

               
def verifyBulk(x,y):
    """
    Verifies there are y-x zeroes of zeta between g_x and g_y on the critical line
    x must be at least 300 for this to be valid methodology as we only showed Turing's bound for t>168pi
    """
    step = 0.1

    if( turing_bound(x, step) == True and turing_bound(y, step)== True):
        print("Turing bound achieved, finding points of failure for Gram's law")
    else:
        if( turing_bound(x, step/10) == True and turing_bound(y, step/10)== True):
            print("Turing bound achieved, finding points of failure for Gram's law")
        else:
            print("Turing bound failed, try a smaller step")
            return 

    difficulties = law([z for z in range(x,y+1)])

    groupsize = 200
    groups = np.array_split(difficulties, np.rint(difficulties.size / groupsize))

    return groups


###############################################################
##############                    #############################
##############  Activation code   #############################
##############                    #############################
###############################################################


start_time = time.perf_counter()

if __name__ == '__main__':

    groups = verifyBulk(300,99999)
    groups2 = []

    for z in groups:
        groups2.append(list(z))

    print(len(groups2), " Lehman's groups to process")

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [ executor.submit(Bulklehmans, z) for z in groups]

        for f in concurrent.futures.as_completed(results):
            print( f.result() )
        

print (time.perf_counter() - start_time, "seconds")


