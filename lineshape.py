import numpy as np

### Defining the single-peak Gaussian function

def Single_Gaussian(params, x):
    x0          = params[0]
    intensity   = params[1]
    FWHM        = params[2]
    bkgnd       = params[3]
    return float(intensity)*np.exp(- 0.5*((x0-x)/(FWHM/2.355))**2) + bkgnd

def Single_Gaussian_integrand(x, x0, intensity, FWHM, bkgnd):
    return float(intensity)*np.exp(- 0.5*((x0-x)/(FWHM/2.355))**2) #+ bkgnd


### Defining the double-peak Gaussian function

def Double_Gaussian(params, x):
    x1          = params[0]
    intensity1  = params[1]
    FWHM1       = params[2]
    bkgnd       = params[3]
    x2          = params[4]
    intensity2  = params[5]
    FWHM2       = params[6]
    return float(intensity1)*np.exp(- 0.5*((x1-x)/(FWHM1/2.355))**2) + float(intensity2)*np.exp(- 0.5*((x2-x)/(FWHM2/2.355))**2) + bkgnd


### Defining the single-peak Crystalball function

def Crystalball(x_array, x0, sigma, alpha, n):
    x = (x_array-x0)/sigma*np.sign(alpha)

    def bigger(x):
        return np.exp(-0.5 * x*x)

    def smaller(x):
        alph = np.abs(alpha)
        b = n / np.abs(alpha) - np.abs(alpha)
        a = ((n / alph)**n) * np.exp(-0.5*alph*alph)
        return a / (b - x) ** n

    y = np.piecewise(x, x >= -np.abs(alpha), [bigger, smaller])

    return y

def Single_Crystalball(params, x_array):
    x0      = params[0]
    N       = params[1]
    sigma   = params[2]
    bkgnd   = params[3]
    alpha   = params[4]
    n       = params[5]
    y = N * Crystalball(x_array, x0, sigma, alpha, n) + bkgnd
    return y

def Single_Crystalball_integrand(x, x0, N, sigma, alpha, n, bkgnd):
    t = (x-x0)/sigma
    if (alpha < 0):
        t = -t
    if (t >= -abs(alpha)):
        y =  np.exp(-0.5*t*t)
    else:
        a =  ((n/abs(alpha))**n)*np.exp(-0.5*abs(alpha)*abs(alpha))
        b = n/abs(alpha) - abs(alpha)
        y = a/(b - t)**n
    return N*y


### Defining the double-peak Crystalball function

def Double_Crystalball(params, x_array):
    x1      = params[0]
    N1      = params[1]
    sigma1  = params[2]
    bkgnd   = params[3]
    x2      = params[4]
    N2      = params[5]
    sigma2  = params[6]
    alpha   = params[7]
    n       = params[8]

    y1 = N1 * Crystalball(x_array, x1, sigma1, alpha, n)
    y2 = N2 * Crystalball(x_array, x2, sigma2, alpha, n)
    y = y1 + y2 + bkgnd
    return y

