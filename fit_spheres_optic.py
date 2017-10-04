#!/usr/bin/python
# -*- coding: utf-8 -*-
#-----------------------------------------------------#
#                                                     #
# This code is a part of T-matrix fitting project     #
# Contributors:                                       #
#  L. Avakyan <laavakyan@sfedu.ru>                    #
#                                                     #
#-----------------------------------------------------#
"""
  Fitting of particle aggreagate T-matrix spectrum
  to experimental SPR spectrum.
"""
#TODO: rewrite in object-oriented way

import os
from mstm_spectrum import SPR, ExplicitSpheres, Background, FilmBackground
import numpy as np
from scipy import interpolate
import scipy.optimize as so

import matplotlib.pyplot as plt


MAX_FILESIZE = 1024*1024*1024  # 1 Gb limit

MATRIX_MATERIAL = 'AIR'

def getfloats(txt):
    words = [w.strip() for w in txt.replace(',', ' ').split()]
    try:
        return [float(w) for w in words]
    except:
        return None

def read_ascii(fname, sort=False, sort_column=0):
    """ read a column ascii column file, returning a group containing the data from the file.
        (original code: M.Newville)
    Example:
    >>> read_ascii(filename, sort=False, sort_column=0)
    """
    if not os.path.isfile(fname):
        raise OSError("File not found: '%s'" % fname)
    if os.stat(fname).st_size > MAX_FILESIZE:
        raise OSError("File '%s' too big for read_ascii()" % fname)

    finp = open(fname, 'r')
    text = finp.readlines()
    finp.close()

    data = []
    ncol = None

    text.reverse()
    section = 'FOOTER'

    for line in text:
        line = line[:-1].strip()
        if len(line) < 1:
            continue
        # look for section transitions (going from bottom to top)
        if section == 'FOOTER' and getfloats(line) is not None:
            section = 'DATA'
        elif section == 'DATA' and getfloats(line) is None:
            section = 'HEADER'
        # Leon: ignoring footer and header...
        if section == 'DATA':
            rowdat  = getfloats(line)
            if ncol is None:
                ncol = len(rowdat)
            if ncol == len(rowdat):
                data.append(rowdat)
    # reverse header, footer, data, convert to arrays
    data.reverse()
    data = np.array(data).transpose()

    ncols, nrow = data.shape

    if sort and sort_column >= 0 and sort_column < nrow:
         data = data[:,np.argsort(data[sort_column])]

    return data

def rebin(xmin, xmax, N, x, y):
    """
    returns rebinned data
    """
    f = interpolate.interp1d(x, y)
    xnew = np.linspace(xmin, xmax, N)
    ynew = f(xnew)
    return xnew, ynew


calculated_extinction = 0 # for plotting

wavelengths = []  # some globals for fitting
exp = []

chisq = -1

background = Background([])

def get_spectrum(wavelengths, values, n_medium=1.66):
    """
    get spectrum of T-matrix spheres aglomerates
    """
    global calculated_extinction
    #print('Current scale: %f bkg: %f' % (values[0], values[1]))
    spr = SPR( wavelengths )
    spr.environment_material = MATRIX_MATERIAL
    result = np.zeros(len(wavelengths))
    try:
        N = ( len(values)-1-background.number_of_params() )/4
        if N > 0:
            #print 'N=', N
            spr.set_spheres(  ExplicitSpheres( N, np.array(values[-N*4:]), [], [], [], mat_filename='etaGold.txt')  )
            _, extinction = spr.simulate('exct.dat')
            result = np.array(extinction)
        else:
            result = np.zeros(len(wavelengths))
    except Exception as e:
        print e
        result = np.zeros(len(wavelengths))
    calculated_extinction = values[0]*result + background.get_bkg(values[1:]) # scale * extinction + bkg
    return calculated_extinction

def target_func(values, x_dat, y_dat):
    print( 'Scale: %f Bkg: %f'% (values[0], values[1]) )
    y_fit =  get_spectrum( wavelengths, values, MATRIX_MATERIAL )
    global chisq
    #chisq = np.sum( (y_fit**3 - y_dat**3)**2 / y_dat**3 )
    chisq = np.sum( (y_fit - y_dat)**2 * y_dat**3 ) * 1E3
    #print( chisq )
    return chisq

def prepare_fit(wavelengths_, exp_, background_object=''):
    global wavelengths # = wavelengths
    wavelengths = wavelengths_ # this is ugly, mabbe make a class?
    global exp
    exp = exp_
    plt.ion()
    global fig
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(wavelengths, exp, 'ro')
    global line1
    line1, = ax.plot(wavelengths, exp, 'b-')
    fig.canvas.draw()
    global background
    if isinstance(background_object, Background):
        background = background_object
    else:
        background = Background(wavelengths)


def cbplot( values ):
    print( 'Scale: %0.3f Bkg: %0.3f ChiSq: %.8f'% (values[0], values[1], chisq) )
    line1.set_ydata(calculated_extinction)
    fig.canvas.draw()


if __name__ == '__main__':
    data = read_ascii('optic_sample22.dat', True, 0) # read and sort by 0th column

    if max(data[0,:]) < 10:
        print('Data X column is in microns, will rescale to nm.')
        data[0,:] = data[0,:] * 1000; # to nanometers

    #wavelengths, exp = rebin(400, 800, 51, data[0,:], data[1,:])
    wavelengths, exp = rebin(300, 800, 51, data[0,:], data[1,:])

    # prepare plot #
    prepare_fit(wavelengths, exp)

    ### SET INITIAL VALUES ###
    values = [0.02, 0.01] # scale, bkg
    A = 200 # 400
    a = 20
    d = 100
    x = -(A/2.0)
    while x < (A/2.0):
        y = -(A/2.0)
        while y < (A/2.0):
            z = -(A/2.0)
            while z < (A/2.0):
                if (x*x+y*y+z*z < A*A/4.0):
                    values.append(x)
                    values.append(y)
                    values.append(z)
                    values.append(a)
                    #print x, y, z
                z = z + (2*a+d)
            y = y + (2*a+d)
        x = x + (2*a+d)
    N = (len(values)-2)/4
    print 'Number of spheres: ', N
    print 'Number of degrees of freedom: ', len(values)
    MATRIX_MATERIAL = 'Glass'
    print 'Matrix material : ', MATRIX_MATERIAL
    raw_input('Press enter')

    ### OPTIMIZE (FIT) VALUES ###
    result = so.fmin( func=target_func, x0=values, callback=cbplot, xtol=0.0001, ftol=0.001, maxiter=1000, full_output=True, disp=True, args=(wavelengths, exp) )

    ### DEAL WITH RESULTS ###
    print(result)
    values = result[0]
    print values

    #~ y_fit = get_spectrum( wavelengths, values )
    #~ plt.plot( wavelengths, exp, wavelengths, y_fit )
    #~ #plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])
    #~ plt.xlabel('Wavelength, nm')
    #~ plt.ylabel('Exctinction, a.u.')
    #~ plt.show()
    raw_input('press enter')
    print 'It is over.'
