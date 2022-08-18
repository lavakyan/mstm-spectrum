# -*- coding: utf-8 -*-
#-----------------------------------------------------#
#                                                     #
# This code is a part of T-matrix fitting project     #
# Contributors:                                       #
#  L. Avakyan <laavakyan@sfedu.ru>                    #
#  A. Skidanenko <ann.skidanenko@ya.ru>               #
#                                                     #
#-----------------------------------------------------#
"""
  Fitting of particle aggreagate T-matrix spectrum
  to experimental SPR spectrum.
"""
from __future__ import print_function
import os
from mstm_studio.mstm_spectrum import SPR, ExplicitSpheres, SpheresOverlapError
from mstm_studio.contributions import ConstantBackground
import numpy as np
from scipy import interpolate
import scipy.optimize as so
from datetime import datetime
import threading

try:
    import matplotlib.pyplot as plt
except:
    pass

# use input in both python2 and python3
try:
   input = raw_input
except NameError:
   pass
# use xrange in both python2 and python3
try:
    xrange
except NameError:
    xrange = range

class Parameter(object):
    """
    Class for parameter object used for storage of
    parameter's name, value and variation limits.

    Parameter naming conventions:

        `scale` - outer common multiplier

        `ext%i` - extra parameter, like background, peaks or Mie contributions

        `a%i` - sphere radius

        `x%i`, `y%i`, `z%i` - coordinates of sphere center

    where `%i` is a number (0, 1, 2, ...)
    """
    def __init__(self, name, value=1, min=None, max=None, internal_loop=False):
        """
        Parameters:

            name: string
                name of parameter used for constraints etc

            value: float
                initial value of parameter

            min, max: float
                bounds for parameter variation (optional)

            internal_loop : bool
                if `True` the parameter will be allowed to vary in internal
                (fast) loop, which does not require MSTM recalculation.
                Note: this flag will be removed in future.

            varied: bool
                if `True` -- will be changed during fit
        """
        self.name = name
        self.value = self.ini_value = value
        self.min = min
        self.max = max
        self.internal_loop = internal_loop
        self.varied = True
        # ...

    def __str__(self):
        return '%s : %f' % (self.name, self.value)


class Constraint(object):
    """
    Abstract constraint class. All other should inherit from it.
    """
    def apply(self, params):
        """
        Modify the params dict
        according to given constranint algorithm.

        Note: Abstract method!
        """
        pass


class FixConstraint(Constraint):
    def __init__(self, prm, value=None):
        """
        Fix value of parameter with name `prm` to `value`.

        Parameters:

            prm: string
                parameter name

            value: float
                if `None` than initial value will be used.
        """
        self.prm = prm.lower()
        self.value = value

    def apply(self, params):
        """ Apply fix constraint """
        assert self.prm in params
        if self.value is not None:
            params[self.prm].value = self.value
        params[self.prm].varied = False


class EqualityConstraint(Constraint):
    def __init__(self, prm1, prm2):
        """
        Fix two parameters with names `prm1` and `prm2` being equal
        """
        self.prm1 = prm1.lower()
        self.prm2 = prm2.lower()

    def apply(self, params):
        """ Apply equality constraint """
        assert self.prm1 in params
        assert self.prm2 in params
        params[self.prm2].value = params[self.prm1].value
        params[self.prm2].varied = False


class ConcentricConstraint(Constraint):
    def __init__(self, i1, i2):
        """
        Two spheres with common centers.

        `i1` and `i2` -- indexes of spheres
        """
        self.constraints = [EqualityConstraint('x%02i'%i1, 'x%02i'%i2),
                            EqualityConstraint('y%02i'%i1, 'y%02i'%i2),
                            EqualityConstraint('z%02i'%i1, 'z%02i'%i2)]

    def apply(self, params):
        """ Apply concentric constraint """
        for c in self.constraints:
            c.apply(params)


class RatioConstraint(Constraint):
    def __init__(self, prm1, prm2, ratio=1):
        """
        Maintain ratio of two variables, `prm1`/`prm2` = `ratio`
        """
        self.prm1 = prm1.lower()
        self.prm2 = prm2.lower()
        self.set_ratio(ratio)

    def apply(self, params):
        """ Apply Ratio constraint """
        assert self.prm1 in params
        assert self.prm2 in params
        params[self.prm2].value = params[self.prm1].value / self.ratio
        params[self.prm2].varied = False

    def set_ratio(self, ratio):
        """
        Set ratio of :math:`prm1/prm2 = ratio`.
        """
        assert np.abs(ratio) > 1e-10
        self.ratio = ratio


class Fitter(threading.Thread):
    """
    Class to perform fit of experimental Exctinction spectrum

    Field:

        tolerance: float
            stopping criterion, default is 1e-4
    """

    tolerance = 1e-4  # stopping criterion

    def __init__(self, exp_filename, wl_min=300, wl_max=800, wl_npoints=51,
                 extra_contributions=None, plot_progress=False):
        """
        Parameters:

            exp_filename: str
                name of file with experimental data

            wl_min, wl_max: float
                wavelength bounds for fitting (in nm).

            wl_npoints: int
                number of wavelengths where spectra will be calcualted and compared.

            extra_contributions: list of Contribution objects
                If `None`, then ConstantBackground will be used.
                Assuming that first element is a background.
                If you don't want any extra contribution, set to empty list `[]`.

            plot_progress: bool
                Show fitting progress using matplotlib.
                Should be turned off when run on parallel cluster without gui.
        """
        super(Fitter, self).__init__()
        self._stop_event = threading.Event()  # to be able to stop outside

        self.exp_filename = exp_filename
        data = np.loadtxt(self.exp_filename)    # load data
        data = data[np.argsort(data[:,0]),:]    # sort by 0th column
        if np.max(data[:,0]) < 10:      # if values are really low
            print('WARNING: Data X column is probably in mum, automatilcally rescaling to nm.')
            data[:,0] = data[:,0] * 1000

        self.wl_min = max(wl_min, data[ 0, 0])
        self.wl_max = min(wl_max, data[-1, 0])
        self.wl_npoints = wl_npoints
        print('Wavelength limits are setted to: %f < wl < %f'% (self.wl_min, self.wl_max))
        self.wls, self.exp = self._rebin(self.wl_min, self.wl_max, self.wl_npoints,
                                         data[:,0], data[:,1])
        self.params = {}             # dictionaty of parameters objects
        self.spheres = None          # object of Spheres
        self.constraints = []        # list of Constraint objects
        self.calc = np.zeros_like(self.wls)  # calculated spectrum
        self.chisq = -1              # chi square (squared residual)
        # set scale as default
        self.set_scale()
        # add extra contributions
        self.extra_contributions = []
        self.set_extra_contributions(extra_contributions)
        # set matrix material as default
        self.set_matrix()
        # plot, if specified
        self.plot_progress = plot_progress
        if self.plot_progress:
            plt.ion()
            self.fig = plt.figure()
            ax = self.fig.add_subplot(111)
            ax.plot(self.wls, self.exp, 'ro')
            self.line1, = ax.plot(self.wls, self.calc, 'b-')
            self.fig.canvas.draw()
            #~ self.lock = threading.Lock()  # used to sync with main thread where plot
        # callback function supplied outside
        self._cbuser = None

    def _rebin(self, xmin, xmax, N, x, y):
        """
        hidden method used to rebin data to uniform scale
        """
        f = interpolate.interp1d(x, y)
        xnew = np.linspace(xmin, xmax, N)
        ynew = f(xnew)
        return xnew, ynew

    def _print_params(self):
        for key in sorted(self.params):
            print('params[ %s ] \t %s ' % (key, self.params[key]))

    def set_matrix(self, material='AIR'):
        """
        set refraction index of matrix material

        material : {'AIR'|'WATER'|'GLASS'} or float
            the name of material or
            refraction index value.
        """
        self.MATRIX_MATERIAL = material

    def set_scale(self, value=1):
        if 'scale' not in self.params:
            self.params['scale'] = Parameter('scale', value=value, internal_loop=True)
        else:
            self.params['scale'].value = value
            self.params['scale'].ini_value = value

    def set_extra_contributions(self, contributions, initial_values=None):
        """
        Add extra contributions and initialize corresponding params.

        Parameters:

            contributions: list of Contribution objests

            initial_values: float array
        """
        # remove old parameters
        i_tot = 0
        for contribution in self.extra_contributions:
            if contribution is not None:
                n = contribution.number_of_params
                for _ in range(n):
                    self.params.pop('ext%02i' % i_tot)
                    i_tot += 1

        if contributions is None:
            contributions = [ConstantBackground(self.wls, 'ConstBkg')]
        self.extra_contributions = contributions[:]

        # create new parameter objects
        n_tot = 0
        i_tot = 0
        for contribution in self.extra_contributions:
            n = contribution.number_of_params
            for _ in range(n):
                self.params['ext%02i' % i_tot] = Parameter('ext%02i' % i_tot, value=0.1, internal_loop=True)
                i_tot += 1
            n_tot += n
        self.extra_contrib_params_count = n_tot
        if initial_values is not None:
            assert len(initial_values) == n_tot
            for i in range(n_tot):
                if initial_values[i] is not None:
                    self.params['ext%02i' % i].value = initial_values[i]
                    self.params['ext%02i' % i].ini_value = initial_values[i]
        # print(self.params)

    def set_spheres(self, spheres):
        """
        Specify the spheres to be fit.

        Paramerer:

            spheres: list of mstm_spectrum.Sphere objects
                If `None` then MSTM will not be run.
        """
        if self.spheres is not None:  # remove parameters of old spheres
            for i in xrange(self.spheres.N):
                self.params.pop('a%02i' % i)
                self.params.pop('x%02i' % i)
                self.params.pop('y%02i' % i)
                self.params.pop('z%02i' % i)
        if spheres is not None:
            self.spheres = spheres
            for i in xrange(self.spheres.N):
                self.params['a%02i' % i] = Parameter('a%02i' % i, self.spheres.a[i])
                self.params['x%02i' % i] = Parameter('x%02i' % i, self.spheres.x[i])
                self.params['y%02i' % i] = Parameter('y%02i' % i, self.spheres.y[i])
                self.params['z%02i' % i] = Parameter('z%02i' % i, self.spheres.z[i])
        else:
            self.set_spheres(ExplicitSpheres())  # empty spheres object

    def _update_spheres(self):
        """
        Set spheres radii and positions to values from params dict
        """
        assert self.spheres is not None
        for i in xrange(len(self.spheres)):
            self.spheres.a[i] = self.params['a%02i' % i].value
            self.spheres.x[i] = self.params['x%02i' % i].value
            self.spheres.y[i] = self.params['y%02i' % i].value
            self.spheres.z[i] = self.params['z%02i' % i].value

    def _update_params(self, values, internal=False):
        """
        Put values from optimized to params

        internal : bool
            if True than internal variables will be updated (scale, bkg, ..)
        """
        try:  # if not iterable (single value in values)
            len(values)
        except:
            print('WARNING: values is not a list')
            values = [values]
        if internal:  # internal (fast) loop parameters
            self.params['scale'].value = values[0]
            for i in range(self.extra_contrib_params_count):
                self.params['ext%02i' % i].value = values[i+1]  # 0th is scale
            print('inner: ', (self.extra_contrib_params_count+1), values)
        else:
            # apply constraints, -- up to now works only for MSTM
            for c in self.constraints:  # apply before
                c.apply(self.params)
            # update params
            i_tot = 0
            for i in range(len(self.spheres)):
                for key in ('a%02i'%i,'x%02i'%i,'y%02i'%i,'z%02i'%i):
                    if self.params[key].varied:
                        self.params[key].value = values[i_tot]
                        i_tot += 1
            assert i_tot == len(values)
            for c in self.constraints:  # and apply after
                c.apply(self.params)
            self.report_result(msg='[%s] Scale: %.3f Bkg: %.2f\n' % (str(datetime.now()),
                self.params['scale'].value, self.params['ext00'].value))  # may be verbous!

    def add_constraint(self, cs):
        """
        Adds constraints on the parameters.
        Usefull for the case of core-shell and layered structures.

        Parameter:

            cs: Contraint object or list of Contraint objects
        """
        try:
            _ = iter(cs)
        except TypeError:
            cs = [cs]
        for c in cs:
            self.constraints.append(c)

    def _get_spectrum(self):
        """
        Calculate the spectrum of agglomerates using mstm_spectrum module.
        """
        if self.stopped():
            raise Exception('Fitting interrupted')

        #~ self._apply_constraints
        if len(self.spheres) > 0:
            spr = SPR(self.wls)
            spr.environment_material = self.MATRIX_MATERIAL

            self._update_spheres()
            spr.set_spheres(self.spheres)
            #~ self.lock.acquire()
            try:
                _, extinction = spr.simulate()
                self.result = np.array(extinction)
            except SpheresOverlapError as e:
                self.chisq = 666  # big evil value
                return np.zeros_like(self.wls)
            except Exception as e:
                print(e)  # let User decide
                raise e
            #~ finally:
                #~ self.lock.release()
        else:  # emty spheres list
            self.result = np.zeros_like(self.wls)

        # perform fast fit over internal variables (scale, bkg, ..)
        values_internal = []
        values_internal.append(self.params['scale'].value)
        for i in range(self.extra_contrib_params_count):
            values_internal.append(self.params['ext%02i' % i].value)

        def _target_func_int(values):
            """ target function for internal fit (fast loop) """
            self._update_params(values, internal=True)

            y_dat = self.exp
            assert self.params['scale'].value == values[0]
            y_fit = values[0] * self.result
            n_tot = 1  # scale is values[0]
            for contribution in self.extra_contributions:
                n = contribution.number_of_params
                y_fit += contribution.calculate(values[n_tot:n_tot+n])
                n_tot += n
            self.chisq = np.sum((y_fit - y_dat)**2)
            #~ self.chisq = np.sum((y_fit - y_dat)**2 * (y_dat/np.max(y_dat)+0.001)) / np.sum((y_dat/np.max(y_dat)+0.001))
            #~ self.chisq = np.sum( (y_fit - y_dat)**2 * y_dat**3 ) * 1E3
            print(self.chisq)
            return self.chisq

        #~ print('/ Internal fit loop /')
        result_int = so.minimize(fun=_target_func_int, x0=values_internal, method='BFGS', tol=self.tolerance,
                                 options={'maxiter':100, 'disp':False})
        values_internal = result_int.x

        self._update_params(values_internal, internal=True)

        self.calc = self.params['scale'].value * self.result
        n_tot = 1  # 0th is scale
        for contribution in self.extra_contributions:
            n = contribution.number_of_params
            self.calc += contribution.calculate(values_internal[n_tot:n_tot+n])
            n_tot += n
        return self.calc

    def get_extra_contributions(self):
        '''
        Return a list of current extra contributions to the spectrum
        '''
        result = []
        values_internal = []
        values_internal.append(self.params['scale'].value)
        for i in range(self.extra_contrib_params_count):
            values_internal.append(self.params['ext%02i' % i].value)
        n_tot = 1  # scale is values[0]
        for contribution in self.extra_contributions:
            n = contribution.number_of_params
            result.append(contribution.calculate(values_internal[n_tot:n_tot+n]))
            n_tot += n
        return result

    def _target_func(self, values):
        """ main target function """
        self._update_params(values)

        y_dat = self.exp
        y_fit = self._get_spectrum()
        self.chisq = np.sum((y_fit - y_dat)**2)
        #~ self.chisq = np.sum((y_fit - y_dat)**2 * (y_dat/np.max(y_dat)+0.001)) / np.sum((y_dat/np.max(y_dat)+0.001))
        #~ self.chisq = np.sum( (y_fit - y_dat)**2 * (y_dat + np.max(y_dat*0.001))**3)
        #~ self.chisq = np.sum((y_fit - y_dat)**2 * y_dat**3) * 1E3
        #print(chisq)
        return self.chisq

    def set_callback(self, func):
        """
        Set callback function which will be called on
        each step of outer optimization loop.

        Parameter:

            func: function(values)
                where values -- list of values passed from optimization routine
        """
        self._cbuser = func

    def _cbplot(self, values):
        """
        callback function
        """
        #~ self.lock.acquire()  # will wait here
        #~ try:
        #print('Scale: %0.3f Bkg: %0.3f ChiSq: %.8f'% (self.params['scale'].value,
        #      self.params['bkg0'].value, self.chisq) )
        if self._cbuser is not None:  # call user-supplied function
            self._cbuser(self, values)
        if self.plot_progress:
            self.line1.set_ydata(self.calc)
            self.fig.canvas.draw()
            #~ plt.pause(0.05)  # this lead of grabbing of the focus by the plot window
            self.fig.canvas.start_event_loop(0.05)
            #from:
            #https://stackoverflow.com/questions/45729092/make-interactive-matplotlib-window-not-pop-to-front-on-each-update-windows-7
            #~ backend = plt.rcParams['backend']
            #~ import matplotlib
            #~ if backend in matplotlib.rcsetup.interactive_bk:
                #~ figManager = matplotlib._pylab_helpers.Gcf.get_active()
                #~ if figManager is not None:
                    #~ canvas = figManager.canvas
                    #~ if canvas.figure.stale:
                        #~ canvas.draw()
                    #~ canvas.start_event_loop(0.05)
        #~ finally:
            #~ self.lock.release()
        #~ input('pe')

    def _apply_constraints(self):
        for c in self.constraints:
            c.apply(self.params)

    def run(self, maxsteps=400):
        """
        Start fitting.

        Parameters:
            maxsteps: int
                limits number of steps performed
        """
        self._apply_constraints()
        # pack parameters to values
        values = []
        for i in range(len(self.spheres)):
            for key in ('a%02i'%i,'x%02i'%i,'y%02i'%i,'z%02i'%i):
                if self.params[key].varied:
                    values.append(self.params[key].value)
        # run optimizer
        result = so.minimize(fun=self._target_func, x0=values, method='Powell', tol=self.tolerance,
                             options={'maxiter':maxsteps, 'disp':True}, callback=self._cbplot)
        self._update_params(result.x)

    def stop(self):
        # https://stackoverflow.com/questions/323972/is-there-any-way-to-kill-a-thread-in-python
        self._stop_event.set()

    def stopped(self):
        return self._stop_event.is_set()

    def report_freedom(self):
        """
        Returns string with short summary before fitting
        """
        self._apply_constraints()
        N = len(self.spheres)
        s = 'Number of spheres:\t%i\n' % N
        n_tot = 0
        for contribution in self.extra_contributions:
            n = contribution.number_of_params
            s += 'Extra contrib. %s with %i params\n' % (contribution.name, n)
            n_tot += n
        s += 'Total number of extra params:\t%i\n' % n_tot
        n_fast = n_slow = n_fix = 0
        for key in self.params:
            if self.params[key].varied:
                if self.params[key].internal_loop:
                    n_fast += 1
                else:
                    n_slow += 1
            else: # not varied
                n_fix += 1
        assert n_fast + n_slow + n_fix == 1 + 4*N + n_tot
        s += 'Degrees of freedom\n'
        s += '\tfast loop:\t%i\n' % n_fast
        s += '\tslow loop:\t%i\n' % n_slow
        print(s)
        return s

    def report_result(self, msg=None):
        """
        Returns string with short summary of fitting results
        """
        s = 'ChiSq:\t%f\n' % self.chisq
        if msg is None:
            s += 'Optimal parameters'
        else:
            s += msg
        for key in sorted(self.params):
            s += '\n\t%s:\t%f\t(Varied:%s)' % (key, self.params[key].value, str(self.params[key].varied))
        print(s)
        return s

if __name__ == '__main__':
    fitter = Fitter('../example/experiment.dat')
    # test Mie fit
    from mstm_studio.contributions import LinearBackground, MieSingleSphere, MieLognormSpheresCached
    from mstm_studio.alloy_AuAg import AlloyAuAg
    fitter.set_extra_contributions([LinearBackground(fitter.wls, 'lin bkg'),
                                    MieLognormSpheresCached(fitter.wls, 'LN Mie')],
                                    [0.02, -0.001,
                                    0.1, 1.5, 0.5])
                                    #~ MieSingleSphere(fitter.wls, 'Mie')],
                                    #~ [0.02, -0.001,
                                    #~ 0.1, 10])
    fitter.extra_contributions[1].set_material(AlloyAuAg(1.), 1.66)
    fitter.extra_contributions[1].plot([0.1, 1.5, 0.5])
    fitter.extra_contributions[1].plot_distrib([0.1, 1.5, 0.5])
    #~ fitter.extra_contributions[1].plot([0.1, 10])
    fitter.set_spheres(None)  # no spheres, no mstm runs
    fitter.report_freedom()
    input('Press enter to run peak fitting')
    fitter.run()
    fitter.report_result()
    contribs = fitter.get_extra_contributions()
    print(contribs)
    input('Press enter to continue')
    # test peak fit
    from contributions import LinearBackground, LorentzPeak
    fitter.set_extra_contributions([LinearBackground(fitter.wls, 'lin bkg'),
                                    LorentzPeak(fitter.wls, 'lorentz peak')],
                                    [0.02, -0.001,
                                     100, 550, 50])
    # fitter.extra_contributions[1].plot([100, 550, 50])
    fitter.set_spheres(None)  # no spheres, no mstm runs
    fitter.report_freedom()
    input('Press enter to run peak fitting')
    fitter.run()
    fitter.report_result()
    input('Press enter to continue')
    # test MSTM fit
    fitter.set_matrix('glass')
    fitter.set_extra_contributions([LinearBackground(fitter.wls, 'lin bkg')], [0.02, -0.001])
    #                         N    X      Y      Z    radius    materials
    spheres = ExplicitSpheres(2, [-1,1], [-2,2], [-3,3], [14,20], [AlloyAuAg(1.), AlloyAuAg(0.)])
    fitter.set_spheres(spheres)
    fitter.add_constraint(ConcentricConstraint(0, 1))
    fitter.add_constraint(FixConstraint('x00', 0))
    fitter.add_constraint(FixConstraint('y00', 0))
    fitter.add_constraint(FixConstraint('z00', 0))
    fitter.add_constraint(RatioConstraint('a00', 'a01', spheres.a[0]/spheres.a[1]))
    fitter.report_freedom()
    input('Press enter to run MSTM fitting')

    fitter.run()
    #~ fitter.start()  # thread method
    #~ fitter.join()   # wait till end
    fitter.report_result()
    #fitter.plot_result()
    #~ y_fit = _get_spectrum( wavelengths, values )
    #~ plt.plot( wavelengths, exp, wavelengths, y_fit )
    #~ #plt.axis([0, 1, 1.1*np.amin(s), 2*np.amax(s)])
    #~ plt.xlabel('Wavelength, nm')
    #~ plt.ylabel('Exctinction, a.u.')
    #~ plt.show()
    input('Press enter to finish')
    print('It is over.')
