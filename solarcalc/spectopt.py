# Copyright (c) 2018 Kevin Murray <foss@kdmurray.id.au>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as np
from scipy.interpolate import interp1d, UnivariateSpline
from scipy import optimize
from scipy.spatial import distance
import pandas as pd


class Spectrum(object):
    """Holds a spectral density curve, also integrate and interpolate"""

    def __init__(self, wavelengths, values):
        self.wavelengths = wavelengths
        self.values = values

    def interpolated(self, method='linear'):
        '''Create a function to interpolate spectrum to new wavelengths'''
        return interp1d(self.wavelengths, self.values, kind=method)

    def spline(self, k=1, s=0):
        '''Create a function to interpolate spectrum to new wavelengths with a spline'''
        return UnivariateSpline(self.wavelengths, self.values, k=k, s=0)

    def one_nm(self, method='linear'):
        '''Gives an interpolated version of self with 1nm bands'''
        wl = np.arange(min(self.wavelengths), max(self.wavelengths)+1, dtype=int)
        return Spectrum(wl, self.interpolated(method)(wl))

    def banded_integral(self, bands):
        '''Gives the integral within each band, where bands are (lhs, rhs)'''
        s = self.spline()
        integr = []
        for band in bands:
            integr.append(s.integral(*band))
        return integr


class Light(object):
    """
    Loads a CSV describing a light source, interpolates the output, and
    allows extraction to custom domains and optimisation to other spectra
    """

    def __init__(self, filename, interpolation='linear'):
        df = pd.read_csv(filename)
        self.channels = {}
        wavelengths = df.values[:, 0]
        for chan in df.columns[1:]:
            values = df[chan]
            self.channels[chan] = Spectrum(wavelengths, values).interpolated()

    def __len__(self):
        return len(self.channels)

    def fit_wavelengths(self, wavelengths=None, minwl=None, maxwl=None, step=1):
        '''Interpolate each channel to be defined over certian wavelengths'''
        if wavelengths is None and minwl is None and maxwl is None:
            raise ValueError("Either wavelengths or maxwl and minwl must be given")
        if wavelengths is None:
            wavelengths = np.arange(minwl, maxwl, step)
        return np.vstack([chanfunc(wavelengths)
                          for chanfunc in self.channels.values()]).T

    def light_output(self, channel_settings, wavelengths=None, minwl=None, maxwl=None, step=1):
        '''Gives the total output of the light unit over wavelengths given channel_settings'''
        outputs = self.fit_wavelengths(wavelengths, minwl, maxwl, step)
        return Spectrum(wavelengths, np.dot(outputs, channel_settings))


    def optimise_settings(desired, cost_function=SimpleCost()):
        '''Optimise channel settings aiming for desired, using cost_function'''
        initial = np.ones(len(self.channels))
        bounds = [(0.,1.)]*len(initial)
        opt = optimize.minimize(cost_function, initial, (self, desired),
                                options={"maxiter": 10000}, bounds=bounds)
        return opt.x


class CostFunc(object):
    '''Base class of cost functions, designed to be used with Light.optimise_settings'''

    def __call__(self, weights, light, desired):
        got = light.light_output(weights, wavelengths=desired.wavelengths)
        cost = self._cost(got, desired)
        print(type(self).__name__, weights, cost)
        return cost

class SimpleCost(CostFunc):
    '''Simple sum of distances across all '''
    def __init__(self, cost=None):
        self.cost = cost

    def _cost(self, x, y):
        x = x.one_nm()
        y = y.one_nm()
        if self.cost is None:
            dist =  np.sum(np.abs(x.values - y.values))
        else:
            dist =  np.sum(np.abs(x.values - y.values) * cost.values)
        return dist

class BandCost(CostFunc):
    '''Cost defined by the differnce in banded integrals along wavelengths, possibly weighted'''
    def __init__(self, bands=None, weights=None):
        if bands is None:
            bands = [(x, x+100) for x in range(300, 700, 100)]
        if weights is None:
            weights = np.array([1 for x in bands])
        self.bands = bands
        self.weights = weights

    def _cost(self, x, y):
        x = np.array(x.banded_integral(bands=self.bands))
        y = np.array(y.banded_integral(bands=self.bands))
        dist =  np.sum(np.abs(x - y) * self.weights)
        return dist
