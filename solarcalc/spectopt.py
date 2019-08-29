# Copyright (c) 2018 Kevin Murray <foss@kdmurray.id.au> This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as np
from scipy.interpolate import interp1d, UnivariateSpline
from scipy import optimize
from scipy.spatial import distance
from sys import stderr
import pandas as pd


class Spectrum(object):
    """Holds a spectral density curve, in Watts per m2 per nm, also integrate and interpolate"""

    def __init__(self, wavelengths, values):
        self.wavelengths = np.array(wavelengths).astype(float)
        self.values = np.array(values).astype(float)

    def interpolated(self, method='linear'):
        '''Create a function to interpolate spectrum to new wavelengths'''
        def interpfn(x):
            minwl = min(self.wavelengths)
            maxwl = max(self.wavelengths)
            wl = self.wavelengths
            vals = self.values
            if any(x < minwl):
                newwl = np.arange(min(x), minwl)
                wl = np.append(newwl, wl)
                vals = np.append(np.zeros_like(newwl), vals)
            if any(x > maxwl):
                newwl = np.arange(maxwl, max(x)+1)
                wl = np.append(wl, newwl)
                vals = np.append(vals, np.zeros_like(newwl))
            interp = interp1d(wl, vals, kind=method)
            return interp(x)
        return interpfn

    def spline(self, k=1, s=0):
        '''Create a function to interpolate spectrum to new wavelengths with a spline'''
        return UnivariateSpline(self.wavelengths, self.values, k=k, s=0)

    def one_nm(self, method='linear'):
        '''Gives an interpolated version of self with 1nm bands'''
        wl = np.arange(min(self.wavelengths), max(self.wavelengths)+1, dtype=int)
        return Spectrum(wl, self.interpolated(method)(wl))

    def par(self, left: int=400, right: int=700):
        """Calculates PAR as integrated uE between `left` and `right` (default 400-700nm)"""
        from scipy.constants import c, h, Avogadro
        left = max(400, min(self.wavelengths))
        right = min(700, max(self.wavelengths))
        wl = np.arange(left, right+1, dtype=int)
        wsqm = self.interpolated()(wl)
        ep = h * (c / (wl/1e9)) # e = hcL^-1; w/ L = wl in meters 
        n = wsqm / ep
        par = np.sum(n/Avogadro)
        return par * 1e6  #  1e6 is to make it uE, not E

    def banded_integral(self, bands):
        '''Gives the integral within each band, where bands are (lhs, rhs)'''
        s = self.spline()
        integr = []
        for band in bands:
            integr.append(s.integral(*band))
        return integr

    def __getitem__(self, item):
        if isinstance(item, slice):
            keep = np.logical_and(self.wavelengths >= item.start, self.wavelengths <= item.stop)
            return Spectrum(self.wavelengths[keep], self.values[keep])


class CostFunc(object):
    '''Base class of cost functions, designed to be used with Light.optimise_settings'''

    def __call__(self, weights, light, desired):
        got = light.light_output(weights, wavelengths=desired.wavelengths)
        cost = self._cost(got, desired)
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
            dist =  np.sum(np.abs(x.values - y.values) * cost.one_nm().values)
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


class Light(object):
    """
    Loads a CSV describing a light source, interpolates the output, and
    allows extraction to custom domains and optimisation to other spectra
    """

    def __init__(self, filename, interpolation='linear'):
        df = pd.read_csv(filename)
        self.channels = {}
        self.wavelengths = df.values[:, 0]
        for chan in df.columns[1:]:
            values = df[chan]
            self.channels[chan] = Spectrum(self.wavelengths, values)

    def __len__(self):
        return len(self.channels)

    def fit_wavelengths(self, wavelengths=None, minwl=None, maxwl=None, step=1):
        '''Interpolate each channel to be defined over certian wavelengths'''
        if wavelengths is None and minwl is None and maxwl is None:
            #raise ValueError("Either wavelengths or maxwl and minwl must be given")
            wavelengths = self.wavelengths
        if wavelengths is None:
            wavelengths = np.arange(minwl, maxwl, step)
        return np.vstack([ch.interpolated()(wavelengths)
                          for ch in self.channels.values()]).T

    def light_output(self, channel_settings, wavelengths=None):
        '''Gives the total output of the light unit over wavelengths given channel_settings'''
        if wavelengths is None:
            wavelengths = self.wavelengths
        outputs = self.fit_wavelengths(wavelengths)
        return Spectrum(wavelengths, np.dot(outputs, channel_settings))

    def max_output(self):
        return self.light_output(np.ones(len(self.channels)))

    def optimise_settings(self, desired, cost_function=SimpleCost()):
        '''Optimise channel settings aiming for desired, using cost_function'''
        initial = np.ones(len(self.channels))/2
        bounds = [(0.,1.)]*len(initial)
        opt = optimize.minimize(cost_function, initial, (self, desired),
                                options={"maxiter": 10000}, bounds=bounds)
        return opt.x
