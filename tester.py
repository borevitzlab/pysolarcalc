import numpy as np
import matplotlib.pyplot as plt

from solarcalc.spectopt import Light, Spectrum, BandCost, SimpleCost

l  = Light("./data/spectra/pretend_rgb.csv", interpolation='linear')

dom = np.linspace(300, 700, 1000)

x = l.fit_wavelengths(dom)


want = Spectrum(dom, np.ones_like(dom)*20 + np.sin(dom))
opt = l.optimise_settings(want)

cf = BandCost(bands=  [(300, 450), (450, 580), (580, 700)],
              weights=[   1.0,        0.1,        1.0])
opt2 = l.optimise_settings(want, cost_function=cf)


from scipy.interpolate import UnivariateSpline

o1 = l.light_output(opt, dom)
o2 = l.light_output(opt2, dom)
l1 = plt.plot(dom, want.values, label='want')
l2 = plt.plot(dom, o1.values, label='simple')
l3 = plt.plot(dom, o2.values, label='banded')
plt.legend()
plt.show()
