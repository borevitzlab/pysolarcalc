import numpy as np
import datetime
from pysolar import solar, util, radiation, elevation
from scipy.interpolate import CubicSpline
from sys import stderr

from .Solarpos import Solarpos

# np.set_printoptions(precision=12)
np.set_printoptions(suppress=True)
SPO = 1360.0

DEGREES_TO_RADIANS = np.pi / 180.0
RADIANS_TO_DEGREES = 180.0 / np.pi
ONEMINUSEPSILON = 1 - np.finfo(np.float64).eps
WAVELENGTH_MICRONS = np.array(
    [0.3, 0.305, 0.31, 0.315, 0.32, 0.325, 0.33, 0.335, 0.34, 0.345, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42,
     0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.57, 0.593, 0.61, 0.63, 0.656,
     0.6676, 0.69, 0.71, 0.718, 0.7244, 0.74, 0.7525, 0.7575, 0.7625, 0.7675, 0.78, 0.8, 0.816, 0.8237, 0.8315, 0.84,
     0.86, 0.88, 0.905, 0.915, 0.925, 0.93, 0.937, 0.948, 0.965, 0.98, 0.9935, 1.04, 1.07, 1.1, 1.12, 1.13, 1.145,
     1.161, 1.17, 1.2, 1.24, 1.27, 1.29, 1.32, 1.35, 1.395, 1.4425, 1.4625, 1.477, 1.497, 1.52, 1.539, 1.558, 1.578,
     1.592, 1.61, 1.63, 1.646, 1.678, 1.74, 1.8, 1.86, 1.92, 1.96, 1.985, 2.005, 2.035, 2.065, 2.1, 2.148, 2.198, 2.27,
     2.36, 2.45, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0])
ETR_SPECTRUM = np.array(
    [535.9, 558.3, 622.0, 692.7, 715.1, 832.9, 961.9, 931.9, 900.6, 911.3, 975.5, 975.9, 1119.9, 1103.8, 1033.8, 1479.1,
     1701.3, 1740.4, 1587.2, 1837.0, 2005.0, 2043.0, 1987.0, 2027.0, 1896.0, 1909.0, 1927.0, 1831.0, 1891.0, 1898.0,
     1892.0, 1840.0, 1768.0, 1728.0, 1658.0, 1524.0, 1531.0, 1420.0, 1399.0, 1374.0, 1373.0, 1298.0, 1269.0, 1245.0,
     1223.0, 1205.0, 1183.0, 1148.0, 1091.0, 1062.0, 1038.0, 1022.0, 998.7, 947.2, 893.2, 868.2, 829.7, 830.3, 814.0,
     786.9, 768.3, 767.0, 757.6, 688.1, 640.7, 606.2, 585.9, 570.2, 564.1, 544.2, 533.4, 501.6, 477.5, 442.7, 440.0,
     416.8, 391.4, 358.9, 327.5, 317.5, 307.3, 300.4, 292.8, 275.5, 272.1, 259.3, 246.9, 244.0, 243.5, 234.8, 220.5,
     190.8, 171.1, 144.5, 135.7, 123.0, 123.8, 113.0, 108.5, 97.5, 92.4, 82.4, 74.6, 68.3, 63.8, 49.5, 48.5, 38.6, 36.6,
     32.0, 28.1, 24.8, 22.1, 19.6, 17.5, 15.7, 14.1, 12.7, 11.5, 10.4, 9.5, 8.6])
WATER_VAPOR_COEFF = np.array(
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.075, 0.0, 0.0, 0.0, 0.0, 0.016, 0.0125, 1.8, 2.5, 0.061, 0.0008,
     0.0001, 0.00001, 0.00001, 0.0006, 0.036, 1.6, 2.5, 0.5, 0.155, 0.00001, 0.0026, 7.0, 5.0, 5.0, 27.0, 55.0, 45.0,
     4.0, 1.48, 0.1, 0.00001, 0.001, 3.2, 115.0, 70.0, 75.0, 10.0, 5.0, 2.0, 0.002, 0.002, 0.1, 4.0, 200.0, 1000.0,
     185.0, 80.0, 80.0, 12.0, 0.16, 0.002, 0.0005, 0.0001, 0.00001, 0.0001, 0.001, 0.01, 0.036, 1.1, 130.0, 1000.0,
     500.0, 100.0, 4.0, 2.9, 1.0, 0.4, 0.22, 0.25, 0.33, 0.5, 4.0, 80.0, 310.0, 15000.0, 22000.0, 8000.0, 650.0, 240.0,
     230.0, 100.0, 120.0, 19.5, 3.6, 3.1, 2.5, 1.4, 0.17, 0.0045])
OZONE_ABSORBTION_COEFF = np.array(
    [10.0, 4.8, 2.7, 1.35, 0.8, 0.38, 0.16, 0.075, 0.04, 0.019, 0.007, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.003, 0.006, 0.009, 0.01400, 0.021, 0.03, 0.04, 0.048, 0.063, 0.075, 0.085, 0.12, 0.119, 0.12, 0.09, 0.065, 0.051,
     0.028, 0.018, 0.015, 0.012, 0.01, 0.008, 0.007, 0.006, 0.005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
UNIFORMLY_MIXED_GAS_ABSORBTION_COEFF = np.array(
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0,
     0.35, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.05, 0.3, 0.02, 0.0002, 0.00011, 0.00001, 0.05, 0.011, 0.005, 0.0006, 0.0, 0.005, 0.13, 0.04,
     0.06, 0.13, 0.001, 0.0014, 0.0001, 0.00001, 0.00001, 0.0001, 0.001, 4.3, 0.2, 21.0, 0.13, 1.0, 0.08, 0.001,
     0.00038, 0.001, 0.0005, 0.00015, 0.00014, 0.00066, 100.0, 150.0, 0.13, 0.0095, 0.001, 0.8, 1.9, 1.3, 0.075, 0.01,
     0.00195, 0.004, 0.29, 0.025])
TRANSMISSITY_COEFF = np.array(
    [0.3, 0.305, 0.31, 0.315, 0.32, 0.325, 0.33, 0.335, 0.34,
     0.345, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
     0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.57, 0.593, 0.61, 0.63, 0.656, 0.668, 0.69, 0.71, 0.718, 0.724,
     0.74, 0.753, 0.758, 0.763, 0.768, 0.78, 0.80, 0.816, 0.824, 0.832, 0.84, 0.86, 0.88, 0.905, 0.915, 0.925,
     0.93, 0.937, 0.948, 0.965, 0.98, 0.994, 1.04, 1.07, 1.10, 1.12, 1.13, 1.145, 1.161, 1.17, 1.20, 1.24,
     1.27, 1.29, 1.32, 1.35, 1.395, 1.443, 1.463, 1.477, 1.497, 1.52, 1.539, 1.558, 1.578, 1.592, 1.61, 1.63,
     1.646, 1.678, 1.74, 1.80, 1.86, 1.92, 1.96, 1.985, 2.005, 2.035, 2.065, 2.10, 2.148, 2.198, 2.27, 2.36,
     2.45, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90,
     4.000])

v = np.array([ETR_SPECTRUM, WATER_VAPOR_COEFF, OZONE_ABSORBTION_COEFF, UNIFORMLY_MIXED_GAS_ABSORBTION_COEFF])

v = np.swapaxes(v, 0, 1)
coeff_spline = CubicSpline(WAVELENGTH_MICRONS, v)
C = 2.9979244e14
CONS = 5.0340365e14
EVOLT = 1.6021891e-19
H = 6.6261762e-34
OMEG = 0.945
OMEGP = 0.095
E = H * C / EVOLT


class Spectra(object):
    def __init__(self, latitude, longitude,
                 low: float = WAVELENGTH_MICRONS[0],
                 high: float = WAVELENGTH_MICRONS[-1],
                 n_wvls: int = 200,
                 wavelengths=[],
                 elevation=0.0,
                 pressure=None
                 ):

        if low > high:
            low, high = high, low
        if low < WAVELENGTH_MICRONS[0]:
            print("Low wavelength is under data, unpredictable results ahead", file=stderr)
        if high > WAVELENGTH_MICRONS[-1]:
            print("high wl is over source data, unpredictable results ahead", file=stderr)

        self.low, self.high = low, high
        if not len(wavelengths):
            self.wavelengths = np.linspace(WAVELENGTH_MICRONS[0], WAVELENGTH_MICRONS[-1], n_wvls)
        else:
            self.wavelengths = sorted(wavelengths)
        self.n_wvls = len(wavelengths)
        # default values
        self.aspect = 180.0
        self.elevation = elevation
        if pressure is None:
            pressure = elevation.get_pressure_with_elevation(self.elevation)
        self.pressure = pressure
        # default to 1998 because numbers.
        self.longitude = longitude
        self.latitude = latitude

        # Input reflectivities
        self.spcrfl = [0.2] * 6
        # Input reflectivity wavelengths
        self.spcwvr = [0.3, 0.7, 0.8, 1.3, 2.5, 4.0]
        # Power on Angstrom turbidity
        self.alpha = 1.14
        # Aerosol assymetry factor (rural assumed)
        self.assym = 0.65

        # Aerosol optical depth at 0.5 microns, base e
        # self.tau500 = None  # 0 < tau500 <= 1.0 >0

        # Precipitable water vapor (cm)
        self.watvap = None  # > 0

    def calc_ozone(self, doy):
        """
        calculates ozone density

        :param: doy
        :rtype: float
        """

        c1 = 100.0
        c2 = 1.5
        c3 = 30.0
        c4 = 152.625
        c5 = 2.0
        c6 = -75.0
        if self.latitude >= 0.0:
            c1 = 150.0
            c2 = 1.28
            c3 = 40.0
            c4 = -30.0
            c5 = 3.0
            c6 = 0.0
            if self.longitude > 0.0:
                c6 = 20.0
        s1 = np.sin(np.radians(0.9865 * (doy + c4)))
        s2 = np.sin(np.radians(c5 * (self.longitude + c6)))
        s3 = np.sin(np.radians(c2 * self.latitude))
        return 0.235 + (c1 + c3 * s1 + 20.0 * s2) * (s3 * s3) / 1000.0

    def prec_h2o_estimate(self, temp, rh):
        tempr = (300.0 / (temp + 273.3))
        v1 = (2.409E12 * rh * pow(tempr, 4.0) * np.exp(-22.64 * tempr))
        return (v1 / 30.0 / (temp + 273.3))

    def calc_vis_spectral_from_array(self, when: datetime.datetime, temperature, relative_humidity, tau500):

        solar_altitude_deg = solar.get_altitude(self.latitude, self.longitude, when,
                                                elevation=self.elevation,
                                                pressure=self.pressure,
                                                temperature=temperature)

        if solar_altitude_deg < 0:
            # print("NIGHTTIME : {}".format(solar_altitude_deg))
            return np.zeros(41, dtype=np.float64)

        zenref = 90 - abs(solar_altitude_deg)

        if zenref > 90:
            print("Zenref greater than 90 wtf!!! : {}".format(zenref), file=stderr)
            return np.zeros(41, dtype=np.float64)

        # tricky! the way that extraterrestrial_irrad is calculated is by getting the erv.
        etr_total = util.extraterrestrial_irrad(self.latitude, self.longitude, when)
        erv = etr_total / 1367.0

        # these are almost identical. opt for pysolar
        amass = radiation.get_air_mass_ratio(altitude_deg=solar_altitude_deg)
        # amass = 1.0 / np.cos(np.radians(zenref)) + 0.50572 * pow(96.07995 - zenref, -1.6364)

        ampress = amass * self.pressure / 1013000
        #if amass < 0:
        #    print(amass, solar_altitude_deg)

        O3 = self.calc_ozone(int(when.timetuple().tm_yday))

        watvap = self.prec_h2o_estimate(temperature, relative_humidity)

        # spectra array
        spec = np.zeros((122, 5))
        # this array should function as the following 5
        # wavelength
        # delta
        # direct
        # integrated direct
        # etr

        # copy wavelengths
        spec[:, 0] = WAVELENGTH_MICRONS[:]

        # delta between wavlengths
        for idx, (wl1, wl2) in enumerate(zip(spec[:-1, 0], spec[1:, 0])):
            spec[idx + 1, 1] = 0.5 * abs(wl2 - wl1)

        cz = np.cos(np.radians(zenref))

        # Ozone mass
        ozone_mass = 1.003454 / np.sqrt((cz * cz) + 0.006908)

        for idx, (wvl, delta, direct, direct_integrated, etr) in enumerate(spec):
            etr = ETR_SPECTRUM[idx] * erv
            spec[idx, -1] = etr
            watvap_coeff = WATER_VAPOR_COEFF[idx]
            ozone_absorb_coeff = OZONE_ABSORBTION_COEFF[idx]
            unif_mix_gas_ab_coeff = UNIFORMLY_MIXED_GAS_ABSORBTION_COEFF[idx]

            # omegl = OMEG * np.exp(-OMEGP * (np.log(wvl / 0.4) * np.log(wvl / 0.4)))
            c1 = tau500 * pow(wvl * 2.0, -self.alpha)
            # Equation 2-4
            Tr = np.exp(-ampress / ((wvl * wvl * wvl * wvl) * (115.6406 - 1.3366 / (wvl * wvl))))
            # Equation 2-9
            To = np.exp(-ozone_absorb_coeff * O3 * ozone_mass)
            # Equation 2-8
            tw_v = float(1.0 + 20.07 * watvap_coeff * watvap * ampress)

            Tw = np.exp(-0.2385 * watvap_coeff * watvap * ampress / pow(tw_v, 0.45))
            # print(tw_v, Tw)

            # Equation 2-11
            # cast to float makes this value real enough for numpy
            Tu = np.exp(
                -1.41 * unif_mix_gas_ab_coeff * ampress / pow(1.0 + 118.3 * unif_mix_gas_ab_coeff * ampress, 0.45))
            # Equation 2-6, sort of
            Ta = np.exp(-c1 * ampress)
            # ........... Direct energy .............
            # Temporary variable
            c2 = etr * To * Tw * Tu
            # Equation 2-1
            direct = c2 * Tr * Ta
            spec[idx, 2] = direct

            # direct integration
            direct_integrated = delta * (direct + spec[max(idx - 1, 0), 2])

            spec[idx, 3] = direct_integrated
            # prevdirect = direct

        # these probably arent needed.
        # totdirect = np.sum(spec[:, 3])
        # totvis = np.sum(np.clip(spec[14:55, 3], 0.0, np.inf))
        # visible_integration = spec[14:55, 3]
        # visble_wavelength_microns = WAVELENGTH_MICRONS[14:55]
        # trad = np.sum(np.clip(spec[:, 2], 0.0, np.inf))
        # # spec[:, 2] /= (solar_irradiance * trad)
        # vrad = spec[14:55, 2]

        # etr = util.extraterrestrial_irrad(self.latitude, self.longitude, when)
        # only return direct and visible
        return np.append([etr_total]+spec[14:55, 3])

    def calc_vis_spectral(self,
                          when: datetime.datetime,
                          temperature: float,
                          relative_humidity: float,
                          tau500: float):

        if temperature < 1:
            temperature = 1
    

        solar_altitude_deg = solar.get_altitude(self.latitude, self.longitude, when,
                                                elevation=self.elevation,
                                                pressure=self.pressure,
                                                temperature=temperature)

        zenref = 90 - abs(solar_altitude_deg)

        if zenref > 90:
            print("Zenref greater than 90 wtf!!! : {}".format(zenref), file=stderr)
            return np.zeros(self.n_wvls+1, dtype=np.float64)
        # tricky! the way that extraterrestrial_irrad is calculated is by getting the erv.
        etr_total = util.extraterrestrial_irrad(self.latitude, self.longitude, when=when)
        #if when.hour == 12 and when.minute == 0:
        #    print(when,"\t", temperature, "\t",relative_humidity, "\t",etr_total)
        erv = etr_total / 1367.0
        # these are almost identical. DO NOT OPT FOR PYSOLAR UNTIL THEY HAVE IMPLEMENTED
        # Airmass Kasten, F. and Young, A. 1989. Revised optical air mass tables
        # amass = radiation.get_air_mass_ratio(altitude_deg=solar_altitude_deg)

        amass = 1.0 / np.cos(np.radians(zenref)) + 0.50572 * pow(96.07995 - zenref, -1.6364)

        ampress = amass * (self.pressure / 1.013e6)

        O3 = self.calc_ozone(int(when.timetuple().tm_yday))

        watvap = self.prec_h2o_estimate(temperature, relative_humidity)

        cz = np.cos(np.radians(zenref))

        # Ozone mass
        ozone_mass = 1.003454 / np.sqrt((pow(cz, 2)) + 0.006908)


        if not 0 <= abs(solar_altitude_deg) <= 90:
            # print("NIGHTTIME : {}".format(solar_altitude_deg))
            # print("Solar alt out of range!!! : {}".format(solar_altitude_deg))
            # print(when, self.pressure, temperature, solar_altitude_deg, amass)
            return np.zeros(self.n_wvls+1, dtype=np.float64)

        if not 0 <= amass < 1000:
            print("Solar alt out of range!!! : {}".format(solar_altitude_deg), file=stderr)
            return np.zeros(self.n_wvls + 1, dtype=np.float64)

        def calc(wvl):
            etr_spec, watvap_coeff, ozone_absorb_coeff, unif_mix_gas_ab_coeff = coeff_spline(wvl)
            watvap_coeff = max(watvap_coeff, 0.0)
            ozone_absorb_coeff = max(ozone_absorb_coeff, 0.0)
            unif_mix_gas_ab_coeff = max(unif_mix_gas_ab_coeff, 0.0)
            etr_spec = max(etr_spec, 0.0)
            etr = etr_spec * erv

            # omegl = OMEG * np.exp(-OMEGP * (np.log(wvl / 0.4) * np.log(wvl / 0.4)))
            c1 = tau500 * pow(wvl * 2.0, -self.alpha)
            # Equation 2-4
            Tr = np.exp(-ampress / (pow(wvl, 4) * (115.6406 - 1.3366 / pow(wvl, 2))))
            # Equation 2-9
            To = np.exp(-ozone_absorb_coeff * O3 * ozone_mass)
            # Equation 2-8

            tw_v = 1.0 + 20.07 * watvap_coeff * watvap * ampress

            Tw = np.exp(-0.2385 * watvap_coeff * watvap * ampress / pow(tw_v, 0.45))

            # Equation 2-11
            # cast to float makes this value real enough for numpy
            Tu = np.exp(
                -1.41 * unif_mix_gas_ab_coeff * ampress / pow(1.0 + 118.3 * unif_mix_gas_ab_coeff * ampress, 0.45))
            # Equation 2-6, sort of
            Ta = np.exp(-c1 * ampress)
            # ........... Direct energy .............
            # Temporary variable
            c2 = etr * To * Tw * Tu
            # Equation 2-1
            direct = c2 * Tr * Ta
            return direct

        vectorized_calc = np.vectorize(calc)

        direct_spec = vectorized_calc(self.wavelengths)

        def integrate(wl_delta, direct, prev_direct):
            return wl_delta * (direct + prev_direct)

        vectorized_integrate = np.vectorize(integrate)
        integrated_d = vectorized_integrate(np.append([0], np.diff(self.wavelengths)),
                                            direct_spec,
                                            np.append([0], direct_spec[:-1]))
        #return np.append([etr_total], integrated_d)
        # KDM testing
        return np.append([etr_total], direct_spec)
