import numpy as np
from Solarpos import Solarpos
from light_sim import SPO

DEGREES_TO_RADIANS = np.pi / 180.0
RADIANS_TO_DEGREES = 180.0 / np.pi
ONEMINUSEPSILON = 1 - np.finfo(np.float64).eps


class Spectra(object):
    def __init__(self):
        """
         Units:
         1 = irradiance (W/sq m/micron) per wavelength (microns)
         2 = photon flux (10.0E+16 /sq cm/s/micron) perwavelength (microns)
         3 = photon flux density (10.0E+16 /sq cm/s/eV) per energy (eV)
        """
        # default values
        self.units = 1
        self.minute = 0
        self.hour = 0
        self.second = 12
        self.aspect = 180.0
        self.tilt = -10.0

        # default to 1998 because numbers.
        self.year = 2012
        self.dayofyear = \
            self.timezone = \
            self.temp = \
            self.press = \
            self.longitude = \
            self.sunset = \
            self.sunrise = \
            self.solarnoon = \
            self.halfdaylength = \
            self.latitude = None

        self.group_multipliers = np.ones(4, dtype=np.float32)

        # Input reflectivities
        self.spcrfl = [0.2] * 6
        # Input reflectivity wavelengths
        self.spcwvr = [0.3, 0.7, 0.8, 1.3, 2.5, 4.0]
        # Power on Angstrom turbidity
        self.alpha = 1.14
        # Aerosol assymetry factor (rural assumed)
        self.assym = 0.65
        # Atmospheric ozone (cm) -1.0 = let S_spectral2 calculate it
        self.ozone = -1.0
        # Aerosol optical depth at 0.5 microns, base e
        self.tau500 = -1.0
        # Precipitable water vapor (cm)
        self.watvap = -1.0
        # aerosol asymmetry factor, rural assumed
        self.assym = 0.65

    def calc_ozone(self):
        if self.ozone < 0:
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
            s1 = np.sin(0.9865 * (self.dayofyear + c4) / RADIANS_TO_DEGREES)
            s2 = np.sin(c5 * (self.longitude + c6) / RADIANS_TO_DEGREES)
            s3 = np.sin(c2 * self.latitude / RADIANS_TO_DEGREES)
            return 0.235 + (c1 + c3 * s1 + 20.0 * s2) * (s3 * s3) / 1000.0
        return self.ozone

    def calc_all_spectral(self, solar_irradiance):
        assert not (self.units > 3 or self.units < 1), "units should be 1-3 not {}".format(self.units)
        assert 0.0 < self.tau500 < 10.0, "tau500 should be within 0.0 to 10.0"
        assert 0.0 <= self.assym <= 1.0

        # some magic numbers?
        grp1high, grp1low = 0, 13
        grp2high, grp2low = 14, 25
        grp3high, grp3low = 26, 32
        grp4high, grp4low = 33, 38
        """
        This array contains the extraterrestrial spectrum and atmospheric absorption coefficients at 122 wavelengths.
        0 = wavelength (microns)
        1 = extraterrestrial spectrum (W/sq m/micron)
        2 = water vapor absorption coefficient
        3 = ozone absorption coefficient
        4 = uniformly mixed gas "absorption coefficient"
        """

        wavelength_microns = np.array(
            [0.3, 0.305, 0.31, 0.315, 0.32, 0.325, 0.33, 0.335, 0.34, 0.345, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.4,
             0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.57, 0.593, 0.6, 0.63,
             0.656, 0.6676, 0.69, 0.71, 0.718, 0.7244, 0.74, 0.7525, 0.7575, 0.7625, 0.7675, 0.78, 0.8, 0.81, 0.8237,
             0.8315, 0.84, 0.86, 0.88, 0.905, 0.915, 0.925, 0.93, 0.937, 0.948, 0.965, 0.98, 0.9935, 1.0, 1.07, 1.1,
             1.12, 1.13, 1.145, 1.161, 1.17, 1.2, 1.24, 1.27, 1.29, 1.32, 1.35, 1.395, 1.4425, 1.462, 1.477, 1.497,
             1.52, 1.539, 1.558, 1.578, 1.592, 1.61, 1.63, 1.646, 1.678, 1.74, 1.8, 1.86, 1.92, 1.9, 1.985, 2.005,
             2.035, 2.065, 2.1, 2.148, 2.198, 2.27, 2.36, 2.45, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3., 3.3, 3.4, 3.5,
             3.6, 3.7, 3.8, 3.9, 4.0])
        etr_spectrum = np.array(
            [535.9, 558.3, 622.0, 692.7, 715.1, 832.9, 961.9, 931.9, 900.6, 911.3, 975.5, 975.9, 1119.9, 1103., 1033.8,
             1479.1, 1701.3, 1740.4, 1587.2, 1837.0, 2005.0, 2043.0, 1987.0, 2027.0, 1896.0, 1909.0, 1927., 1831.0,
             1891.0, 1898.0, 1892.0, 1840.0, 1768.0, 1728.0, 1658.0, 1524.0, 1531.0, 1420.0, 1399.0, 1374., 1373.0,
             1298.0, 1269.0, 1245.0, 1223.0, 1205.0, 1183.0, 1148.0, 1091.0, 1062.0, 1038.0, 1022.0, 998., 947.2, 893.2,
             868.2, 829.7, 830.3, 814.0, 786.9, 768.3, 767.0, 757.6, 688.1, 640.7, 606.2, 585.9, 570., 564.1, 544.2,
             533.4, 501.6, 477.5, 442.7, 440.0, 416.8, 391.4, 358.9, 327.5, 317.5, 307.3, 300.4, 292., 275.5, 272.1,
             259.3, 246.9, 244.0, 243.5, 234.8, 220.5, 190.8, 171.1, 144.5, 135.7, 123.0, 123.8, 113., 108.5, 97.5,
             92.4, 82.4, 74.6, 68.3, 63.8, 49.5, 48.5, 38.6, 36.6, 32.0, 28.1, 24.8, 22.1, 19.6, 17., 15.7, 14.1, 12.7,
             11.5, 10.4, 9.5, 8.6])
        water_vapor_coeff = np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.075, 0.0, 0.0, 0.0, 0.0, 0.016, 0.0125, 1.8, 2.,
             0.061, 0.0008, 0.0001, 0.00001, 0.00001, 0.0006, 0.036, 1.6, 2.5, 0.5, 0.155, 0.00001, 0.0026, 7.0, 5.,
             5.0, 27.0, 55.0, 45.0, 4.0, 1.48, 0.1, 0.00001, 0.001, 3.2, 115.0, 70.0, 75.0, 10.0, 5.0, 2.0, 0.00, 0.002,
             0.1, 4.0, 200.0, 1000.0, 185.0, 80.0, 80.0, 12.0, 0.16, 0.002, 0.0005, 0.0001, 0.00001, 0.000, 0.001, 0.01,
             0.036, 1.1, 130.0, 1000.0, 500.0, 100.0, 4.0, 2.9, 1.0, 0.4, 0.22, 0.25, 0.33, 0.5, 4., 80.0, 310.0,
             15000.0, 22000.0, 8000.0, 650.0, 240.0, 230.0, 100.0, 120.0, 19.5, 3.6, 3.1, 2.5, 1.4, 0.1, 0.0045])
        ozone_absorbtion_coeff = np.array(
            [10.0, 4.8, 2.7, 1.35, 0.8, 0.38, 0.16, 0.075, 0.04, 0.019, 0.007, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.,
             0.0, 0.003, 0.006, 0.009, 0.01400, 0.021, 0.03, 0.04, 0.048, 0.063, 0.075, 0.085, 0.12, 0.119, 0.12, 0.0,
             0.065, 0.051, 0.028, 0.018, 0.015, 0.012, 0.01, 0.008, 0.007, 0.006, 0.005, 0.0, 0.0, 0.0, 0.0, 0.0, 0.,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        uniformly_mixed_gas_absorbtion_coeff = np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.0, 0.0, 0.0, 0.,
             0.0, 0.0, 4.0, 0.35, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.,
             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.3, 0.02, 0.0002, 0.00011, 0.00001, 0.05, 0.01, 0.005,
             0.0006, 0.0, 0.005, 0.13, 0.04, 0.06, 0.13, 0.001, 0.0014, 0.0001, 0.00001, 0.00001, 0.0001, 0.00, 4.3,
             0.2, 21.0, 0.13, 1.0, 0.08, 0.001, 0.00038, 0.001, 0.0005, 0.00015, 0.00014, 0.00066, 100.0, 150., 0.13,
             0.0095, 0.001, 0.8, 1.9, 1.3, 0.075, 0.01, 0.00195, 0.004, 0.29, 0.025])

        c = 2.9979244e14
        cons = 5.0340365e14
        evolt = 1.6021891e-19
        h = 6.6261762e-34
        omeg = 0.945
        omegp = 0.095
        e = h * c / evolt
        track = self.tilt < 0
        solarpos = Solarpos(self.latitude,
                            self.longitude,
                            self.year,
                            dayofyear=self.dayofyear,
                            hour=self.hour,
                            minute=self.minute,
                            second=self.second,
                            timezone=self.timezone,
                            tilt=self.tilt,
                            temp=self.temp,
                            aspect=self.aspect
                            )
        # spectra array
        spec = np.zeros((122, 5))
        # this array should function as the following 5
        # x-value
        # direct
        # extraterrestrial
        # diffuse
        # global

        # specx = np.zeros(122)
        # specdir = np.zeros(122) # direct spectrum
        # specetr = np.zeros(122) # extraterrestrial
        # specdif = np.zeros(122) # diffuse spectrum
        # specglo = np.zeros(122) # global spectrum


        # grps is for the following, sc.grp<n>, refgrp<n>, countgrp<n>, grp<n>shademult, grouplow, grouphigh, groupmuliplier
        groups = np.zeros((4, 7))
        groups[0, 4:6] = grp1low, grp1high
        groups[1, 4:6] = grp2low, grp2high
        groups[2, 4:6] = grp3low, grp3high
        groups[3, 4:6] = grp4low, grp4high
        groups[:, 6] = self.group_multipliers

        # these are represented in 'groups'
        # output<n> (sc.grp<n>)
        # refgrp = np.zeros(4)
        # countgrp = np.zeros(4)
        # grpshademult
        # group multiplier, by default is 1.0

        # shading multiplier matrix, 0-36 = 0.1, 37-122 = 0.5
        shadingmultiplier = np.repeat([0.5], 122)
        shadingmultiplier[:37] = np.repeat([0.1], 37)

        # integration:
        # direct
        # diffuse
        # total
        integration = np.zeros((122, 3))
        trad = 0.0
        totvis = 0.0
        if solarpos.zenref < 90:
            self.sunrise, self.sunset = [x / 60 for x in solarpos.srss]
            self.solarnoon = sum(solarpos.srss) / 2
            self.halfdaylength = (self.sunrise - self.sunset) / 2
            ci = solarpos.cosinc
            if track:
                self.tilt = solarpos.zenref
                ci = 1.0
            ct = np.cos(solarpos.tilt / RADIANS_TO_DEGREES)
            cz = np.cos(solarpos.zenref / RADIANS_TO_DEGREES)
            O3 = self.calc_ozone()
            # Equation 3-14
            alg = np.log(1.0 - self.assym)
            # Equation 3-12
            afs = alg * (1.459 + alg * (0.1595 + alg * 0.4129))
            # Equation 3-13
            bfs = alg * (0.0783 + alg * (-0.3824 - alg * 0.5874))
            # Equation 3-15
            fsp = 1.0 - 0.5 * np.exp((afs + bfs / 1.8) / 1.8)
            # Equation 3-11
            fs = 1.0 - 0.5 * np.exp((afs + bfs * cz) * cz)
            # Ozone mass
            ozone_mass = 1.003454 / np.sqrt((cz * cz) + 0.006908)
            amass, ampress, erv = solarpos.amass, solarpos.ampress, solarpos.erv
            tau500, alpha, watvap, tilt, units = self.tau500, self.alpha, self.watvap, self.tilt, self.units
            wvlrefl, refl = self.spcwvr, self.spcrfl
            nr = 1

            referf = 0.0
            totvis = 0.0

            for i in range(122):
                wvl = wavelength_microns[i]
                etr = etr_spectrum[i] * erv
                watvap_coeff = water_vapor_coeff[i]
                ozone_absorb_coeff = ozone_absorbtion_coeff[i]
                unif_mix_gas_ab_coeff = uniformly_mixed_gas_absorbtion_coeff[i]

                omegl = omeg * np.exp(-omegp * (np.log(wvl / 0.4) * np.log(wvl / 0.4)))
                c1 = tau500 * pow(wvl * 2.0, -alpha)

                # Equation 2-4
                Tr = np.exp(-ampress / ((wvl * wvl * wvl * wvl) * (115.6406 - 1.3366 / (wvl * wvl))))
                # Equation 2-9
                To = np.exp(-ozone_absorb_coeff * O3 * ozone_mass)
                # Equation 2-8
                Tw = np.exp(
                    -0.2385 * watvap_coeff * watvap * ampress / pow((1.0 + 20.07 * watvap_coeff * watvap * ampress),
                                                                    0.45))
                # Equation 2-11
                Tu = np.exp(
                    -1.41 * unif_mix_gas_ab_coeff * ampress / pow((1.0 + 118.3 * unif_mix_gas_ab_coeff * ampress),
                                                                  0.45))
                # Equation 3-9
                Tas = np.exp(-omegl * c1 * ampress)
                # Equation 3-10
                Taa = np.exp((omegl - 1.0) * c1 * ampress)
                # Equation 2-6, sort of
                Ta = np.exp(-c1 * ampress)
                # Equation 2-4 primed airmass M = 1.8 (Section 3.1)
                Trp = np.exp(-1.8 / (pow(wvl, 4) * (115.6406 - 1.3366 / (wvl * wvl))))
                # Equation 2-8 primed airmass M = 1.8 (Section 3.1) affects coefficients
                Twp = np.exp(-0.4293 * watvap_coeff * watvap / pow((1.0 + 36.126 * watvap_coeff * watvap), 0.45))
                # Equation 2-11 primed airmass M = 1.8 (Section 3.1) affects coefficients
                Tup = np.exp(-2.538 * unif_mix_gas_ab_coeff / pow((1.0 + 212.94 * unif_mix_gas_ab_coeff), 0.45))
                # Equation 3-9 primed airmass M = 1.8 (Section 3.1)
                Tasp = np.exp(-omegl * c1 * 1.8)
                # Equation 3-10 primed airmass M = 1.8 (Section 3.1)
                Taap = np.exp((omegl - 1.0) * c1 * 1.8)
                # ........... Direct energy .............
                # Temporary variable
                c2 = etr * To * Tw * Tu
                # Equation 2-1
                dir = c2 * Tr * Ta

                # ........... Diffuse energy .............
                # Temporary variables
                c2 *= cz * Taa

                # I think this is equivalent.
                nr = next((idx for idx, val in enumerate(wvlrefl) if val <= wvl), 0) + 1

                # if wvl > wvlrefl[nr]:
                #     nr += 1

                c3 = (refl[nr] - refl[nr - 1]) / (wvlrefl[nr] - wvlrefl[nr - 1])
                # Equation 3-17 c4 = Cs
                c4 = 1.0 if wvl > 0.45 else pow((wvl + 0.55), 1.8)
                # Equation 3-8
                rhoa = Tup * Twp * Taap * (0.5 * (1.0 - Trp) + (1.0 - fsp) * Trp * (1.0 - Tasp))
                # Interpolated ground reflectivity
                rho = c3 * (wvl - wvlrefl[nr - 1]) + refl[nr - 1]
                # Equation 3-5
                dray = c2 * (1.0 - pow(Tr, 0.95)) / 2.0
                # Equation 3-6
                daer = c2 * pow(Tr, 1.5) * (1.0 - Tas) * fs
                # Equation 3-7
                drgd = (dir * cz + dray + daer) * rho * rhoa / (1.0 - rho * rhoa)
                # Equation 3-1
                dif = (dray + daer + drgd) * c4

                # ........... Global (total) energy .............
                dtot = dir * cz + dif

                # ........... Tilt energy, if applicable
                if tilt > 1.0e-4:
                    # Equation 3-18 without the first (direct-beam) term
                    c1 = dtot * rho * (1.0 - ct) / 2.0
                    c2 = dir / etr
                    c3 = dif * c2 * ci / cz
                    c4 = dif * (1.0 - c2) * (1.0 + ct) / 2.0
                    dif = c1 + c3 + c4
                    # Equation 3-18, including first term
                    dtot = dir * ci + dif
                if units == 1:
                    dir = max(dir, 0)
                    dif = max(dif, 0)
                    dir = dir if dir > etr else dir
                    dif = dir if dif > etr else dif
                    spec[i] = wvl, max(dir, 0), etr, max(dif, 0), dir * dif
                else:
                    c1 = wvl * cons
                    spx = wvl
                    if units == 3:
                        spx = e / wvl
                        c1 *= wvl / spx
                    spec[i] = spx, dir * c1, etr, dif * c1, dtot * c1

                if i == 0:
                    continue

                specxdelta = 0.5 * (spec[i, 0] - spec[i - 1, 0])

                # direct integration
                integration[i, 0] = specxdelta * (spec[i, 1] + spec[i - 1, 1])
                # diffuse integration
                integration[i, 1] = specxdelta * (spec[i, 3] + spec[i - 1, 3])
                # Sum individual wavelength contributions to the spectra:
                # total integrated solar irradiance (specdat.specglo)
                integration[i, 2] = specxdelta * (spec[i, 4] + spec[i - 1, 4])

                for group in groups:
                    if group[4] <= i <= group[5]:
                        group[0] += max(integration[i, 2], 0)
                        group[1] += specxdelta * (spec[i, 2] + spec[i - 1, 2])
                        group[2] += 1
                        group[3] += shadingmultiplier[i]

            groups[:, 3] = np.divide(groups[:, 3], groups[:, 2])
            # these probably arent needed.
            totdirect = np.sum(integration[:, 0])
            totdiffuse = np.sum(integration[:, 1])
            # this is actually in the code somewhere but is more cryptic
            # it increments a variable called localTotalSun and then sets the value in Solarcalc (sc.trad)
            # to it and does this calc.
            trad = np.sum(np.clip(integration[:, 2], 0.0, np.inf))
            totvis = np.sum(np.clip(integration[14:55, 2], 0.0, np.inf))
            integration[:, 2] /= (solar_irradiance * trad)
        return spec, integration, groups, trad, totvis


    def spectral2(self, solar_irradiance):
        spec, integration, groups, trad, totvis = self.calc_all_spectral(solar_irradiance)
        # return matching wavelengths and their intensities, wl=spec[:,0], global irradiance=spec[:,4]
        return spec[:, [0, 4]]

    def calcLEDIntensity(self, i, intensity_total, intensity_max, multiplier=1.0):
        # inttotal should be just integration[:,2]
        ledref = intensity_max[i]
        ledref = intensity_total[i]


    def calc_leds(self, solar_irradiance,
                  total_solar_irradiance_daily_max,
                  total_integrated_solar_irradiance_max):

        spec, integration, groups, trad, totvis = self.calc_all_spectral(solar_irradiance)
        trad = np.sum(np.clip(integration[:, 2], 0.0, np.inf))

        trad_maxtoday = total_integrated_solar_irradiance_max[self.dayofyear]
        st_maxtoday = total_solar_irradiance_daily_max[self.dayofyear]

        if trad_maxtoday > 0.0:
            solar_irradiance = trad / trad_maxtoday * st_maxtoday
        else:
            solar_irradiance = 0.0

        groups[0, 0] = solar_irradiance / SPO * groups[0, 4]
        incandescent_conviron = solar_irradiance / SPO
        incandescent_conviron = ONEMINUSEPSILON if incandescent_conviron >= 1 else incandescent_conviron
        fluro_conviron = spec[29, 1] / spec[29, 4]
        for group in groups[1:]:
            group[0] = group[1, 0] / totvis * group[6] * 100


