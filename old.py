# def calc_all_spectral(self, solar_irradiance):
# """
# calculates the spectrum, integrated irradiation, spectral groups, total radiation, and total visible light
#
# BROKEN, spectra are NAN
#
# :param solar_irradiance: solar irradiance as calculated in light_sim.py
# :return: tuple of spectrum, integrated irradiation, spectral groups, total radiation, and total visible light
# :rtype: tuple[np.ndarray, np.ndarray, np.ndarray, float, float]
# """
# assert not (self.units > 3 or self.units < 1), "units should be 1-3 not {}".format(self.units)
# assert (0.0 <= self.tau500 <= 10.0), "tau500 should be within 0.0 to 10.0"
# assert (0.0 <= self.assym <= 1.0), "assym should be between 0 and 1"
#
# """
# This array contains the extraterrestrial spectrum and atmospheric absorption coefficients at 122 wavelengths.
# 0 = wavelength (microns)
# 1 = extraterrestrial spectrum (W/sq m/micron)
# 2 = water vapor absorption coefficient
# 3 = ozone absorption coefficient
# 4 = uniformly mixed gas "absorption coefficient"
# """
# # spline these
#
# track = self.tilt < 0
# solarpos = Solarpos(self.latitude,
#                     self.longitude,
#                     self.year,
#                     dayofyear=self.dayofyear,
#                     hour=self.hour,
#                     minute=self.minute,
#                     second=self.second,
#                     timezone=self.timezone,
#                     tilt=self.tilt,
#                     temp=self.temp,
#                     aspect=self.aspect
#                     )
# # spectra array
# spec = np.zeros((122, 5))
# # this array should function as the following 5
# # x-value
# # direct
# # extraterrestrial
# # diffuse
# # global
#
# # specx = np.zeros(122)
# # specdir = np.zeros(122) # direct spectrum
# # specetr = np.zeros(122) # extraterrestrial
# # specdif = np.zeros(122) # diffuse spectrum
# # specglo = np.zeros(122) # global spectrum
#
# # these are represented in 'groups'
# # output<n> (sc.grp<n>)
# # refgrp = np.zeros(4)
# # countgrp = np.zeros(4)
# # grpshademult
# # group multiplier, by default is 1.0
#
# # shading multiplier matrix, 0-36 = 0.1, 37-122 = 0.5
# shadingmultiplier = np.repeat([0.5], 122)
# shadingmultiplier[:37] = np.repeat([0.1], 37)
#
# # integration:
# # direct
# # diffuse
# # total
# integration = np.zeros((122, 2))
# trad = 0.0
# totvis = 0.0
# if solarpos.zenref < 90:
#     self.sunrise, self.sunset = [x / 60 for x in solarpos.srss]
#     self.solarnoon = sum(solarpos.srss) / 2
#     self.halfdaylength = (self.sunrise - self.sunset) / 2
#     ci = solarpos.cosinc
#     if track:
#         self.tilt = solarpos.zenref
#         ci = 1.0
#     ct = np.cos(np.radians(solarpos.tilt))
#     cz = np.cos(np.radians(solarpos.zenref))
#     O3 = self.calc_ozone()
#     # Equation 3-14
#     alg = np.log(1.0 - self.assym)
#     # Equation 3-12
#     afs = alg * (1.459 + alg * (0.1595 + alg * 0.4129))
#     # Equation 3-13
#     bfs = alg * (0.0783 + alg * (-0.3824 - alg * 0.5874))
#     # Equation 3-15
#     fsp = 1.0 - 0.5 * np.exp((afs + bfs / 1.8) / 1.8)
#     # Equation 3-11
#     fs = 1.0 - 0.5 * np.exp((afs + bfs * cz) * cz)
#     # Ozone mass
#     ozone_mass = 1.003454 / np.sqrt((cz * cz) + 0.006908)
#     amass, ampress, erv = float(solarpos.amass), float(solarpos.ampress), float(solarpos.erv)
#     tau500, alpha, watvap, tilt, units = self.tau500, self.alpha, self.watvap, self.tilt, self.units
#     wvlrefl, refl = self.spcwvr, self.spcrfl
#     nr = 1
#
#     referf = 0.0
#     totvis = 0.0
#
#     for i in range(122):
#         wvl = WAVELENGTH_MICRONS[i]
#         etr = ETR_SPECTRUM[i] * erv
#         watvap_coeff = WATER_VAPOR_COEFF[i]
#         ozone_absorb_coeff = OZONE_ABSORBTION_COEFF[i]
#         unif_mix_gas_ab_coeff = UNIFORMLY_MIXED_GAS_ABSORBTION_COEFF[i]
#
#         omegl = OMEG * np.exp(-OMEGP * (np.log(wvl / 0.4) * np.log(wvl / 0.4)))
#         c1 = tau500 * np.power(wvl * 2.0, -alpha)
#
#         # Equation 2-4
#         Tr = np.exp(-ampress / ((wvl * wvl * wvl * wvl) * (115.6406 - 1.3366 / (wvl * wvl))))
#         # Equation 2-9
#         To = np.exp(-ozone_absorb_coeff * O3 * ozone_mass)
#         # Equation 2-8
#         Tw = np.exp(-0.2385 * watvap_coeff * watvap * ampress / np.power(
#             (1.0 + 20.07 * watvap_coeff * watvap * ampress), 0.45))
#         # print(watvap, ampress, Tw)
#         # Equation 2-11
#         Tu = np.exp(
#             -1.41 * unif_mix_gas_ab_coeff * ampress / np.power(1.0 + 118.3 * unif_mix_gas_ab_coeff * ampress,
#                                                                0.45))
#         # Equation 3-9
#         Tas = np.exp(-omegl * c1 * ampress)
#         # Equation 3-10
#         Taa = np.exp((omegl - 1.0) * c1 * ampress)
#         # Equation 2-6, sort of
#         Ta = np.exp(-c1 * ampress)
#         # Equation 2-4 primed airmass M = 1.8 (Section 3.1)
#         Trp = np.exp(-1.8 / (np.power(wvl, 4) * (115.6406 - 1.3366 / (wvl * wvl))))
#         # Equation 2-8 primed airmass M = 1.8 (Section 3.1) affects coefficients
#         Twp = np.exp(-0.4293 * watvap_coeff * watvap / np.power((1.0 + 36.126 * watvap_coeff * watvap), 0.45))
#         # Equation 2-11 primed airmass M = 1.8 (Section 3.1) affects coefficients
#         Tup = np.exp(-2.538 * unif_mix_gas_ab_coeff / np.power((1.0 + 212.94 * unif_mix_gas_ab_coeff), 0.45))
#         # Equation 3-9 primed airmass M = 1.8 (Section 3.1)
#         Tasp = np.exp(-omegl * c1 * 1.8)
#         # Equation 3-10 primed airmass M = 1.8 (Section 3.1)
#         Taap = np.exp((omegl - 1.0) * c1 * 1.8)
#         # ........... Direct energy .............
#         # Temporary variable
#         c2 = etr * To * Tw * Tu
#         # Equation 2-1
#         direct = c2 * Tr * Ta
#
#         # ........... Diffuse energy .............
#         # Temporary variables
#         c2 *= cz * Taa
#
#         # I think this is equivalent.
#         # nr = next((idx for idx, val in enumerate(wvlrefl) if val <= wvl), 0) + 1
#         # to this
#         if wvl > wvlrefl[nr]:
#             nr += 1
#
#         c3 = (refl[nr] - refl[nr - 1]) / (wvlrefl[nr] - wvlrefl[nr - 1])
#         # Equation 3-17 c4 = Cs
#         c4 = 1.0 if wvl > 0.45 else np.power((wvl + 0.55), 1.8)
#         # Equation 3-8
#         rhoa = Tup * Twp * Taap * (0.5 * (1.0 - Trp) + (1.0 - fsp) * Trp * (1.0 - Tasp))
#         # Interpolated ground reflectivity
#         rho = c3 * (wvl - wvlrefl[nr - 1]) + refl[nr - 1]
#         # Equation 3-5
#         dray = c2 * (1.0 - np.power(Tr, 0.95)) / 2.0
#         # Equation 3-6
#         daer = c2 * np.power(Tr, 1.5) * (1.0 - Tas) * fs
#         # Equation 3-7
#         drgd = (direct * cz + dray + daer) * rho * rhoa / (1.0 - rho * rhoa)
#         # Equation 3-1
#         diffuse = (dray + daer + drgd) * c4
#
#         # ........... Global (total) energy .............
#         dtot = direct * cz + diffuse
#
#         # ........... Tilt energy, if applicable
#         if tilt > 1.0e-4:
#             # Equation 3-18 without the first (direct-beam) term
#             c1 = dtot * rho * (1.0 - ct) / 2.0
#             c2 = direct / etr
#             c3 = diffuse * c2 * ci / cz
#             c4 = diffuse * (1.0 - c2) * (1.0 + ct) / 2.0
#             diffuse = c1 + c3 + c4
#             # Equation 3-18, including first term
#             dtot = direct * ci + diffuse
#         if units == 1:
#             direct = max(direct, 0)
#             diffuse = max(diffuse, 0)
#             direct = direct if direct > etr else direct
#             diffuse = direct if diffuse > etr else diffuse
#             spec[i] = wvl, max(direct, 0), etr, max(diffuse, 0), direct * diffuse
#         else:
#             c1 = wvl * CONS
#             spx = wvl
#             if units == 3:
#                 spx = E / wvl
#                 c1 *= wvl / spx
#             spec[i] = spx, direct * c1, etr, diffuse * c1, dtot * c1
#
#         if i == 0:
#             continue
#
#         specxdelta = 0.5 * (spec[i, 0] - spec[i - 1, 0])
#
#         # direct integration
#         integration[i, 0] = specxdelta * (spec[i, 1] + spec[i - 1, 1])
#         # Sum individual wavelength contributions to the spectra:
#         # total integrated solar irradiance (specdat.specglo)
#         integration[i, 1] = specxdelta * (spec[i, 4] + spec[i - 1, 4])
#
#     # these probably arent needed.
#     totdirect = np.sum(integration[:, 0])
#
#     # this is actually in the code somewhere but is more cryptic
#     # it increments a variable called localTotalSun and then sets the value in Solarcalc (sc.trad)
#     # to it and does this calc.
#     trad = np.sum(np.clip(integration[:, 1], 0.0, np.inf))
#     totvis = np.sum(np.clip(integration[14:55, 1], 0.0, np.inf))
#     integration[:, 1] /= (solar_irradiance * trad)
#
# # integration  (W/sq m/micron)
# return spec, integration, trad, totvis

# def spectral2(self, solar_irradiance):
#     """
#     Should operate in the same way as SolarCalc/spectral.java:L25: S_spectral2
#
#     :param solar_irradiance: solar irradiance to drive calculation from
#     :return: numpy array of spectrum and integrated spectral radiation
#     """
#
#     spec, integration, groups, trad, totvis = self.calc_all_spectral(solar_irradiance)
#     # return matching wavelengths and their intensities, wl=spec[:,0], global irradiance=spec[:,4]
#     return spec[:, [0, 4]]


# def calc_maxes(self):
#     yearly_solarmax = 0.0
#     daily_solarmax = np.zeros(365)
#     daily_tradmax = np.zeros(365)
#
#     intmax = np.zeros(122)
#
#     spectral29_max = 0.0
#     spectral_intensity_max = [0] * 122
#     lc = self.longitudal_correction()
#     temp_rh = list()
#     tempspline = self.tempsim.get_daily_temp_spline()
#     deltaspline = self.tempsim.get_daily_deltat_spline()
#     pressure = self.pressure_kPa()
#     spectra = Spectra(latitude=self.latitude,
#                       longitude=self.longitude)
#
#     for day in range(1, 366):
#         temp = tempspline(day)
#         deltat = deltaspline(day) * 0.5
#         mintemp, maxtemp = temp - deltat, temp + deltat
#         for hour in range(8, 16):
#             solarnoon = 12 - lc - self.equation_of_time_correction(day)
#             solar_declination = self.solar_declination(day)
#             zenith_angle = self.zenith_angle(solar_declination, hour, solarnoon)
#             irradiance = self.diffuse_sky_irradiance(pressure, zenith_angle)
#             yearly_solarmax = max(yearly_solarmax, irradiance)
#             daily_solarmax[day - 1] = max(daily_solarmax[day - 1], irradiance)
#
#             halfdaylength = self.calc_half_day_length(solar_declination)
#             sunrise, sunset = solarnoon - halfdaylength, solarnoon + halfdaylength
#
#             # approximatation of sp
#             # tou is atmospheric tramsmission :
#             # overcast = 0.4 --> from Liu and Jordan (1960)
#             # clear = 0.70 --> as given in Gates (1980)
#
#             if halfdaylength < 10.5:
#                 d1 = 0 if day > 365 else day + 1
#                 if self.use_ces:
#                     # tomorrow_mintemp = temps[d1][0]
#                     temp_rh.extend(self.fit_temp_rh(mintemp, maxtemp, tomorrow_mintemp, sunrise, sunset))
#                 else:
#                     temp_rh.extend(self.fit_sine_temp_rh(mintemp, maxtemp))
#             else:
#                 temp_rh.extend(self.fit_sine_temp_rh(mintemp, maxtemp))
#
#             air_temp = (mintemp + maxtemp) / 2
#
#             spec = spectra.calc_all_spectral()
#             spectra.longitude = self.longitude
#             spectra.latitude = self.latitude
#             spectra.units = 1
#             spectra.watvap = 1
#             spectra.press = pressure
#             spectra.year = 2012
#             spectra.dayofyear = day
#             spectra.tau500 = self.calc_tao(day)
#             spectra.temp = temp_rh[-1][0]
#             spectra.timezone = -11
#
#             spectra.hour = hour
#             spec, integration, groups, trad, totvis = spectra.calc_all_spectral(irradiance)
#             daily_tradmax[day - 1] = max(daily_tradmax[day - 1], trad)
#             np.maximum(intmax, integration[:, 2], intmax)
#     # return daily_solarmax, daily_tradmax, intmax
