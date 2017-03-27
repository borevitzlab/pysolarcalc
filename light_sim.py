import numpy as np
import datetime
from world_temp_sim import TempSim
from functools import wraps
from Spectra import Spectra
iso8601 = "%Y-%m-%dT%H:%M:%S"

DEGREES_TO_RADIANS = np.pi / 180.0
RADIANS_TO_DEGREES = 180.0 / np.pi
# solar constant?
SPO = 1360.0

class LightSim(object):
    def __init__(self,
                 latitude: float = 0.0, longitude: float = 90.0, elevation: float = 0.0,
                 number_of_wls: int = 10,
                 cesfit: bool = False):
        self.output_list = ["datetime", "simtime", "temp", "relativehumidity",
                            *["LED{}".format(x) for x in range(number_of_wls)]]
        self.latitude = latitude
        self.longitude = longitude
        self.latrad = self.latitude * DEGREES_TO_RADIANS
        self.lonrad = self.longitude * DEGREES_TO_RADIANS
        self.altitude = elevation
        self.tempsim = TempSim()
        self.temps = self.tempsim.get_daily_temps(self.latitude, self.longitude)

        self.rain = np.array([])

        # this is for the climate change fitting from climateAPI
        # TODO: write the climateapi and solve that problem.
        self.use_ces = cesfit

    def read_weather_file(self, fn):
        """
        fillout to read csv of rain, minairtemp, maxairtemp

        assign to the tempsims variables.

        :param fn: filename
        :return:
        """
        pass

    def main_calc(self, start, end):
        try:
            start = datetime.datetime.strptime(start, iso8601)
            end = datetime.datetime.strptime(end, iso8601)
        except Exception as e:
            print("main_calc, given non-iso8601 start/end dates, format: {}".format(iso8601))

    def calc_minutes(self):
        pass

    def calc_maxes(self):
        yearly_solarmax = 0.0
        daily_solarmax = np.zeros(365)
        daily_tradmax = np.zeros(365)

        intmax = np.zeros(122)

        spectral29_max = 0.0
        spectral_intensity_max = [0] * 122
        lc = self.longitudal_correction()
        temp_rh = list()
        for day in range(1, 366):
            for hour in range(8, 16):
                solarnoon = 12 - lc - self.equation_of_time_correction(day)
                solar_declination = self.solar_declination(day)
                zenith_angle = self.zenith_angle(solar_declination, hour, solarnoon)
                pressure = self.pressure()
                irradiance = self.diffuse_sky_irradiance(pressure, zenith_angle)
                yearly_solarmax = max(yearly_solarmax, irradiance)
                daily_solarmax[day-1] = max(daily_solarmax[day-1], irradiance)

                halfdaylength = self.calc_half_day_length(solar_declination)
                sunrise, sunset = solarnoon - halfdaylength, solarnoon + halfdaylength

                # approximatation of sp
                # tou is atmospheric tramsmission :
                # overcast = 0.4 --> from Liu and Jordan (1960)
                # clear = 0.70 --> as given in Gates (1980)
                mintemp, maxtemp = self.tempsim.avgdailyt[day-1]
                if halfdaylength < 10.5:
                    d1 = 0 if day > 365 else day + 1
                    if self.use_ces:
                        tomorrow_mintemp = self.tempsim.avgdailyt[d1][0]
                        temp_rh.extend(self.fit_temp_rh(mintemp, maxtemp, tomorrow_mintemp, sunrise, sunset))
                    else:
                        temp_rh.extend(self.fit_sine_temp_rh(mintemp, maxtemp))
                else:
                    temp_rh.extend(self.fit_sine_temp_rh(mintemp, maxtemp))
                tao = self.calc_tao(day)
                air_temp = (mintemp + maxtemp) / 2

                spectra = Spectra()
                spectra.longitude = self.longitude
                spectra.latitude = self.latitude
                spectra.units = 1
                spectra.watvap = 1
                spectra.press = pressure

                spectra.year = 2012
                spectra.dayofyear = day
                spectra.temp = temp_rh[-1][0]

                spectra.hour = hour
                spec, integration, groups, trad, totvis = spectra.calc_all_spectral(irradiance)
                daily_tradmax[day-1] = max(daily_tradmax[day-1], trad)
                np.maximum(intmax, integration[:, 2], intmax)
        return daily_solarmax, daily_tradmax, intmax


    def fit_sine_temp_rh(self, maxtemp: float, mintemp: float) -> list:
        """
        does the same as fit_temp_rh, but uses a sine wave rather than the solar days.

        :param mintemp: minimum temperature for the day
        :param maxtemp: maximum temperature for the day
        :return:
        """
        avg = (mintemp + maxtemp) / 2
        amplitude = (maxtemp - mintemp) / 2.0
        temp_rh = list()
        for t in range(24 * 60 + 1):
            # assume daily fluctuation mimics sinewave...period 24 h
            # freq = 2 * pi / 86400 sec = 7.27 E-5 sec-1
            temp = avg - amplitude * np.sin(t * 7.27E-5 * 60.0)
            rh = self.relative_humidity(mintemp, temp)
            temp_rh.append((temp, rh))
        return temp_rh

    def fit_temp_rh(self, mintemp: float, maxtemp: float, tomorrow_min: float, sunrise: float, sunset: float) -> list:
        """
        fits temperature and relative humidiy for a day given the maxtemp, mintemp and mintemp2
        also needs sunset and sunrise
        From Cesaraccio et al. 2001 Int.J. Biometeorol. 161-169

        :param mintemp: minimum temperature for the day
        :param maxtemp: maximum temperature for the day
        :param tomorrow_min: tomorrows minimum temperature
        :param sunset: sunset time (solar time?)
        :param sunrise: sunrise time (solar time?)
        :return: list of tuples with temperature,relative_humidity
        :rtype: list(tuple[float,float])
        """
        alpha = maxtemp - mintemp
        to = maxtemp - 0.39 * (maxtemp - tomorrow_min)
        r = maxtemp - to
        hx = sunset - 4.0
        hp = sunrise + 24.0
        b = (tomorrow_min - to) / np.sqrt(hp - sunset)
        t = 0.0
        temp_rh = list()
        for hr in range(0, 25):
            for min in range(0, 60):
                time = hr + min / 60.0
                if hx >= time > sunrise:
                    t = mintemp + alpha * np.sin((time - sunrise) / (hx - sunrise) * np.pi / 2.0)
                if hx < time < sunset:
                    t = to + r * np.sin(np.pi / 2.0 + (time - hx) / 8.0 * np.pi)
                if sunset < time <= hp:
                    t = to + b * np.sqrt(time - sunset)
                if time <= sunrise:
                    t = to + b * np.sqrt(time + 24.0 - sunset)
                r = self.relative_humidity(mintemp, t)
                temp_rh.append((t, r))
        return temp_rh

    def calc_tao(self, dayofyear: int) -> float:
        """
        calculates tao from cloud cover (rain) and deltaT values.

        :param dayofyear: day of the year to calculate for.
        :return: tao value
        :rtype: float
        """
        tao = 0.70
        cur_deltat = self.tempsim.deltat[dayofyear]
        # only do cloudcover if
        if len(self.rain):
            cur = self.rain[dayofyear]
            # assume raining is overcast
            if cur:
                tao = 0.4

            # avoid going out of bounds
            if dayofyear != 0:
                # if its been raining for two days then denser cloud cover
                prev = self.rain[dayofyear - 1]
                if cur and prev:
                    tao = 0.3

                # assign pre-rain days to 80% of tao value ?
                # todo: what does this mean? should this be done better?
                if (not cur) and prev:
                    tao = 0.6

        # if air temperature rise is less than 10 --> lower tao value
        # unless near poles.
        if abs(self.latitude / np.pi * 180) < 60:
            if cur_deltat <= 10 and cur_deltat != 0:
                tao /= 11.0 - cur_deltat

        return tao

    @staticmethod
    def relative_humidity(mintemp: float, temp: float) -> float:
        """
        calculates saturated water vapor concentration with a minimum temperature for the day and temperature.
        from FAO Allen, 1985 ??

        :param mintemp: minumum temperature (C) for the day
        :param temp: current temperature (C)
        :return: estimated relative humidity
        :rtype: float
        """

        # these values shouldnt be under 0
        eo = max(0.6108 * np.exp((17.27 * temp) / (temp + 273.3)), 0.0)
        eos = max(0.6108 * np.exp((17.27 * mintemp) / (mintemp + 273.3)), 0.0)
        # rh shouldnt be over 1.0
        rh = min(eos / eo, 1.0)
        return rh * 100.0

    def pressure(self):
        return 101.0 * np.exp(-1 * self.altitude / 8200.0)

    def diffuse_sky_irradiance(self, pressure: float, zenith_angle: float) -> float:
        """
        Diffuse sky irradiance on horizontal plane (Sd)
        calc using tao, m , and zangle.
        formula given in Campbell and Norman (1998)

        :param zenith_angle: zenith angle
        :param pressure: air pressure
        :return: diffuse sky irradiance
        """
        tao = 0.7
        m = pressure / 101.3 / np.cos(zenith_angle)
        sp = SPO * pow(tao, m)
        sd = 0.3 * (1.0 - pow(tao, m)) * np.cos(zenith_angle) * SPO
        # beam irradiance on a horizontal surface
        sb = sp * np.cos(zenith_angle)
        return sb + sd

    def calc_half_day_length(self, solar_declination: float) -> float:
        """
        Calculates 1/2 solar day length

        :param solar_declination: solar declination
        :return: half a solar day
        :rtype: float
        """
        v0 = (np.cos(90.0 * np.pi / 180.0) - np.sin(self.latrad) * np.sin(solar_declination))
        v1 = (np.cos(solar_declination) * np.cos(self.latrad))
        return (np.arccos(v0 / v1) * 180.0 / np.pi) / 15.0

    def zenith_angle(self, solar_declination: float, time: float, solarnoon: float) -> float:
        """
        calculates the zenith angle with regards to the solar declination, time and solar noon

        :param solar_declination:
        :param time: time for zenith
        :param solarnoon: solar noon value
        :return: zenith angle
        :rtype: float
        """

        v0 = np.sin(self.latrad) * np.sin(solar_declination)
        v1 = np.cos(self.latrad) * np.cos(solar_declination)
        v2 = np.cos(15 * (time - solarnoon) * DEGREES_TO_RADIANS)

        return np.arccos(v0 + v1 * v2)

    def longitudal_correction(self):
        """
        gets the longitudal correction
        assumes the self.longrad is in decimal format
        translates to 4 minutes/degree

        :return: longitudal correction value
        :rtype: float
        """
        return self.lonrad / 360 * 24

    @staticmethod
    def equation_of_time_correction(dayofyear: int) -> float:
        """
        calculates equation of time correction typically a 15-20 minute correct depending on calendar day inputs

        :param dayofyear: day of the year
        :type: int
        :return: eot correction
        :rtype: float
        """
        et_calc = (279.575 + 0.9856 * dayofyear) * DEGREES_TO_RADIANS
        v0 = -104.7 * np.sin(et_calc) + 596.2 * np.sin(et_calc * 2) + 4.3
        v1 = np.sin(3.0 * et_calc) + -12.7
        v2 = np.sin(4.0 * et_calc) + -429.3 * np.cos(et_calc)
        v3 = -2.0 * np.cos(2 * et_calc) + 19.3 * np.cos(3.0 * et_calc)
        return (v0 * v1 * v2 + v3) / 3600.0

    @staticmethod
    def solar_declination(dayofyear: int) -> float:
        """
        Calculate solar declination angle...
        formula from Campbell and Norman, 1998 [Eq. 11.2]
        corrected for day of year (Jan. 1 = 1, etc.)

        :param dayofyear: day of the year
        :type: int
        :return: solar declination angle
        :rtype: float
        """
        v1 = np.sin((356.6 + 0.9856 * dayofyear) * DEGREES_TO_RADIANS)
        v2 = np.sin((278.97 + 0.9856 * dayofyear + 1.9165 * v1) * DEGREES_TO_RADIANS)
        return np.arcsin(0.39785 * v2)

    @staticmethod
    def stefan_boltzman_radiation_watts(airtemp: float) -> float:
        """
        Stefan-Boltzman law returns black body radiation in units of Watts/m2 emitted from body.

        :param airtemp: air temperature
        :return: watts/m2 radiation value
        :rtype: float
        """
        return 5.67E-08 * pow((airtemp + 273.16), 4)
