import numpy as np
import datetime
from pysolar import solar, util
from timezonefinder import TimezoneFinder
import pytz
from world_temp_sim import TempSim
from Spectra import Spectra

iso8601 = "%Y-%m-%dT%H:%M:%S"
from scipy.interpolate import CubicSpline
DEGREES_TO_RADIANS = np.pi / 180.0
RADIANS_TO_DEGREES = 180.0 / np.pi
# solar constant?
SPO = 1360.0

import multiprocessing

# daily sine freq
daily_freq = 2 * np.pi / 1440.0


WAVELENGTH_MICRONS = np.array(
    [0.3, 0.305, 0.31, 0.315, 0.32, 0.325, 0.33, 0.335, 0.34, 0.345, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42,
     0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.57, 0.593, 0.61, 0.63, 0.656,
     0.6676, 0.69, 0.71, 0.718, 0.7244, 0.74, 0.7525, 0.7575, 0.7625, 0.7675, 0.78, 0.8, 0.816, 0.8237, 0.8315, 0.84,
     0.86, 0.88, 0.905, 0.915, 0.925, 0.93, 0.937, 0.948, 0.965, 0.98, 0.9935, 1.04, 1.07, 1.1, 1.12, 1.13, 1.145,
     1.161, 1.17, 1.2, 1.24, 1.27, 1.29, 1.32, 1.35, 1.395, 1.4425, 1.4625, 1.477, 1.497, 1.52, 1.539, 1.558, 1.578,
     1.592, 1.61, 1.63, 1.646, 1.678, 1.74, 1.8, 1.86, 1.92, 1.96, 1.985, 2.005, 2.035, 2.065, 2.1, 2.148, 2.198, 2.27,
     2.36, 2.45, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0])
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

transmissivity_spline = CubicSpline(WAVELENGTH_MICRONS, TRANSMISSITY_COEFF)


def day_both_calc(dayvalues):
    day, temp_avg, temp_amp, temp_min, temp_delta, latitude, longitude, elevation, pressure, wavelengths = dayvalues
    spectra = Spectra(latitude=latitude,
                      longitude=longitude,
                      elevation=elevation,
                      pressure=pressure,
                      wavelengths=wavelengths)
    rvals = []
    sunrise_dt, sunset_dt = util.get_sunrise_sunset(latitude, longitude, day)
    # offset the sine so that the coldest time of the day is 30m after sunrise
    # see http://cliffmass.blogspot.com.au/2011/01/what-is-coldest-time-of-day.html
    sr_offset = 360 - (sunrise_dt.hour * 60.0 + sunrise_dt.minute + 30.0)

    for d in daterange(day, minutes=1):
        # this was wrong....
        minute = d.minute + (d.hour * 60.0)

        # temp = tavgamp * np.sin(minute * 7.27220521664304e-05 * 60.0)
        temp = temp_avg - temp_amp * np.sin((minute + sr_offset) * daily_freq)
        rh = 100 * np.exp((17.625 * temp_min) / (243.04 + temp_min)) / \
                        np.exp((17.625 * temp) / (243.04 + temp))

        tao = 0.7
        if abs(latitude / np.pi * 180) < 60:
            if temp_delta <= 10 and temp_delta != 0:
                tao /= 11.0 - temp_delta
        th = np.array([temp, rh], dtype=np.float64)
        if sunrise_dt <= d <= sunset_dt:
            rvals.append(np.append(th, spectra.calc_vis_spectral(d, temp, rh, tao)))
        else:
            rvals.append(np.append(th, np.append(0, np.zeros_like(wavelengths))))
    return rvals


def daterange(start_date: datetime.datetime, end_date: datetime.datetime = None, **kwargs):
    if not len(kwargs):
        kwargs['days'] = 1
    end_date = end_date if end_date else start_date + datetime.timedelta(days=1)
    while start_date < end_date:
        start_date += datetime.timedelta(**kwargs)
        yield start_date


visble_wavelengths = np.linspace(0.4, 1.0, 20)


class LightSim(object):
    """
    Light simulation.
    Main simulation class for the time being.

    can take a list of wavelengths in microns (float) or nanometers.
    if any of the wavelengths are greater than 50, it assumes they are in nanomenters not microns.

    """

    def __init__(self,
                 start: datetime.datetime, end: datetime.datetime,
                 latitude: float = 0.0, longitude: float = 90.0, elevation: float = 0.0,
                 wavelengths: list = visble_wavelengths,
                 max_light_intensity: float = 1000.0,
                 cesfit: bool = False):

        self.latitude = latitude
        self.longitude = longitude
        # can use nm or microns...
        if any(i > 50 for i in wavelengths):
            wavelengths = np.array(wavelengths) / 1000
        self.wavelengths = wavelengths
        self.output_list = ["datetime", "temp", "relativehumidity", 'total_etr',
                            *["{}nm".format(int(x * 1000)) for x in self.wavelengths]]


        self.max_light_intensity = max_light_intensity

        for x in range(60):
            tz = TimezoneFinder().closest_timezone_at(lat=latitude, lng=longitude, delta_degree=x)
            if tz is not None:
                print("Timezone: {}".format(tz))
                self.tz = pytz.timezone(tz)
                self.start = start.replace(tzinfo=self.tz)
                self.end = end.replace(tzinfo=self.tz)
                break
        else:
            print("No Timezone, using utc")
            self.start = start
            self.end = end

        self.latrad = self.latitude * DEGREES_TO_RADIANS
        self.lonrad = self.longitude * DEGREES_TO_RADIANS
        self.elevation = elevation

        self.tempsim = TempSim(latitude, longitude)
        self.year = 2012
        self.rain = np.array([])

        # this is for the climate change fitting from climateAPI
        # TODO: write the climateapi and solve that problem.
        self.use_ces = cesfit
        self.monthly_temperature_spline = self.tempsim.get_splines()
        # self.temp_hum_spline = self.calc_temp_humidity_spline()
        # self.spectra_spline = self.calc_spectra_spline()
        self.maxes = None
        self.combined_spline = self.main_calc()


    def read_weather_file(self, fn):
        """
        fillout to read csv of rain, minairtemp, maxairtemp

        assign to the tempsims variables.

        :param fn: filename
        :return:
        """
        pass

    def ces_temp_fit(self, m, mintemp, maxtemp, tomorrow_min, sunrise, sunset):
        alpha = maxtemp - mintemp
        to = maxtemp - 0.39 * (maxtemp - tomorrow_min)
        r = maxtemp - to
        hx = sunset - 4.0
        hp = sunrise + 24.0
        b = (tomorrow_min - to) / np.sqrt(hp - sunset)
        t = 0.0

        time = float(m) / 24.0
        if hx >= time > sunrise:
            return mintemp + alpha * np.sin((time - sunrise) / (hx - sunrise) * np.pi / 2.0)
        if hx < time < sunset:
            return to + r * np.sin(np.pi / 2.0 + (time - hx) / 8.0 * np.pi)
        if sunset < time <= hp:
            return to + b * np.sqrt(time - sunset)
        if time <= sunrise:
            return to + b * np.sqrt(time + 24.0 - sunset)
        return t

    def write_file(self, fn, **kwargs):
        float_formatter = lambda x: "%.2f" % x

        def getval(d):
            doyf = d.timetuple().tm_yday + d.minute / 1440.0 + d.hour / 24.0
            # np.append(np.around(v[:2], 1), np.around(v[2:], 0))
            return self.combined_spline(doyf)

        vgetval = np.vectorize(getval)

        dts = np.array(list(daterange(self.start, self.end, **kwargs)))
        # r = vgetval(dts)

        # np.savetxt(fn, r, delimiter=",")
        with open(fn, 'w+') as f:
            f.write(",".join(self.output_list) + "\n")
            for d in dts:
                vals = getval(d)
                th = np.around(vals[:2], 1)
                wvl = vals[2:].astype(int)
                f.write(d.isoformat() + "," + ",".join(map(str, th))+ "," + ",".join(map(str, wvl)) + "\n")

    def calc_temp_humidity_spline(self):
        """
        calculates a temperature humidity spline with deltat.

        returns a spline of [temp[:], humidity[:], deltat[:]]
        """

        d = self.start
        xx = list()
        ff = [[], []]

        while d < self.end:
            # increment time by one day
            d += datetime.timedelta(days=1)
            # calculate sr ss
            sunrise_dt, sunset_dt = util.get_sunrise_sunset(self.latitude, self.longitude, d)
            sunrise = solar.get_solar_time(self.longitude, sunrise_dt)
            sunset = solar.get_solar_time(self.longitude, sunset_dt)
            doy = d.timetuple().tm_yday

            temp, deltat = self.monthly_temperature_spline(doy)
            deltat = abs(deltat)
            mintemp, maxtemp = temp - deltat, temp + deltat
            t_amplitude = (maxtemp - mintemp) / 2.0
            t_avg = (mintemp + maxtemp) / 2.0

            if self.use_ces:
                tomorrow = d + datetime.timedelta(hours=24)
                tdoy = tomorrow.timetuple().tm_yday
                tomorrow_temp, tomorrow_deltat = self.monthly_temperature_spline(tdoy)
                tomorrow_mintemp = tomorrow_temp - abs(tomorrow_deltat)
                for m in range(1440):
                    t = self.ces_temp_fit(m, mintemp, maxtemp, tomorrow_mintemp, sunrise, sunset)
                    ff[0].append(t)
                    ff[1].append(self.relative_humidity(mintemp, t))
                    xx.append(doy + (m / 1440))
            else:
                for m in range(1440):
                    # assume daily fluctuation mimics sinewave...period 24 h
                    # freq = 2 * pi / 86400 sec = 7.27 E-5 sec-1
                    t = t_avg - t_amplitude * np.sin(m * 7.27220521664304e-05 * 60.0)
                    ff[0].append(t)
                    ff[1].append(self.relative_humidity(mintemp, t))
                    xx.append(doy + (m / 1440))

        xx.append(xx[-1] + (1 / 1440))
        ff[0].append(ff[0][0])
        ff[1].append(ff[1][0])
        ff = np.array(ff)
        ff = np.swapaxes(ff, 0, 1)
        xx = np.array(xx)
        return CubicSpline(xx, ff)

    def main_calc(self):
        d = self.start
        days = []
        xx = [x.timetuple().tm_yday + x.minute / 1440.0 + x.hour / 24.0 for x in
              daterange(self.start, self.end, minutes=1)]
        xx = np.array(xx)

        print("Calculating daily temp,hum,tao...")
        while d < self.end:
            d += datetime.timedelta(days=1)
            doy = d.timetuple().tm_yday

            temp, deltat = self.monthly_temperature_spline(doy)
            deltat = abs(deltat)

            mintemp, maxtemp = temp - deltat, temp + deltat
            t_amplitude = (maxtemp - mintemp) / 2.0
            t_avg = (mintemp + maxtemp) / 2.0
            # need to clip temps
            mintemp = min(max(1, mintemp), maxtemp)
            maxtemp = max(maxtemp, mintemp)

            days.append((d, t_avg, t_amplitude, mintemp, deltat,
                         self.latitude,
                         self.longitude,
                         self.elevation,
                         self.pressure,
                         np.array(self.wavelengths)))

        pool = multiprocessing.Pool(processes=4)

        print("Calculating spectra...")
        ff = pool.map(day_both_calc, days)
        ff = [item for sublist in ff for item in sublist]
        ff = np.array(ff)
        self.maxes = ff[np.argmax(ff, axis=0), np.arange(len(ff[0]))]

        # clip temperature and humidty
        ff[:, 0] = np.clip(ff[:, 0], 1.0, 50.0)
        ff[:, 1] = np.clip(ff[:, 1], 0, 100.0)

        # for i, wvl in enumerate(self.wavelengths):
        #     if self.maxes[3+i] != 0.0:
        #         ff[:, 3+i] = ff[:, 3+i] / self.maxes[3+i] * self.max_light_intensity * transmissivity_spline(wvl)
        # ff[:, 3:] = np.clip(ff[:, 3:], 0.0, self.max_light_intensity)

        return CubicSpline(xx, ff)

    def fit_sine_temp_rh(self, maxtemp: float, mintemp: float) -> list:
        """
        does the same as fit_temp_rh, but uses a sine wave rather than the solar days.

        :param mintemp: minimum temperature for the day
        :param maxtemp: maximum temperature for the day
        :return:
        """
        avg = (mintemp + maxtemp) / 2
        amplitude = (maxtemp - mintemp) / 2.0
        # assume daily fluctuation mimics sinewave...period 24 h
        # freq = 2 * pi / 86400 sec = 7.27 E-5 sec-1
        return [avg - amplitude * np.sin(t * 7.27E-5 * 60.0) for t in range(24 * 60 + 1)]

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
        temps = list()
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
                temps.append(t)
        return temps

    def calc_tao(self, dayofyear: int, deltat: float) -> float:
        """
        calculates tao from cloud cover (rain) and deltaT values.

        :param dayofyear: day of the year to calculate for.
        :return: tao value
        :rtype: float
        """
        tao = 0.70
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
            if deltat <= 10 and deltat != 0:
                tao /= 11.0 - deltat

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

        # # these values shouldnt be under 0
        # eo = max(0.6108 * np.exp((17.27 * temp) / (temp + 273.3)), 0.0)
        #
        # eos = max(0.6108 * np.exp((17.27 * mintemp) / (mintemp + 273.3)), 0.0)
        # # rh shouldnt be over 1.0
        # rh = min(eos / eo, 1.0)
        rh = np.exp((17.625 * mintemp) / (243.04 + mintemp)) / np.exp((17.625 * temp) / (243.04 + temp))
        return rh * 100.0

    @property
    def pressure_kPa(self):
        """
        calculates the pressure in kPa

        :rtype: float
        """
        return 0.1 * ((44331.514 - self.elevation) / 11880.516) ** (1 / 0.1902632)
        # return 101.0 * np.exp(-1 * self.altitude / 8200.0)

    @property
    def pressure(self):
        """
        calculates the pressure in Pa

        :rtype: float
        """
        return self.pressure_kPa * 1000.0

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
        sp = SPO * np.power(tao, m)
        sd = 0.3 * (1.0 - np.power(tao, m)) * np.cos(zenith_angle) * SPO
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
        v2 = np.cos(np.radians(15 * (time - solarnoon)))

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
        et_calc = np.radians(279.575 + 0.9856 * dayofyear)
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
        v1 = np.sin(np.radians(356.6 + 0.9856 * dayofyear))
        v2 = np.sin(np.radians(278.97 + 0.9856 * dayofyear + 1.9165 * v1))
        return np.arcsin(0.39785 * v2)

    @staticmethod
    def solar_declination2(dayofyear: int) -> float:
        """
        Calculate solar declination angle...
        formula from Campbell and Norman, 1998 [Eq. 11.2]
        corrected for day of year (Jan. 1 = 1, etc.)

        :param dayofyear: day of the year
        :type: int
        :return: solar declination angle
        :rtype: float
        """
        return 23.45 * np.sin((2 * np.pi / 365.0) * (dayofyear - 81))

    @staticmethod
    def stefan_boltzman_radiation_watts(airtemp: float) -> float:
        """
        Stefan-Boltzman law returns black body radiation in units of Watts/m2 emitted from body.

        :param airtemp: air temperature
        :return: watts/m2 radiation value
        :rtype: float
        """
        return 5.67E-08 * pow((airtemp + 273.16), 4)
