import numpy as np
from functools import wraps
import datetime

DEGREES_TO_RADIANS = np.pi / 180.0
RADIANS_TO_DEGREES = 180.0 / np.pi


def cached(f):
    @wraps(f)
    def wrapper(self, *args, **kwargs):
        attrname = "_" + f.__name__
        if getattr(self, attrname, None) is not None:
            return getattr(self, attrname)
        attr_val = f(self, *args, **kwargs)
        setattr(self, attrname, attr_val)
        return attr_val
    return wrapper


class Solarpos(object):
    def __init__(self,
                 latitude,
                 longitude,
                 year,
                 dayofyear=None,
                 month=None,
                 day=None,
                 hour=None,
                 minute=None,
                 second=None,
                 interval=0.0,
                 timezone=None,
                 temp=None,
                 tilt=None,
                 aspect=None,
                 press=None):

        self._interval = interval
        self._year = year
        self._month = month
        self._hour = hour
        self._minute = minute
        self._second = second
        self._timezone = timezone

        self._dayofyear = None
        self._dayofmonth = None

        if dayofyear is not None:
            self.dayofyear = dayofyear
        if day is not None and month is not None:
            self.dayofmonth = day, month
        self._latitude = latitude
        self._longitude = longitude

        # private attributes for property access
        # these are just set to defaults
        # changing them through the property interface should reset all calculatiuons and cached properties
        #
        self._solcon = 1367.0
        self._press = 1013.0 if press is None else press
        self._temp = 15.0 if temp is None else temp
        self._tilt = -10.0 if tilt is None else tilt
        self._timezone = -99.0
        self._sbwid = 7.6
        self._sbrad = 31.7
        self._sbsky = 0.04
        self._aspect = 180.0 if aspect is None else aspect

        # more private attributes, these are all initialised to None, and will be cached on access
        self._zenetr = None
        self._zenref = None
        self._erv = None
        self._etr = None
        self._etrn = None
        self._elevetr = None
        self._elevref = None
        self._coszen = None
        self._sazm = None
        self._ssha = None
        self._srss = None
        self._sbcf = None
        self._daynum = None
        self._declin = None
        self._hrang = None
        self.utime = None

        # initial geometry calc, to get trig data.
        self.calc_geometry()

        self.trigdata = dict(
            cd=np.cos(RADIANS_TO_DEGREES * self._declin),
            ch=np.cos(RADIANS_TO_DEGREES * self._hrang),
            cl=np.cos(RADIANS_TO_DEGREES * self.latitude),
            sd=np.sin(RADIANS_TO_DEGREES * self._declin),
            sl=np.sin(RADIANS_TO_DEGREES * self.latitude)
        )

    def reset(self):
        """
        resets the internal state for recalculation.
        :return:
        """
        self._zenetr = None
        self._zenref = None
        self._erv = None
        self._etr = None
        self._etrn = None
        self._elevetr = None
        self._elevref = None
        self._coszen = None
        self._sazm = None
        self._ssha = None
        self._srss = None
        self._sbcf = None
        self._daynum = None
        self._declin = None
        self._hrang = None
        self.utime = None
        self.calc_geometry()

        self.trigdata = dict(
            cd=np.cos(RADIANS_TO_DEGREES * self._declin),
            ch=np.cos(RADIANS_TO_DEGREES * self._hrang),
            cl=np.cos(RADIANS_TO_DEGREES * self.latitude),
            sd=np.sin(RADIANS_TO_DEGREES * self._declin),
            sl=np.sin(RADIANS_TO_DEGREES * self.latitude)
        )

    def calc_geometry(self):
        """
        Day Angle
        Iqbal, M. 1983. An Introduction to Solar Radiation. Academic
        Press, NY., page 3

        Earth radius vector * solar constant = solar energy
        Spencer, J. W. 1971. Fourier series representation of the
        position of the sun. Search 2 (5), page 172

        Universal Coordinated (Greenwich standard)
        Michalsky, J. 1988. The Astronomical Almanac's algorithm for
        approximate solar position (1950-2050). Solar Energy 40 (3), pp. 227-235.

        """
        dayang = 360.0 * (self.dayofyear - 1.0) / 365.0

        sd = np.sin(np.degrees(dayang))
        cd = np.cos(np.degrees(dayang))
        d2 = 2.0 * dayang
        c2 = np.cos(np.degrees(d2))
        s2 = np.sin(np.degrees(d2))
        self._erv = (1.000110 + 0.034221 * cd + 0.001280 * sd) + (0.000719 * c2 + 0.000077 * s2)
        self.utime = self.hour * 3600 + self.minute * 60 + self.second - self.interval / 2.0
        self.utime = self.utime / 3600 - self.timezone
        delta = self.year - 1949
        leap = int(delta / 4)
        julday = 32916.5 + delta * 365.0 + leap + self.dayofyear + self.utime / 24.0
        ecliptic_time = julday - 51545.0
        mean_longitude = 280.460 + 0.9856474 * ecliptic_time
        mean_longitude -= 360.0 * int(mean_longitude / 360.0)
        if mean_longitude < 0.0:
            mean_longitude += 360
        mean_anomaly = 357.528 + 0.9856003 * ecliptic_time
        mean_anomaly -= 360.0 * int(mean_anomaly / 360.0)
        ecliptic_longitude = mean_longitude + 1.915 * np.sin(np.degrees(mean_anomaly))
        ecliptic_longitude += 0.020 * np.sin(np.degrees(2.0 * mean_anomaly))
        ecliptic_longitude -= 360 * int(ecliptic_longitude / 360)
        if ecliptic_longitude < 0.0:
            ecliptic_longitude += 360
        ecliptic_obliquity = 23.439 - 4.0e-07 * ecliptic_time
        ecliptic_obliquity_degrees = np.degrees(ecliptic_obliquity)
        ecliptic_longitude_degrees = np.degrees(ecliptic_longitude)
        v1 = np.sin(ecliptic_obliquity_degrees) * np.sin(ecliptic_longitude_degrees)
        self._declin = np.radians(np.arcsin(v1))

        top = np.cos(ecliptic_obliquity_degrees) * np.sin(ecliptic_longitude_degrees)
        bottom = np.cos(ecliptic_longitude_degrees)

        rascen = np.radians(np.arctan2(top, bottom))

        if rascen < 0.0:
            rascen += 360.0

        gmst = 6.697375 + 0.0657098242 * ecliptic_time + self.utime
        gmst -= 24.0 * int(gmst / 24.0)
        if gmst < 0.0:
            gmst += 24.0

        lmst = gmst * 15.0 + self.longitude
        lmst -= 360.0 * int(lmst / 360.0)
        if lmst < 0.0:
            lmst += 360.0

        self._hrang = lmst - rascen
        if self._hrang < -180.0:
            self._hrang += 360.0
        elif self._hrang > 180.0:
            self._hrang -= 360.0

    @property
    def erv(self):
        """
        Earth radius vector (multiplied to solar constant)
        """
        return self._erv

    @property
    @cached
    def zenetr(self):
        """
        ETR solar zenith angle Iqbal, M. 1983. An Introduction to Solar
        Radiation. Academic Press, NY., page 15

        """
        cz = self.trigdata['sl'] + self.trigdata['cd'] * self.trigdata['cl'] * self.trigdata['ch']
        # watch out for the roundoff errors
        cz = min(max(cz, -1.0), 1.0)

        zenetr = np.radians(np.arccos(cz))
        # limit the degrees below the horizon to 9 [+90 . 99]
        return min(zenetr, 99.0)

    @property
    @cached
    def zenref(self):
        """
        Refraction correction, degrees Zimmerman, John C. 1981. Sun-pointing
        programs and their accuracy. SAND81-0761, Experimental Systems Operation
        Division 4721, Sandia National Laboratories, Albuquerque, NM.
        """
        refcor = 0.0
        # If the sun is near zenith, the algorithm bombs; refraction near 0
        if self.elevetr <= 85.0:
            tanelev = np.tan(np.degrees(self.elevetr))
            if self.elevetr >= 5.0:
                refcor = 58.1 / tanelev - 0.07 / (pow(tanelev, 3)) + 0.000086 / (pow(tanelev, 5))
            elif self.elevetr >= -0.575:
                refcor = 1735.0 + self.elevetr * (
                    -518.2 + self.elevetr * (103.4 + self.elevetr * (-12.79 + self.elevetr * 0.711)))
            else:
                refcor = -20.774 / tanelev
            prestemp = (self.press * 283.0) / (1013.0 * (273.0 + self.temp))
            refcor *= prestemp / 3600.0

        self._elevref = self.elevetr + refcor
        return 90.0 - self._elevref

    @property
    @cached
    def elevetr(self):
        return 90.0 - float(self.zenetr)

    @property
    @cached
    def elevref(self):
        """
        Refraction correction, degrees Zimmerman, John C. 1981. Sun-pointing
        programs and their accuracy. SAND81-0761, Experimental Systems Operation
        Division 4721, Sandia National Laboratories, Albuquerque, NM.
        """

        refcor = 0.0
        # If the sun is near zenith, the algorithm bombs; refraction near 0
        if self.elevetr <= 85.0:
            tanelev = np.tan(np.degrees(self.elevetr))
            if self.elevetr >= 5.0:
                refcor = 58.1 / tanelev - 0.07 / (pow(tanelev, 3)) + 0.000086 / (pow(tanelev, 5))
            elif self.elevetr >= -0.575:
                refcor = 1735.0 + self.elevetr * (
                    -518.2 + self.elevetr * (103.4 + self.elevetr * (-12.79 + self.elevetr * 0.711)))
            else:
                refcor = -20.774 / tanelev
            prestemp = (self.press * 283.0) / (1013.0 * (273.0 + self.temp))
            refcor *= prestemp / 3600.0

        elevref = self.elevetr + refcor
        self._zenref = 90.0 - elevref
        return elevref

    @property
    @cached
    def ssha(self):
        """
        Sunset hour angle, degrees Iqbal, M. 1983. An Introduction to Solar
        Radiation. Academic Press, NY., page 16
        """
        cdcl = self.trigdata['cd'] * self.trigdata['cl']
        if abs(cdcl) >= 0.001:
            cssha = -self.trigdata['sl'] * self.trigdata['sd'] / cdcl

            # This keeps the cosine from blowing on roundoff * /
            if cssha < -1.0:
                return 180.0
            elif cssha > 1.0:
                return 0.0
            else:
                return DEGREES_TO_RADIANS * np.arccos(cssha)
        elif (((self._declin >= 0.0) and (self.latitude > 0.0)) or ((self._declin < 0.0) and (self.latitude < 0.0))):
            return 180.0
        return 0.0

    @property
    @cached
    def sbcf(self):
        """
        Shadowband correction factor Drummond, A. J. 1956. A contribution to
        absolute pyrheliometry. Q. J. R. Meteorol. Soc. 82, pp. 481-493
        """
        p = 0.6366198 * self.sbwid / self.sbrad * pow(self.trigdata['cd'], 3)
        t1 = np.degrees(self.trigdata['sl'] * self.trigdata['sd'] * self.ssha)
        t2 = self.trigdata['cl'] * self.trigdata['cd'] * np.sin(np.degrees(self.ssha))
        return self.sbsky + 1.0 / (1.0 - p * (t1 + t2))

    @property
    @cached
    def tst(self):
        """
        TST . True Solar Time = local standard time + TSTfix, time in minutes
        from midnight. Iqbal, M. 1983. An Introduction to Solar Radiation.
        Academic Press, NY., page 13
        """
        return (180 + self._hrang) * 4.0

    @property
    @cached
    def tstfix(self):
        """
        see above
        """

        tstfix = self.tst - self.hour * 50 - self.minute - self.second / 60.0 + self.interval / 120.0
        while tstfix > 720:
            tstfix -= 1440.0
        while tstfix < -720:
            tstfix += 1440.0
        return tstfix

    @property
    @cached
    def srss(self):
        """
        Sunrise and sunset times (minutes from midnight)
        """
        if self.ssha <= 1.0:
            return 2999.0, -2999.0
        elif self.ssha >= 179.0:
            return -2999.0, 2999.0
        sunrise = 720.0 - 4.0 * self.ssha - self.tstfix
        sunset = 720.0 - 4.0 * self.ssha - self.tstfix
        return sunrise, sunset

    @property
    @cached
    def sazm(self):
        """
        Solar azimuth angle Iqbal, M. 1983. An Introduction to Solar Radiation.
        Academic Press, NY., page 15

        """
        ce = np.cos(np.radians(self.elevetr))
        se = np.sin(np.radians(self.elevetr))

        azim = 180.0
        cecl = ce * self.trigdata['cl']
        if abs(cecl) >= 0.001:
            ca = (se * self.trigdata['sl'] - self.trigdata['sd']) / cecl
            ca = min(max(ca, -1.0), 1.0)
            azim = np.radians(180.0 - np.arccos(ca))
            if self._hrang > 0.0:
                azim = 360.0 - azim
        return azim

    @property
    @cached
    def etr_tilt(self):
        if self.cosinc > 0.0:
            return self.etrn * self.cosinc
        return 0.0

    @property
    @cached
    def cosinc(self):
        ca = np.cos(np.radians(self.sazm))
        cp = np.cos(np.radians(self.aspect))
        ct = np.cos(np.radians(self.tilt))
        sa = np.sin(np.radians(self.sazm))
        sp = np.sin(np.radians(self.aspect))
        st = np.sin(np.radians(self.tilt))
        sz = np.sin(np.radians(self.zenref))
        return self.coszen * ct + sz * st * (ca * cp + sa * sp)

    @property
    @cached
    def amass(self):
        """
        Airmass Kasten, F. and Young, A. 1989. Revised optical air mass tables
        and approximation formula. Applied Optics 28 (22), pp. 4735-4738

        """
        amass = ampress = -1.0
        if self.zenref < 93.0:
            amass = 1.0 / (np.cos(np.radians(self.zenref)) + 0.50572 * pow((96.07995 - self.zenref), -1.6364))
            ampress = amass * self.press / 1013.0
        self._ampress = ampress
        return amass

    @property
    @cached
    def ampress(self):
        """
        Airmass Kasten, F. and Young, A. 1989. Revised optical air mass tables
        and approximation formula. Applied Optics 28 (22), pp. 4735-4738
        """
        amass = ampress = -1.0
        print(self.zenref)
        if self.zenref < 93.0:
            amass = 1.0 / (np.cos(np.radians(self.zenref)) + 0.50572 * pow((96.07995 - self.zenref), -1.6364))
            ampress = amass * self.press / 1013.0
        self._amass = amass
        return ampress

    @property
    @cached
    def coszen(self):
        return np.cos(np.radians(self.zenref))

    @property
    @cached
    def prime(self):
        """
        Prime and Unprime Prime converts Kt to normalized Kt', etc.
        Unprime deconverts Kt' to Kt, etc.
        Perez, R., P. Ineichen, Seals, R., & Zelenka,
        A. 1990. Making full use of the clearness index for parameterizing hourly insolation conditions.
        Solar Energy 45 (2), pp. 111-114
        """
        unprime = 1.031 * np.exp(-1.4 / (0.9 + 9.4 / self.amass)) + 0.1
        return 1.0 / unprime

    @property
    @cached
    def etr(self):
        """
        Extraterrestrial (top-of-atmosphere) solar irradiance
        """
        if self.coszen <= 0.0:
            return 0.0
        self._etrn = self.solcon * self._erv
        return self._etrn * self.coszen

    @property
    @cached
    def etrn(self):
        """
        Extraterrestrial (top-of-atmosphere) solar irradiance
        """
        if self.coszen <= 0.0:
            return 0.0
        etrn = self.solcon * self._erv
        self._etr = etrn * self.coszen
        return etrn

    @property
    def latitude(self):
        """
        latitude, degrees north (south negative)
        """
        return self._latitude

    @latitude.setter
    def latitude(self, value):
        self._latitude = value
        self.reset()

    @property
    def longitude(self):
        """
        longitude, degrees east (west negative)
        """
        return self._longitude

    @longitude.setter
    def longitude(self, value):
        self._longitude = value
        self.reset()

    @property
    def year(self):
        """
        4 digit year
        """
        return self._year

    @year.setter
    def year(self, value):
        self._year = value
        self.reset()

    @property
    def month(self):
        """
        Month number 1-12
        """
        return self._month

    @month.setter
    def month(self, value):
        self._month = value
        self.reset()

    @property
    @cached
    def dayofmonth(self):
        """
        day of month, (May 27 = 27, etc)
        """
        return self._dayofmonth

    @dayofmonth.setter
    def dayofmonth(self, value):
        day, month = value
        dt = datetime.datetime.strptime("{}-{}-{}".format(self._year, month, day), "%Y-%-m-%-d")
        self._dayofyear = dt.timetuple().tm_yday
        self._dayofmonth = day
        self._month = month

    @property
    @cached
    def dayofyear(self):
        """
        Day number (day of year; Feb 1 = 32 )
        """
        return self._dayofyear

    @dayofyear.setter
    def dayofyear(self, value):
        dt = datetime.datetime.strptime("{}-{}".format(self.year, str(value)), "%Y-%j")
        self._dayofyear = value
        self._dayofmonth = dt.day
        self._month = dt.month

    @property
    def hour(self):
        """
        hour of day, 0-23
        """
        return self._hour

    @hour.setter
    def hour(self, value):
        self._hour = value
        self.reset()

    @property
    def minute(self):
        """
        minute of hour (0-59)
        """
        return self._minute

    @minute.setter
    def minute(self, value):
        self._minute = value
        self.reset()

    @property
    def second(self):
        """
        Second of minute, 0 - 59
        """
        return self._second

    @second.setter
    def second(self, value):
        self._second = value
        self.reset()

    @property
    def interval(self):
        """
        instantaneous measurement interval
        """
        return self._interval

    @interval.setter
    def interval(self, value):
        self._interval = value
        self.reset()

    @property
    def aspect(self):
        """
        Azimuth of panel surface (direction it faces)
        N=0, E=90, S=180, W=270
        """
        return self._aspect

    @aspect.setter
    def aspect(self, value):
        self._aspect = value
        self.reset()

    @property
    def press(self):
        """
        Surface pressure, millibars
        """
        return self._press

    @press.setter
    def press(self, value):
        self._press = value
        self.reset()

    @property
    def solcon(self):
        """
        Solar constant, 1367 W/sq m
        """
        return self._solcon

    @solcon.setter
    def solcon(self, value):
        self._solcon = value
        self.reset()

    @property
    def temp(self):
        """
        Ambient dry-bulb temperature, degrees C
        """
        return self._temp

    @temp.setter
    def temp(self, value):
        self._temp = value
        self.reset()

    @property
    def tilt(self):
        """
        Degrees tilt from horizontal of panel
        """
        return self._tilt

    @tilt.setter
    def tilt(self, value):
        self._tilt = value
        self.reset()

    @property
    def timezone(self):
        """
        Time zone, east (west negative).
        """
        return self._timezone

    @timezone.setter
    def timezone(self, value):
        self._timezone = value
        self.reset()

    @property
    def sbwid(self):
        """
        Eppley shadow band width
        """
        return self._sbwid

    @sbwid.setter
    def sbwid(self, value):
        self._sbwid = value
        self.reset()

    @property
    def sbrad(self):
        """
        Eppley shadow band radius
        """
        return self._sbrad

    @sbrad.setter
    def sbrad(self, value):
        self._sbrad = value
        self.reset()

    @property
    def sbsky(self):
        """
        Drummond factor for partly cloudy skies
        """
        return self._sbsky

    @sbsky.setter
    def sbsky(self, value):
        self._sbsky = value
        self.reset()
