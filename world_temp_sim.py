from scipy.interpolate import CubicSpline
import gzip
import numpy as np
import datetime


def get_doy(m):
    """
    gets the day of the year of the middle of the month m except in special circumstances:

     returns 0 if m is 0
     returns the doy/m + 31 * number of months past 12

    todo: clean this so that it handles out of range values more elegantly.

    :param m: month number, 0 - 13
    """
    midmonth = 14.5
    leapday = 0.25
    if m == 0:
        return 0
    mo = min(max(m, 0), 12)
    dy = datetime.datetime.strptime("{}1".format(mo), "%m%d").timetuple().tm_yday
    dy += midmonth + (leapday if not mo == 2 else 0)
    dy += 31 * (m % 12) * (m // 12)
    return dy


class TempSim(object):
    def __init__(self, latitude, longitude, climate_change: bool = False, temperature_offsets: list = [0.0] * 12):
        self._climate_change = False
        self.tempoff = temperature_offsets
        self.lat = latitude
        self.lon = longitude
        self.avgtemp = np.array([])
        self.avgdailyt = np.array([])
        self.diurnaltemp = TempSim.read_deltat()

        self.year_avg = \
            self._min_temp = \
            self._max_temp = None

    @staticmethod
    def read_deltat():
        """
        reads temperature deltas from the "deltat.gz" data file.
        temperatures are in degrees centigrade * 10.0 apparently.
        the data file has headers that should be length 8 and contain the following info headers:
        grid_sz, xmin, ymin, xmax, ymax, n_cols, n_rows, n_months

        the most important of these are n_months, n_rows, and n_cols

        :return: array of shape (nrows, ncols, n_months) of temperature deltas
        :rtype: np.array
        """
        fname = "data/deltat.gz"
        diurnaltemp = list()

        ar1 = np.genfromtxt(fname, delimiter=5, skip_header=2)
        return ar1.reshape((12, 360, 720)).swapaxes(0, 2)
        # with gzip.open(fname, "rb") as f:
        #     lines = f.read().decode("utf-8").split("\n")[1:]
        #     header_data = lines.pop(0).split()
        #     assert len(header_data) == 8, "Header not 8 long: {}".format(header_data)
        #     grid_sz, xmin, ymin, xmax, ymax, n_cols, n_rows, n_months = header_data
        #     lness = list()
        #     for mth in range(0, int(n_months)):
        #         mthlist = list()
        #         for hh in range(0, int(n_rows)):
        #             line = lines[mth * hh]
        #             mthlist.append([float(line[i:i + 5]) for i in range(0, len(line), 5)])
        #             assert len(mthlist[-1]) == int(n_cols), str(mthlist[-1])
        #         lness.append(mthlist)
        #     # clear the list...
        #     for x in range(int(n_rows)):
        #         diurnaltemp.append(list())
        #         for y in range(int(n_cols)):
        #             mlist = [lness[m][x][y] for m in range(12)]
        #             diurnaltemp[-1].append(mlist)

        # return np.array(diurnaltemp)


    @staticmethod
    def read_avgt(lat2: float, lon2: float):
        """
        takes a latitude and a longitude and looks up the monthly average temperature for a year from the
        "avgtemp.gz" data file.
        The temp file should contain lat,long,jan,feb,mar,apr,may,jun,jul,sep,oct,nov,dec,yearlyavg

        :param lat2: latitude
        :param lon2: longitude
        :return: array containing a 12 length array of temperatures for jan-dec and average temperature for the year
        :rtype: np.array, float
        """
        fname = "data/avgtemp.gz"
        ar1 = np.genfromtxt(fname)

        def find_nearest_vector(array, value):
            idx = np.array([np.linalg.norm(x + y) for (x, y) in np.abs(array[:, :2] - value)]).argmin()
            return array[idx]

        v = find_nearest_vector(ar1, [lon2, lat2])

        return v[2:-1], v[-1]
        #
        # with gzip.open(fname, "rb") as f:
        #     data = [[float(z) for z in x.split()] for x in f.read().decode("utf-8").split("\n")]
        #     data.pop()
        #     data = np.array(data)
        #
        #     def match(v):
        #         return abs(float(v[0]) - lon2) < 0.5 and abs(float(v[1]) - lat2) < 0.5
        #
        #     matching = next((x[2:] for x in data if match(x)), None)
        #     if matching is None:
        #         raise ValueError("{},{} not found".format(lat2, lon2))
        #
        # return np.array(matching[:-1]), float(matching[-1])

    def get_daily_temps(self, orglon: float = None):
        """
        gets the daily temperature ranges for a year.

        :param orglon: -180 to 180 longitude
        :return: 366,2 shape array of temperature min,max temperatures for each day.
        """

        lat2 = int((-2.0 * self.lat) + 180.5)
        lon2 = int((self.lon * 2) + 0.25)
        lt = (-0.5 * lat2) + 90.25
        ln = (0.5 * lon2) - 0.25

        self.deltat = self._calc_spline_delta_t(lat2, lon2)
        # not sure why we use different.
        self.avgdailyt = self._calc_spline_avg_t(lt, ln)

        temps = np.zeros((366, 2))
        for day in range(temps.shape[1]):
            temps[day][0] = self.avgdailyt[day] - 0.5 * self.deltat[day]
            temps[day][1] = self.avgdailyt[day] + 0.5 * self.deltat[day]

        return temps

    def get_splines(self, bc_type='not-a-knot'):
        """
        A combination fo the spline calculation methods. 
        
        Returns a spline along the day of the year for the temperature with temperature values and half deltat.
        
        :param orglon: 
        :return: temperature/deltat spline for the year
        :rtype: scipy.interpolate.CubicSpline
        """
        xx = [(m, datetime.datetime(2012, m+1, 1).timetuple().tm_yday) for m in range(0, 12)]
        ff = [[], []]

        lat2 = int((-2.0 * self.lat) + 180.5) % 360
        lon2 = int((self.lon * 2) + 0.25) % 720
        lt = (-0.5 * lat2) + 90.25
        ln = (0.5 * lon2) - 0.25
        north = self.lat >= 0

        self.avgtemp, self.year_avg = TempSim.read_avgt(self.lat, self.lon)

        for month, doy in xx:
            if float(self.diurnaltemp[lon2, lat2, month]) == -9999:
                print("Invalid Delta T -- approximating with Geerts (2002)")
                # ll = abs(lat2 + 0.001)
                ll = abs(self.lat + 0.001)
                tropics = ll < 40

                wint = 12.8 - (0.114 * ll)
                the_sum = 12.34 - (3.8 * np.cos(np.pi * ll / 35.0))

                if tropics:
                    wint = 14.31 - (3.41 * np.cos(np.pi * ll / 20.0))
                    the_sum = 12.34 - (3.80 * np.cos(np.pi * ll / 35.0))

                xxx = [0.0, 1.0, 3.5, 7.0, 13.0]
                fff = [the_sum, the_sum, (wint + the_sum) / 2.0, wint, the_sum]
                if north:
                    # northern hemisphere
                    fff = [wint, wint, (wint + the_sum) / 2.0, the_sum, wint]
                fff = list(map(abs, fff))
                monthly_delta = CubicSpline(xxx, fff, bc_type=bc_type)
                for month, doy in xx:
                    t = self.avgtemp[month % len(self.avgtemp)]
                    if self._climate_change:
                        t += self.tempoff[month % len(self.tempoff)]
                    ff[0].append(t)
                    ff[1].append(monthly_delta(month)/2)
                break
        else:
            print("Have data!!! Calulating temp.")
            for month, doy in xx:
                t = self.avgtemp[month % len(self.avgtemp)]
                if self._climate_change:
                    t += self.tempoff[month % len(self.tempoff)]
                ff[0].append(t)
                # the values in the datafile are integers multiplied by 10 (bloody fortranners)
                deltat = self.diurnaltemp[lon2, lat2, month % 12] / 20.0
                ff[1].append(abs(deltat))

        # ensure ends are identical.
        xx = list(np.array(xx)[:, 1])

        # ff = np.array(ff)
        # ff = np.swapaxes(ff, 0, 1)
        # return CubicSpline(xx, ff)
        if bc_type == 'periodic':
            xx.append(xx[-1] + 31)
            ff[0].append(ff[0][0])
            ff[1].append(ff[1][0])

        xx = np.array(xx)
        ff = np.array(ff)
        ff = np.swapaxes(ff, 0, 1)

        return CubicSpline(xx, ff, bc_type=bc_type)

    def get_daily_temp_spline(self, orglon: float = None):
        """
        calculates the average temperature for a year from the internal attributes avgtemp and tempoff
        attributes avgtemp must be set, and tempoff must be set in order to use climate change.

        :param lt: latitude
        :param ln: longitude
        :return: spline for the temperature for that year
        :rtype: scipy.interpolation.CubicSpline
        """
        xx = list()
        ff = list()
        cc = list()

        if orglon and orglon < 0:
            orglon += 360
        self.lon = self.lon if not orglon else orglon
        lat2 = int((-2.0 * self.lat) + 180.5)
        lon2 = int((self.lon * 2) + 0.25)
        lt = (-0.5 * lat2) + 90.25
        ln = (0.5 * lon2) - 0.25

        self.avgtemp, self.year_avg = TempSim.read_avgt(lt, ln)

        # we have already added month 0 and apparently we will add an extra for some reason.

        for month in range(0, 12):
            doy = datetime.datetime(2012, month+1, 1).timetuple().tm_yday
            xx.append(doy)
            ff.append(self.avgtemp[month % len(self.avgtemp)])
            cc.append(self.tempoff[month % len(self.tempoff)])
        xx.append(xx[-1]+31)
        ff.append(ff[0])
        cc.append(cc[0])
        if self._climate_change:
            return CubicSpline(xx, cc)
        return CubicSpline(xx, ff, bc_type='periodic')

    def get_daily_deltat_spline(self, orglon: float = None):
        """
        calculates the temperature delta per day for a year using the temperature delta data provided by `func:read_deltat`
        Also assigns values to those unknown if there are unknown values for the delta in the data file.

        :orglon: 
        :return: Cubic Spline of the deltat over a year, periodic boundary
        :rtype: scipy.interpolation.CubicSpline
        """
        if orglon and orglon < 0:
            orglon += 360
        self.lon = self.lon if not orglon else orglon
        lat2 = int((-2.0 * self.lat) + 180.5)
        lon2 = int((self.lon * 2) + 0.25)

        north = lat2 >= 0

        # get proper naming here....
        xx = list()
        ff = list()

        for m in range(12):
            # TODO: possible semantic bug! should lon2 and lat2 be swapped here?
            if float(self.diurnaltemp[lat2, lon2, m]) / 10.0 == -999.90:
                print("Invalid Delta T -- approximating with Geerts (2002)")
                ll = abs(lat2 + 0.001)
                tropics = ll < 40

                wint = 12.8 - (0.114 * ll)
                the_sum = 12.34 - (3.80 * np.cos(np.pi * ll / 35.0))

                if tropics:
                    wint = 14.31 - (3.41 * np.cos(np.pi * ll / 20.0))
                    the_sum = 12.34 - (3.80 * np.cos(np.pi * ll / 35.0))

                xxx = [0.0, 1.0, 3.5, 7.0, 13.0]
                fff = [the_sum, the_sum, (wint + the_sum) / 2.0, wint, the_sum]
                if north:
                    # northern hemisphere
                    fff = [wint, wint, (wint + the_sum) / 2.0, the_sum, wint]

                monthly_delta = CubicSpline(xxx, fff)
                for month in range(12):
                    doy = datetime.datetime(2012, month + 1, 1).timetuple().tm_yday
                    xx.append(doy)
                    ff.append(monthly_delta(month))
                break
        else:
            for month in range(0, 12):
                doy = datetime.datetime(2012, month+1, 1).timetuple().tm_yday
                xx.append(doy)
                # the values in the datafile are integers multiplied by 10 (bloody fortranners)
                ff.append(self.diurnaltemp[lon2, lat2, month % 12] / 10.0)

        ff.append(ff[0])
        xx.append(xx[-1]+31)
        # create the spline
        return CubicSpline(xx, ff, bc_type='periodic')

    def get_daily_temps_with_climate_change(self, lat: float, lon: float = None, orglon: float = None):
        """
        identical to above but with climate change.

        :param lat: latitcude
        :param lon: 0-360 longitude
        :param orglon: -180 to 180 longitude
        :return: 366,2 shape array of temperature min,max temperatures for each day.
        """
        if not orglon and not lon:
            print("no longitude or orglon")
            return
        if orglon and orglon < 0:
            orglon += 360
        lon = lon if lon is not None else orglon

        self._climate_change = True
        return self.get_daily_temps(lat, lon)

    def _calc_spline_avg_t(self, lt, ln):
        """
        calculates the average temperature for a year from the internal attributes avgtemp and tempoff
        attributes avgtemp must be set, and tempoff must be set in order to use climate change.

        :param lt: latitude
        :param ln: longitude
        :return: 366 length array of temperatures in degrees centigrade
        :rtype: np.ndarray
        """
        xx = list()
        ff = list()
        cc = list()

        self.avgtemp, self.year_avg = TempSim.read_avgt(lt, ln)

        if not self.tempoff:  # giant state machine oldness. making sure that tempoff indexes are defined
            self.tempoff = [0.0] * len(self.avgtemp)

        # we have already added month 0 and apparently we will add an extra for some reason.
        for month in range(0, 14):
            xx.append(get_doy(month))

            idx = month - 1
            idx = len(self.avgtemp) - 1 if month == 0 else idx

            ff.append(self.avgtemp[idx % len(self.avgtemp)])
            cc.append(self.tempoff[idx % len(self.tempoff)])

        avgt_spline = CubicSpline(xx, ff)
        climate_change_offset_spline = CubicSpline(xx, cc)

        return np.array(
            [avgt_spline(d) + (climate_change_offset_spline(d) if self._climate_change else 0.0) for d in range(366)])

    def _calc_spline_delta_t(self, lat2: int, lon2: int):
        """
        calculates the temperature delta per day for a year using the temperature delta data provided by `func:read_deltat`
        Also assigns values to those unknown if there are unknown values for the delta in the data file.

        :param lat2: latitude
        :param lon2: longitude
        :return: 366 length array of deltas in degrees centigrade
        :rtype: np.ndarray
        """
        invalid = False
        north = lat2 >= 0

        for m in range(12):
            # TODO: possible semantic bug! should lon2 and lat2 be swapped here?
            if float(self.diurnaltemp[lat2, lon2, m]) / 10.0 == -999.90:
                invalid = True
                print("Invalid Delta T -- approximating with Geerts(2002)")
                break

        if invalid:
            ll = abs(lat2 + 0.001)
            tropics = ll < 40

            wint = 12.8 - (0.114 * ll)
            the_sum = 12.34 - (3.80 * np.cos(np.pi * ll / 35.0))

            if tropics:
                wint = 14.31 - (3.41 * np.cos(np.pi * ll / 20.0))
                the_sum = 12.34 - (3.80 * np.cos(np.pi * ll / 35.0))

            xxx = [0.0, 1.0, 3.5, 7.0, 13.0]
            fff = [the_sum, the_sum, (wint + the_sum) / 2.0, wint, the_sum]
            if north:
                # northern hemisphere
                fff = [wint, wint, (wint + the_sum) / 2.0, the_sum, wint]

            monthly_delta = CubicSpline(xxx, fff)
            for x in range(12):
                y = monthly_delta(x)
                self.diurnaltemp[lon2, lat2, x] = y * 10.0

        # get proper naming here....
        xx = list()
        ff = list()

        for month in range(0, 14):
            xx.append(get_doy(month))
            ff.append(self.diurnaltemp[lon2, lat2, month % 12] / 10.0)

        # create the spline
        deltat_spline = CubicSpline(xx, ff)
        # return an array for the year.
        return np.array([deltat_spline(i) for i in range(366)])
