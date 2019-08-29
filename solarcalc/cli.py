import argparse
import numpy as np
import datetime
import os
from dateutil import parser

from solarcalc.light_sim import LightSim, daterange
from solarcalc.spectopt import Light, Spectrum, BandCost, SimpleCost

def valid_date(s):
    try:
        if isinstance(s, datetime.datetime) or isinstance(s, datetime.date):
            return s
        return datetime.datetime.strptime(s, "%Y-%m-%d")
    except ValueError:
        msg = "Not a valid date in Y-m-d form: '{0}'.".format(s)
        raise argparse.ArgumentTypeError(msg)

def main():

    p = argparse.ArgumentParser()
    p.add_argument("-s", "--start", type=valid_date, default="2000-01-01",
                   help="Simulation start date")
    p.add_argument("-e", "--end", type=valid_date, default="2000-12-31",
                   help="Simulation end date")
    p.add_argument("-i", "--interval", type=int, default=5,
                   help="Simulation interval (minutes)")
    p.add_argument("-S", "--chamber-start", type=valid_date, default="2020-01-01",
                   help="Real date on which conditions will be run in a chamber")
    p.add_argument("-r", "--scale-factor", type=float, default=0.5,
                   help="Proportion of expected total irradiance to replicate in chambers. Controls for the fact that chambers can't get to full intensity. Set to approx max(chamber) / max(sunshine)")
    p.add_argument("-p", "--place", type=float, required=True, nargs=3,
                   help="Location of simulation as lat, long, elevation in decimal degrees")
    p.add_argument("-l", "--light", type=str, required=True,
                   help="Light definition file")
    p.add_argument("-o", "--output", type=argparse.FileType('w'), default='-',
                   help="Output file")
    args = p.parse_args()

    lat, lon, elev = args.place
    wavelengths= np.arange(350, 800, 5, dtype=float)
    sim = LightSim(args.start, args.end, lat, lon, elev, wavelengths=wavelengths)

    light = Light(args.light)

    print("datetime", "model_datetime", "temp", "humidity", "total_solar", *light.channels, sep="\t", file=args.output)
    date_offset = args.chamber_start - args.start
    for d in daterange(args.start, args.end, minutes=args.interval):
        doyf = d.timetuple().tm_yday + d.minute / 1440.0 + d.hour / 24.0
        theory = sim.combined_spline(doyf)
        temp, rh, totsrad, *spectra = theory
        spectra = Spectrum(wavelengths=wavelengths,
                           values=np.array(spectra) / 1000 * args.scale_factor)
        opt = light.optimise_settings(spectra)
        reald = d + date_offset
        f2 = lambda x: "{:0.2f}".format(np.round(x, 2))
        print(reald.isoformat(timespec='seconds'), d.isoformat(timespec='seconds'),
              f2(temp), f2(rh), f2(totsrad), *map(f2, opt), sep="\t", file=args.output)
    

if __name__ == "__main__":
    main()
