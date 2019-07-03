#!/bin/python3
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import datetime
from light_sim import LightSim
import os
from dateutil import parser

start = parser.parse(os.environ.get("START_DATE", "2012-01-01"))
end = parser.parse(os.environ.get("END_DATE", "2012-01-02"))



wavelengths= np.arange(int(os.environ.get("START_WAVELENGTH", 350)) , int(os.environ.get("END_WAVELENGTH", 750)), int(os.environ.get("WAVELENGTH_STEP", 10)))
# wavelengths=[420, 530, 660 ]

lat = os.environ.get("LATITUDE", -35.278454)
lon = os.environ.get("LONGITUDE", 149.121483)

if __name__ == "__main__":

    
    elevation = 577

    l = LightSim(start, end, lat, lon, elevation, wavelengths=wavelengths)

    l.write_file("output.csv", minutes=10)
