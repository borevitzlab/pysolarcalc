#!/bin/python3
from light_sim import LightSim
import datetime


if __name__ == "__main__":

    start = datetime.datetime.strptime("2015-02-02", "%Y-%m-%d")
    end = start + datetime.timedelta(weeks=4)
    latitude = -35.278454
    longitude = 149.121483
    elevation = 577
    l = LightSim(start, end, latitude, longitude, elevation)
    l.write_file("output.csv", output_start=datetime.datetime.strptime("2019-01-01", "%Y-%m-%d") ,minutes=10)

