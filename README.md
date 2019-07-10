# pysolarcalc

A python port of SolarCalc by Kurt Spokas

Follows formulae largely from [Spokas and Forcella (2006)](https://pubag.nal.usda.gov/catalog/1910) and [this NREL tech report](https://rredc.nrel.gov/solar/pubs/spectral/model/section3.html)

### Usage:


```bash
# pip install -r requirements.txt && ./setup.py install OR
# conda env create -f dev-environment.yml && ./setup.py install


pysolarcalc --place $lat $lon $elev \
	--light data/spectra/growtainer-wm2.csv  \
	--start 2019-01-01 \
	--end 2019-12-31 \
	--chamber-start 2020-01-01 \
	--scale-factor 0.25  # see below
```


The main params are obvious (lat/long/elevation/start & end dates).
`--scale-factor` is a multiplier that accounts for the fact that a chamber
cannot produce the quantity of light the sun does, and allows one to optimise
spectra assuming that the sun was `scale_factor` times as bright. Defaults to
0.5, which is about right for our high light chambers. It should be set to
approximately `max(chamber par) / max(sun)`.
