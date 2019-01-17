import argparse as ap

from solarcalc.spectopt import Light, Spectrum, BandCost, SimpleCost

def parse_theory(tfile):
    with open(tfile) as fh:
        header = list(map(lambda x: x.strip(), fh.readline().split(",")))
        domain = [int(x.rstrip("nm")) for x in header[5:]]
        for line in fh:
            cells = list(map(lambda x: x.strip(), line.split(",")))
            date, mdate, temp, hum, tot, *wls = cells
            wls = list(map(float, wls))
            yield (date, mdate, temp, hum, tot, Spectrum(domain, wls))


def optimise_for_light():
    p = ap.ArgumentParser()
    p.add_argument("-l", "--light", type=str, required=True,
                   help="Light definition file")
    p.add_argument("theoretical", type=str,
                   help="Theoretical model output as a CSV")
    args = p.parse_args()

    light = Light(args.light, interpolation='linear')

    cf = BandCost(bands=[(400, 500), (500, 600), (600, 700)], weights=[1, 0.1, 1])
    light = Light("./data/spectra/growtainer-actual.csv", interpolation='linear')
    print("datetime", "model_datetime", "temp", "humidity", "total_solar", *light.channels, sep="\t")
    for line in parse_theory("./data/output.csv"):
        want = line[5][min(light.wavelengths):max(light.wavelengths)]
        opt = light.optimise_settings(want, cost_function=cf)
        print(*line[0:5], *opt, sep="\t")


if __name__ == "__main__":
    optimise_for_light()
