import scipy as sp
import numpy as np
from scipy import optimize
from scipy.spatial import distance

def spectfun(weights, desired, machinespec):
    return distance.cityblock(desired, np.dot(weights.T, machinespec))

def find_weights(desired, machinespec):
    opt = optimize.minimize(spectfun, np.zeros(4), (want, testmachine),
                            options={"maxiter": 10000}, bounds=[(0.,1.)]*4)
    return opt.x

testmachine = np.eye(4)
testmachine[0,1] = 0.2
want = np.array([0.3, 0.2, 0.2, 0.4])
find_weights(want, testmachine)


