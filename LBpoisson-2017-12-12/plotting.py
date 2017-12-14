#!/usr/bin/python3
# Python script to plot results

import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("rmsError.dat")

sizes = []
errors = []

for i in data:
    sizes.append(i[0])
    errors.append(i[1])

plt.loglog(sizes, errors)
#plt.plot(sizes, errors)
plt.show()