#!/usr/bin/python3

import sys
import pandas as pd
import matplotlib.pyplot as plt

if (len(sys.argv) < 2):
	rawdata = pd.read_csv(sys.stdin)
else:
	rawdata = pd.read_csv(sys.argv[1])

time = rawdata['time'].values
ke = rawdata['kinetic'].values
pe = rawdata['potential'].values
te = rawdata['total'].values
p = rawdata['momentum'].values

# energy plot
plt.figure(1)

plt.plot(time, ke, label='Kinetic')
plt.plot(time, pe, label='Potential')
plt.plot(time, te, label='Total')

plt.xlabel('Time')
plt.ylabel('Energy')
plt.axes().set_ylim(0, 1.25*max(te))
plt.legend()

# momentum plot
plt.figure(2)

plt.plot(time, p)

plt.xlabel('Time')
plt.ylabel('Momentum')
plt.axes().set_ylim(1.25*min(0,min(p)), 1.25*max(0,max(p)))

plt.show()
