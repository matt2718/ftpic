#!/usr/bin/python3

import sys
import pandas as pd
import matplotlib.pyplot as plt

wp = -1 # plasma frequency specified on command line

i = 1
fname = ''
while i < len(sys.argv):
	if sys.argv[i] == '-w':
		i += 1
		if i < len(sys.argv):
			wp = float(sys.argv[i])
	else:
		fname = argv[i]
	i += 1

if fname:
	rawdata = pd.read_csv(fname)
else:
	rawdata = pd.read_csv(sys.stdin)

time = rawdata['time'].values
ke = rawdata['kinetic'].values
pe = rawdata['potential'].values
te = rawdata['total'].values
p = rawdata['momentum'].values

# energy plot
plt.figure(1)

if wp > 0: time = time * wp

plt.plot(time, ke, label='Kinetic')
plt.plot(time, pe, label='Potential')
plt.plot(time, te, label='Total')

if wp > 0:
	plt.xlabel(r'$\omega_p t$', fontsize=20)
else:
	plt.xlabel('Time', fontsize=14)

plt.ylabel('Energy', fontsize=14)
plt.axes().set_ylim(0, 1.3*te[0])
plt.legend()

# momentum plot
plt.figure(2)

plt.plot(time, p)

if wp > 0:
	plt.xlabel(r'$\omega_p t$', fontsize=20)
else:
	plt.xlabel('Time', fontsize=14)

plt.ylabel('Momentum', fontsize=14)
plt.axes().set_ylim(1.25*min(0,min(p)), 1.25*max(0,max(p)))

plt.show()
