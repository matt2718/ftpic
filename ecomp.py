#!/usr/bin/python3

import sys
import pandas as pd
import matplotlib.pyplot as plt

wp = -1 # plasma frequency specified on command line

i = 1
names = []
fnames = []
while i < len(sys.argv):
	arg = sys.argv[i]
	if arg == '-w':
		i += 1
		if i < len(sys.argv):
			wp = float(sys.argv[i])
	else:
		# arguments take the form name,file
		split = arg.find(',')
		if split > 0:
			name = arg[:split]
			fname = arg[split+1:]
			names.append(name)
			fnames.append(fname)
	i += 1

print(names)
print(fnames)
print(wp)
times = []
energies = []
for fname in fnames:
	rawdata = pd.read_csv(fname)
	times.append(rawdata['time'].values)
	energies.append(rawdata['total'].values)

# energy plot
plt.figure(1, figsize=(8,4))

norm = energies[0][0]
for i in range(len(energies)):
	# normalize
	energies[i] /= norm
	if wp > 0: times[i] *= wp

	plt.plot(times[i], energies[i], label=names[i])

if wp > 0:
	plt.xlabel(r'$\omega_p t$', fontsize=20)
else:
	plt.xlabel('Time', fontsize=14)

plt.ylabel('Total Energy', fontsize=14)
plt.axes().set_ylim(0.98, 1.1)
plt.legend()

plt.tight_layout()
plt.show()
