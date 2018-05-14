#!/usr/bin/python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# parse arguments
decay = None
files = {}
plots = []

i = 1
while i < len(sys.argv):
	arg = sys.argv[i]
	if arg == '-d':
		# theoretical decay rate
		i += 1
		if i == len(sys.argv): quit(2)
		decay = float(sys.argv[i])
	else:
		# arguments take the form mode,file
		split = arg.find(',')
		if split > 0:
			mode = arg[:split]
			fname = arg[split+1:]
			plots.append((mode,fname))
			files[fname] = None
	i += 1

# read files
for fname in files:
	files[fname] = pd.read_csv(fname)

plt.figure()

for modename,fname in plots:
	rawdata = files[fname]

	time = rawdata['time'].values
	mode = rawdata['m' + modename].values

	plt.plot(time, np.log10(mode), label = modename + ',' + fname)

if decay != None:
	y0 = np.log10(mode[0]) - 0.2
	plt.plot([0, 5], [y0, y0 - 5/np.log(10) * decay])

plt.xlabel('Time')
plt.ylabel('log(E_' + modename + ')')

plt.legend(loc=4)

plt.show()
