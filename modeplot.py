#!/usr/bin/python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# parse arguments
decay = None
wp = None
offset = 0

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

	elif arg == '-o':
		# decay line offset
		i += 1
		if i == len(sys.argv): quit(2)
		offset = float(sys.argv[i])

	elif arg == '-w':
		# plasma frequency
		i += 1
		if i == len(sys.argv): quit(2)
		wp = float(sys.argv[i])

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
	if wp != None: time = time * wp
	mode = rawdata['m' + modename].values

	tmpname = fname

	## hack for specific files
	#if 'asdff' in fname: tmpname = 'PIF'
	#if 'asdfo' in fname: tmpname = 'PIC'

	plt.plot(time, np.log10(mode), label = 'Mode ' + modename + ', ' + tmpname)
	#plt.plot(time, np.log10(mode), label = 'mode ' + modename)

if decay != None:
	y0 = np.log10(mode[0]) + offset
	if wp != None:
		plt.plot([0, 3 * wp], [y0, y0 - 3 / np.log(10) * decay])
	else:
		plt.plot([0, 3], [y0, y0 - 3/np.log(10) * decay])

if wp != None:
	plt.xlabel(r'$\omega_p t$', fontsize=20)
else:
	plt.xlabel('Time', fontsize=16)
	
plt.ylabel('log(E)', fontsize=16)

plt.legend(loc=4)

plt.show()
