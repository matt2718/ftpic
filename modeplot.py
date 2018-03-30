#!/usr/bin/python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 3: exit()

modename = sys.argv[1]

plt.figure(num=modename)

for myfile in sys.argv[2:]:
	rawdata = pd.read_csv(myfile)

	time = rawdata['time'].values
	mode = rawdata[modename].values

	plt.plot(time, np.log10(mode), label=myfile)

plt.xlabel('Time')
plt.ylabel('log(E_' + modename + ')')

plt.legend(loc=4)

plt.show()
