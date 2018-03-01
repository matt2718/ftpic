#!/usr/bin/python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.figure(1)

if len(sys.argv) < 3: exit()

modename = sys.argv[1]

for myfile in sys.argv[2:]:
	rawdata = pd.read_csv(myfile)

	time = rawdata['time'].values
	mode = rawdata[modename].values

	plt.plot(time, np.log10(mode), label=myfile)

plt.xlabel('Time')
plt.ylabel('log(E_' + modename + ')')

plt.legend()

plt.show()
