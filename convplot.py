#!/usr/bin/python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

dirname = 'dtscale/'
dt_ar = []
error_o = []
error_f = []

for i in range(10):
	dt = 0.0001 * (1 << i)
	dt_ar.append(dt)
	
	rawdata = pd.read_csv(dirname + 'o-' + str(dt) + '.csv')
	etotal = rawdata['total'].values
	error_o.append(max(abs(etotal - etotal[0])) / etotal[0])

	rawdata = pd.read_csv(dirname + 'f-' + str(dt) + '.csv')
	etotal = rawdata['total'].values
	error_f.append(max(abs(etotal - etotal[0])) / etotal[0])

plt.figure(1)
plt.loglog(dt_ar, error_o, label='o')
plt.loglog(dt_ar, error_f, label='f')

plt.xlabel('DT')
plt.ylabel('max error')

plt.legend()
plt.grid()

plt.show()
