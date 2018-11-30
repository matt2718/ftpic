#!/usr/bin/python3

import sys
import pandas as pd
import numpy as np

rawdata = pd.read_csv(sys.argv[1])

time = rawdata['time'].values
te = rawdata['total'].values

te = te / te[0]

print(max(np.abs(np.diff(te))))
