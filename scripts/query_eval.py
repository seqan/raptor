#!/usr/bin/env python

import pandas as pd
import sys

fin = sys.argv[1]
initial = int(sys.argv[2])
increment = int(sys.argv[3])

batch = int(fin.split("/")[-1][:-5])
last_ref = initial - 1 + (increment * batch)

data = pd.read_csv(fin, sep="\t", header=None, usecols=[0,1], names=["read", "ref"], converters={0: lambda v: int(v[4:][:-5])}).fillna(-1)
data.ref = data.ref.replace(to_replace=['bin_', '.fa'], value='', regex=True).astype(int)

data["TP"] = (data.ref <= last_ref) & (data.ref == data.read)
data["FP"] = (data.ref <= last_ref) & (data.ref != data.read) & (data.ref >= 0)
data["TN"] = (data.read > last_ref) & (data.ref == -1)
data["FN"] = (data.read <= last_ref) & (data.ref == -1)

data[["TP", "FP", "TN", "FN"]].sum().to_csv(sys.stdout, sep="\t", header=None)