import sys, os
import argparse as args
import numpy as np

parser = args.ArgumentParser()
parser.add_argument('--file_dir',default='within_6A_of_888/')
args = parser.parse_args()

files = [os.path.join(args.file_dir,file) for file in os.listdir(args.file_dir)]
cluster_range = sorted([int(line.split(' ')[0]) for line in open(files[0])])
# Removes redundant cluster indices.
cluster_range = [cluster_range[i] for i in range(len(cluster_range)) if cluster_range[i] != cluster_range[i-1]]

nrmvs = []
for cluster_number in cluster_range:
  cluster_nrmvs = []
  for file in files:

    with open(file) as fn:
      data = [(int(line.split(' ')[0]),float(line.split(' ')[1].strip('\n'))) for line in fn]
      for d in data:
        if d[0] == cluster_number:
          cluster_nrmvs.append(d[1])
  average = np.average(cluster_nrmvs)
  stdev = np.std(cluster_nrmvs)
  nrmvs.append((cluster_number,average,stdev))

from netCDF4 import Dataset
new_fn = args.file_dir.strip('/')+'-averageScores.nc'
new_fn = Dataset(new_fn,'w',format='NETCDF4')
new_dim = new_fn.createDimension('data_set')
new_dim = new_fn.createDimension('data')
new_var = new_fn.createVariable('scores','f8',('data_set','data'))
new_var[:] = nrmvs
