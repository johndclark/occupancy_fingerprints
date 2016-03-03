## John D. Clark (Written approximately Sept. 2015) -- Minh Lab (Second Stage
## Screening group) -- Illinois Institute of Technology, Chicago, IL

## NRMV Calculation script:

## The overall function of this script is to calculate the NRMV of a particular
## distance cutoff given a parameter file and an assignments file. Both files 
## must contain a list of values whose indices reflect that of the original .pdb
## file. The output of this script can be plotted using the provided scripts.

## The NRMV is calculated simply by taking the mean of the variance within each
## cluster. Here's a quick example: if an assignments file with 3 frames whose
## values in the parameter file are 1.0, 2.0, and 3.0 (respectively), and the
## assignments file looks like 1, 1, and 2, then the first cluster is [1.0,2.0]
## and the second cluster is just [3.0]. Since there is only one value in the
## second cluster, it equals 0. The first cluster's variance is 0.25. The mean
## of the two values is 0.125. Now take the root, equals ~0.354. Finally,
## we normalize this value against the standard deviation of the whole dataset
## (since we took the root of the mean variance). The standard deviation is
## ~0.816, so the final answer is ~0.433.

def calcNRMV(parameter_array,assignments_array):
  if len(parameter_array) == len(assignments_array):
    import numpy as np
    
    parameter_values = [parameter_array[i] for i in range(len(parameter_array)) if parameter_array[i] != 0]
    total_stdev = np.std(parameter_values)
    
    cluster_variances = []
    
    if max(assignments_array) == 1:
      cluster_variances = np.var([parameter_array[frame] for frame in range(len(assignments_array)) if parameter_array[frame] != 0])
      
    else:
      for cluster_index in range(max(assignments_array)+1):
        parameter_bin = []
        for frame in range(len(assignments_array)):
          if assignments_array[frame] == cluster_index:
            if parameter_array[frame] != 0:
              parameter_bin.append(parameter_array[frame])
        if len(parameter_bin) > 1:
          cluster_variance = np.var(parameter_bin)
          cluster_variances.append(cluster_variance)
        elif len(parameter_bin) == 1: 
          cluster_variance = 0
          cluster_variances.append(cluster_variance)
        else:
          pass
#          cluster_variance = 0
#          cluster_variances.append(cluster_variance)
    average = np.mean(cluster_variances)
    rmv = np.sqrt(average)
    nrmv = rmv/float(total_stdev)
    return nrmv
  else: 
    print 'Assignments array is not the same length as parameter array'

## Import packages and init arg parser.
import sys, os 
import argparse as args
import fileLoader

parser = args.ArgumentParser()
parser.add_argument('--parameter_fileDir',
       help='Path to directory containing parameter files',
       default='../../pocket_volume/volumes_s1000.nc'
       )

parser.add_argument('--assignments_fileDir',
       help='Path to directory containing assignments files',
       default='split_assignments/'
       )
args = parser.parse_args()

## Checks for parameter files in a directory. Parameter file should be a .nc
## file, which contains a list of values whose indices reflect the index of
## the original .pdb file you clustered. If you are missing a value in the 
## parameter file, set every value that you are missing to 0.
if os.path.isdir(args.parameter_fileDir):
  parameter_files = [os.path.join(args.parameter_fileDir,file) for file in os.listdir(args.parameter_fileDir) if os.path.isfile(os.path.join(args.parameter_fileDir,file))]
else: parameter_files = [args.parameter_fileDir]

assignment_files = [os.path.join(args.assignments_fileDir,file) for file in os.listdir(args.assignments_fileDir)]

## Calculate NRMV and save files.
from netCDF4 import Dataset

for parameter_file in parameter_files:
  
  parameter_array = fileLoader.load_nc(parameter_file)
  for assignment_file in assignment_files:
    assignment_array = fileLoader.load_dat(assignment_file)
    noc = max(assignment_array)
    nrmv = calcNRMV(parameter_array,assignment_array)
    if not os.path.exists('nrmvs'): os.mkdir('nrmvs')
    nrmv_fn = 'nrmvs/nrmv-%s-%s'%(
parameter_file.split('/')[len(parameter_file.split('/'))-1].strip('.nc'),
assignment_file.split('-')[len(assignment_file.split('-'))-1]
)
    with open(nrmv_fn,'a') as fn:
      nrmv_line = str(noc) + ' ' + str(nrmv) + '\n'
      fn.write(nrmv_line)
      fn.close()
    
