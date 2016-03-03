## John D. Clark (written approx Sept. 2015)
## Complete-linkage hierarchical clustering script 

## The overall function of this script is to cluster a .pdb trajectory using 
## SciPy's complete-linkage hierarchical clustering algorithm over a range of 
## cutoff distances. It will produce assignments files 

## The purpose of this script is to compare RMSD as a distance metric to 
## occupancy fingerprints. For RMSD based clustering, you only need a .pdb file,
## an atom indices file (0 based), and parameters for clustering, as the 
## distance matrix is computed in this script. A precalculated flattened 
## distance matrix must be supplied (this can be done in the 
## occupancy_fingerprint-project repository).

## An alternate function of this script is to detect whether there are outlier
## data points that are clustered incorrectly. This is done by breaking the
## distance matrix by stride (e.g. for stride 2, 2 matrices will be produced:
## one will have indices 0,2,..,i+2,n; the other will be 1,3,...,i+2,n.
## By calculating the NRMV plots across cutoff distance, one can inter whether 
## there is a particular set of points that result in an outlier cluster.

def calcDistMatrix(coordsets):
  from pyRMSD.matrixHandler import MatrixHandler
  matrix = MatrixHandler().createMatrix(coordsets,'NOSUP_SERIAL_CALCULATOR')
  return matrix.get_data()

def clusterData(matrix, method='complete'):
  from scipy.cluster import hierarchy
  cluster = hierarchy.linkage(matrix,method=method)
  return cluster

def distRange(start,stop,step):
  values = []
  value = start
  while value <= stop:
    values.append(value)
    value = value + step
  return values

def createAssignments(dist_range,cluster_object,file_name,f_num=-1):
  from scipy.cluster import hierarchy
  for dist in dist_range:
    fcluster = hierarchy.fcluster(cluster_object,dist,criterion='distance')
    if f_num > -1:  
      fn = 'split_assignments/%s_assignments-%s-%s.dat'%(dist,file_name,f_num)
    else: fn = 'assignments_files/%s_assignments-%s.dat'%(dist,file_name)
    with open(fn,'w') as fn:
      for assignment in fcluster:
        fn.write((str(assignment))+'\n')
      fn.close()

def arraySplitter(array,num_arrays,original_size=(2880, 2880)):
  import numpy as np
  array = np.array(array)

  ## Converts condensed matrix to 2D triangular matrix of size (M x M)
  ## Matrix is upper triangular. Transpose and diagonal == 0.
  new_array = np.zeros(original_size)
  
  k = 0
  for i in range(len(new_array)):
    for j in range(i+1,len(new_array)):
      new_array[i,j] = array[k]
      k = k+1
  
  ## Splits upper triangular matrix into n submatrices:
  matrix_set = []
  for s in range(num_arrays):
    matrix = []
    for i in range(s,len(new_array),num_arrays):
      for j in range(i+num_arrays,len(new_array),num_arrays):
        matrix.append(new_array[i,j])
    matrix_set.append(matrix)

  return matrix_set

## Import packages and init arg parser.
import sys, os
import fileLoader
import argparse as args
import numpy as np
parser = args.ArgumentParser()
parser.add_argument('--pdb_file',default='../all_trajs-s1000.pdb')
parser.add_argument('--atom_indices',default=
'../../msmbuilder-hybrid/CA/AtomIndices-raw.dat \
../../msmbuilder-hybrid/within_6A_of_888/AtomIndices-raw.dat'
)
parser.add_argument('--fingerprint_matrix',default=
'../../pocket_fingerprint/fpMatrix_s1000.nc'
)
parser.add_argument('--split_array', default = 1,
help = 'An option to split the arrays to detect spurious cluster assignments.\
        If the value is greater than 1, the output will be a set of arrays\
        of NRMVs for each distance cutoff'
)
parser.add_argument('--distance_range', default = [1.0, 7.0, 0.2],
help = 'Distance range used for clustering using RMSD. Units are in Angstroms.\
        Format for input is *start_value end_value interval* '
)
args = parser.parse_args()


## Load .pdb file:
pdb_file = fileLoader.load_pdb(args.pdb_file)

## Loads atom indices to be used for alignment and RMSD calculation:
atom_indices = []
for file in args.atom_indices.split(' '):
  indices_name = file.split('/')[len(file.split('/'))-2]
  lines = fileLoader.load_dat(file)
  atom_indices.append((lines,indices_name))

## Distance range parameters for RMSD clustering:
start_dist, end_dist, interval = args.distance_range[0:]


## Loads precalculated fingerprint matrix.
fingerprints_matrix = fileLoader.load_nc(args.fingerprint_matrix)

from netCDF4 import Dataset
## user-defined distance cutoffs.
from prody import alignCoordsets
for atom_inds in atom_indices:
  selection = pdb_file[atom_inds[0]]
  align_coords = alignCoordsets(selection)
  aligned_coords = align_coords.getCoordsets()
  distance_matrix = calcDistMatrix(aligned_coords)
  
#  dmat_file = Dataset('%s-distanceMatrix.nc'%atom_inds[1],mode = 'w',format='NETCDF4')
#  dmat_dim = dmat_file.createDimension( 'distance_matrix')
#  dmat_var = dmat_file.createVariable('distance_matrix','f',('distance_matrix',))
#  dmat_var[:] = distance_matrix
#  continue
#  break
  if args.split_array > 1:
    if not os.path.exists('split_assignments'): os.mkdir('split_assignments')    
    split_matrices = arraySplitter(distance_matrix,args.split_array)
    for i in range(len(split_matrices)):
      cluster = clusterData(split_matrices[i])
      dist_range = distRange(start_dist,end_dist,interval)
      
      createAssignments(dist_range,cluster,atom_inds[1],f_num=i)
      print 'Finished creating assignments for submatrix %s for RMSD of %s' % (str(i),atom_inds[1])
  else:
    cluster = clusterData(distance_matrix)
    dist_range = distRange(start_dist,end_dist,interval)
     
    if not os.path.exists(os.path.join(os.getcwd(),'assignments_files')):
      os.mkdir(os.path.join(os.getcwd(),'assignments_files'))
  
    createAssignments(dist_range,cluster,atom_inds[1])


if os.path.exists(args.fingerprint_matrix):
  if args.split_array > 1:
    split_matrices = arraySplitter(fingerprints_matrix,args.split_array)
    for i in range(args.split_array):
      cluster = clusterData(np.array(split_matrices[i]))
      dist_range = distRange(0.001,0.45,0.01)
      createAssignments(dist_range,cluster,'fingerprints',f_num=i)
      print 'Finished creating assignments for submatrix %s of fingerprints'%str(i)
  else:
    cluster = clusterData(fingerprints_matrix)
    dist_range = distRange(0.001,0.45,0.01)
    createAssignments(dist_range,cluster,'fingerprints')

