## John D. Clark (written approx. Oct 2015) -- Minh Lab (Second Stage Screening
## group) -- Illinois Institute of Technology, Chicago, IL

## Expectation value calculation script:

## This script will calculate a particular type of expected value (see help
## menu, python calculateExpectationValue.py --help, for the different types
## that can be calculated), given an assignments file and parameter file of the
## same length whose indices reflect the original indices of the MD simulation
## frames of the .pdb file. The output is stored in a directory called
## 'expectation_values'. The purpose of this script is to calculate the expected
## value of a given parameter of biophysical interest by sampling from clusters
## of a cluster assignments file.

## The expected value of a particular parameter can be considered to be the long
## -run average of repetitions of an experiment. In our case, the experiment is
## to do runs of sampling from each cluster (one from each cluster, over
## multiple repititions), and measure how quickly (lowest # of samples) one can
## achieve the solution. We can assume based on the law of large numbers
## that the arithmetic mean of the values obtained from the experiment will
## almost surely converge upon the expected value.

## Now, given a sufficient number of samples, you should see on the plots of the
## expectation value. This is good, now the next step would be to measure the 
## variance as a function of the number of samples. From either plot, you can
## infer which distance metric you used in clustering gives you the best results
## If the former case is not true, you might want to turn up the number of
## samples, or you can see which distance metric converges most quickly by
## directly measuring the expectation value as a function of the number of
## samples.

## Note that there is a function to calculate the analytical solution (this is
## literally every nonredundant combination of possible selections from each
## cluster), and as you can imagine, it isn't going to work well for anything
## in the 1000+ frame regime because of the exploding number of permutations
## that tend to arise well before then. So, the exhaustiveSearch function is
## effectively broken. There is also a random cluster assignment generator that
## you can use to establish a baseline, if you feel so inclined.

### Returns an array of integer values whose indices correspond to MD frames.
def assignmentsFileReader(file):
  distance = file.split('_')[0]
  type = file.split('-')[1].strip('.dat')
  return distance, type

### Splits an array into a specified number of subarrays.
def arraySplitter(array,stride=3):
  split_arrays = []
  for i in range(stride):
    split_array = [array[j] for j in range(i,len(array),stride)]
    split_arrays.append(split_array)
  return split_arrays

### Calculates the analytical value of the NRMV.
def exhaustiveSearch(assignments,parameter_arr,stride=100):
  import itertools
  import numpy as np
  
  std_of_parameter_set = np.std(parameter_arr)
  cluster_range = list(set(assignments))
  cluster_values = []

  for cluster_index in cluster_range:
      
    values = [parameter_arr[frame] for frame in range(len(assignments)) if assignments[frame] == cluster_index]
    values = [values[i] for i in range(0,len(values),stride)]
    print values
    if len(values) > 0:
      cluster_values.append(values)
    else: pass
  combinations = list(itertools.product(*cluster_values))  
  sqr_avg_of_variances =np.sqrt(np.average([np.var(combo) for combo in combinations]))/float(std_of_parameter_set)

### Returns an expectation value and standard deviation for a number of
### types of estimators (weighted, exponential, maximum, minimum).
def getExpectation(assignments,param_arr,num_samples,calculate_exponential=False,calculate_maximum=False,calculate_minimum=False):
  cluster_samples = []

  num_samples = int(num_samples)
  cluster_range = list(set(assignments))
  for cluster_index in cluster_range:
    frames = [parameter_arr[frame] for frame in range(len(assignments)) if assignments[frame] == cluster_index and parameter_arr[frame] != 0]
    weight = len(frames)/float(len(assignments))

    if len(frames)>1:
      samples = []
      for sample in range(num_samples):
        
        new_sample = np.random.choice(frames)
        samples.append(new_sample)
      if calculate_exponential == True:
        samples = [np.exp(-sample) for sample in samples]
    elif len(frames) == 1:
      samples = param_arr[frames[0]]*np.ones((num_samples,))
      if calculate_exponential == True:
        samples = [np.exp(-sample) for sample in samples]
    ## If we are taking max or min, we don't want statistical weight.
    else: pass

    if calculate_maximum == True or calculate_minimum == True:
      cluster_samples.append(samples)
    else:
      cluster_samples.append([sample*weight for sample in samples])

  ## Takes the average and stdev of the maximum value from each column.
  if calculate_maximum == True:
    average = np.average(np.max(np.array(cluster_samples),axis=0))
    stdev = np.std(np.max(np.array(cluster_samples),axis=0))
  
  ## Takes the average and stdev of the minimum value from each column.
  elif calculate_minimum == True:
    average = np.average(np.min(np.array(cluster_samples),axis=0))
    stdev = np.std(np.min(np.array(cluster_samples),axis=0))

  ## Else, sums the columns.
  else:
    expectations = np.sum(np.array(cluster_samples),axis=0)
    average,stdev = np.average(expectations),np.std(expectations)

    ## Take log of both average and stdev
    if calculate_exponential == True:
      average,stdev = -np.log(average),-np.log(stdev)
  return average,stdev

### Generates an array of random values at a specified length and maximum value.
def randomCluster(length,noc):
  import numpy as np
  import random
  
  array = np.zeros(length)
  for i in range(len(array)):
    if noc == 0:
      pass
    else:
      array[i] = random.randint(0,noc)
  return array

### Import packages and init argument parser. ##################################
import sys, os
import numpy as np
import random
import fileLoader
import argparse as args

parser = args.ArgumentParser()
parser.add_argument('--assignments_directory',
                    default='assignments_files/')
parser.add_argument('--parameter_file',
                    default='parameter_files/volumes_s1000.nc')
parser.add_argument('--number_samples',
                    default=30000)
parser.add_argument('--exponential',
                    default=False)
parser.add_argument('--maximum',
                    default=True)
parser.add_argument('--minimum',
                    default=False)
parser.add_argument('--exhaustive_search',
                    default=False)
parser.add_argument('--random_cluster',
                    default=False)
parser.add_argument('--get_expectation',
                     default=True)
parser.add_argument('--split_assignments',
                     default=False)
args = parser.parse_args()

### Options: ###################################################################
num_samples = args.number_samples
exponential = args.exponential
maximum = args.maximum
minimum = args.minimum
exhaustive_search = args.exhaustive_search
random_cluster = args.random_cluster
get_expectation = args.get_expectation
split_assignments = args.split_assignments


### Gathers assignments files and loads parameter file. ########################
assignments_files = [file for file in os.listdir(args.assignments_directory)]
parameter_arr = fileLoader.load_nc(args.parameter_file)

if split_assignments == True: 
  num_splits = [assignments_files[i].split('-')[len(assignments_files[i].split('-'))-1] for i in range(len(assignments_files))]
  num_splits = len(set([int(num_splits[i].strip('.dat')) for i in range(len(num_splits))]))
  
  parameter_arrs = arraySplitter(parameter_arr,stride=num_splits)

### Calculates expectation value of cluster assignments. #######################
if get_expectation == True:
  data = []
  for assignments_file in assignments_files:
    print 'Processing '+assignments_file
    assignments = fileLoader.load_dat(os.path.join(args.assignments_directory,assignments_file))
    assignments = [int(assignment) for assignment in assignments]
    assignment_info = assignmentsFileReader(assignments_file)
    if split_assignments == True:
      submatrix_id = int(assignments_file.split('-')[len(assignments_file.split('-'))-1].strip('.dat'))
      
      expectation_value = getExpectation(
                                     assignments,
                                     parameter_arrs[submatrix_id],
                                     num_samples,
                                     calculate_exponential=exponential,
                                     calculate_maximum=maximum,
                                     calculate_minimum=minimum
                                     )
      data.append((assignment_info,expectation_value,max(assignments),submatrix_id))
    else:
      expectation_value = getExpectation(
                                     assignments,
                                     parameter_arr,
                                     num_samples,
                                     calculate_exponential=exponential,
                                     calculate_maximum=maximum,
                                     calculate_minimum=minimum
                                     )
      data.append((assignment_info,expectation_value,max(assignments)))

  for dataset in sorted(data):
    if not os.path.exists('expectation_values'):
      os.mkdir('expectation_values')
    fn = args.parameter_file.split('/')[len(args.parameter_file.split('/'))-1].strip('.nc') + '-' + dataset[0][1]+'.dat'
    noc = str(dataset[2])
    stdev = str(dataset[1][1])
    avg = str(dataset[1][0])
    distance = str(dataset[0][0])

### File writer for split assignment expectation values.
    if split_assignments == True:
      if not os.path.exists('split_expectations'):
        os.mkdir('split_expectations')
      submatrix_id = str(dataset[3])

      distline = 'Distance: %s\n'%distance
      submatline = 'Submatrix ID: %s\n'%submatrix_id
      nocline = 'NOC: %s\n'%noc
      evline = 'Expectation Value: %s\n'%avg
      stdline = 'Standard Deviation: %s\n'%stdev
      with open('split_expectations/%s'%fn,'a') as fn:
        fn.write(distline)
        fn.write(submatline)
        fn.write(nocline)
        fn.write(evline)
        fn.write(stdline)
        fn.write('\n')

### File writer for expectation values.
    else:
      nocline = 'NOC: %s\n'%noc
      evline = 'Expectation Value: %s\n'%avg
      stdline = 'Standard Deviation: %s\n'%stdev
      with open('expectation_values/%s'%fn,'a') as fn:
        fn.write(nocline)
        fn.write(evline)
        fn.write(stdline)
        fn.write('\n')

### Calculates analytical solution for a particular assignments file. ##########
if exhaustive_search == True:
  data = []
  for assignments_file in assignments_files:
    print 'Processing '+assignments_file
    assignments = fileLoader.load_dat(os.path.join(args.assignments_directory,assignments_file))
    assignments = [int(assignment) for assignment in assignments]
    assignment_info = assignmentsFileReader(assignments_file)
    exhaustive_solution = exhaustiveSearch(assignments,parameter_arr,stride=10)
    print 'Finished '+ assignments_file     
    data.append((max(assignments),exhaustive_solution))

  with open('analytical_solution/volume_solution.dat','a') as fn: 
    for dataset in sorted(data):
      noc = dataset[0]
      solution = dataset[1]
      line = str(noc)+' '+str(solution)+'\n'
      fn.write(line)

### Generation of random assignments and expectation. ##########################
if random_cluster == True:

  data = []
  for num_clusters in range(0,len(parameter_arr),100):
    data_set = []
    for sample in range(num_samples):
    
      random_arr = randomCluster(len(parameter_arr),num_clusters)
      expectation_value = getExpectation( 
                                     random_arr,
                                     parameter_arr,
                                     num_samples,
                                     calculate_exponential=exponential,
                                     calculate_maximum=maximum,
                                     calculate_minimum=minimum
                                     )

      data_set.append(expectation_value)
    
    average = np.average(data_set)
    stdev = np.std(data_set)
    data.append((num_clusters,average,stdev))
  with open('expectation_values/random-volumes.dat','a') as fn:
  
    for data in sorted(data):
      noc = data[0]
      avg = data[1]
      stdev = data[2]

      line = str(noc) + ' ' + str(avg) + ' ' + str(stdev) + '\n' 
      fn.write(line)

