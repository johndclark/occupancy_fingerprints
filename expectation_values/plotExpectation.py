## John D. Clark (written approx. Nov 2015) -- Minh Lab (Second Stage Screening
## group) -- Illinois Institute of Technology, Chicago, IL

## Expectation value plotting script:

## Just a basic plotting script for data obtained from calculateExpectationValue
## .py. Make sure the data is organized into this directory format: 
## '/*X*_samples/parameter_type/' where *X* is the number of samples.
## The plot_type argument should reflect the type of expectation value you've
## calculated. The plot_var argument will plot the variance instead of the
## expected value. The truncate_clusters argument allows you to focus in on a 
## range of number of clusters from 1 cluster to X clusters.


import os, sys
import argparse as args

parser = args.ArgumentParser()

parser.add_argument('--file_dir',default='30k_samples/volume/')
parser.add_argument('--plot_type',default='maximum')
parser.add_argument('--plot_var',default='False')
parser.add_argument('--truncate_clusters',default='None')
args = parser.parse_args()

files = [os.path.join(args.file_dir,f) for f in os.listdir(args.file_dir) if args.plot_type in f]

####
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

colors = ['red','green','blue']
color_index = 0
patches = []
for f in files:

  file_info = f.split('/')[len(f.split('/'))-1].split('-')
  ligand = file_info[0]
  metric = file_info[1]
  expectation_type = file_info[len(file_info)-1].strip('.dat')

  with open(f,'r') as fn:
    x_vals = [int(line.split(' ')[len(line.split(' '))-1]) for line in fn if 'NOC' in line]
  with open(f,'r') as fn:
    y_vals = [float(line.split(' ')[len(line.split(' '))-1]) for line in fn if 'Expectation' in line]
  with open(f,'r') as fn:
    stdevs = [float(line.split(' ')[len(line.split(' '))-1]) for line in fn if 'Standard' in line]

  if args.truncate_clusters != 'None':
    x_vals = [x[0] for x in x_vals if x[0] < int(args.truncate_clusters)]
    y_vals = [y[1] for y in y_vals if y[0] < int(args.truncate_clusters)]
    stdevs = [z[2] for z in stdevs if z[0] < int(args.truncate_clusters)]
  if args.plot_var == 'True':
    line = plt.plot(x_vals,[np.log(stdev*stdev) for stdev in stdevs],marker='o',color=colors[color_index])

  else:
    line = plt.plot(x_vals,y_vals,marker='o',color=colors[color_index])
    err = plt.errorbar(x_vals,y_vals,stdevs,linestyle="None",marker='o',color=colors[color_index])
  patches.append(mpatches.Patch(color=colors[color_index],label=metric))
  color_index = color_index+1

print len(patches)
plt.legend(handles=[patches[0],patches[1],patches[2]],loc=0)
if args.plot_var == 'True':
  plt.title('Log variance of %s expectation for %s'%(args.plot_type,str(args.file_dir).split('/')[len(args.file_dir.split('/'))-2]),size=12)
  plt.xlabel('Number of clusters (NOC)')
  plt.ylabel('log($\sigma^2$)')
else:
  plt.title('Expectation of %s dock score for %s'%(args.plot_type,str(args.file_dir).split('/')[len(args.file_dir.split('/'))-2]),size=12)
  plt.xlabel('Number of clusters (NOC)')
  plt.ylabel('<Dock Score>')
plt.show()
