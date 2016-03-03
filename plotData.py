

import sys, os
import argparse as args
import fileLoader

parser = args.ArgumentParser()
parser.add_argument('--data_dir',
                   help='Path to a directory containing data files in the form:\
                    x y.',
                   default='nrmvs/dock_score'
                   )
args = parser.parse_args()


from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

plot_type = args.data_dir.split('/')[len(args.data_dir.split('/'))-1]

files = [os.path.join(args.data_dir,file) for file in os.listdir(args.data_dir) if '.nc' in file]
markers = ['D','o','^']
colors = ['red','green','blue']

for i in range(len(files)):
  fn = fileLoader.load_nc(files[i])
  
  if 'CA' in files[i]:
    label = 'CA'
  if 'within' in files[i]:
    label = 'within 6A of ligand CoM'
  if 'fingerprint' in files[i]:
    label = 'occupancy fingerprint'

  if len(fn[0]) == 3:
    exec(
'line%s = plt.errorbar([x[0] for x in fn],[y[1] for y in fn],[z[2] for z in fn],linestyle="None",color=colors[%i],marker=markers[%i])'%(i,int(i),int(i)))
  else:
    exec(
'line%s = plt.plot([x[0] for x in fn],[y[1] for y in fn],linestyle="None",marker=markers[%i])'%(i,int(i)))
  exec('patch%s = mpatches.Patch(color=colors[%i],label="%s")'%(i,int(i),label))

legend = plt.legend(handles=[patch0,patch1,patch2])
plt.title('Comparison of selected distance metrics by NRMV of %s'%plot_type,size=13)
#plt.text(1700,0.7,'$\sigma$ = 1538.9667 $\AA^3$')
plt.xlabel('Number of Clusters (NOC)')
plt.ylabel('NRMV of %s'%plot_type)
plt.show()
