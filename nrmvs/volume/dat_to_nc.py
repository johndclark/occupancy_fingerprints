import sys, os
from netCDF4 import Dataset

files = [fn for fn in os.listdir(os.getcwd()) if '.dat' in fn]

for f in files:
  data = []
  with open(f) as fn:
    for line in fn:
      x = int(line.split(' ')[0])
      y = float(line.split(' ')[1])

      data.append((x,y))

  new_filename = f.split('/')[len(f.split('/'))-1].strip('.dat')+'.nc'
  new_fn = Dataset(new_filename,'w',format='NETCDF4')
  new_dim = new_fn.createDimension('volume_data_x')
  new_dim = new_fn.createDimension('volume_data_y')
  new_var = new_fn.createVariable('volume_data','f8',('volume_data_x','volume_data_y'))
  new_var[:] = data


