import sys, os

def load_gz(file):
  from gzip import open
  gz = open(file)
  return gz

def load_nc(file):
  from netCDF4 import Dataset
  nc = Dataset(file)
  key = nc.variables.keys()[0]
  nc = nc.variables[key][:]
  return nc

def load_h5(file):
  from msmbuilder import io
  h5 = io.loadh(file)
  return h5 

def load_pdb(file):
  from prody import parsePDB
  pdb = parsePDB(file)
  return pdb

## also loads other text files.
def load_dat(file):
  try:
    lines = [int(line.strip('\n')) for line in open(file,'r')]
    return lines
  except:
    lines = [int(line) for line in open(file,'r')]
    return lines
