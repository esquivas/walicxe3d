#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
from constants import RJUP, AMH, KB

# Estas 3 lineas permiten editar el contenido en walicxe3d_utils y que los
# cambios surtan efecto sin tener que reiniciar la sesion de python
import importlib
import walicxe3d_utils
importlib.reload(walicxe3d_utils)

from walicxe3d_utils import *

try:
  nout = int( sys.argv[1] )
except:
  # default output # if not given as argument
  nout = 0
try:
  model_name = str( sys.argv[2] )
except:
  # default model name if not given as argument
  model_name = 'M1'
try:
  fnum = int( sys.argv[3] )
except:
  # default figure number if not given as argument
  fnum = 0

################################################################################
#  Simulation parameters
#path = '/Users/esquivel/Desktop/diable_storage/storage5/esquivel/wexo/' \
#     + model_name+'/output/'
path = '/storage5/esquivel/wexo/'+ model_name+'/output/'
neq     =  0
neqtot  =  7
maxlev  =  5
ncellsB = 16
nbrootx =  4
nbrooty =  1
nbrootz =  4
nprocs  = 40
CV      = 1.5   # 1/(gamma-1)
verbose = True  # optional


dens = read_3d_eq(path, nout, neq, neqtot, maxlev, ncellsB, \
                  nbrootx, nbrooty, nbrootz, nprocs, CV, verbose = verbose )

plt.ion()

plt.figure(2) ; plt.clf()
plt.imshow(dens[:,127 ,:].T, origin='lower', norm=LogNorm(vmin=100, vmax=1e5), \
           cmap = 'cividis' )

plt.colorbar()
