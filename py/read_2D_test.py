#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
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
path = '/Users/esquivel/Desktop/diable_storage/storage5/esquivel/wexo/' \
     + model_name+'/output/'
#path = '/storage5/esquivel/wexo/'+ model_name+'/output/'
neq     =  6
neqtot  =  7
#cut     = 128
maxlev  =  5
ncellsB = 16
nbrootx =  4
nbrooty =  1
nbrootz =  4
nprocs  = 32
#axis    = 'Y'
CV      = 1.5   # 1/(gamma-1)
verbose = False  # optional

plt.ion()

#-------------------------------------------------------------------------------
#  YZ plane
cut  = 512
axis = 'X'
dens_cutX = read_2d_cut(path, nout, neq, neqtot, cut, maxlev, ncellsB, \
                       nbrootx, nbrooty, nbrootz, nprocs, CV, axis=axis, \
                       verbose = verbose )

mesh_cutX = read_mesh_2D_cut(path, nout, neqtot, cut, maxlev, ncellsB, nbrootx,\
                             nbrooty, nbrootz, nprocs, axis, verbose=verbose)

fig = plt.figure(num = fnum)  ; plt.clf()
ax  = fig.add_subplot(111)

im = ax.imshow(dens_cutX, origin='lower', norm=LogNorm(vmin=100, vmax=1e5), \
           cmap = 'cividis' )
cb = plt.colorbar(im)

for nb in range( mesh_cutX.shape[0] ):
  rect = patches.Rectangle((mesh_cutX[nb,0], mesh_cutX[nb,1]), mesh_cutX[nb,2],
                            mesh_cutX[nb,3], linewidth=0.5,
                            edgecolor='lightgray', facecolor='none', alpha=0.25 )
  ax.add_patch(rect)

  ax.set_xlabel(f'Y [$R_{{p}}$]')````
ax.set_ylabel(f'Z [$R_{{p}}$]')
#-------------------------------------------------------------------------------
#  XZ plane
cut  = 128
axis = 'Y'
dens_cutY = read_2d_cut(path, nout, neq, neqtot, cut, maxlev, ncellsB, \
                       nbrootx, nbrooty, nbrootz, nprocs, CV, axis=axis, \
                       verbose = verbose )

mesh_cutY = read_mesh_2D_cut(path, nout, neqtot, cut, maxlev, ncellsB, nbrootx,\
                             nbrooty, nbrootz, nprocs, axis, verbose=verbose)

fig = plt.figure(num = fnum+1) ; plt.clf()
ax  = fig.add_subplot(111)

im = ax.imshow(dens_cutY, origin='lower', norm=LogNorm(vmin=100, vmax=1e5), \
           cmap = 'cividis' )
cb = plt.colorbar(im)

for nb in range( mesh_cutY.shape[0] ):
  rect = patches.Rectangle((mesh_cutY[nb,0], mesh_cutY[nb,1]), mesh_cutY[nb,2],
                            mesh_cutY[nb,3], linewidth=0.5,
                            edgecolor='lightgray', facecolor='none', alpha=0.25 )
  ax.add_patch(rect)

ax.set_xlabel(f'X [$R_{{p}}$]')
ax.set_ylabel(f'Z [$R_{{p}}$]')

#-------------------------------------------------------------------------------
#  XY plane
cut  = 512
axis = 'Z'
dens_cutZ = read_2d_cut(path, nout, neq, neqtot, cut, maxlev, ncellsB, \
                       nbrootx, nbrooty, nbrootz, nprocs, CV, axis=axis, \
                       verbose = verbose )

mesh_cutZ = read_mesh_2D_cut(path, nout, neqtot, cut, maxlev, ncellsB, nbrootx,\
                             nbrooty, nbrootz, nprocs, axis, verbose=verbose)

fig = plt.figure(num = fnum+2) ; plt.clf()
ax  = fig.add_subplot(111)

im = ax.imshow(dens_cutZ.T, origin='lower', norm=LogNorm(vmin=100, vmax=1e5), \
           cmap = 'cividis' )
cb = plt.colorbar(im)

for nb in range( mesh_cutZ.shape[0] ):
  rect = patches.Rectangle((mesh_cutZ[nb,0], mesh_cutZ[nb,1]), mesh_cutZ[nb,2],
                            mesh_cutZ[nb,3], linewidth=0.5,
                            edgecolor='lightgray', facecolor='none', alpha=0.25 )
  ax.add_patch(rect)


ax.set_xlabel(f'X [$R_{{p}}$]')
ax.set_ylabel(f'Y [$R_{{p}}$]')
