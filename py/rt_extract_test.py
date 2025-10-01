#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
import sys
from constants import RJUP, AMH, KB

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
  fnum = 10
################################################################################
neq = 7
nx  = 1024
ny  =  256
nz  = 1024

#path = '/Users/esquivel/Desktop/diable_storage/storage5/esquivel/wexo/' \
#     + model_name+'/output/'

path = '/storage5/esquivel/wexo/'+ model_name+'/output/'

plt.ion()
################################################################################
#  scalings
l_sc = 1.98*RJUP   # length scale (cm)
d_sc = 1.0*AMH     # density scale (g cm^-3)
v_sc = 1.0e5       # velocity scale (cm s^-1)
pas_sc = d_sc      # passive scalars scale
#  Derived scalings
p_sc = d_sc*v_sc**2
e_sc = p_sc
t_sc = l_sc/v_sc
B_sc = np.sqrt(4.0*np.pi*p_sc)
################################################################################
file_name = path + f'CutY.{nout:04d}.bin'
data = read_cut( file_name, neq, nx, ny, nz, axis='Y')

# Unpack data
# Ojo el extract escribe las primitivas
rho = data[:,:,0]
vx  = data[:,:,1]
vy  = data[:,:,2]
vz  = data[:,:,3]
Pg  = data[:,:,4]
nh0 = data[:,:,5]
Tau = data[:,:,6]

#  Ahora escalamos (a unidades de codigo)
rho = data[:,:,0] / d_sc
vx  = data[:,:,1] / v_sc
vy  = data[:,:,2] / v_sc
vz  = data[:,:,3] / v_sc
Pg  = data[:,:,4] / p_sc
nh0 = data[:,:,5]
#  set floor to y0
nh0[ nh0 <= 0 ] = 1e-3
y0  = nh0/rho

#Temp = Pg/rho*(1.3*AMH*p_sc/d_sc/KB)
Temp = Pg/(2.0*rho-nh0)*(1.3*AMH*p_sc/d_sc/KB)

vmag = np.sqrt( vx**2 + vy**2 + vz**2 )

print('Pressure ', end=': ')
minmax(Pg)
print('Temp     ', end=': ')
minmax(Temp)

#  Read, scale and center mesh
filemesh = path + f'MeshCutY.{nout:04d}.bin'
dx = 50.0/nx
mesh = read_mesh(filemesh, dx, center = True, xc=25.0, yc=25.0)

################################################################################
#  simulate rad transfer here
a0   = 6.3e-18     #  Fotoionization cross section
#S0st = 2.e33      #  Photon rate [ s^-1 ]
F0   = 2.e13       #  Ioniizng flux [ cm^-2 s^-1]
E0   = 3.8e-12     #  energy gain per ionization [ erg ]
dz   = 50*RJUP*1.98/float(nz)
dV   = dz**3
Tau2 = np.zeros((nz,nx))
phi  = np.zeros((nz,nx))
phi2 = np.zeros((nz,nx))
psi  = np.zeros((nz,nx))
psi2 = np.zeros((nz,nx))
flux = np.zeros((nz,nx))

#  emulate Rad Transfer post-processing
for i in range(nx):
  for k in range(nz):
    if (k==0):
      dtau = a0 * dz * nh0[k,i]
      flux[k,i] = F0
      phi [k,i] =     F0*(1-np.exp(-dtau)) / dz / nh0[k,i]
      psi [k,i] = E0* F0*(1-np.exp(-dtau)) / dz
      Tau2[k,i] = dtau
    else:
      dtau = a0 * dz * nh0[k,i]
      phi [k,i] =      flux[k-1,i]*(1.0-np.exp(-dtau)) / dz / nh0[k,i]
      psi [k,i] = E0 * flux[k-1,i]*(1.0-np.exp(-dtau)) / dz
      Tau2[k,i] = Tau2[k-1,i] + dtau
      flux[k,i] = flux[k-1,i]*np.exp(-dtau)
################################################################################
for i in range(nx):
  for k in range(nz):
    if (k==0):
      dtau = a0 * dz * nh0[k,i]
      phi2[k,i] =      F0 * ( 1-np.exp(-dtau) ) /dz /nh0[k,i]
      psi2[k,i] = E0 * F0 * ( 1-np.exp(-dtau) ) /dz
      #Tau[k,i]  = dtau
    else :
      dtau = a0 * dz * nh0[k,i]
      phi2[k,i] =    F0*np.exp(-Tau[k,i]) * (1.0-np.exp(-dtau)) /dz /nh0[k,i]
      psi2[k,i] = E0*F0*np.exp(-Tau[k,i]) * (1.0-np.exp(-dtau)) /dz

print('phi  ', end=': ')
minmax(phi)
print('phi2 ', end=': ')
minmax(phi2)

print('psi  ', end=': ')
minmax(psi)
print('psi2 ', end=': ')
minmax(psi2)

################################################################################
#  Plots
#plt.figure(num=fnum, figsize=(16,11) ) ; plt.clf()
fig, axs = plt.subplots(3, 4, layout='tight', num=fnum, \
                        sharex=True, sharey=True, figsize=(16,12) )

left   = 0.045    # the left side of the subplots of the figure
right  = 0.99    # the right side of the subplots of the figure
bottom = 0.02    # the bottom of the subplots of the figure
top    = 0.95    # the top of the subplots of the figure
wspace = 0.05   # the amount of width reserved for blank space between subplots
hspace = 0.0   # the amount of height reserved for white space between subplots

plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, \
                    wspace=wspace, hspace=hspace)


fig.suptitle(f'Model '+model_name+f' @ time: {nout*0.1:.2f} days')
#-------------------------------------------------------------------------------
ax = axs[0,0]
vmin = 1e2
vmax = 1e5
im = ax.imshow(rho, cmap='cividis', extent = [-25,25,-25,25], \
           norm=LogNorm(vmin=vmin,vmax=vmax), origin='lower' )
cb = plt.colorbar(im, ax=ax, extend=colorbar_ext(rho, vmin, vmax), shrink=0.7)
ax.set_title(r'Total Density [ cm$^{-3}$]')
ax.set_ylabel(f'Z [$R_{{p}}$]')

for nb in range( mesh.shape[0] ):
  rect = patches.Rectangle((mesh[nb,0], mesh[nb,1]), mesh[nb,2], mesh[nb,3],
                         linewidth=0.5, edgecolor='lightgray', facecolor='none',
                         alpha=0.25 )
  ax.add_patch(rect)

#-------------------------------------------------------------------------------
ax=axs[0,1]
vmin = 0 #1e-4
vmax = 1
im = ax.imshow( y0 , cmap='Spectral', extent = [-25,25,-25,25],
               vmin=vmin, vmax=vmax, origin='lower')
#           norm=LogNorm(vmin=vmin, vmax=vmax ), origin='lower' )
ax.set_title(f'Neutral fraction')
cb = plt.colorbar(im, ax=ax, extend=colorbar_ext(y0, vmin, vmax), shrink=0.7)

#-------------------------------------------------------------------------------
ax=axs[0,2]
vmin = 1e6
vmax = 1e8
im = ax.imshow(Pg, cmap='magma', extent = [-25,25,-25,25],
           norm=LogNorm(vmin=vmin,vmax=vmax), origin='lower' )
ax.set_title(r'Pressure [ erg cm$^{-3}$]')
cb = plt.colorbar(im, ax=ax, extend=colorbar_ext(Pg, vmin, vmax), shrink=0.7)

X,Z = np.meshgrid( np.linspace(-25.0,25.0,num=Pg.shape[0]),
                   np.linspace(-25.0,25.0,num=Pg.shape[1]) )
vec_step = 43
vec = ax.quiver( X[::vec_step,::vec_step],  Z[::vec_step,::vec_step],
                 vx[::vec_step,::vec_step], vz[::vec_step,::vec_step],
                 color='gray', alpha=0.5, width=0.005, scale =3000,
                 pivot='mid', cmap='inferno_r')

#-------------------------------------------------------------------------------
ax=axs[0,3]
vmin = 1e4
vmax = 1e6
im = ax.imshow(Temp, cmap='gist_heat', extent = [-25,25,-25,25],
           norm=LogNorm(vmin=vmin,vmax=vmax), origin='lower' )
ax.set_title(f'Temperature [ K ]')
plt.colorbar(im, ax=ax, extend=colorbar_ext(Temp, vmin, vmax), shrink=0.7)

#-------------------------------------------------------------------------------
ax = axs[1,0]
vmin = 1e-3
vmax = 10
im = ax.imshow(Tau, cmap='viridis', extent = [-25,25,-25,25], \
           norm=LogNorm(vmin=vmin,vmax=vmax), origin='lower' )
cb = plt.colorbar(im, ax=ax ,extend=colorbar_ext(Tau, vmin, vmax), shrink=0.7)
ax.set_title(r'$\tau$ (from file)')
ax.set_ylabel(f'Z [$R_{{p}}$]')

#-------------------------------------------------------------------------------
ax = axs[1,1]
vmin= 1e-11  #phi2.min()
vmax= 1e-8   #phi2.max()
im = ax.imshow(phi2, cmap='cividis', extent = [-25,25,-25,25], \
           norm=LogNorm(vmin=vmin, vmax=vmax, clip=True) , origin='lower' )
cb = plt.colorbar(im, ax=ax, extend=colorbar_ext(phi2, vmin, vmax), shrink = 0.7 )
ax.set_title(r'$\phi$ (post-process) [ s$^{-1}$ ]')

#-------------------------------------------------------------------------------
ax = axs[1,2]
vmin= 2e10
vmax= 2e13
F02 = F0 * np.exp(-Tau)
im = ax.imshow(F02  , cmap='inferno', extent = [-25,25,-25,25], \
           norm=LogNorm(vmin=vmin, vmax=vmax, clip=True), origin='lower' )
cb= plt.colorbar(im, ax=ax, extend=colorbar_ext(F02, vmin, vmax), shrink=0.7)
ax.set_title(r'$F$ (post-process) [ cm$^{-2}$ s$^{-1}$ ]')

#-------------------------------------------------------------------------------
ax = axs[1,3]
vmin= 1e-14    #psi2.max()/10
vmax= 1e-10    #psi2.max()
im = ax.imshow(psi2, cmap='inferno', extent = [-25,25,-25,25], \
           norm=LogNorm(vmin=vmin, vmax=vmax, clip=True), origin='lower' )
cb = plt.colorbar(im, ax=ax, extend=colorbar_ext(psi2, vmin, vmax), shrink=0.7)
ax.set_title(r'$\psi$ (post-process) [ erg cm$^{-3}$ s$^{-1}$]')

#-------------------------------------------------------------------------------
ax = axs[2,0]
vmin = 1e-3
vmax = 10
im = ax.imshow(Tau2, cmap='viridis', extent = [-25,25,-25,25], \
           norm=LogNorm(vmin=vmin, vmax=vmax), origin='lower' )
cb = plt.colorbar(im, ax=ax, extend=colorbar_ext(Tau2, vmin, vmax),shrink=0.7 )
ax.set_title(r'$\tau$ (post-process (2))')
ax.set_ylabel(f'Z [$R_{{p}}$]')

#-------------------------------------------------------------------------------
ax = axs[2,1]
vmin= 1e-11
vmax= 1e-8
im = ax.imshow(phi, cmap='cividis', extent = [-25,25,-25,25], \
           norm=LogNorm(vmin=vmin, vmax=vmax, clip=True), origin='lower' )
cb = plt.colorbar(im, ax=ax, extend=colorbar_ext(phi, vmin, vmax), shrink=0.7)
ax.set_title(r'$\phi$ (post-process 2) [ s$^{-1}$ ]')
ax.set_xlabel(f'X [$R_{{p}}$]')

#-------------------------------------------------------------------------------
ax = axs[2,2]
vmin= 2e10
vmax= 2e13
im = ax.imshow(flux, cmap='inferno', extent = [-25,25,-25,25], \
           norm=LogNorm(vmin=vmin, vmax=vmax, clip=True), origin='lower' )
cb = plt.colorbar(im, ax=ax, extend=colorbar_ext(flux, vmin, vmax), shrink=0.7)
ax.set_title(r'$F$ (post-process 2) [ cm$^{-2}$ s$^{-1}$ ]')
ax.set_xlabel(f'X [$R_{{p}}$]')

#-------------------------------------------------------------------------------
ax = axs[2,3]
vmin= 1e-14     #psi.max()/10
vmax= 1e-10     #psi.max()
im = ax.imshow(psi, cmap='inferno', extent = [-25,25,-25,25], \
           norm=LogNorm(vmin=vmin, vmax=vmax, clip=True), origin='lower' )
cb = plt.colorbar(im, ax=ax, extend=colorbar_ext(psi, vmin, vmax), shrink=0.7)
ax.set_title(r'$\psi$ (post-process 2) [ erg cm$^{-3}$ s$^{-1}$]')
ax.set_xlabel(f'X [$R_{{p}}$]')

#-------------------------------------------------------------------------------
#plt.show()
