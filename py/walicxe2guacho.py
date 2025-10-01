#!/usr/bin/python
import numpy as np
import struct
from constants import RJUP, AMH, KB
from walicxe3d_utils import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from guacho_utils import *


################################################################################
nout       =   2    # Output number
model_name = 'M2'   # Model na
# Simulations paramweters (Walicxe)
path_in = '/storage5/esquivel/wexo/'+ model_name+'/output/'
nxtot   = 1024   # Number of cells @highest resolution in x axis
nytot   =  256   # Number of cells @highest resolution in y axis
nztot   = 1024   # Number of cells @highest resolution in z axis
neqtot  = 7      # Total number of equations
neqdyn  = 5      # Number of dynamical equations
nghost  = 2      # Number of host cells

nbrootx_w = 4    # Number of MPI blocks in x axis
nbrooty_w = 1    # Number of MPI blocks in y axis
nbrootz_w = 4    # Number of MPI blocks in z axis
maxlev  = 5      # Maximum level of refinement
ncellsB = 16     # Number of cells per side, per block

nprocs_w = 32    # Number of procs used

xsize = 50.0     # physical size in code units in x axis
ysize = 12.5     # physical size in code units in y axis
zsize = 50.0     # physical size in code units in z axis

#  scalings
rsc = 1.98*RJUP  # length scale (cm)
vsc = 1.0e5      # velocity scale (cm s^-1)
rhosc = 1.0*AMH  # density scale (g cm^-3)

Cv=  1.5         # 1/(gamma-1)
#-------------------------------------------------------------------------------
# Path for guacho files
path_out = '/storage5/esquivel/wexo/'+ model_name+'/BIN/'

# Fill the info that is contained in the header for guacho
mpi_x_g = 4  # number of mpi blocks in x axis
mpi_y_g = 1  # number of mpi blocks in y axis
mpi_z_g = 4  # number of mpi blocks in z axis

nprocs_g = mpi_x_g * mpi_y_g * mpi_z_g  #  number of blocks (files) for guacho

#  sizes (in cells) and spacings (in code units)
nx = nxtot // mpi_x_g ; ny = nytot // mpi_y_g ; nz = nztot // mpi_z_g
dx = xsize / nxtot    ; dy = ysize / nytot    ; dz = zsize / nztot
f_kind = 'd'  # default in both codes is double precision

# Pack the info into a header holder
header = ( 0, f_kind,(nx,ny,nz),(dx,dy,dz),[0,0,0], \
           (nbrootx_w, nbrooty_w, nbrootz_w, ), neqtot ,neqdyn,nghost,\
           (rsc,vsc,rhosc), Cv )

################################################################################
# Load ALL conserved variables interpolated to highest resolution

U_big=np.zeros(shape=(nztot+2*nghost ,nytot+2*nghost, nxtot+2*nghost, neqtot))
for neq in range(neqtot):

    data_in = read_3d_eq(path_in, nout, neq, neqtot, maxlev, ncellsB, \
                         nbrootx_w, nbrooty_w, nbrootz_w, nprocs_w, Cv, \
                         verbose = True, conserved=True )
    U_big[nghost:-nghost,nghost:-nghost,nghost:-nghost,neq]= data_in[:,:,:]

plt.ion()
plt.figure(1) ; plt.clf()
# qvoid plotting ghost cells for big array
plt.imshow(U_big[nghost:-nghost,128,nghost:-nghost,0].T, \
           norm=LogNorm(vmin=100, vmax=1e5), origin='lower', cmap = 'cividis' )
plt.colorbar()

################################################################################

def distributed_write(U_big, nout, header, path_out=path_out ,\
                      base_out='points',verbose=True):

    proc = 0
    # unpack some header info
    nx,    ny,    nz    = header[2]
    mpi_x, mpi_y, mpi_z = header[5]
    neqs                = header[6]

    for ip in range(mpi_x) :
      for jp in range(mpi_y) :
        for kp in range(mpi_z) :

          file_out = path_out+base_out+str(proc).zfill(3)+\
                     '.'+str(nout).zfill(3)+'.bin'
          if verbose: print (f'Writing {file_out}')

          # adjust position of mpi blocks
          header[4][0] = nx * ip
          header[4][1] = ny * jp
          header[4][2] = nz * kp
          x0, y0, z0   = header[4]

          data_out = np.zeros(shape=(nz+2*nghost,ny+2*nghost,nx+2*nghost,neqs) )
          data_out[:,:,:,:]=\
                 U_big[z0:z0+nz+2*nghost,y0:y0+ny+2*nghost,x0:x0+nx+2*nghost,:]

          copy_header(header , file_out, verbose=verbose)

          f = open(file_out,"ab")
          f.write(data_out.tobytes() )
          f.close()
          proc += 1
    return
################################################################################

# replicate the original headers
def copy_header(header, file_out, verbose=True):

    f = open(file_out,"w")
    s = "**************** Output for Guacho v1.3****************\n"
    s += "Dimensions    : "+ str(header[2][0])+" "+str(header[2][1])+" "\
                           + str(header[2][2])+"\n"
    s += "Spacings      : "+ str(header[3][0])+" "+str(header[3][1])+" "\
                           + str(header[3][2])+"\n"
    s += "Block Origin, cells    : " + str(header[4][0])+" "+str(header[4][1])\
                                 +" "+ str(header[4][2])+"\n"
    s += "MPI blocks (X, Y, Z)   : " + str(header[5][0])+" "+str(header[5][1])+\
                                  " "+ str(header[5][2])+"\n"
    s += "Number of Equations/dynamical ones  " + str(header[6]) + "/" \
                                                + str(header[7])+"\n"
    s += "Number of Ghost Cells  " + str(header[8])+"\n"
    s += "Scalings " + str(header[9][0])+" "+str(header[9][1])+" " \
                     + str(header[9][2])+"\n"
    s += "Specfic heat at constant volume Cv: " + str(header[10])+"\n"
    if(header[1]=='d'):
        s += "Double precision 8 byte floats"+"\n"
    else:
        s += "single precision 4 byte floats"+"\n"
    s += "*******************************************************"+"\n"
    f.write(s)
    f.close()
    f = open(file_out,"ab")
    f.write(b"\xff")
    f.write(b"\n")
    f.write(b'd')
    chunks = []
    chunks += [struct.pack('<i', int(i))   for i in header[2]]
    chunks += [struct.pack('<d', float(x)) for x in header[3]]
    chunks += [struct.pack('<i', int(i))   for i in header[4]]
    chunks += [struct.pack('<i', int(i))   for i in header[5]]
    chunks += [struct.pack('<i', int(header[6])),
               struct.pack('<i', int(header[7])),
               struct.pack('<i', int(header[8]))]
    chunks += [struct.pack('<d', float(x)) for x in header[9]]
    chunks += [struct.pack('<d', float(header[10]))]

    binary_data = b''.join(chunks)

    f.write(binary_data)
    f.close()
    return

################################################################################

# write points for guacho
distributed_write(U_big,nout,header,path_out=path_out,verbose=False)

# free memory
U_big = []

################################################################################

# verify reading the blocks created``
dens = readbin3d_all(nout,0,path=path_out, scale=False )

plt.figure(2) ; plt.clf()
# qvoid plotting ghost cells for big array
plt.imshow(dens[:,128,:].T, norm=LogNorm(vmin=100, vmax=1e5), \
           origin='lower', cmap = 'cividis' )
plt.colorbar()