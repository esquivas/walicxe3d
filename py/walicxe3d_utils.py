#!/usr/bin/python
# walicxe3d_utils.py
################################################################################
# This file contains a series of utilities to read and analyze the output
# from WALICXE-3D
################################################################################
import numpy as np
import os
from scipy.ndimage import zoom
from numpy.typing import NDArray


#-------------------------------------------------------------------------------
# print min and max value of an array (q)
def minmax(q):
  print(f"min = {q.min():10.5g},  max = {q.max():10.5g}")

#-------------------------------------------------------------------------------
# Returns the 'correct' extend to use in colorbar
def colorbar_ext(field, vmin, vmax):
  low  = True if( field.min() < vmin) else False
  high = True if( field.max() > vmax) else False

  if (low and high):
    extend = 'both'
  elif (low) :
    extend = 'min'
  elif (high) :
    extend = 'max'
  else :
    extend = 'neither'
  return extend

#-------------------------------------------------------------------------------
#  Reads 2D cuts written by extract
def read_cut(filename: str, neq:int, nx: int, ny: int, nz:int, axis: str='Z' ):
  #  filename : name of file to read (including path)
  #  neq      : number of equations stored by extract utility
  #  nx       : number of cells (max resolution) in x axis
  #  ny       : number of cells (max resolution) in y axis
  #  nz       : number of cells (max resolution) in z axis
  # axis      : 'X', 'Y', 'Z' axis perpendicular to cut
  f = open(filename,'rb')
  print('Reading file: ',filename)
  print(' ')
  match axis:
    case 'X':
      nxp = nx
      nyp = ny
    case 'Y':
      nxp = nx
      nyp = nz
    case 'Z':
      nxp = nx
      nyp = ny
    case _:
      print('Invalid axis')
  data = np.fromfile( f, dtype = 'd',count=(neq*nxp*nyp) ).\
    reshape(neq, nxp , nyp , order='F').T
  f.close()
  return data

#-------------------------------------------------------------------------------
# Read, scale and center mesh
def read_mesh(filemesh: str, dx: float, center: bool = False, \
              xc: float = 0, yc: float = 0):
  # filemesh : name of file to read (including path)
  # dx       : cell spacing
  # center   : (bool) whether to center mesh
  # xc       : x position of center from corner
  # yc       : y position of center from corner
  # determine how many blocks are stored in file
  n_blocks = os.path.getsize(filemesh) // (4 * np.dtype(np.int32).itemsize )
  # Read data
  mesh = np.fromfile(filemesh, dtype=np.int32).reshape(n_blocks, 4)
  # scale
  mesh = mesh * dx
  # center
  if (center) :
    mesh[:,0] -= xc
    mesh[:,1] -= yc
  return mesh

#-------------------------------------------------------------------------------
# returns the level of block bID
def get_level(bID: int, maxlev: int, nbrootx: int, nbrooty:int, nbrootz: int):
  # bID    : block ID
  # maxlev : maximum level
  # nbrootx : number of root blocks in x
  # nbrooty : number of root blocks in y
  # nbrootz : number of root blocks in z
  maxID = 0
  for level in range(1,maxlev+1):
    minID = maxID + 1
    maxID = maxID + nbrootx*nbrooty*nbrootz*8**(level-1)
    if (bID >= minID and bID <= maxID):
      return level

#-------------------------------------------------------------------------------
# Returns the offset associated with respect of level
# the first block on level l will have bID = offset(l) + 1
def levelOffset(level:int, nbrootx:int, nbrooty:int, nbrootz:int):
  # level   : refinement level
  # nbrootx : number of root blocks in x
  # nbrooty : number of root blocks in y
  # nbrootz : number of root blocks in z
  bOffset = 0
  for ilev in range(1, level):
    bOffset += nbrootx*nbrooty*nbrootz*8**(ilev -1)
  return bOffset

#-------------------------------------------------------------------------------
# Returns the block coordinates of a block at the block's mesh level
def bCoords(bID: int, maxlev:int, nbrootx:int, nbrooty:int, nbrootz:int):
  # bID     : Block (absolute) ID
  # maxlev  : maximum level
  # nbrootx : number of root blocks in x
  # nbrooty : number of root blocks in y
  # nbrootz : number of root blocks in z
  level  = get_level(bID, maxlev, nbrootx, nbrooty, nbrootz)
  nx = nbrootx*2**(level-1)
  ny = nbrooty*2**(level-1)
  nz = nbrootz*2**(level-1)
  localID = bID - levelOffset(level, nbrootx, nbrooty, nbrootz)
  ix =  (localID - 1)               % nx + 1
  iy = ((localID - 1) // nx)        % ny + 1
  iz = ((localID - 1) // (nx * ny)) % nz + 1
  return ix, iy, iz

#-------------------------------------------------------------------------------
# Extract equation from data, and convert to primitives if needed
def u2prim(Ublock: NDArray[np.float64], neq: int, CV: float, \
           conserved: bool = False):
  # Ublock: 4 dimensional arrayt containing a complete block of Us
  # neq   : Number of equation to retrieve (-1 relative to the value in fortran)
  # CV    : Specific heat at constant volume

  if conserved :
    data = Ublock[neq] # missing axes are implicitly :
  else:
    Pblock = Ublock
    if (neq == 0):
      # density
      data = Pblock[neq]
    elif (neq >=1 and neq<=3):
      # velocity components
      data = Ublock[neq]/Ublock[0]
    elif (neq == 4):
      # thermal pressure
      Pblock[1:4,:] = Ublock[1:4,:]/Ublock[0,:]
      data = ( Ublock[4] - 0.5*Pblock[0]* \
               ( Pblock[1]**2 + Pblock[2]**2 + Pblock[3]**2 ) ) / CV
    else:
      data = Pblock[neq]

  return data

#-------------------------------------------------------------------------------
# Returns one variable ***** at this point the variables are unscaled ****
def read_3d_eq(path: str, nout: int, neq: int, neqtot: int, maxlev: int, \
               ncellsB : int, nbrootx : int, nbrooty: int, nbrootz : int,\
               nprocs  : int, CV: float, conserved: bool = False,        \
               verbose : bool = True):
  # path    : full path to output directory
  # nout    : number of oputput to read
  # neq     : number of equation to retrieve (-1 relative to fortran)
  # neqtot  : total number of equations
  # maxlev  : maximum level of refinement
  # ncellsB : number of cels per side per block
  # nbrootx : number of root blocks in x
  # nbrooty : number of root blocks in y
  # nbrootz : number of root blocks in z
  # nprocs  : number of processors in the sim
  # CV    : Specific heat at constant volume
  # conserved : returns conserved variables (if false -> returns primitives)
  # verbose : set to False to inhibit screen output

  if verbose: print(f'Reading eqn. {neq}, interpolated to maximum resolution')

  nx = nbrootx*ncellsB*2**(maxlev-1)
  ny = nbrooty*ncellsB*2**(maxlev-1)
  nz = nbrootz*ncellsB*2**(maxlev-1)

  # Output array
  field = np.zeros(shape=(nz,ny,nx), dtype = 'd')

  # Loop over each procesor's output
  for proc in range(nprocs):

    filename = path + f'Blocks{proc:03d}.{nout:04d}.bin'
    if verbose: print('Processing:', filename)

    f = open(filename,'rb')
    nblocks = np.fromfile(f, dtype=np.int32, count=1)[0]
    blocksize = neqtot * ncellsB**3

    #  loops over all blocks
    for nb in range(nblocks):
        # read data and retrieve equation neq
        bID     = np.fromfile(f, dtype=np.int32, count=1)[0]
        block   = np.fromfile(f, dtype='d', count = blocksize ).\
          reshape(neqtot, ncellsB, ncellsB, ncellsB , order='F')
        #data    = block[neq,:,:,:]
        data = u2prim( block[:,:,:,:], neq, CV, conserved=conserved )

        # Get level and zoom factor
        level = get_level(bID, maxlev, nbrootx, nbrooty, nbrootz)
        ZoomFact = 2**(maxlev-level)

        # Zoom if needed
        if (level < maxlev):
          data = zoom( data, ZoomFact, order=0 )

        #  Determine target bounds
        coords = bCoords(bID, maxlev, nbrootx, nbrooty, nbrootz)
        i0 = (coords[0]-1) * ncellsB*ZoomFact ; i1 = i0 +ZoomFact*ncellsB
        j0 = (coords[1]-1) * ncellsB*ZoomFact ; j1 = j0 +ZoomFact*ncellsB
        k0 = (coords[2]-1) * ncellsB*ZoomFact ; k1 = k0 +ZoomFact*ncellsB

        #  Copy data
        field[i0:i1, j0:j1, k0:k1] = data[:,:,:]
    f.close()
  if verbose: print('')
  return field

#-------------------------------------------------------------------------------
# Returns a 2D cut of one variable ** unscaled **
def read_2d_cut(path: str, nout: int, neq: int, neqtot: int, cut: int,\
                maxlev: int, ncellsB : int, nbrootx : int, nbrooty: int, \
                nbrootz : int, nprocs  : int, CV: float, axis: str='Z',\
                conserved: bool = False, verbose : bool = True):
  # path    : full path to output directory
  # nout    : number of oputput to read
  # neq     : number of equation to retrieve (-1 relative to the fortran)
  # neqtot  : total number of equations
  # cut     : cell (integer number @ highest res.) of 2D cut along 'axis'
  # maxlev  : maximum level of refinement
  # ncellsB : number of cels per side per block
  # nbrootx : number of root blocks in x
  # nbrooty : number of root blocks in y
  # nbrootz : number of root blocks in z
  # nprocs  : number of processors in the sim
  # CV      : Specific heat at constant volume
  # axis    : axis:  'X'-->YZ plane, 'Y'-->XZ plane, 'Z'-->XY plane
  # conserved : returns conserved variables (default false -> returns primitives)
  # verbose : set to False to inhibit screen output

  if verbose:
    print(f'Reading 2D of eqn {neq} cut nterpolated to maximum resolution')
    print(f'Position: {cut:0d} along the {axis} axis')

  nx = nbrootx*ncellsB*2**(maxlev-1)
  ny = nbrooty*ncellsB*2**(maxlev-1)
  nz = nbrootz*ncellsB*2**(maxlev-1)
  match axis:
    case 'X':
      nxp = nx
      nyp = ny
    case 'Y':
      nxp = nx
      nyp = nz
    case 'Z':
      nxp = nx
      nyp = ny
    case _:
      print('Invalid axis')

  # Output array
  field = np.zeros(shape=(nyp,nxp), dtype = 'd')

    # Loop over each procesor's output
  for proc in range(nprocs):

    filename = path + f'Blocks{proc:03d}.{nout:04d}.bin'
    if verbose: print('Processing:', filename)

    f = open(filename,'rb')
    nblocks = np.fromfile(f, dtype=np.int32, count=1)[0]
    blocksize = neqtot * ncellsB**3

    #  loops over all blocks
    for nb in range(nblocks):
      # read data
      bID     = np.fromfile(f, dtype=np.int32, count=1)[0]
      block   = np.fromfile(f, dtype='d', count = blocksize ).\
          reshape(neqtot, ncellsB, ncellsB, ncellsB , order='F')

      level = get_level(bID, maxlev, nbrootx, nbrooty, nbrootz)
      ZoomFact = 2**(maxlev-level)

      #  Determine target bounds
      coords = bCoords(bID, maxlev, nbrootx, nbrooty, nbrootz)
      i0 = (coords[0]-1) * ncellsB*ZoomFact ; i1 = i0 +ZoomFact*ncellsB
      j0 = (coords[1]-1) * ncellsB*ZoomFact ; j1 = j0 +ZoomFact*ncellsB
      k0 = (coords[2]-1) * ncellsB*ZoomFact ; k1 = k0 +ZoomFact*ncellsB

      useBlock = False

      match axis:
        case 'X':
          if ( cut >= i0 and cut < i1):
            useBlock = True
            cutloc = (cut - i0) // ZoomFact
            slice = block[:,cutloc,:,:]

        case 'Y':
          if ( cut >= j0 and cut < j1):
            useBlock = True
            cutloc = (cut - j0 ) // ZoomFact
            slice = block[:,:,cutloc,:]

        case 'Z':
          if ( cut >= k0 and cut < k1):
            useBlock = True
            cutloc = (cut - k0) // ZoomFact
            slice = block[:,:,:,cutloc]

        case _:
          print('Invalid axis')

      if (useBlock) :

        data = u2prim( slice, neq, CV, conserved=conserved )
        # Zoom if needed
        if (level < maxlev):
          data = zoom( data, ZoomFact, order=0 )

        #  Copy data
        match axis:
          case 'X':
            #  Returns YZ plane
            field[ j0:j1, k0:k1 ] = data[:,:]
          case 'Y':
            #  Returns XZ plane
            field[ i0:i1, k0:k1 ] = data[:,:]
          case 'Z':
            #  Returns YX.T = XY plane
            data = data.T
            field[ j0:j1, i0:i1 ] = data[:,:]
            #field = field.T
          case _:
            print('Invalid axis')
  if verbose: print('')
  return field.T

#-------------------------------------------------------------------------------
# Returns the (blocks) mesh of a 2D cut
def read_mesh_2D_cut(path: str, nout: int, neqtot: int, cut: int, maxlev: int, \
                     ncellsB : int, nbrootx : int, nbrooty: int, nbrootz : int,\
                     nprocs  : int, axis: str='Z', verbose : bool = True):
  # path    : full path to output directory
  # nout    : number of oputput to read
  # neqtot  : total number of equations
  # cut     : cell (integer number @ highest res.) of 2D cut along 'axis'
  # maxlev  : maximum level of refinement
  # ncellsB : number of cels per side per block
  # nbrootx : number of root blocks in x
  # nbrooty : number of root blocks in y
  # nbrootz : number of root blocks in z
  # nprocs  : number of processors in the sim
  # axis    : axis: 'X'-->YZ plane, 'Y'-->XZ plane, 'Z'-->XY plane
  # verbose : set to False to inhibit screen output

  if verbose:
    print(f'Reading mesh in a 2D cut ')
    print(f'Position: {cut:0d} along the {axis} axis ')

  mesh = []

    # Loop over each procesor's output
  for proc in range(nprocs):

    filename = path + f'Blocks{proc:03d}.{nout:04d}.bin'
    if verbose: print('Processing:', filename)

    f = open(filename,'rb')
    nblocks = np.fromfile(f, dtype=np.int32, count=1)[0]
    blocksize = neqtot * ncellsB**3

    #  loops over all blocks
    for nb in range(nblocks):
      # read data
      bID     = np.fromfile(f, dtype=np.int32, count=1)[0]
      block   = np.fromfile(f, dtype='d', count = blocksize )

      level = get_level(bID, maxlev, nbrootx, nbrooty, nbrootz)
      ZoomFact = 2**(maxlev-level)

      #  Determine target bounds
      coords = bCoords(bID, maxlev, nbrootx, nbrooty, nbrootz)
      i0 = (coords[0]-1) * ncellsB*ZoomFact ; i1 = i0 +ZoomFact*ncellsB
      j0 = (coords[1]-1) * ncellsB*ZoomFact ; j1 = j0 +ZoomFact*ncellsB
      k0 = (coords[2]-1) * ncellsB*ZoomFact ; k1 = k0 +ZoomFact*ncellsB

      match axis:
        case 'X':
          if ( cut >= i0 and cut < i1):
           b_size = ncellsB*ZoomFact
           mesh.append([j0, k0, b_size, b_size])

        case 'Y':
          if ( cut >= j0 and cut < j1):
            b_size = ncellsB*ZoomFact
            mesh.append([i0, k0, b_size, b_size])

        case 'Z':
          if ( cut >= k0 and cut < k1):
            b_size = ncellsB*ZoomFact
            # notice transposiition
            mesh.append([i0, j0, b_size, b_size])

        case _:
          print('Invalid axis')
  if verbose: print()
  return np.asarray(mesh)
