#!/usr/bin/env python3
# walicxe3d_utils.py
################################################################################
# This file contains a series of utilities to read and analyze the output
# from WALICXE-3D
################################################################################
import os
import math

import numpy as np


#-------------------------------------------------------------------------------
# Internal helpers

def _validate_axis(axis: str) -> str:
  axis = axis.upper()
  if axis not in ('X', 'Y', 'Z'):
    raise ValueError(f"Invalid axis {axis!r}; expected 'X', 'Y' or 'Z'")
  return axis


def _grid_shape(maxlev: int, ncellsB: int, nbrootx: int, nbrooty: int,
                nbrootz: int):
  if maxlev < 1:
    raise ValueError('maxlev must be >= 1')
  if ncellsB <= 0:
    raise ValueError('ncellsB must be > 0')
  if nbrootx <= 0 or nbrooty <= 0 or nbrootz <= 0:
    raise ValueError('nbrootx, nbrooty and nbrootz must be > 0')

  scale = 2**(maxlev - 1)
  nx = nbrootx * ncellsB * scale
  ny = nbrooty * ncellsB * scale
  nz = nbrootz * ncellsB * scale
  return nx, ny, nz


def _validate_equation(neq: int, neqtot: int):
  if neqtot <= 0:
    raise ValueError('neqtot must be > 0')
  if neq < 0 or neq >= neqtot:
    raise ValueError(f'neq={neq} is outside the valid range [0, {neqtot - 1}]')


def _read_exact(f, dtype, count: int, what: str):
  data = np.fromfile(f, dtype=dtype, count=count)
  if data.size != count:
    raise EOFError(
      f'Unexpected end of file while reading {what}: '
      f'expected {count} values, found {data.size}'
    )
  return data


def _prolong_piecewise_constant(data: np.ndarray, factor: int):
  """Prolong cell-centered data by exact piecewise-constant replication."""
  if factor < 1:
    raise ValueError('prolongation factor must be >= 1')
  if factor == 1:
    return data

  out = data
  for axis in range(data.ndim):
    out = np.repeat(out, factor, axis=axis)
  return out


def _skip_exact(f, nbytes: int, what: str):
  """Advance exactly nbytes, rejecting truncated files."""
  if nbytes < 0:
    raise ValueError('nbytes must be >= 0')
  start = f.tell()
  end = start + nbytes
  filesize = os.fstat(f.fileno()).st_size
  if end > filesize:
    raise EOFError(
      f'Unexpected end of file while skipping {what}: '
      f'needed {nbytes} bytes, only {filesize - start} remain'
    )
  f.seek(nbytes, os.SEEK_CUR)


#-------------------------------------------------------------------------------
# print min and max value of an array (q)
def minmax(q):
  print(f"min = {q.min():10.5g},  max = {q.max():10.5g}")


#-------------------------------------------------------------------------------
# Returns the 'correct' extend to use in colorbar
def colorbar_ext(field, vmin, vmax):
  low  = field.min() < vmin
  high = field.max() > vmax

  if low and high:
    extend = 'both'
  elif low:
    extend = 'min'
  elif high:
    extend = 'max'
  else:
    extend = 'neither'
  return extend


#-------------------------------------------------------------------------------
# Reads 2D cuts written by extract
#
# Returned array convention:
#   axis='X' -> (z, y, equation)
#   axis='Y' -> (z, x, equation)
#   axis='Z' -> (y, x, equation)
def read_cut(filename: str, neq: int, nx: int, ny: int, nz: int,
             axis: str = 'Z'):
  # filename : name of file to read (including path)
  # neq      : number of equations stored by extract utility
  # nx       : number of cells (max resolution) in x axis
  # ny       : number of cells (max resolution) in y axis
  # nz       : number of cells (max resolution) in z axis
  # axis     : 'X', 'Y', 'Z' axis perpendicular to cut
  axis = _validate_axis(axis)
  if neq <= 0:
    raise ValueError('neq must be > 0')
  if nx <= 0 or ny <= 0 or nz <= 0:
    raise ValueError('nx, ny and nz must be > 0')

  # horizontal and vertical sizes in the returned plane
  if axis == 'X':
    nh = ny
    nv = nz
  elif axis == 'Y':
    nh = nx
    nv = nz
  else:  # axis == 'Z'
    nh = nx
    nv = ny

  print('Reading file: ', filename)
  print(' ')

  with open(filename, 'rb') as f:
    raw = _read_exact(f, dtype='d', count=neq * nh * nv,
                      what='2D cut data')

  # File ordering is (equation, horizontal, vertical) in Fortran order.
  # Return conventional image ordering (vertical, horizontal, equation).
  data = raw.reshape(neq, nh, nv, order='F').transpose(2, 1, 0)
  return data


#-------------------------------------------------------------------------------
# Read, scale and center mesh
def read_mesh(filemesh: str, dx: float, center: bool = False,
              xc: float = 0, yc: float = 0):
  # filemesh : name of file to read (including path)
  # dx       : cell spacing
  # center   : (bool) whether to center mesh
  # xc       : x position of center from corner
  # yc       : y position of center from corner
  itemsize = np.dtype(np.int32).itemsize
  filesize = os.path.getsize(filemesh)
  record_size = 4 * itemsize
  if filesize % record_size != 0:
    raise ValueError(
      f'Mesh file size ({filesize} bytes) is not a multiple of '
      f'{record_size} bytes per block record'
    )

  n_blocks = filesize // record_size
  mesh = np.fromfile(filemesh, dtype=np.int32).reshape(n_blocks, 4)

  # scale; multiplication by float intentionally promotes to floating point
  mesh = mesh * dx

  if center:
    mesh[:, 0] -= xc
    mesh[:, 1] -= yc
  return mesh


#-------------------------------------------------------------------------------
# returns the level of block bID
def get_level(bID: int, maxlev: int, nbrootx: int, nbrooty: int, nbrootz: int):
  # bID     : block ID
  # maxlev  : maximum level
  # nbrootx : number of root blocks in x
  # nbrooty : number of root blocks in y
  # nbrootz : number of root blocks in z
  if bID < 1:
    raise ValueError('bID must be >= 1')
  if maxlev < 1:
    raise ValueError('maxlev must be >= 1')
  if nbrootx <= 0 or nbrooty <= 0 or nbrootz <= 0:
    raise ValueError('nbrootx, nbrooty and nbrootz must be > 0')

  maxID = 0
  for level in range(1, maxlev + 1):
    minID = maxID + 1
    maxID = maxID + nbrootx * nbrooty * nbrootz * 8**(level - 1)
    if minID <= bID <= maxID:
      return level

  raise ValueError(f'bID={bID} does not exist for maxlev={maxlev}')


#-------------------------------------------------------------------------------
# Returns the offset associated with respect of level
# the first block on level l will have bID = offset(l) + 1
def levelOffset(level: int, nbrootx: int, nbrooty: int, nbrootz: int):
  # level   : refinement level
  # nbrootx : number of root blocks in x
  # nbrooty : number of root blocks in y
  # nbrootz : number of root blocks in z
  if level < 1:
    raise ValueError('level must be >= 1')
  if nbrootx <= 0 or nbrooty <= 0 or nbrootz <= 0:
    raise ValueError('nbrootx, nbrooty and nbrootz must be > 0')

  bOffset = 0
  for ilev in range(1, level):
    bOffset += nbrootx * nbrooty * nbrootz * 8**(ilev - 1)
  return bOffset


#-------------------------------------------------------------------------------
# Returns the block coordinates of a block at the block's mesh level
def bCoords(bID: int, maxlev: int, nbrootx: int, nbrooty: int, nbrootz: int):
  # bID     : Block (absolute) ID
  # maxlev  : maximum level
  # nbrootx : number of root blocks in x
  # nbrooty : number of root blocks in y
  # nbrootz : number of root blocks in z
  level = get_level(bID, maxlev, nbrootx, nbrooty, nbrootz)
  nx = nbrootx * 2**(level - 1)
  ny = nbrooty * 2**(level - 1)

  localID = bID - levelOffset(level, nbrootx, nbrooty, nbrootz)
  ix =  (localID - 1)               % nx + 1
  iy = ((localID - 1) // nx)        % ny + 1
  iz = ((localID - 1) // (nx * ny)) + 1
  return ix, iy, iz


#-------------------------------------------------------------------------------
# Extract equation from data, and convert to primitives if needed
def u2prim(Ublock: np.ndarray, neq: int, CV: float,
           conserved: bool = False):
  # Ublock: array containing a complete block/slice of conserved variables.
  #         Axis 0 is the equation index; all following axes are spatial.
  # neq   : Number of equation to retrieve (-1 relative to the value in Fortran)
  # CV    : Specific heat at constant volume
  #
  # NOTE: For neq == 4 this preserves the original HD conversion:
  #       thermal energy = total energy - kinetic energy.
  #       No magnetic-energy subtraction is introduced here because the
  #       magnetic variable ordering/normalization is not specified by this
  #       utility.
  if neq < 0 or neq >= Ublock.shape[0]:
    raise ValueError(
      f'neq={neq} is outside the valid range [0, {Ublock.shape[0] - 1}]'
    )

  if conserved:
    return Ublock[neq]

  if neq == 0:
    # density
    return Ublock[0]

  if 1 <= neq <= 3:
    # velocity components
    return Ublock[neq] / Ublock[0]

  if neq == 4:
    # thermal pressure, following the original HD convention.
    # Compute without modifying Ublock.
    rho = Ublock[0]
    vx = Ublock[1] / rho
    vy = Ublock[2] / rho
    vz = Ublock[3] / rho
    kinetic = 0.5 * rho * (vx**2 + vy**2 + vz**2)
    return (Ublock[4] - kinetic) / CV

  # Additional variables are returned unchanged, as in the original code.
  return Ublock[neq]


#-------------------------------------------------------------------------------
# Returns one variable ***** at this point the variables are unscaled ****
#
# Returned array convention: field[z, y, x], shape = (nz, ny, nx)
def read_3d_eq(path: str, nout: int, neq: int, neqtot: int, maxlev: int,
               ncellsB: int, nbrootx: int, nbrooty: int, nbrootz: int,
               nprocs: int, CV: float, conserved: bool = False,
               verbose: bool = True):
  # path      : full path to output directory
  # nout      : number of output to read
  # neq       : number of equation to retrieve (-1 relative to Fortran)
  # neqtot    : total number of equations
  # maxlev    : maximum level of refinement
  # ncellsB   : number of cells per side per block
  # nbrootx   : number of root blocks in x
  # nbrooty   : number of root blocks in y
  # nbrootz   : number of root blocks in z
  # nprocs    : number of processors in the sim
  # CV        : specific heat at constant volume
  # conserved : returns conserved variables (if False -> returns primitives)
  # verbose   : set to False to inhibit screen output
  _validate_equation(neq, neqtot)
  nx, ny, nz = _grid_shape(maxlev, ncellsB, nbrootx, nbrooty, nbrootz)
  if nprocs <= 0:
    raise ValueError('nprocs must be > 0')

  if verbose:
    print(f'Reading eqn. {neq}, prolonged to maximum resolution')

  field = np.zeros(shape=(nz, ny, nx), dtype='d')
  blocksize = neqtot * ncellsB**3

  for proc in range(nprocs):
    filename = os.path.join(path, f'Blocks{proc:03d}.{nout:04d}.bin')
    if verbose:
      print('Processing:', filename)

    with open(filename, 'rb') as f:
      nblocks = int(_read_exact(f, np.int32, 1, 'number of blocks')[0])
      if nblocks < 0:
        raise ValueError(f'Invalid negative block count ({nblocks}) in {filename}')

      for nb in range(nblocks):
        bID = int(_read_exact(f, np.int32, 1,
                              f'block ID {nb} in {filename}')[0])
        block = _read_exact(
          f, 'd', blocksize, f'block {bID} data in {filename}'
        ).reshape(neqtot, ncellsB, ncellsB, ncellsB, order='F')

        data = u2prim(block, neq, CV, conserved=conserved)

        level = get_level(bID, maxlev, nbrootx, nbrooty, nbrootz)
        ZoomFact = 2**(maxlev - level)

        if level < maxlev:
          data = _prolong_piecewise_constant(data, ZoomFact)

        ix, iy, iz = bCoords(bID, maxlev, nbrootx, nbrooty, nbrootz)
        i0 = (ix - 1) * ncellsB * ZoomFact
        i1 = i0 + ZoomFact * ncellsB
        j0 = (iy - 1) * ncellsB * ZoomFact
        j1 = j0 + ZoomFact * ncellsB
        k0 = (iz - 1) * ncellsB * ZoomFact
        k1 = k0 + ZoomFact * ncellsB

        # block/data ordering is (x, y, z); field ordering is (z, y, x)
        field[k0:k1, j0:j1, i0:i1] = data.transpose(2, 1, 0)

  if verbose:
    print('')
  return field


#-------------------------------------------------------------------------------
# Returns a 2D cut of one variable ** unscaled **
#
# Returned array convention:
#   axis='X' -> field[z, y], shape = (nz, ny)  (YZ plane)
#   axis='Y' -> field[z, x], shape = (nz, nx)  (XZ plane)
#   axis='Z' -> field[y, x], shape = (ny, nx)  (XY plane)
def read_2d_cut(path: str, nout: int, neq: int, neqtot: int, cut: int,
                maxlev: int, ncellsB: int, nbrootx: int, nbrooty: int,
                nbrootz: int, nprocs: int, CV: float, axis: str = 'Z',
                conserved: bool = False, verbose: bool = True,
                return_scale: bool = False):
  # path      : full path to output directory
  # nout      : number of output to read
  # neq       : number of equation to retrieve (-1 relative to Fortran)
  # neqtot    : total number of equations
  # cut       : cell (integer number @ highest res.) of 2D cut along 'axis'
  # maxlev    : maximum level of refinement
  # ncellsB   : number of cells per side per block
  # nbrootx   : number of root blocks in x
  # nbrooty   : number of root blocks in y
  # nbrootz   : number of root blocks in z
  # nprocs    : number of processors in the sim
  # CV        : specific heat at constant volume
  # axis      : 'X'->YZ plane, 'Y'->XZ plane, 'Z'->XY plane
  # conserved : returns conserved variables (default False -> primitives)
  # verbose   : set to False to inhibit screen output
  # return_scale : if True, also return a map containing the local AMR
  #                coarse-cell size in finest-grid cells:
  #                1 at maxlev, 2 one level coarser, 4 two levels coarser, ...
  axis = _validate_axis(axis)
  _validate_equation(neq, neqtot)
  nx, ny, nz = _grid_shape(maxlev, ncellsB, nbrootx, nbrooty, nbrootz)
  if nprocs <= 0:
    raise ValueError('nprocs must be > 0')

  if axis == 'X':
    if cut < 0 or cut >= nx:
      raise ValueError(f'cut={cut} is outside X range [0, {nx - 1}]')
    field = np.zeros(shape=(nz, ny), dtype='d')
  elif axis == 'Y':
    if cut < 0 or cut >= ny:
      raise ValueError(f'cut={cut} is outside Y range [0, {ny - 1}]')
    field = np.zeros(shape=(nz, nx), dtype='d')
  else:  # axis == 'Z'
    if cut < 0 or cut >= nz:
      raise ValueError(f'cut={cut} is outside Z range [0, {nz - 1}]')
    field = np.zeros(shape=(ny, nx), dtype='d')

  scale_map = np.zeros(field.shape, dtype=np.int32) if return_scale else None

  if verbose:
    print(f'Reading 2D cut of eqn. {neq}, prolonged to maximum resolution')
    print(f'Position: {cut:d} along the {axis} axis')

  blocksize = neqtot * ncellsB**3

  for proc in range(nprocs):
    filename = os.path.join(path, f'Blocks{proc:03d}.{nout:04d}.bin')
    if verbose:
      print('Processing:', filename)

    with open(filename, 'rb') as f:
      nblocks = int(_read_exact(f, np.int32, 1, 'number of blocks')[0])
      if nblocks < 0:
        raise ValueError(f'Invalid negative block count ({nblocks}) in {filename}')

      for nb in range(nblocks):
        bID = int(_read_exact(f, np.int32, 1,
                              f'block ID {nb} in {filename}')[0])
        block = _read_exact(
          f, 'd', blocksize, f'block {bID} data in {filename}'
        ).reshape(neqtot, ncellsB, ncellsB, ncellsB, order='F')

        level = get_level(bID, maxlev, nbrootx, nbrooty, nbrootz)
        ZoomFact = 2**(maxlev - level)

        ix, iy, iz = bCoords(bID, maxlev, nbrootx, nbrooty, nbrootz)
        i0 = (ix - 1) * ncellsB * ZoomFact
        i1 = i0 + ZoomFact * ncellsB
        j0 = (iy - 1) * ncellsB * ZoomFact
        j1 = j0 + ZoomFact * ncellsB
        k0 = (iz - 1) * ncellsB * ZoomFact
        k1 = k0 + ZoomFact * ncellsB

        useBlock = False

        if axis == 'X':
          if i0 <= cut < i1:
            useBlock = True
            cutloc = (cut - i0) // ZoomFact
            block_slice = block[:, cutloc, :, :]
        elif axis == 'Y':
          if j0 <= cut < j1:
            useBlock = True
            cutloc = (cut - j0) // ZoomFact
            block_slice = block[:, :, cutloc, :]
        else:  # axis == 'Z'
          if k0 <= cut < k1:
            useBlock = True
            cutloc = (cut - k0) // ZoomFact
            block_slice = block[:, :, :, cutloc]

        if useBlock:
          data = u2prim(block_slice, neq, CV, conserved=conserved)

          if level < maxlev:
            data = _prolong_piecewise_constant(data, ZoomFact)

          # Slice data order before transpose:
          #   X cut -> (y, z)
          #   Y cut -> (x, z)
          #   Z cut -> (x, y)
          # Output uses conventional image order (vertical, horizontal).
          if axis == 'X':
            field[k0:k1, j0:j1] = data.T
            if return_scale:
              scale_map[k0:k1, j0:j1] = ZoomFact
          elif axis == 'Y':
            field[k0:k1, i0:i1] = data.T
            if return_scale:
              scale_map[k0:k1, i0:i1] = ZoomFact
          else:  # axis == 'Z'
            field[j0:j1, i0:i1] = data.T
            if return_scale:
              scale_map[j0:j1, i0:i1] = ZoomFact

  if verbose:
    print('')
  if return_scale:
    return field, scale_map
  return field


#-------------------------------------------------------------------------------
# Returns the (blocks) mesh of a 2D cut
def read_mesh_2D_cut(path: str, nout: int, neqtot: int, cut: int, maxlev: int,
                     ncellsB: int, nbrootx: int, nbrooty: int, nbrootz: int,
                     nprocs: int, axis: str = 'Z', verbose: bool = True):
  # path    : full path to output directory
  # nout    : number of output to read
  # neqtot  : total number of equations
  # cut     : cell (integer number @ highest res.) of 2D cut along 'axis'
  # maxlev  : maximum level of refinement
  # ncellsB : number of cells per side per block
  # nbrootx : number of root blocks in x
  # nbrooty : number of root blocks in y
  # nbrootz : number of root blocks in z
  # nprocs  : number of processors in the sim
  # axis    : 'X'->YZ plane, 'Y'->XZ plane, 'Z'->XY plane
  # verbose : set to False to inhibit screen output
  axis = _validate_axis(axis)
  if neqtot <= 0:
    raise ValueError('neqtot must be > 0')
  nx, ny, nz = _grid_shape(maxlev, ncellsB, nbrootx, nbrooty, nbrootz)
  if nprocs <= 0:
    raise ValueError('nprocs must be > 0')

  if axis == 'X':
    if cut < 0 or cut >= nx:
      raise ValueError(f'cut={cut} is outside X range [0, {nx - 1}]')
  elif axis == 'Y':
    if cut < 0 or cut >= ny:
      raise ValueError(f'cut={cut} is outside Y range [0, {ny - 1}]')
  else:  # axis == 'Z'
    if cut < 0 or cut >= nz:
      raise ValueError(f'cut={cut} is outside Z range [0, {nz - 1}]')

  if verbose:
    print('Reading mesh in a 2D cut')
    print(f'Position: {cut:d} along the {axis} axis')

  mesh = []
  blocksize = neqtot * ncellsB**3
  blockbytes = blocksize * np.dtype('d').itemsize

  for proc in range(nprocs):
    filename = os.path.join(path, f'Blocks{proc:03d}.{nout:04d}.bin')
    if verbose:
      print('Processing:', filename)

    with open(filename, 'rb') as f:
      nblocks = int(_read_exact(f, np.int32, 1, 'number of blocks')[0])
      if nblocks < 0:
        raise ValueError(f'Invalid negative block count ({nblocks}) in {filename}')

      for nb in range(nblocks):
        bID = int(_read_exact(f, np.int32, 1,
                              f'block ID {nb} in {filename}')[0])

        # The mesh routine needs only the block ID; skip the payload, but
        # still detect truncated files.
        _skip_exact(f, blockbytes, f'block {bID} data in {filename}')

        level = get_level(bID, maxlev, nbrootx, nbrooty, nbrootz)
        ZoomFact = 2**(maxlev - level)

        ix, iy, iz = bCoords(bID, maxlev, nbrootx, nbrooty, nbrootz)
        i0 = (ix - 1) * ncellsB * ZoomFact
        i1 = i0 + ZoomFact * ncellsB
        j0 = (iy - 1) * ncellsB * ZoomFact
        j1 = j0 + ZoomFact * ncellsB
        k0 = (iz - 1) * ncellsB * ZoomFact
        k1 = k0 + ZoomFact * ncellsB
        b_size = ncellsB * ZoomFact

        # Mesh records are [horizontal, vertical, width, height], suitable for
        # matplotlib.patches.Rectangle in the corresponding cut plane.
        if axis == 'X':
          if i0 <= cut < i1:
            mesh.append([j0, k0, b_size, b_size])
        elif axis == 'Y':
          if j0 <= cut < j1:
            mesh.append([i0, k0, b_size, b_size])
        else:  # axis == 'Z'
          if k0 <= cut < k1:
            mesh.append([i0, j0, b_size, b_size])

  if verbose:
    print('')
  return np.asarray(mesh).reshape(-1, 4)

#-------------------------------------------------------------------------------
# Finite-difference utilities
#
# These routines generate centered finite-difference stencils automatically.
# They are intended to replace hand-written gradient/Laplacian expressions.
#
# Conventions for 2D arrays:
#   axis 0 -> vertical coordinate (y or z, depending on the cut)
#   axis 1 -> horizontal coordinate (x or y, depending on the cut)

def fd_weights(derivative: int = 1, accuracy: int = 2):
  """Return centered finite-difference offsets and coefficients.

  Parameters
  ----------
  derivative : int
      Derivative order. Currently 1 or 2.
  accuracy : int
      Even formal accuracy order: 2, 4, 6, ...

  Returns
  -------
  offsets : ndarray of int
      Integer stencil offsets.
  coeffs : ndarray of float
      Coefficients for unit grid spacing. Divide by h**derivative.
  """
  if derivative not in (1, 2):
    raise ValueError('derivative must be 1 or 2')
  if accuracy < 2 or accuracy % 2 != 0:
    raise ValueError('accuracy must be a positive even integer >= 2')

  radius = accuracy // 2
  if derivative > 2 * radius:
    raise ValueError('stencil is too small for requested derivative')

  offsets = np.arange(-radius, radius + 1, dtype=float)

  # Taylor constraints:
  #   sum_j c_j x_j^m = m! for m == derivative, otherwise 0.
  matrix = np.empty((2 * radius + 1, 2 * radius + 1), dtype=float)
  rhs = np.zeros(2 * radius + 1, dtype=float)
  for m in range(2 * radius + 1):
    matrix[m, :] = offsets**m
  rhs[derivative] = float(math.factorial(derivative))
  coeffs = np.linalg.solve(matrix, rhs)

  return offsets.astype(np.int32), coeffs


def _fd_axis_2d(field: np.ndarray, axis: int, derivative: int = 1,
                accuracy: int = 2, step: int = 1, spacing: float = 1.0):
  """Centered finite difference along one axis of a 2D field.

  The result has the same shape as ``field``. Points where the full stencil
  does not fit are NaN. ``step`` is measured in array cells; the physical
  stencil spacing is ``step * spacing``.
  """
  field = np.asarray(field, dtype=float)
  if field.ndim != 2:
    raise ValueError('field must be a 2D array')
  if axis not in (0, 1):
    raise ValueError('axis must be 0 or 1')
  if step < 1:
    raise ValueError('step must be >= 1')
  if spacing <= 0:
    raise ValueError('spacing must be > 0')

  offsets, coeffs = fd_weights(derivative, accuracy)
  radius = int(np.max(np.abs(offsets)))
  margin = radius * step

  out = np.zeros(field.shape, dtype=float)
  for offset, coeff in zip(offsets, coeffs):
    out += coeff * np.roll(field, -int(offset) * step, axis=axis)

  out /= (step * spacing)**derivative

  # np.roll wraps around; explicitly invalidate those edge points.
  if margin > 0:
    if axis == 0:
      out[:margin, :] = np.nan
      out[-margin:, :] = np.nan
    else:
      out[:, :margin] = np.nan
      out[:, -margin:] = np.nan

  return out


def gradient_2d(field: np.ndarray, spacing=(1.0, 1.0), accuracy: int = 2,
                step: int = 1, return_components: bool = False):
  """Centered gradient of a 2D field.

  ``spacing`` is (vertical_spacing, horizontal_spacing). ``step`` is the
  stencil separation in array cells.
  """
  dy, dx = spacing
  gy = _fd_axis_2d(field, 0, derivative=1, accuracy=accuracy,
                   step=step, spacing=dy)
  gx = _fd_axis_2d(field, 1, derivative=1, accuracy=accuracy,
                   step=step, spacing=dx)
  magnitude = np.sqrt(gx**2 + gy**2)
  if return_components:
    return magnitude, gx, gy
  return magnitude


def laplacian_2d(field: np.ndarray, spacing=(1.0, 1.0), accuracy: int = 2,
                 step: int = 1):
  """Centered Cartesian Laplacian of a 2D field.

  ``accuracy=2`` gives the standard 5-point cross. ``accuracy=4`` gives the
  axial fourth-order stencil. Higher even orders are generated automatically.
  """
  dy, dx = spacing
  d2y = _fd_axis_2d(field, 0, derivative=2, accuracy=accuracy,
                    step=step, spacing=dy)
  d2x = _fd_axis_2d(field, 1, derivative=2, accuracy=accuracy,
                    step=step, spacing=dx)
  return d2x + d2y


def _same_scale_stencil_mask(scale_map: np.ndarray, axis: int, step: int,
                             accuracy: int):
  """Mask points whose full centered stencil stays at one AMR scale."""
  scale_map = np.asarray(scale_map)
  offsets, _ = fd_weights(1, accuracy)
  radius = int(np.max(np.abs(offsets)))
  margin = radius * step

  mask = (scale_map == step)
  for offset in offsets:
    mask &= (np.roll(scale_map, -int(offset) * step, axis=axis) == step)

  if margin > 0:
    if axis == 0:
      mask[:margin, :] = False
      mask[-margin:, :] = False
    else:
      mask[:, :margin] = False
      mask[:, -margin:] = False
  return mask


def amr_gradient_2d(field: np.ndarray, scale_map: np.ndarray,
                    spacing=(1.0, 1.0), accuracy: int = 2,
                    mask_level_interfaces: bool = True,
                    return_components: bool = False):
  """Gradient on a max-resolution AMR cut without differentiating replicas.

  ``scale_map`` must come from ``read_2d_cut(..., return_scale=True)``.
  Its value is the local coarse-cell width in finest-grid cells (1, 2, 4, ...).

  At scale ``s`` the derivative uses samples separated by ``s`` finest-grid
  cells instead of adjacent replicated samples.

  If ``mask_level_interfaces`` is True, points whose stencil crosses a
  refinement-level interface are returned as NaN.
  """
  field = np.asarray(field, dtype=float)
  scale_map = np.asarray(scale_map)
  if field.ndim != 2 or scale_map.shape != field.shape:
    raise ValueError('field and scale_map must be 2D arrays with equal shape')

  gy = np.full(field.shape, np.nan, dtype=float)
  gx = np.full(field.shape, np.nan, dtype=float)

  scales = np.unique(scale_map)
  scales = scales[scales > 0]

  for scale_value in scales:
    step = int(scale_value)
    if step != scale_value or step < 1:
      raise ValueError('scale_map values must be positive integer cell widths')

    gy_s = _fd_axis_2d(field, 0, derivative=1, accuracy=accuracy,
                       step=step, spacing=spacing[0])
    gx_s = _fd_axis_2d(field, 1, derivative=1, accuracy=accuracy,
                       step=step, spacing=spacing[1])

    mask = (scale_map == step)
    if mask_level_interfaces:
      mask_y = mask & _same_scale_stencil_mask(scale_map, 0, step, accuracy)
      mask_x = mask & _same_scale_stencil_mask(scale_map, 1, step, accuracy)
    else:
      mask_y = mask
      mask_x = mask

    gy[mask_y] = gy_s[mask_y]
    gx[mask_x] = gx_s[mask_x]

  magnitude = np.sqrt(gx**2 + gy**2)
  if return_components:
    return magnitude, gx, gy
  return magnitude


def amr_laplacian_2d(field: np.ndarray, scale_map: np.ndarray,
                     spacing=(1.0, 1.0), accuracy: int = 2,
                     mask_level_interfaces: bool = True):
  """Cartesian Laplacian on a max-resolution AMR cut.

  The local stencil spacing follows ``scale_map``. If
  ``mask_level_interfaces`` is True, points whose axial stencil crosses a
  refinement-level interface are returned as NaN.
  """
  field = np.asarray(field, dtype=float)
  scale_map = np.asarray(scale_map)
  if field.ndim != 2 or scale_map.shape != field.shape:
    raise ValueError('field and scale_map must be 2D arrays with equal shape')

  d2y = np.full(field.shape, np.nan, dtype=float)
  d2x = np.full(field.shape, np.nan, dtype=float)

  scales = np.unique(scale_map)
  scales = scales[scales > 0]

  for scale_value in scales:
    step = int(scale_value)
    if step != scale_value or step < 1:
      raise ValueError('scale_map values must be positive integer cell widths')

    d2y_s = _fd_axis_2d(field, 0, derivative=2, accuracy=accuracy,
                        step=step, spacing=spacing[0])
    d2x_s = _fd_axis_2d(field, 1, derivative=2, accuracy=accuracy,
                        step=step, spacing=spacing[1])

    mask = (scale_map == step)
    if mask_level_interfaces:
      mask_y = mask & _same_scale_stencil_mask(scale_map, 0, step, accuracy)
      mask_x = mask & _same_scale_stencil_mask(scale_map, 1, step, accuracy)
    else:
      mask_y = mask
      mask_x = mask

    d2y[mask_y] = d2y_s[mask_y]
    d2x[mask_x] = d2x_s[mask_x]

  return d2x + d2y
