!===============================================================================
!> @file admesh.f90
!> @brief Adaptive mesh subroutines
!> @author Juan C. Toledo
!> @date 2/Dic/2011

! Copyright (c) 2014 Juan C. Toledo and Alejandro Esquivel
!
! This file is part of Walicxe3D.
!
! Walicxe3D is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.

!===============================================================================

!> @brief Adaptive mesh subroutines
!> @details The module contains the subroutines and utilities that implement the
!> adaptive mesh grid

module amr

  use utils,      only : find, put, pop

  implicit none

  !! The following is a solution to refresh the ghost cells and faces avoiding a
  !! circular dependency (need to call boundaries that is later in the Makefile)

  !> Interface of the ghost-refresh hook
  abstract interface
    subroutine ghost_refresher ()
    end subroutine ghost_refresher
  end interface

  !> Refreshes one layer of ghost cells/faces (and primitives) of all local
  !! blocks. Set by the main program (ghosts module): amr cannot use
  !! boundaries or flux_ct directly, since they depend on amr. If not set,
  !! refineZone works without refreshing (as before).
  procedure(ghost_refresher), pointer :: refresh_ghosts => null()

contains

!===============================================================================

!> @brief High-level routine that triggers mesh refinements
!> @details This routine is called in the main program, and is in charge of
!! deciding and triggering mesh refinement and coarsening, as needed based
!! on the physical data contained by the grid.
subroutine admesh ()

  use parameters
  use globals
  use tictoc
  implicit none

  integer :: nb, bID, flag
  integer :: refcount, crscount
  integer :: localFlags (nbMaxProc)
  integer :: mark
  logical :: test(nbMaxGlobal)  ! DEBUG
  logical :: debug = .false.
  integer :: nf
  integer, allocatable :: fathers(:)

  if (verbosity > 3) call tic(mark)
  if (verbosity > 1) then
    write(logu,*) ""
    write(logu,'(1x,a)') "============================================"
    write(logu,'(1x,a)') " Updating Adaptive Mesh ..."
    write(logu,'(1x,a)') "============================================"
    write(logu,*) ""
  end if

  ! Each process inspects its own active blocks and flags those that meet
  ! the physical criteria for refinement or coarsening
  if (verbosity > 2) write(logu,'(1x,a)')                                      &
                               "> Flagging blocks due to physical gradients ..."
  refcount = 0
  crscount = 0
  localFlags(:) = FLAG_NONE
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then
      call markByPhysical (nb, bID, flag)
      localFlags(nb) = flag
      if (flag.eq.FLAG_REFINE) then
        refcount = refcount + 1
      else if (flag.eq.FLAG_COARSE) then
        crscount = crscount + 1
      end if
    end if
  end do
  if (verbosity > 2)  then
    write(logu,'(1x,i8,a)') refcount, " local blocks marked for refinement"
    write(logu,'(1x,i8,a)') crscount, " local blocks marked for coarsening"
  end if

  ! Syncronize refinement flags
  call mpi_allgather (localFlags, nbMaxProc, mpi_integer, &
    refineFlags, nbMaxProc, mpi_integer, mpi_comm_world, ierr)

  if (debug) then
    test = refineFlags.eq.FLAG_REFINE
    write(logu,'(1x,i8,a)') count(test), " global blocks marked for refinement"
    test = refineFlags.eq.FLAG_COARSE
    write(logu,'(1x,i8,a)') count(test), " global blocks marked for coarsening"
  end if

  ! Sweep all levels except the finest looking to:
  ! a) flag aditional blocks for refinement, and/or
  ! b) inhibit block coarsening
  ! when proximity to a block marked for refinement requires such change
  if (verbosity > 2) write(logu,'(1x,a)') "> Checking proximity criterion ..."
  call checkProximity (maxlev-1)

  ! Blocks about to be refined (needed by the face-field fix-up)
  nf = count(refineFlags == FLAG_REFINE)
  allocate (fathers(nf))
  fathers = pack(globalBlocks, refineFlags == FLAG_REFINE)

  ! Refine local blocks marked for refinement
  if (verbosity > 2) write(logu,'(1x,a)') "> Refining local blocks ..."
  do nb=nbmin,nbmax
    bID = globalBlocks(nb)
    flag = refineFlags(nb)
    if (flag.eq.FLAG_REFINE) then
        if (verbosity > 1) write(logu,'(1x,a,i8)') "Refining block", bID
      call refineBlock(bID)
    end if
  end do

  ! Update global block list
  call syncBlockLists ()

#ifdef BFIELD
  if (enable_flux_ct) call fix_new_faces (fathers, nf)
#endif
  deallocate (fathers)

  ! Do block coarsening
  if (verbosity > 2) write(logu,'(1x,a)') "> Coarsening block families ..."
  call coarseBlocks ()

  ! Update global block list
  call syncBlockLists ()

  if (verbosity > 3) then
    write(logu,*) ""
    write(logu,'(1x,a,a)') "> Mesh refined/coarsened in ", nicetoc(mark)
  end if

end subroutine admesh

!===============================================================================

!> @brief Marks a block for refinement or coarsening if it meets physical criteria
!> @details Note: this requires updated primitives
!> @param locInd Local index of the block to be checked
!> @param bID bID of the block
!> @param flag Refinement flag
subroutine markByPhysical (locIndx, bID, flag)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: locIndx
  integer, intent(in) :: bID
  integer, intent(out) :: flag
  integer :: i, j, k, ilev
  real :: grad, maxgrad, gradx, grady, gradz
  logical, parameter  :: debug = .false.


  call meshlevel (bID, ilev)
  grad = 0.0
  maxgrad = 0.0

  do i=1,ncells_x
    do j=1,ncells_y
      do k=1,ncells_z

        ! Pressure gradient check

        gradx = abs(PRIM(locIndx,5,i+1,j,k)-PRIM(locIndx,5,i-1,j,k)) / &
                PRIM(locIndx,5,i,j,k) / 2.0  !(2*dx(ilev))  use relative change
        grady = abs(PRIM(locIndx,5,i,j+1,k)-PRIM(locIndx,5,i,j-1,k)) / &
                PRIM(locIndx,5,i,j,k) / 2.0  !(2*dx(ilev))
        gradz = abs(PRIM(locIndx,5,i,j,k+1)-PRIM(locIndx,5,i,j,k-1)) / &
                PRIM(locIndx,5,i,j,k) / 2.0  !(2*dx(ilev))
        grad = max( gradx, grady, gradz )
        !   DEBUG ##############################################################
        !grad = 0.
        !   DEBUG ##############################################################

        maxgrad = max( maxgrad, grad )

        ! If the gradient is larger than the refinement threshold, mark the block
        ! unless it is already at the maximum level; exit subroutine regardless
        if (grad.ge.refineThres) then
          if (ilev.eq.maxlev) then
            if (debug) then
              write(logu,'(a,i8,a)') "Block ", bID, " can't be refined any further - you might wanna increase the number of levels"
            end if
            flag = FLAG_NONE
            return
          else if (ilev.lt.maxlev) then
            flag = FLAG_REFINE
            return
          end if
        end if

        ! Density gradient check

        gradx = abs(PRIM(locIndx,1,i+1,j,k)-PRIM(locIndx,1,i-1,j,k)) / &
                PRIM(locIndx,1,i,j,k) / 2.0  !(2*dx(ilev))
        grady = abs(PRIM(locIndx,1,i,j+1,k)-PRIM(locIndx,1,i,j-1,k)) / &
                PRIM(locIndx,1,i,j,k) / 2.0  !(2*dx(ilev))
        gradz = abs(PRIM(locIndx,1,i,j,k+1)-PRIM(locIndx,1,i,j,k-1)) / &
                PRIM(locIndx,1,i,j,k) / 2.0  !(2*dx(ilev))
        grad = max( gradx, grady, gradz )
        maxgrad = max( maxgrad, grad )

        ! If the gradient is larger than the refinement threshold, mark the block
        ! unless it is already at the maximum level; exit subroutine regardless
        if (grad.ge.refineThres) then
          if (ilev.eq.maxlev) then
            if (debug) then
              write(logu,'(a,i8,a)') "Block ", bID, " can't be refined any further - you might wanna increase the number of levels"
            end if
            flag = FLAG_NONE
            return
          else if (ilev.lt.maxlev) then
            flag = FLAG_REFINE
            return
          end if
        end if

        ! Passive scalar check (absolute jumps)
        if (refine_passive .and. (npassive >= 1)) then
          gradx = abs(PRIM(locIndx,firstpas,i+1,j,k)-PRIM(locIndx,firstpas,i-1,j,k)) / 2.0
          grady = abs(PRIM(locIndx,firstpas,i,j+1,k)-PRIM(locIndx,firstpas,i,j-1,k)) / 2.0
          gradz = abs(PRIM(locIndx,firstpas,i,j,k+1)-PRIM(locIndx,firstpas,i,j,k-1)) / 2.0
          grad = max( gradx, grady, gradz )
          maxgrad = max( maxgrad, grad )
          if (grad.ge.refineThres) then
            if (ilev.lt.maxlev) then
              flag = FLAG_REFINE
            else
              flag = FLAG_NONE
            end if
            return
          end if
        end if


      end do
    end do
  end do

  ! If the maximum gradient is smaller than the coarsening threshold, mark
  ! it for coarsening, unless it is a root block
  if ((maxgrad.le.coarseThres).and.(ilev.gt.1)) then
    flag = FLAG_COARSE
  else
    flag = FLAG_NONE
  end if

  return

end subroutine markByPhysical

!===============================================================================

! NO LONGER NEEDED
!
!> @brief Marks a block's *neighbors* for refinement by proximity
!> @details The neighbors are checked.
!subroutine markByProximity (bID)

!  use parameters
!  use globals
!  implicit none

!  integer, intent(in) :: bID

!  integer :: dir, nID, fID, fowner, nbf
!  integer :: neighType, neighList(4)

!!  write(logu,*) "PROXIMITY CHECK for bID", bID

!  do dir=1,6
!    call neighbors(bID, dir, neighType, neighList)
!    ! DEBUG
!!    select case (dir)
!!    case (1)
!!      write(logu,*) "Direction:   LEFT"
!!    case (2)
!!      write(logu,*) "Direction:   RIGHT"
!!    case (3)
!!      write(logu,*) "Direction:   FRONT"
!!    case (4)
!!      write(logu,*) "Direction:   BACK"
!!    case (5)
!!      write(logu,*) "Direction:   BOTTOM"
!!    case (6)
!!      write(logu,*) "Direction:   TOP"
!!    end select
!!    write(logu,*) "Neighbor Type:", neighType
!!    write(logu,*) "Neighbor ID(s):", neighList
!    ! DEBUG
!  end do

!end subroutine markByProximity

!===============================================================================

!> @brief Refines a block into 8 child blocks
!> @details This requires the following steps:
!! 1) Obtaining the bIDs of all the child blocks
!! 2) Registering the children in the local block registry
!! 3) Interpolating (limited linear) the father's data into the children
!! 4) Unregistering the father from the local block registry
!> @param bID The father block's (absolute) ID
subroutine refineBlock (fatherID)

  use parameters
  use globals
  use clean_quit, only : clean_abort
  use hydro_core, only : flow2prim, limited_slope, calcPrimsBlock

  implicit none

  integer, intent(in) :: fatherID

  integer :: ilev, nb, fatherIndex, childIndex, childID
  integer :: c, i, j, k, ip, jp, kp, sx, sy, sz
  integer :: ic, jc, kc, di, dj, dk, istat
  integer :: childList(8)
#ifdef BFIELD
  integer :: badcells
#endif
  real    :: c0(neqtot), sl(neqtot,3), fv(neqtot)
  logical :: ok


  ! TODO: THIS CHECK MIGHT BE ELIMINATED FOR EFFICIENCY
  call find (fatherID, localBlocks, nbMaxProc, fatherIndex)
  if (fatherIndex.eq.-1) then
      write(logu,*) ""
      write(logu,'(a,i5,a)') "Block ", fatherID, " is not a local block; can't refine it!"
      write(logu,'(1x,a)') "***ABORTING***"
    call clean_abort (ERROR_LOCAL_BID_NOT_FOUND)
  end if

  ! Don't refine past the highest mesh level
  ! TODO: THIS CHECK MIGHT BE REDUNDANT
  call meshlevel(fatherID, ilev)
  if (ilev.eq.maxlev) then
    write(logu,*) ""
    write(logu,'(a,i5,a)') "Trying to refine block ", fatherID, " past max mesh level! Aborting!"
    write(logu,'(a)') "***ABORTING***"
    call clean_abort (ERROR_ALREADY_MAX_LEV)
  end if

  ! Check if there are enough local block slots
  if (nbLocal+8>nbMaxProc) then
    write(logu,*) ""
    write(logu,'(a,i4)') "Not enough memory space for child blocks of bID ", fatherID
    write(logu,'(a)') "You might want to increase nbMaxProc."
    write(logu,'(a)') "***ABORTING***"
    call clean_abort (ERROR_INSUFICIENT_NBMAXPROC)
  end if

  ! If previous checks passed, go ahead

  ! Get list of children bIDs
  call children (fatherID, childList)

  ! Create children: assign children's bIDs to empty slots in localBlocks
  ! and copy corresponding father's data in children's data space
  do c=1,8

    childID = childList(c)

    ! Assign child's bID to empty local block slot
    call put (childID, localBlocks, nbMaxProc, childIndex)
    if (childIndex.eq.-1) then
      write(logu,*) ""
      write(logu,'(a)') "Child block couldn't be assigned to an empty slot!"
      write(logu,'(a)') "***ABORTING***"
      call clean_abort (ERROR_REGISTER_CHILD)
    !    else
    !      write(logu,'(a,i8,a,i5)') " Child block with bID ", childID, &
    !                                " assigned to local slot ", childIndex
    end if

    ! Determine child's sibling coordinates (xs, ys, zs)
    call siblingCoords (childID, sx, sy, sz)

    ! Fill the child by limited linear interpolation of the father's data.
    ! Each father cell (ip,jp,kp) gives the 2x2x2 child cells (i,j,k), whose
    ! average equals the father value. Slopes use father neighbors. including
    ! its refreshed ghost cells (in main program).  If any of the 8 children
    ! fails flow2prim, all 8 get the father value (injection), which keeps the
    ! refinement conservative.
    do kc=1,ncells_z/2
      do jc=1,ncells_y/2
        do ic=1,ncells_x/2

          ip = (ncells_x/2)*sx + ic
          jp = (ncells_y/2)*sy + jc
          kp = (ncells_z/2)*sz + kc
          c0(:) = U(fatherIndex,:,ip,jp,kp)

          sl(:,:) = 0.0
          call limited_slope ( c0 - U(fatherIndex,:,ip-1,jp,kp),               &
                                    U(fatherIndex,:,ip+1,jp,kp) - c0, sl(:,1))
          call limited_slope ( c0 - U(fatherIndex,:,ip,jp-1,kp),               &
                                    U(fatherIndex,:,ip,jp+1,kp) - c0, sl(:,2))
          call limited_slope ( c0 - U(fatherIndex,:,ip,jp,kp-1),               &
                                    U(fatherIndex,:,ip,jp,kp+1) - c0, sl(:,3))

          ok = .true.
          do dk=0,1
            do dj=0,1
              do di=0,1
                i = 2*ic - 1 + di
                j = 2*jc - 1 + dj
                k = 2*kc - 1 + dk
                fv(:) = c0(:) + 0.5*(di-0.5)*sl(:,1) + 0.5*(dj-0.5)*sl(:,2)   &
                              + 0.5*(dk-0.5)*sl(:,3)
                U(childIndex,:,i,j,k) = fv(:)
                call flow2prim (fv, PRIM(childIndex,:,i,j,k), istat)
                if (istat /= 0) ok = .false.
              end do
            end do
          end do

          ! Positivity safeguard: inject the father value into all 8 children
          if (.not.ok) then
            do dk=0,1
              do dj=0,1
                do di=0,1
                  i = 2*ic - 1 + di
                  j = 2*jc - 1 + dj
                  k = 2*kc - 1 + dk
                  U(childIndex,:,i,j,k)    = c0(:)
                  PRIM(childIndex,:,i,j,k) = PRIM(fatherIndex,:,ip,jp,kp)
                end do
              end do
            end do
          end if
#ifdef BFIELD
          ! Divergence-free prolongation of the face field (Toth & Roe 2002)
          if (enable_flux_ct) call prolong_cell_faces (fatherIndex, childIndex,&
                                       ip, jp, kp, 2*ic-1, 2*jc-1, 2*kc-1)
#endif
        end do
      end do
    end do

#ifdef BFIELD
    ! Cell-centred B from the prolonged faces, then primitives
    if (enable_flux_ct) then
      call faces_to_cells (childIndex)
      call calcPrimsBlock (U, PRIM, childIndex, CELLS_PHYS, badcells)
    end if
#endif

  end do

  ! Free father's slot
  !  do nb=1,nbMaxProc
  !    if (localBlocks(nb).eq.fatherID) then
  !      localBlocks(nb) = -1
  !!      write(logu,'(a,i8)') " Freeing slot of father block, bID ", fatherID
  !      exit
  !    end if
  !  end do
  call pop (fatherID, localBlocks, nbMaxProc, nb)

  ! Update local block count
  nbLocal = nbLocal + 7

end subroutine refineBlock

!===============================================================================

!> @brief Synchronizes the global block registry
!> @details Simply merges the lists of local blocks
subroutine syncBlockLists()

  use parameters
  use globals
  implicit none

  integer :: nb

#ifndef MPIP

  ! TODO

#else

  if (verbosity > 2) then
    write(logu,*) ""
    write(logu,'(1x,a)') "Synchronizing block lists ..."
  end if
  ! Reset global block registry
  globalBlocks(:) = -1

  ! Synchronize globalBlocks
  call mpi_allgather (localBlocks, nbMaxProc, mpi_integer, &
      globalBlocks, nbMaxProc, mpi_integer, mpi_comm_world, ierr)

  ! Update nbActive count
  nbActive = 0
  do nb=1,nbMaxGlobal
    if (globalBlocks(nb).ne.-1) then
      nbActive = nbActive + 1
    end if
  end do

  if (verbosity > 2) write(logu,'(1x,a,i0,a)') "There are ", nbActive,         &
                          " active blocks globally."

  ! Everybody stops here
  call mpi_barrier (mpi_comm_world, ierr)

#endif
end subroutine syncBlockLists

!===============================================================================

!> @brief Returns the neighbors of a block at the same level of refinement
!> @details Returns the bID, regardless of whether it is active or not
!> @param bID The (absolute) block ID of which one needs a neighbor
!> @param dir The direction of inquiry, one of the named constants
!! LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK
subroutine neighborLevel(bID, dir, neighID)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(in) :: dir
  integer, intent(out) :: neighID

  !integer :: offsets
  integer :: x, y, z, ilev, nx, ny, nz
  integer :: xp, yp, zp

  call bcoords(bID, x, y, z)
  call meshlevel(bID, ilev)
  nx = nbx(ilev)
  ny = nby(ilev)
  nz = nbz(ilev)

  neighID = 0
  xp = x
  yp = y
  zp = z

  select case (dir)

  case (LEFT)
    if (x.eq.1) then
      if (bc_left.eq.BC_PERIODIC) then
        xp = nx
      else
        neighID = -1
        return
      end if
    else
      xp = x-1
    end if

  case (RIGHT)
    if (x.eq.nx) then
      if (bc_right.eq.BC_PERIODIC) then
        xp = 1
      else
        neighID = -1
        return
      end if
    else
      xp = x+1
    end if

  case (FRONT)
    if (y.eq.1) then
      if (bc_front.eq.BC_PERIODIC) then
        yp = ny
      else
        neighID = -1
        return
      end if
    else
      yp = y-1
    end if

  case (BACK)
    if (y.eq.ny) then
      if (bc_back.eq.BC_PERIODIC) then
        yp = 1
      else
        neighID = -1
        return
      end if
    else
      yp = y+1
    end if

  case (BOTTOM)
    if (z.eq.1) then
      if (bc_bottom.eq.BC_PERIODIC) then
        zp = nz
      else
        neighID = -1
        return
      end if
    else
      zp = z-1
    end if

  case (TOP)
    if (z.eq.nz) then
      if (bc_top.eq.BC_PERIODIC) then
        zp = 1
      else
        neighID = -1
        return
      end if
    else
      zp = z+1
    end if

  case DEFAULT
    print*, "INVALID neighbor direction:", dir
    return

  end select

  ! Calculate actual neighbor bID
  neighID = 1 + (xp-1) + (yp-1)*nx + (zp-1)*nx*ny + offsets(ilev)

end subroutine neighborLevel

!===============================================================================

!> @brief Returns the bID of actual active neighbors of a given block
!> @details Returns an integer neighType which describes the kind of
!! neighbors:
!! @n NEIGH_SAME     = neighbors is as the same level of refinement
!! @n NEIGH_COARSER  = neighbors is coarser
!! @n NEIGH_FINER    = higher-level (finer) neighbors (four)
!! @n NEIGH_BOUNDARY = simulation box boundary
!> @param bID The (absolute) block ID of which one needs a neighbor
!> @param dir The direction of inquiry, one of the named constants
!! LEFT, RIGHT, TOP, BOTTOM, FRONT, BACK
!> @param neighType The type of neighbor(s)
!! @param neighList List of neighbor(s) bID(s)
subroutine neighbors (bID, dir, neighType, neighList)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(in) :: dir
  integer, intent(out) :: neighType
  integer, intent(out) :: neighList(4)

  integer :: nID, fID, nb, childList(8)

  neighType = -1
  neighList(:) = -1

  ! Invalid bID
  if (bID.eq.-1) then
    neighType = -1
    return
  end if

  ! Get the direct neighbor's bID and its father's bID
  call neighborLevel(bID, dir, nID)
  call father(nID, fID)

  ! Neighbor is grid boundary
  if (nID.eq.-1) then
    neighType = NEIGH_BOUNDARY
    return
  end if

  ! Check which is active, if any
  do nb=1,nbMaxGlobal
    if (globalBlocks(nb).eq.nID) then
      neighType = NEIGH_SAME
      neighList(1) = nID
      return
    else if ((fID.ne.-1).and.(globalBlocks(nb).eq.fID)) then
      neighType = NEIGH_COARSER
      neighList(1) = fID
      return
    end if
  end do

  ! If neither is found, neighbors are four higher-level blocks, the children
  ! of the same-level neighbor
  neighType = NEIGH_FINER
  call children(nID, childList)

  ! Determine which four of the eight children are neighbors
  select case(dir)

    case (LEFT)
      neighList(1) = childList(2)
      neighList(2) = childList(4)
      neighList(3) = childList(6)
      neighList(4) = childList(8)

    case (RIGHT)
      neighList(1) = childList(1)
      neighList(2) = childList(3)
      neighList(3) = childList(5)
      neighList(4) = childList(7)

    case (FRONT)
      neighList(1) = childList(3)
      neighList(2) = childList(4)
      neighList(3) = childList(7)
      neighList(4) = childList(8)

    case (BACK)
      neighList(1) = childList(1)
      neighList(2) = childList(2)
      neighList(3) = childList(5)
      neighList(4) = childList(6)

    case (BOTTOM)
      neighList(1) = childList(5)
      neighList(2) = childList(6)
      neighList(3) = childList(7)
      neighList(4) = childList(8)

    case (TOP)
      neighList(1) = childList(1)
      neighList(2) = childList(2)
      neighList(3) = childList(3)
      neighList(4) = childList(4)

  end select
  return

end subroutine neighbors

!===============================================================================

!> @brief Returns the refinement level of a block
!> @param bID The (absolute) block ID of which one wants the mesh level
subroutine meshlevel(bID, level)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: level

  if (bID.eq.-1) then
    level = -1
    return
  end if

  level = -1
  do level=1,maxlev
    if ((bID.ge.minID(level)).and.(bID.le.maxID(level))) then
      return
    end if
  end do

end subroutine meshlevel

!===============================================================================
!> @brief Returns the number of blocks in all previous levels
!> @param level The mesh level for which one wants the block offset
integer function offsets(level)

  use globals, only : nbrootx, nbrooty, nbrootz

  implicit none
  integer, intent(in) :: level
  integer :: ilev

  offsets = 0
  do ilev=2,level
    offsets = offsets + nbrootx*nbrooty*nbrootz*8**(ilev-2)
  end do
  return

end function offsets

!===============================================================================

!> @brief Returns the bID of a block's firstborn
!> @param bID The (absolute) block ID of the father
!> @param fbID The (absolute) block ID of the block's firstborn child
subroutine firstborn(bID, firstbornID)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: firstbornID

  !integer :: offsets
  integer :: x, y, z, ilev, nx, ny, nz
  integer :: xp, yp, zp

  call bcoords(bID, x, y, z)
  call meshlevel(bID, ilev)
  ! DEBUG
  !write(logu,*) "Calling firstborn..."
  !write(logu,*) ilev
  !if (ilev.gt.maxlev) then
  !  write(logu,*) "Block", bID, "has firstborn beyond the max mesh level!!!"
  !  write(logu,*) "Y U ASKING FOR IT?!"
  !end if
  ! DEBUG
  nx = nbx(ilev+1)
  ny = nby(ilev+1)
  nz = nbz(ilev+1)

  xp = 2*x-1
  yp = 2*y-1
  zp = 2*z-1

  firstbornID = 1 + (xp-1) + (yp-1)*nx + (zp-1)*nx*ny + offsets(ilev+1)

end subroutine firstborn

!===============================================================================

!> @brief Returns a list of bIDs with all 8 siblings
!> @param bID The block ID of one of the siblings
!> @param sibling_list The list of the 8 (absolute) block IDs of the block's
!! siblings, including itself.
!> @details The siblings are returned in "natural order": first x, then y,
!! last z, producing the following sibling ordering:
!!    ---------
!!    | 7 | 8 |
!!    ---------   z = top
!!    | 5 | 6 |
!!    ---------
!!    ---------
!!    | 3 | 4 |
!!    ---------   z = bottom
!!    | 1 | 2 |
!! ^  ---------
!! |
!! y  x -->
subroutine siblings(bID, sibling_list)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: sibling_list(8)

  integer :: fatherID, firstbornID, ilev, nx, ny, nz, i

  call father(bID, fatherID)

  if (fatherID.eq.-1) then
    do i=1,8
      sibling_list(i) = -1
    end do
    return
  end if
  call firstborn(fatherID, firstbornID)
  call meshlevel(bID, ilev)

  nx = nbx(ilev)
  ny = nby(ilev)
  nz = nbz(ilev)

  sibling_list(1) = firstbornID
  sibling_list(2) = firstbornID + 1
  sibling_list(3) = firstbornID + nx
  sibling_list(4) = firstbornID + nx + 1

  sibling_list(5) = firstbornID + nx*ny
  sibling_list(6) = firstbornID + nx*ny + 1
  sibling_list(7) = firstbornID + nx*ny + nx
  sibling_list(8) = firstbornID + nx*ny + nx + 1

end subroutine siblings

!===============================================================================

! >> OBSOLETE
!> @brief Calculates the sibling ID of a block
!> @param bID The (absolute) block ID
!> @param sID The sibling ID of the block, following the "natural" order:
!!    ---------
!!    | 7 | 8 |
!!    ---------   z = top
!!    | 5 | 6 |
!!    ---------
!!    ---------
!!    | 3 | 4 |
!!    ---------   z = bottom
!!    | 1 | 2 |
!! ^  ---------
!! |
!! y  x -->
subroutine siblingID(bID, sID)

  use parameters
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: sID

  integer :: sibling_list(8), i, fatherID

  call father(bID, fatherID)
  !  print*, "Father of", bID, "is", fatherID  ! DEBUG
  if (fatherID.eq.-1) then
    ! return -1 if block is root block (no father or siblings)
    sID = -1
    return
  end if

  call siblings(bID, sibling_list)

  sID = -1
  do i=1,8
    if (bID.eq.sibling_list(i)) then
      sID = i
      exit
    end if
  end do

end subroutine siblingID
!===============================================================================

!> @brief Calculates the sibling position (coords) of a block
!> @param bID The (absolute) block ID
!> @param sx The sibling's x-coordinate (0 or 1)
!> @param sy The sibling's y-coordinate (0 or 1)
!> @param sz The sibling's z-coordinate (0 or 1)
subroutine siblingCoords (bID, sx, sy, sz)

  use parameters
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: sx, sy, sz

  integer :: sibling_list(8), i, fatherID, sID

  call father(bID, fatherID)
  if (fatherID.eq.-1) then
    ! return -1 if block is root block (no father or siblings)
    sx = -1;  sy = -1;  sz = -1;
    return
  end if

  call siblings(bID, sibling_list)

  sID = -1
  do i=1,8
    if (bID.eq.sibling_list(i)) then
      sID = i
      exit
    end if
  end do

  sx = 0
  sy = 0
  sz = 0
  if (mod(sID,2).eq.0) then
    sx = 1
  end if
  if ((sID.eq.3).or.(sID.eq.4).or.(sID.eq.7).or.(sID.eq.8)) then
    sy = 1
  end if
  if (sID.gt.4) then
    sz = 1
  end if

  return

end subroutine siblingCoords

!===============================================================================

!> @brief Returns a list of bIDs with all 8 children of a block
!> @param bID The (absolute) block ID of one of the siblings
!> @param sibling_list The list of the 8 (absolute) block IDs of the block's
!! siblings, including itself. These are given in "natural order": first x,
!! then y, then z
subroutine children(bID, children_list)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: children_list(8)

  integer :: firstbornID, ilev, nx, ny, nz

  call firstborn(bID, firstbornID)
  call meshlevel(bID, ilev)
  nx = nbx(ilev+1)
  ny = nby(ilev+1)
  nz = nbz(ilev+1)

  children_list(1) = firstbornID
  children_list(2) = firstbornID + 1
  children_list(3) = firstbornID + nx
  children_list(4) = firstbornID + nx + 1

  children_list(5) = firstbornID + nx*ny
  children_list(6) = firstbornID + nx*ny + 1
  children_list(7) = firstbornID + nx*ny + nx
  children_list(8) = firstbornID + nx*ny + nx + 1

end subroutine children

!===============================================================================

!> @brief Returns a block's father's bID
!> @param bID The block's (absolute) bID
!> @param fatherbID The block's father's (absolute) bID
subroutine father(bID, fatherID)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: fatherID

  !integer :: offsets
  integer :: x, y, z, ilev, nx, ny, nz
  integer :: xp, yp, zp

  if (bID.eq.-1) then
    fatherID = -1
    return
  end if

  call meshlevel(bID, ilev)

  if (ilev.eq.1) then

    fatherID = -1
    return

  else

    call bcoords(bID, x, y, z)

    xp = floor((x+1)/2.0)
    yp = floor((y+1)/2.0)
    zp = floor((z+1)/2.0)

    nx = nbx(ilev-1)
    ny = nby(ilev-1)
    nz = nbz(ilev-1)

    fatherID = 1 + (xp-1) + (yp-1)*nx + (zp-1)*nx*ny + offsets(ilev-1)

  !print*, "Block", bID, ", x=", x, ", y=", y, ", z=", z  ! DEBUG
  !print*, bID, "'s father coords:", xp, yp, zp  ! DEBUG
  !print*, ''  ! DEBUG

  end if

end subroutine father


!===============================================================================

!> @brief Returns the (x,y,z) local integer coordinates of a block at
!! the block's mesh level
!> @param bID The (absolute) block ID of the block
!> @param x The level-specific integer x-coordinate
!> @param y The level-specific integer y-coordinate
!> @param z The level-specific integer z-coordinate
subroutine bcoords(bID, x, y, z)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: x, y, z

  !integer :: offsets
  integer :: ilev, nx, ny, nz, localID

  call meshlevel(bID, ilev)
  nx = nbx(ilev)
  ny = nby(ilev)
  nz = nbz(ilev)

  localID = bID - offsets(ilev)
  x = mod(localID,nx)
  if (x.eq.0) x=nx
  y = mod(ceiling(localID*1.0/(nx)),ny)
  if (y.eq.0) y=ny
  z = mod(ceiling(localID*1.0/(nx*ny)),nz)
  if (z.eq.0) z=nz

end subroutine bcoords

!===============================================================================

!> @brief Returns the physical coordinates (in code units) of a block's
!! reference corner
!> @param bID The (absolute) block ID of the block
!> @param x The level-specific integer x-coordinate
!> @param y The level-specific integer y-coordinate
!> @param z The level-specific integer z-coordinate
subroutine getRefCorner(bID, xx, yy, zz)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  real, intent(out) :: xx, yy, zz

  integer :: x, y, z, ilev

  call meshlevel(bID, ilev)
  call bcoords(bID, x, y, z)

  xx = (x-1)*ncells_x*dx(ilev)
  yy = (y-1)*ncells_y*dy(ilev)
  zz = (z-1)*ncells_z*dz(ilev)

end subroutine getRefCorner

!===============================================================================

!> @brief Returns the physical position (in code units) of a single cell
!> @param bID The block ID of the host block
!> @param i The cell's integer x-position within the block (1 to ncells_x)
!> @param j The cell's integer y-position within the block (1 to ncells_y)
!> @param k The cell's integer z-position within the block (1 to ncells_z)
!> @param xx The x-position of the cell's physical center, in code units
!> @param yy The y-position of the cell's physical center, in code units
!> @param zz The z-position of the cell's physical center, in code units
!> @param center, if true returns the position respect to the center of
!> the grid, by default returns respect to one corner
subroutine cellPos (bID, i, j, k, x, y, z, center)

  use parameters
  use globals
  implicit none

  integer,           intent(in)  :: bID
  integer,           intent(in)  :: i, j, k
  real,              intent(out) :: x, y, z
  logical, optional, intent(in)  :: center

  integer :: ilev, xb, yb, zb

  if (bID.eq.-1) then
    x = -1.0
    y = -1.0
    z = -1.0
    return
  end if

  call meshlevel(bID, ilev)
  call bcoords(bID, xb, yb, zb)

  x = ((xb-1)*ncells_x + (i-1) + 0.5) * dx(ilev)
  y = ((yb-1)*ncells_y + (j-1) + 0.5) * dy(ilev)
  z = ((zb-1)*ncells_z + (k-1) + 0.5) * dz(ilev)

  if (present(center) .and. center) then
    x = x - xphystot/2.0/l_sc
    y = y - yphystot/2.0/l_sc
    z = z - zphystot/2.0/l_sc
  end if

  return

end subroutine cellPos

!===============================================================================

!> @brief Returns the rank of the process who owns this bID
!> @details Returns MPI_PROC_NULL if the bID is not in globalBlocks
!> @param bID The (absolute) block ID of the block
!> @param process The block owner's rank
subroutine getOwner (bID, owner)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: owner

  integer :: gIndx

  if (bID.eq.-1) then
    owner = MPI_PROC_NULL
    return
  end if

  call find (bID, globalBlocks, nbMaxGlobal, gIndx)

  if (gIndx.ne.-1) then
    owner = (gIndx-1)/nbMaxProc
  else
    owner = MPI_PROC_NULL
  end if
  return

end subroutine getOwner

!===============================================================================

!> @brief Ensure the proximity criterion is on all levels less than or equal
!! to finestLev
!> @details This routine will flag blocks for refinement or un-flag blocks
!! for coarsening if they would violate the level proximity condition. Only
!! level finestLev and any coarser levels are checked. Assumesthe refinement
!! flags were previously calculated and synchronized.
!> @param startLev Initial (finer) level to check
subroutine checkProximity (finestLev)

  use globals
  use parameters
  implicit none

  integer, intent(in) :: finestLev

  integer :: level, nb, nb1, nb2
  integer :: bID, nID, dir, i, bflag, nflag
  integer :: neighType, neighList(4), owner
  integer :: sib_list(8)  ! DEBUG
  logical :: inhibited

  ! DEBUG
  !  write(logu,*) ""
  !  write(logu,*) "Starting Proximity checks ..."
  ! DEBUG

  do level=finestLev,1,-1
    ! For every level from startLev up to root ...

    ! DEBUG
    !    write(logu,*) "Checking blocks in level", level
    !    write(logu,*) "GLOBAL refinement list:"
    !    do nb=1,nbMaxProc
    !      if (gToRefine(nb).ne.-1) then
    !        write(logu,*) gToRefine(nb)
    !      end if
    !    end do
    ! DEBUG

    do nb=1,nbMaxGlobal
      bID = globalBlocks(nb)
      if ((bID.ne.-1).and.(bID.ge.minID(level)).and.(bID.le.maxID(level))) then

        ! This block's current refinement flag
        bflag = refineFlags(nb)

        ! Step 1: apply proximity refinement
        !
        ! For every block in this level currently marked for refinement, check if:
        ! a) an unmarked neighbor must be marked for refinement due to proximity
        ! b) a neighbor marked for coarsening must be inhibited and/or marked
        ! for refinement due to proximity
        ! Note: refinement is never inhibited
        !
        if (bflag.eq.FLAG_REFINE) then
          do dir=1,6
            ! Check all neighbors

            call neighbors (bID, dir, neighType, neighList)
            call getOwner (bID, owner)

            ! Inhibit coarsening of same-level neighbors
            if (neighType.eq.NEIGH_SAME) then

              nID = neighList(1)
              call find (nID, globalBlocks, nbMaxGlobal, nb1)
              nflag = refineFlags(nb1)

              if (nflag.eq.FLAG_COARSE) then
                refineFlags(nb1) = FLAG_NONE
              end if

            ! Force refinement of coarser neighbors
            else if (neighType.eq.NEIGH_COARSER) then

              nID = neighList(1)
              call find (nID, globalBlocks, nbMaxGlobal, nb2)
              nflag = refineFlags(nb2)

              if (nflag.ne.FLAG_REFINE) then
                refineFlags(nb2) = FLAG_REFINE
              end if

            end if

          end do
        end if

        ! Step 2: inhibit coarsening due to proximity
        !
        ! Now check all blocks marked for coarsening. Inhibit coarsening if:
        ! a) any finer neighbor is not marked for coarsening too
        !
        if (bflag.eq.FLAG_COARSE) then

          inhibited = .false.
          ! Check all neighbors
          do dir=1,6

            !            call neighbors(bID, dir, neighType, neighList)
            !            if (neighType.eq.NEIGH_FINER) then
            !              do i=1,4
            !                nID = neighList(i)
            !                call find (nID, globalBlocks, nbMaxGlobal, nb1)
            !                if (nb1.eq.-1) then
            !                  refineFlags(nb) = FLAG_NONE
            !                  inhibited = .true.
            !                  exit
            !                end if
            !              end do
            !            end if

            ! Skip remaining checks if already inhibited
            !if (inhibited) exit

            call neighbors(bID, dir, neighType, neighList)
            if (neighType.eq.NEIGH_FINER) then

              call siblings(neighList(1), sib_list)
              do i=1,8

                nID = sib_list(i)
                call find (nID, globalBlocks, nbMaxGlobal, nb1)
                if (nb1.ne.-1) then
                  nflag = refineFlags(nb1)
                  if (nflag.ne.FLAG_COARSE) then
                    refineFlags(nb) = FLAG_NONE
                    inhibited = .true.
                    exit
                  end if
                end if

              end do
            end if

            ! Skip remaining checks if already inhibited
            if (inhibited) exit

          end do
        end if

      end if

    end do
  end do

!  write(logu,*) "Proximity checks done for all levels"   ! DEBUG


end subroutine checkProximity

!===============================================================================

!> @brief Refines a zone of the grid to the requested level
!> @details The rectangular zone to be refined is defined by the
!! provided zone array, which contains the bounding box of the zone:
!!  zone(1): left edge
!!  zone(2): right edge
!!  zone(3): front edge
!!  zone(4): back edge
!!  zone(5): bottom edge
!!  zone(6): top edge
!! These must be given in cgs units (cm). The second parameter, z_level,
!! specifies the level at which the zone is to be refined.
subroutine refineZone (zone, z_level)

  use parameters
  use globals
  use clean_quit, only : clean_abort

  implicit none

  real, intent(in) :: zone(6)
  integer, intent(in) :: z_level

  integer :: ilev, nb, nb1, nb3, bID
  integer :: flag
  integer :: localFlags(nbMaxProc)
  real :: z_left, z_right, z_front, z_back, z_bottom, z_top
  real :: b_left, b_right, b_back, b_front, b_bottom, b_top
  integer :: nf
  integer, allocatable :: fathers(:)

  ! Verbosity flag (for debugging, or if you're feeling lonely)
  logical, parameter  :: debug = .false.

  ! Sanity check
  if (z_level.gt.maxlev) then
    write(logu,*) "Impossible to refine IC zone to level", z_level, "!"
    write(logu,*) "The maximum refinement level is", maxlev
    write(logu,*) "***ABORTING***"
    call clean_abort (ERROR_REFINING_PAST_MAX)
  end if

  ! Unpack zone bounding box
  z_left   = zone(1)
  z_right  = zone(2)
  z_front  = zone(3)
  z_back   = zone(4)
  z_bottom = zone(5)
  z_top    = zone(6)

  if (verbosity > 2) then
    write(logu,*) ""
    write(logu,'(1x,a,i0)') "> Refining zone to level ", z_level
    write(logu,'(1x,a,es12.5,1x,es12.5)') "x-range: ", z_left, z_right
    write(logu,'(1x,a,es12.5,1x,es12.5)') "y-range: ", z_front, z_back
    write(logu,'(1x,a,es12.5,1x,es12.5)') "z-range: ", z_bottom, z_top
  end if
  ! For all refinement levels up to the requestes level
  do ilev=1,z_level-1

    if (debug) then
      write(logu,*) ""
      write(logu,*) "Checking blocks in level", ilev
      write(logu,*) "Local blocks:"
      do nb=1,nbMaxProc
        if (localBlocks(nb).ne.-1) then
          write(logu,*) localBlocks(nb)
        end if
      end do
    end if

    ! One layer of fresh ghost cells and faces: refineBlock interpolates from
    ! them, and the blocks created in the previous pass have none yet
    if (associated(refresh_ghosts)) call refresh_ghosts ()

    ! Clear local refinement flags
    localFlags(:) = -1

    ! Flag for refinement all *local* blocks that contain the zone
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if ((bID.ne.-1).and.(bID.ge.minID(ilev)).and.(bID.le.maxID(ilev))) then

        ! Get reference corner position and convert to cgs
        call getRefCorner(bID, b_left, b_front, b_bottom)
        b_right = (b_left + dx(ilev)*ncells_x) * l_sc
        b_back = (b_front + dy(ilev)*nCells_y) * l_sc
        b_top = (b_bottom + dz(ilev)*nCells_z) * l_sc
        b_left = b_left * l_sc
        b_front = b_front * l_sc
        b_bottom = b_bottom * l_sc

        if (debug) then
          write(logu,*) ""
          write(logu,*) "Local block", bID, "has bounding volume"
          write(logu,*) "x:", b_left, b_right
          write(logu,*) "y:", b_front, b_back
          write(logu,*) "z:", b_bottom, b_top
        end if

        ! If block within bounding box, flag it for refinement
        if ((z_left.lt.b_right).and.(z_right.gt.b_left).and.&
          (z_front.lt.b_back).and.(z_back.gt.b_front).and.&
          (z_bottom.lt.b_top).and.(z_top.gt.b_bottom)) then
          localFlags(nb) = FLAG_REFINE
          if (debug) write(logu,*) "Block is contained - will be refined"
        else
          localFlags(nb) = FLAG_NONE
          if (debug) write(logu,*) "Block is not contained"
        end if

      end if
    end do

    if (debug) then
      write(logu,*) "localFlags after marking"
      write(logu,*) "         bID         flag"
      do nb3=1,nbMaxProc
        if (localFlags(nb3).ne.-1) then
          write(logu,*) localBlocks(nb3), localFlags(nb3)
        end if
      end do
    end if

    ! Syncronize refinement flags
    call mpi_allgather (localFlags, nbMaxProc, mpi_integer, &
    refineFlags, nbMaxProc, mpi_integer, mpi_comm_world, ierr)

    if (debug) then
      write(logu,*) "refineFlags after synchronization"
      write(logu,*) "         bID         flag"
      do nb3=1,nbMaxGlobal
        if (refineFlags(nb3).ne.-1) then
          write(logu,*) globalBlocks(nb3), refineFlags(nb3)
        end if
      end do
    end if

    ! Check proximity criterion for all blocks in this level and coarser levels
    if (verbosity > 2) write(logu,'(1x,a)') "Doing proximity checks ..."
    call checkProximity(ilev)

    ! Blocks about to be refined (needed by the face-field fix-up)
    nf = count(refineFlags == FLAG_REFINE)
    allocate (fathers(nf))
    fathers = pack(globalBlocks, refineFlags == FLAG_REFINE)

    if (debug) then
      write(logu,*) "         bID        owner        flag"
      do nb3=1,nbMaxGlobal
        if (refineFlags(nb3).ne.-1) then
          call getOwner(globalBlocks(nb3), nb1)
          write(logu,*) globalBlocks(nb3), nb1, refineFlags(nb3)
        end if
      end do
    end if

    ! Refine final list of local blocks
    if (verbosity > 2) write(logu,'(1x,a)') "Refining local blocks ..."
    do nb=nbmin,nbmax
      bID = globalBlocks(nb)
      flag = refineFlags(nb)
      if (flag.eq.FLAG_REFINE) then
      !if ((bID.ne.-1).and.(flag.eq.FLAG_REFINE)) then
        call refineBlock(bID)
      end if
    end do

    ! Update block registry
    call syncBlockLists ()

#ifdef BFIELD
    if (enable_flux_ct) call fix_new_faces (fathers, nf)
#endif
  deallocate (fathers)

    ! Do load balance
    call doBalance ()

    ! Done for this level - continue with next
  end do

  if (verbosity > 2) write(logu,'(1x,a)') "Finished refining Zone..."

end subroutine refineZone

!===============================================================================

!> @brief Coarsens families of blocks into a single father block
!> @details Determines which families (sets of siblings) are all set to be
!! coarsened, then does the coarsening.
subroutine coarseBlocks ()

  use globals
  use parameters
  use clean_quit, only : clean_abort
  use hydro_core, only : calcPrimsBlock

  implicit none

  integer :: nb, nb1, bID, eID, sID, fID, fatherIndex, locIndx
  integer :: i, j, k, ip, jp, kp, i1, i2, j1, j2, k1, k2, sx, sy, sz
  integer :: next, flag, sib, sib_owner, elder_owner, ieq, badcells
  integer :: nData, mpistatus(MPI_STATUS_SIZE)
  integer :: family (nbActive,2), sib_list(8)
  real :: buf(neqtot, ncells_x/2, ncells_y/2, ncells_z/2)
  logical, parameter :: debug = .false.
  logical :: counted
  integer, parameter :: m3 = (ncells_x/2+1)*(ncells_x/2)**2
#ifdef BFIELD
  real :: fbuf(3*m3)
  real :: r1(0:ncells_x/2,  ncells_y/2,    ncells_z/2)
  real :: r2( ncells_x/2,  0:ncells_y/2,   ncells_z/2)
  real :: r3( ncells_x/2,  ncells_y/2,   0:ncells_z/2)
#endif

  ! 1. Build the list of families

  ! Sweep refineFlags and add every block marked for coarseningto its family.
  ! The family array has the following meaning:
  !  family(:,1) = bID of the elder block
  !  family(:,2) = number of member of this family marked for coarsening
  family(:,1) = -1
  family(:,2) = 0
  next = 1

  do nb=1,nbMaxGlobal
    flag = refineFlags(nb)
    bID = globalBlocks(nb)
    if ((bID.ne.-1).and.(flag.eq.FLAG_COARSE)) then

      if (debug) write(logu,*) "Block", bID, "marked for coarsening"
      ! Obtain elder sibling ID
      call siblings (bID, sib_list)
      eID = sib_list(1)

      if (debug) write(logu,*) "Family of", bID, ":", sib_list
      ! Check if this family is already registered. If so, increase its count.
      counted = .false.
      do nb1=1,next-1
        if (family(nb1,1).eq.eID) then
          if (debug) write(logu,*) bID, "'s family exists; counting block"
          family(nb1,2) = family(nb1,2) + 1
          counted = .true.
          exit
        end if
      end do

      ! If it's not, add it
      if (.not.counted) then
        family(next,1) = eID
        family(next,2) = 1
        next = next + 1
        if (debug) write(logu,*) bID, "'s family doesnt exist; adding new family"
      end if

    end if
  end do

  if (debug) write (logu,*) "Done counting block families"

  ! 2. Apply coarsening for any family were all 8 siblings were marked

  do nb=1,nbActive

    ! Proceed if the family has 8 marked members
    eID = family(nb,1)
    if ((eID.ne.-1).and.(family(nb,2).eq.8)) then

      call getOwner (eID, elder_owner)
      call siblings (eID, sib_list)
      call father (eID, fID)
      if (verbosity > 2) write(logu,'(1x,a,i0)') "Coarsening children of block ", fID
      if (debug) then
        write(logu,*) "Elder owner:", elder_owner
        write(logu,*) "family:", sib_list
        write(logu,*) "father's bID:", fID
      end if

      ! The owner of the elder sibling allocates space for the new father block
      if (rank.eq.elder_owner) then

        if (debug) write(logu,*) "Registering new father block", fID
        call put (fID, localBlocks, nbMaxProc, fatherIndex)

        if (fatherIndex.eq.-1) then
          write(logu,*) ""
          write(logu,'(a)') "Father block couldn't be assigned an empty slot!"
          write(logu,'(a)') "***ABORTING***"
          call clean_abort (ERROR_REGISTER_FATHER)
        end if

        if ((fatherIndex.lt.1).or.(fatherIndex.gt.nbMaxProc)) then
          write(logu,*) "Invalid fatherIndex returned!"
          write(logu,*) "fatherIndex:", fatherIndex
          write(logu,*) "***ABORTING***"
          call clean_abort (ERROR_INVALID_LOC_INDEX)
        end if

      end if

      ! Now copy the data of each sibling into the father's data space
      ! The siblings might be local or on another process' memory
      do sib=1,8

        sID = sib_list(sib)
        call getOwner (sID, sib_owner)
        if (debug) write(logu,*) "Sibling", sID, "is owned by", sib_owner

        ! SENDER-side operations

        ! If rank owns this sibling but does not own the elder, average
        ! the sibling's data and send it through MPI to the elder's owner
        if ((rank.ne.elder_owner).and.(rank.eq.sib_owner)) then

          ! Average sibling's data into buffer
          call find (sID, localBlocks, nbMaxProc, locIndx)
          do i=1,ncells_x/2
            do j=1,ncells_y/2
              do k=1,ncells_z/2
                ip = i*2-1
                jp = j*2-1
                kp = k*2-1
                do ieq=1,neqtot
                  buf(ieq,i,j,k) = sum( U(locIndx,ieq,ip:ip+1,jp:jp+1,kp:kp+1) ) / 8.0
                end do
              end do
            end do
          end do

          ! Send average data through MPI
          nData = neqtot*(ncells_x/2)*(ncells_y/2)*(ncells_z/2)
          call MPI_SEND ( buf, nData, mpi_real_kind, elder_owner, sID, &
            mpi_comm_world, ierr)
          if (debug) write(logu,*) "Sent", nData, "values to process", elder_owner
#ifdef BFIELD
          if (enable_flux_ct) then
            call restrict_child_faces (locIndx, r1, r2, r3)
            fbuf(1:m3)          = reshape( r1, [m3] )
            fbuf(m3+1:2*m3)     = reshape( r2, [m3] )
            fbuf(2*m3+1:3*m3)   = reshape( r3, [m3] )
            call MPI_SEND ( fbuf, 3*m3, mpi_real_kind, elder_owner, sID,       &
            mpi_comm_world, ierr)
          end if
#endif
          ! Free this sibling's local slot
          call pop (sID, localBlocks, nbMaxProc, nb1)

        end if

        ! RECEIVER-side operations

        ! The owner of the elder sibling either receives this sibling's data
        ! through MPI ...
        if ((rank.eq.elder_owner).and.(rank.ne.sib_owner)) then

          ! Received MPI data into appropriate memory space
          call siblingCoords (sID, sx, sy, sz)
          i1 = 1 + ncells_x/2*sx
          i2 = ncells_x/2*(1+sx)
          j1 = 1 + ncells_y/2*sy
          j2 = ncells_y/2*(1+sy)
          k1 = 1 + ncells_z/2*sz
          k2 = ncells_z/2*(1+sz)
          nData = neqtot*(i2-i1+1)*(j2-j1+1)*(k2-k1+1)
          call MPI_RECV ( U(fatherIndex, 1:neqtot, i1:i2, j1:j2, k1:k2), nData, &
            mpi_real_kind, sib_owner, sID, mpi_comm_world, mpistatus, ierr)
          if (debug) write(logu,*) "Received", nData, "values into local space", i1,i2,j1,j2,k1,k2
#ifdef BFIELD
          if (enable_flux_ct) then
            call MPI_RECV ( fbuf, 3*m3, mpi_real_kind, sib_owner, sID,         &
                            mpi_comm_world, mpistatus, ierr)
            r1 = reshape( fbuf(1:m3),        [ncells_x/2+1, ncells_x/2, ncells_x/2] )
            r2 = reshape( fbuf(m3+1:2*m3),   [ncells_x/2, ncells_x/2+1, ncells_x/2] )
            r3 = reshape( fbuf(2*m3+1:3*m3), [ncells_x/2, ncells_x/2, ncells_x/2+1] )
            call put_child_faces (fatherIndex, sx, sy, sz, r1, r2, r3)
          end if
#endif

        ! ... or averages and copies it directly from local memory
        else if ((rank.eq.elder_owner).and.(rank.eq.sib_owner)) then

          call siblingCoords (sID, sx, sy, sz)
          call find (sID, localBlocks, nbMaxProc, locIndx)

          do i=1,ncells_x/2
            do j=1,ncells_y/2
              do k=1,ncells_z/2
                i1 = i*2-1
                j1 = j*2-1
                k1 = k*2-1
                ip = i + sx*ncells_x/2
                jp = j + sy*ncells_y/2
                kp = k + sz*ncells_z/2

                do ieq=1,neqtot

                if ((fatherIndex.lt.1).or.(fatherIndex.gt.nbMaxProc)) then
                  write(logu,*) "invalid fatherIndex !!"
                  write(logu,*) "fatherIndex=", fatherIndex
                  call clean_abort (ERROR_INVALID_LOC_INDEX)
                end if
                if (locIndx.eq.0) then
                  write(logu,*) "locIndx = 0 !!"
                  call clean_abort (ERROR_INVALID_LOC_INDEX)
                end if

                  U(fatherIndex,ieq,ip,jp,kp) = sum( U(locIndx,ieq,i1:i1+1,j1:j1+1,k1:k1+1) ) / 8.0

                end do

              end do
            end do
          end do

          if (debug) write(logu,*) "Averaged and copied local data"

#ifdef BFIELD
          if (enable_flux_ct) then
            call restrict_child_faces (locIndx, r1, r2, r3)
            call put_child_faces (fatherIndex, sx, sy, sz, r1, r2, r3)
          end if
#endif
          ! Free this sibling's local slot
          call pop (sID, localBlocks, nbMaxProc, nb1)

        end if

        if (debug) then
          if ((rank.ne.elder_owner).and.(rank.ne.sib_owner)) then
            write(logu,*) "No action required."
          end if
        end if

      end do

      ! Calculate primitives on new father block
      if (rank.eq.elder_owner) then
#ifdef BFIELD
        if (enable_flux_ct) call faces_to_cells (fatherIndex)
#endif
        call calcPrimsBlock (U, PRIM, fatherIndex, CELLS_PHYS, badcells)
      end if


      ! Calculate primitives on new father block
      if (rank.eq.elder_owner) then
        call calcPrimsBlock (U, PRIM, fatherIndex, CELLS_PHYS, badcells)
      end if

    end if
  end do

end subroutine coarseBlocks

!===============================================================================

!> @brief Returns whether the cell coordinates are invalid
!> @details A cell is considered invalid if:
!! a) At least one of its indices exceeds the ranges of physical+ghost cells
!! b) At least two of its indices are in ghost cell range
subroutine validCell (i, j, k, isValid)

  use parameters
  implicit none

  integer, intent(in) :: i
  integer, intent(in) :: j
  integer, intent(in) :: k
  logical, intent(out) :: isValid

  integer :: counter

  if ((i.lt.nxmin).or.(i.gt.nxmax).or.   &
      (j.lt.nymin).or.(j.gt.nymax).or.   &
      (k.lt.nzmin).or.(k.gt.nzmax)) then
    isValid = .false.
    return
  end if

  counter = 0
  if ((i.lt.1).or.(i.gt.ncells_x)) counter = counter + 1
  if ((j.lt.1).or.(j.gt.ncells_y)) counter = counter + 1
  if ((k.lt.1).or.(k.gt.ncells_z)) counter = counter + 1
  if (counter.le.1) then
    isValid = .true.
  else
    isValid = .false.
  end if

end subroutine

!===============================================================================

!> @brief Returns the global cell coords in max-resolution grid
!> @details Given the cell position i,j,k of block bID, returns the
!! ip,jp,kp global coordinates of that cell as if the whole simulation
!! box were at the maximum refinement level. Note that this will use
!! the reference corner (as opposed to the cell center) of the cell
!! when the block is not at the highest refinement level.
!> @param i,j,k The block-level integer coordinates of the cell
!> @param ip,jp,kp The returned absolute integer coordinates of the cell
!! at the maximum refinement level
subroutine absCoords (bID, i, j, k, ip, jp, kp)

  use parameters
  use globals
  implicit none

  integer, intent(in) :: bID
  integer, intent(in) :: i, j, k
  integer, intent(out) :: ip, jp, kp

  integer :: ilev, bx, by, bz

  ! Obtain block coords
  call meshlevel (bID, ilev)
  call bcoords(bID, bx, by, bz)

  ip = (bx-1)*ncells_x*2**(maxlev-ilev) + (i-1)*2**(maxlev-ilev) + 1
  jp = (by-1)*ncells_y*2**(maxlev-ilev) + (j-1)*2**(maxlev-ilev) + 1
  kp = (bz-1)*ncells_z*2**(maxlev-ilev) + (k-1)*2**(maxlev-ilev) + 1

end subroutine absCoords


!===============================================================================

!> @brief Performs load balancing among all processors
!> @details It first builds a loadOrder of the blocks following a Hilbert walk.
!! Then, it builds a loadScheme, which contains a list of all the active
!! blocks and the current and new owner for each block. Finally, it applies
!! the computed loadScheme.
!! skip primit optional  flag allows to use before Uvars or Primit
!! have been set
subroutine doBalance ()

  use parameters
  use globals
  use tictoc
  use clean_quit, only : clean_abort
  use hydro_core, only : calcPrimsBlock
  implicit none
  integer :: loadOrder(nbMaxGlobal)
  integer :: numBlocks(0:nProcs-1)
  integer :: loadScheme(nbActive, 3), tots(nProcs)
  integer :: nb, nbloc, p, c, bID, nData
  integer :: old_owner, new_owner, blocks, badcells
  integer :: mpistatus(MPI_STATUS_SIZE)
  integer :: mark
  integer :: sent, received
  logical, parameter :: debug = .false.

  ! The goal is to build and apply the new load balancing scheme, which is
  ! containted in the loadScheme array. Its structure is:
  !   loadScheme(:,1) = block's bID
  !   loadScheme(:,2) = block's current owner
  !   loadScheme(:,3) = block's new owner
  ! After it's been built, it is used walked one entry at a time to transfer
  ! blocks between processes

  if (verbosity > 3) call tic(mark)

  ! First, build an ordered list of blocks according to its position in
  ! in a Hilbert walk at the finest level
  loadOrder(:) = -1
  call HilbertOrder (loadOrder)

  ! Determine how many blocks each process gets. Process 0 is assigned last.
  numBlocks(:) = floor(real(nbActive)/nProcs)
  do nb=1,mod(nbActive,nProcs)
    p = mod(nb,nProcs)
    numBlocks(p) = numBlocks(p) + 1
  end do

  ! Initialize the loadScheme list using loadOrder and the current owner
  c = 1
  do nb=1,nbMaxGlobal
    if (loadOrder(nb).ne.-1) then
      bID = loadOrder(nb)
      call getOwner (bID, old_owner)
      loadScheme(c,1) = bID
      loadScheme(c,2) = old_owner
      loadScheme(c,3) = -1
      c = c + 1
    end if
  end do

  ! Next, build the new loadScheme
  nb = 1
  do p=0,nProcs-1
    blocks = numBlocks(p)
    do c=1,blocks
      loadScheme(nb,3) = p
      nb = nb + 1
    end do
  end do

  ! DEBUG
  if (debug) then
    ! Report the loadScheme
    write(logu,'(1x,a)') "The new load scheme is:"
    write(logu,'(6x,a)') "bID   old owner  ->   new owner"
    tots(:) = 0
    do nb=1,nbActive
      bID = loadScheme(nb,1)
      old_owner = loadScheme(nb,2)
      new_owner = loadScheme(nb,3)
      write(logu,'(1x,i8,i8,6x,a,i8)') bID, old_owner, "->", new_owner
      tots(new_owner+1) = tots(new_owner+1) + 1
    end do
    write(logu,'(1x,a)') "Summary:"
    do p=0,nProcs-1
      write(logu,'(2x,a,i8,a,i4,a)') "Process ", p, " gets ", tots(p+1), " blocks"
    end do
  end if

  ! Everybody stop here
  call mpi_barrier (mpi_comm_world, ierr)

  ! DEBUG
  if (debug) then
    write(logu,*) "There are", nbActive, "blocks globally"
    write(logu,*) "I have", nbLocal, "local blocks"
  end if
  ! DEBUG

  ! Finally, distribute the blocks via MPI messages
  !
  ! For every block, compare current and new owner.
  ! Then, one of four things will happen:
  ! 1) if the old owner is also the new owner, ignore this code
  !   (this block doesn't change process)
  ! 2) if they are different, but neither is my rank, ignore this code
  !   (I'm not involved in this exchange)
  ! 3) if they are different, and I'm the old owner, MPI_SEND the block to
  !    the new owner
  ! 4) if they are different, and I'm the new owner, MPI_RECV the block from
  !    the old owner
  if (verbosity > 3) call tic(mark)
  sent = 0
  received = 0
  do nb=1,nbActive

    bID = loadScheme(nb,1)
    old_owner = loadScheme(nb,2)
    new_owner = loadScheme(nb,3)

    if (old_owner.ne.new_owner) then
      if (rank.eq.old_owner) then

        ! Transmit this block (flow variables on physical cells)

        call find (bID, localBlocks, nbMaxProc, nbloc)
        if (nbloc.eq.-1) then
          write(logu,*) "Couldn't find bID", bID, "in local list!!"
          write(logu,*) "***ABORTING!***"
          call clean_abort (ERROR_LOADBAL_NO_LOCAL_BID)
        end if

        if (debug) then
          write(logu,'(1x,a,i8,a,i4,a,i4)') "Transmitting block ", bID, &
            " with local index ", nbloc, " to process ", new_owner
        end if
        nData = neqtot*ncells_x*ncells_y*ncells_z
        call MPI_SEND( U(nbloc, 1:neqtot, 1:ncells_x, 1:ncells_y, 1:ncells_z), &
             nData, mpi_real_kind, new_owner, bID, mpi_comm_world, ierr)
#ifdef BFIELD
        if (enable_flux_ct) then
          nData = 3*(ncells_x+1)*(ncells_y+1)*(ncells_z+1)
          call MPI_SEND( BF(nbloc, 1:3, 0:ncells_x, 0:ncells_y, 0:ncells_z), &
               nData, mpi_real_kind, new_owner, bID, mpi_comm_world, ierr)
        end if
#endif
        localBlocks(nbloc) = -1

        sent = sent + 1

      else if (rank.eq.new_owner) then

        ! Receive this block (flow variables on physical cells)

        call find (-1, localBlocks, nbMaxProc, nbloc)
        if (nbloc.eq.-1) then
          write(logu,*) "Couldn't find a free local slot!"
          write(logu,*) "***ABORTING!***"
          call clean_abort (ERROR_LOADBAL_INSUFFICIENT_RAM)
        end if

        if (debug) then
          write(logu,'(1x,a,i8,a,i4,a,i4)') "Receiving block ", bID, &
            " into free local index ", nbloc, " from process ", old_owner
        end if
        nData = neqtot*ncells_x*ncells_y*ncells_z
        call MPI_RECV( U(nbloc, 1:neqtot, 1:ncells_x, 1:ncells_y, 1:ncells_z), &
             nData, mpi_real_kind, old_owner, bID, mpi_comm_world, mpistatus, ierr)
#ifdef BFIELD
        if (enable_flux_ct) then
          nData = 3*(ncells_x+1)*(ncells_y+1)*(ncells_z+1)
          call MPI_RECV( BF(nbloc, 1:3, 0:ncells_x, 0:ncells_y, 0:ncells_z), &
               nData, mpi_real_kind, old_owner, bID, mpi_comm_world, mpistatus, ierr)
        end if
#endif
        localBlocks(nbloc) = bID

        ! Update primitives of received block
         call calcPrimsBlock (U, PRIM, nbloc, CELLS_PHYS, badcells)

         received = received + 1

      end if
    end if

  end do


  ! Re-calculate number of local blocks
  nbLocal = 0
  do nb=1,nbMaxProc
    if (localBlocks(nb).ne.-1) then
      nbLocal = nbLocal + 1
    end if
  end do
  if (verbosity > 2)  then
    ! Report load balance results
    write(logu,'(1x,a,i0,a,i0,a)') "Sent ", sent, " blocks; received ", received, " blocks"
    write(logu,'(1x,a,i0,a)') "Now have ", nbLocal, " local blocks"
  end if
  ! Synchronize global block registry
  call syncBlockLists ()

  ! DEBUG
  if (debug) then
    write(logu,*) "Local blocks after syncBlockLists:"
    do nb=1,nbMaxProc
      if (localBlocks(nb).ne.-1) then
        write(logu,*) nb, localBlocks(nb)
      end if
    end do
  end if

end subroutine doBalance

!===============================================================================

!> @brief Returns a load-balancing ordered list of blocks
!> @details The list is obtained by calculating the index of each block
!! along a Hilbert curve, and then sorting the block list using the indices
!> @param loadOrder The load-balanced ordered list of bIDs
subroutine HilbertOrder (loadOrder)

  use parameters
  use globals
  use utils, only : QuickSort
  implicit none

  integer, intent(out) :: loadOrder(nbMaxGlobal)

  integer :: nb, next, bID, istat
  integer :: hk
  integer, allocatable :: blocks(:)
  integer, allocatable :: keys(:)

  ! Re-calculate nbActive, declare blocks(:) and keys(:), and fill them
  ! with blocks and Hilbert Keys

  nbActive = 0
  do nb=1,nbMaxGlobal
    if (globalBlocks(nb).ne.-1) then
      nbActive = nbActive + 1
    end if
  end do

  allocate( blocks(nbActive), stat=istat)
  allocate( keys(nbActive), stat=istat)
  blocks(:) = -1
  keys(:) = -1

  ! Build lists of blocks with corresponding Hilbert indices as keys
  next = 1
  do nb=1,nbMaxGlobal
    if (globalBlocks(nb).ne.-1) then
      ! Obtain this block's Hilbert key at the finest mesh level
      bID = globalBlocks(nb)
      call getHilbertKey(bID, hk)
      blocks(next) = bID
      keys(next) = hk     ! Implicit conversion to integer here
      next = next + 1
    end if
  end do

  ! Sort 'blocks' list using 'keys' through QuickSort
  call QuickSort (nbActive, blocks, keys, 1, nbActive)

  ! Fill loadOrder list with sorted blocks list
  loadOrder(:) = -1
  do nb=1,nbActive
    loadOrder(nb) = blocks(nb)
  end do

  ! DEBUG
!  if (rank.eq.master) then
!    write(logu,*) "Hilbert indexing:"
!    do nb=1,nbActive
!      bID = loadOrder(nb)
!      call bcoords(bID, x, y, z)
!      write(logu,*) bID, x, y, z, keys(nb)
!    end do
!  end if
  ! DEBUG

  ! Don't forget to deallocate these arrays :)
  deallocate (blocks)
  deallocate (keys)

end subroutine HilbertOrder

!===============================================================================

!> @brief Returns a load-balancing ordered list of blocks (test routine)
!> @details The list is obtained by order in the globalBlocks list
!> @param loadOrder The load-balanced ordered list of bIDs
subroutine NaiveOrder (loadOrder)

  use parameters
  use globals
  implicit none

  integer, intent(out) :: loadOrder(nbMaxGlobal)

  integer :: nb, next

  loadOrder(:) = -1
  next = 1
  do nb=1,nbMaxGlobal
    if (globalBlocks(nb).ne.-1) then
      loadOrder(next) = globalBlocks(nb)
      next = next + 1
    end if
  end do

end subroutine NaiveOrder

!===============================================================================


!> @brief Calculates the Hilbert Key of a bID at the finest mesh level
!> @details If the block's level is not the finest,
!> @param bID The global ID of the block
!> @param hk The block's Hilbert key
subroutine getHilbertKey (bID, hk)

  use parameters
  use globals
  use hilbert, only : HKey3

  implicit none

  integer, intent(in) :: bID
  integer, intent(out) :: hk

  integer :: maxnb, order, x, y, z, ilev, xp, yp, zp

  ! Determine from the base grid geometry the "bounding cube" for the
  ! Hilbert curve, and hence the order of the curve
  maxnb = max(nbx(maxlev), nby(maxlev), nbz(maxlev))
  order = ceiling(log(1.0*maxnb)/log(2.0))

  ! Obtain the block's local coords (own mesh level) and convert
  ! them to coords at the finest mesh level
  call bcoords(bID, x, y, z)
  call meshlevel(bID, ilev)
  xp = (x-1)*2**(maxlev-ilev) + 1
  yp = (y-1)*2**(maxlev-ilev) + 1
  zp = (z-1)*2**(maxlev-ilev) + 1

  ! Get the Hilbert Key from the subroutine in hilbert.f90
  hk = -1
  call HKey3(xp-1, yp-1, zp-1, order, hk)

end subroutine getHilbertKey

!===============================================================================

!> @brief Divergence-free prolongation of the face field of one parent cell
!! (Toth & Roe 2002, sect. 3.2), for cubic cells
!> @details Exterior fine faces (eq. 4): parent face value plus limited
!! transverse slopes (+-1/4 of a parent cell), so the four fine faces of a
!! parent face average to it (eq. 3). Interior fine faces (eqs. 8-12): all
!! eight fine cells get the parent's divergence. With dx=dy=dz all face areas
!! are equal and the formulas hold for field values (the h^2 factors of the
!! third-order terms cancel with the 2h^2 of eq. 12). The parent's ghost
!! faces (depth 1) are used for the slopes at the block edges.
!> @param fInd Local index of the parent block
!> @param cInd Local index of the child block
!> @param ip,jp,kp Parent cell (parent block indices)
!> @param i0,j0,k0 Low fine cell of this parent cell in the child block
subroutine prolong_cell_faces (fInd, cInd, ip, jp, kp, i0, j0, k0)

  use parameters
  use globals,    only : BF
  use hydro_core, only : average
  implicit none

  integer, intent(in) :: fInd, cInd, ip, jp, kp, i0, j0, k0

  ! Fine faces: fx(n,j,k) = u^{n,j,k}, fy(i,n,k) = v^{i,n,k}, fz(i,j,n) = w^{i,j,n}
  ! with n = 0,1,2 for T&R -2,0,+2 and transverse 0,1 for -1,+1
  real    :: fx(0:2,0:1,0:1), fy(0:1,0:2,0:1), fz(0:1,0:1,0:2)
  real    :: bc, s1, s2
  integer :: a, b, c, l

  ! --- Exterior faces (eq. 4); a = 0: low parent face, a = 2: high parent face
  do a=0,2,2
    l = a/2 - 1                ! parent face index offset: -1 (low) or 0 (high)

    ! x-faces: u = U + (j Uy + k Uz)/4
    bc = BF(fInd,1,ip+l,jp,kp)
    s1 = average( bc - BF(fInd,1,ip+l,jp-1,kp), BF(fInd,1,ip+l,jp+1,kp) - bc, limiter_type )
    s2 = average( bc - BF(fInd,1,ip+l,jp,kp-1), BF(fInd,1,ip+l,jp,kp+1) - bc, limiter_type )
    do c=0,1
      do b=0,1
        fx(a,b,c) = bc + 0.25*( (2*b-1)*s1 + (2*c-1)*s2 )
      end do
    end do

    ! y-faces: v = V + (i Vx + k Vz)/4
    bc = BF(fInd,2,ip,jp+l,kp)
    s1 = average( bc - BF(fInd,2,ip-1,jp+l,kp), BF(fInd,2,ip+1,jp+l,kp) - bc, limiter_type )
    s2 = average( bc - BF(fInd,2,ip,jp+l,kp-1), BF(fInd,2,ip,jp+l,kp+1) - bc, limiter_type )
    do c=0,1
      do b=0,1
        fy(b,a,c) = bc + 0.25*( (2*b-1)*s1 + (2*c-1)*s2 )
      end do
    end do

    ! z-faces: w = W + (i Wx + j Wy)/4
    bc = BF(fInd,3,ip,jp,kp+l)
    s1 = average( bc - BF(fInd,3,ip-1,jp,kp+l), BF(fInd,3,ip+1,jp,kp+l) - bc, limiter_type )
    s2 = average( bc - BF(fInd,3,ip,jp-1,kp+l), BF(fInd,3,ip,jp+1,kp+l) - bc, limiter_type )
    do c=0,1
      do b=0,1
        fz(b,c,a) = bc + 0.25*( (2*b-1)*s1 + (2*c-1)*s2 )
      end do
    end do
  end do

  call tr_interior(fx, fy, fz)

  ! --- Store in the child block (face index i = face i+1/2)
  do c=0,1
    do b=0,1
      do a=0,2
        BF(cInd,1,i0-1+a,j0+b,k0+c) = fx(a,b,c)
        BF(cInd,2,i0+b,j0-1+a,k0+c) = fy(b,a,c)
        BF(cInd,3,i0+b,j0+c,k0-1+a) = fz(b,c,a)
      end do
    end do
  end do

end subroutine prolong_cell_faces

!===============================================================================

!> @brief Interior fine faces of one parent cell from its 24 exterior fine
!! faces (Toth & Roe 2002, eqs. 8-12, cubic cells)
!> @param fx,fy,fz Fine faces as in prolong_cell_faces; on input the exterior
!! ones (n = 0, 2), on output also the interior ones (n = 1)
subroutine tr_interior (fx, fy, fz)

  implicit none

  real, intent(inout) :: fx(0:2,0:1,0:1), fy(0:1,0:2,0:1), fz(0:1,0:1,0:2)

  real    :: t1, t2, si, sj, sk, uu, vv, ww
  real    :: uxx, vyy, wzz, uxyz, vxyz, wxyz
  integer :: b, c, ii, jj, kk

  ! --- Second and third derivatives (eqs. 11, 12 and cyclic permutations)
  uxx = 0.0 ; vyy = 0.0 ; wzz = 0.0 ; uxyz = 0.0 ; vxyz = 0.0 ; wxyz = 0.0

  do kk=0,1
    do jj=0,1
      do ii=0,1
        si = real(2*ii-1)
        sj = real(2*jj-1)
        sk = real(2*kk-1)
        uu = fx(2*ii, jj, kk)          ! u^{2i,j,k}
        vv = fy(ii, 2*jj, kk)          ! v^{i,2j,k}
        ww = fz(ii, jj, 2*kk)          ! w^{i,j,2k}
        uxx  = uxx  + si*sj*vv + si*sk*ww
        vyy  = vyy  + sj*sk*ww + sj*si*uu
        wzz  = wzz  + sk*si*uu + sk*sj*vv
        uxyz = uxyz + si*sj*sk*uu
        vxyz = vxyz + si*sj*sk*vv
        wxyz = wxyz + si*sj*sk*ww
      end do
    end do
  end do
  uxx  = 0.125*uxx
  vyy  = 0.125*vyy
  wzz  = 0.125*wzz
  uxyz = 0.0625*uxyz           ! = (dy^2) * U_xyz of eq. 12, for cubic cells
  vxyz = 0.0625*vxyz
  wxyz = 0.0625*wxyz

  ! --- Interior faces (eqs. 8-10); t1, t2 = signs of the two transverse indices
  do c=0,1
    do b=0,1
      t1 = real(2*b-1)
      t2 = real(2*c-1)
      ! u^{0,j,k}: j = t1, k = t2
      fx(1,b,c) = 0.5*(fx(2,b,c) + fx(0,b,c)) + uxx + t2*vxyz + t1*wxyz
      ! v^{i,0,k}: i = t1, k = t2
      fy(b,1,c) = 0.5*(fy(b,2,c) + fy(b,0,c)) + vyy + t1*wxyz + t2*uxyz
      ! w^{i,j,0}: i = t1, j = t2
      fz(b,c,1) = 0.5*(fz(b,c,2) + fz(b,c,0)) + wzz + t2*uxyz + t1*vxyz
    end do
  end do

end subroutine tr_interior

!===============================================================================

!> @brief Recomputes the 12 interior fine faces of one parent cell of a
!! block from the block's own exterior faces (Toth & Roe eqs. 8-12)
!> @param cInd Local index of the (child) block
!> @param i0,j0,k0 Low fine cell of the parent cell
subroutine reset_interior_faces (cInd, i0, j0, k0)

  use parameters
  use globals, only : BF
  implicit none

  integer, intent(in) :: cInd, i0, j0, k0

  real    :: fx(0:2,0:1,0:1), fy(0:1,0:2,0:1), fz(0:1,0:1,0:2)
  integer :: a, b, c

  do c=0,1
    do b=0,1
      do a=0,2,2
        fx(a,b,c) = BF(cInd,1,i0-1+a,j0+b,k0+c)
        fy(b,a,c) = BF(cInd,2,i0+b,j0-1+a,k0+c)
        fz(b,c,a) = BF(cInd,3,i0+b,j0+c,k0-1+a)
      end do
    end do
  end do

  call tr_interior (fx, fy, fz)

  do c=0,1
    do b=0,1
      BF(cInd,1,i0,  j0+b,k0+c) = fx(1,b,c)
      BF(cInd,2,i0+b,j0,  k0+c) = fy(b,1,c)
      BF(cInd,3,i0+b,j0+c,k0  ) = fz(b,c,1)
    end do
  end do

end subroutine reset_interior_faces

!===============================================================================

!> @brief Face-field fix-up after refinement (Toth & Roe 2002, sect. 3.2)
!> @details Blocks created by refineBlock get their exterior faces from the
!! father (eq. 4). Where a new block touches a same-level block that already
!! existed, the shared faces are replaced by that block's (evolved) faces;
!! where two new blocks touch, both copies are averaged. Each group of four
!! fine faces still averages to the father's face, so the interior faces of
!! the new blocks are then recomputed (eqs. 8-12) and every fine cell gets
!! the father's divergence. Must be called by ALL ranks, after refinement
!! and syncBlockLists, before any load balance or coarsening.
!! Where a new block touches finer blocks (created in the same pass from a
!! block one level finer than its father), its face takes the average of
!! their faces, i.e. the old fine block's faces.
!> @param fathers bIDs of the blocks just refined (same list on all ranks)
!> @param nf Number of entries in fathers
subroutine fix_new_faces (fathers, nf)

  use parameters
  use globals
  use hydro_core, only : calcPrimsBlock
  implicit none

  integer, intent(in) :: nf
  integer, intent(in) :: fathers(nf)

  integer :: axis, nb, aID, bID, aOwner, bOwner, aInd, bInd, ntype, neighs(4)
  integer :: ic, jc, kc, nData, badcells
  integer :: mpistatus(MPI_STATUS_SIZE)
  logical :: aNew, bNew
  real    :: fa(ncells_x,ncells_x), fb(ncells_x,ncells_x), dmax, gmax
  integer :: dir, ib, h, lc, lf, sx, sy, sz, oa, ob, a, b
  real    :: fr(ncells_x/2,ncells_x/2)  ! ncells_x = ncells_block


  if (nf == 0) return              ! nf is the same on all ranks

  nData = ncells_x**2
  aInd  = -1
  bInd  = -1
  dmax  = 0.0

  ! 1. Shared faces of same-level pairs with at least one new block
  do axis=1,3              ! RIGHT/BACK/TOP face of a = LEFT/FRONT/BOTTOM of b
    do nb=1,nbMaxGlobal
      aID = globalBlocks(nb)
      if (aID == -1) cycle
      call neighbors (aID, 2*axis, ntype, neighs)
      if (ntype /= NEIGH_SAME) cycle
      bID  = neighs(1)
      aNew = is_new(aID)
      bNew = is_new(bID)
      if (.not.(aNew.or.bNew)) cycle
      call getOwner (aID, aOwner)
      call getOwner (bID, bOwner)
      if ((rank /= aOwner).and.(rank /= bOwner)) cycle

      if (rank == aOwner) then
        call find (aID, localBlocks, nbMaxProc, aInd)
        call get_plane (aInd, axis, ncells_x, fa)          ! a's high face
      end if
      if (rank == bOwner) then
        call find (bID, localBlocks, nbMaxProc, bInd)
        call get_plane (bInd, axis, 0, fb)                 ! b's low face
      end if
      if ((rank == aOwner).and.(rank /= bOwner)) then
        call MPI_SENDRECV (fa, nData, mpi_real_kind, bOwner, aID,             &
                           fb, nData, mpi_real_kind, bOwner, bID,             &
                           mpi_comm_world, mpistatus, ierr)
      else if ((rank == bOwner).and.(rank /= aOwner)) then
        call MPI_SENDRECV (fb, nData, mpi_real_kind, aOwner, bID,             &
                           fa, nData, mpi_real_kind, aOwner, aID,             &
                           mpi_comm_world, mpistatus, ierr)
      end if

      dmax = max( dmax, maxval(abs(fa - fb)) )
      if (aNew.and.bNew) then
        fa = 0.5*(fa + fb)
      else if (aNew) then
        fa = fb                    ! b existed: its faces win
      end if                       ! (only b new: a's faces win, fa as is)
      if ((rank == aOwner).and.aNew) call put_plane (aInd, axis, ncells_x, fa)
      if ((rank == bOwner).and.bNew) call put_plane (bInd, axis, 0, fa)
    end do
  end do

  ! 1b. Faces of new blocks next to finer blocks. Those were created in this
  !     same pass by refining a block one level finer than this block's
  !     father (forced by the 2:1 balance). Their faces average exactly to
  !     the old fine block's faces, which this block must take (Toth & Roe:
  !     copy the fine fluxes) instead of the ones prolonged from its father.
  h = ncells_x/2
  oa = 0
  ob = 0
  do dir=1,6
    axis = (dir+1)/2
    lc   = merge(0, ncells_x, mod(dir,2) == 1)  ! this block's face
    lf   = ncells_x - lc                        ! the finer blocks' face
    do nb=1,nbMaxGlobal
      aID = globalBlocks(nb)
      if (aID == -1) cycle
      if (.not.is_new(aID)) cycle
      call neighbors (aID, dir, ntype, neighs)
      if (ntype /= NEIGH_FINER) cycle
      call getOwner (aID, aOwner)

      do ib=1,4
        bID = neighs(ib)
        call getOwner (bID, bOwner)

        ! Fine side: average of each 2x2 group of faces
        if (rank == bOwner) then
          call find (bID, localBlocks, nbMaxProc, bInd)
          call get_plane (bInd, axis, lf, fb)
          do b=1,h
            do a=1,h
              fr(a,b) = 0.25*( fb(2*a-1,2*b-1) + fb(2*a,2*b-1)                &
                             + fb(2*a-1,2*b  ) + fb(2*a,2*b  ) )
            end do
          end do
          if (rank /= aOwner) call MPI_SEND (fr, h*h, mpi_real_kind, aOwner,  &
                                             bID, mpi_comm_world, ierr)
        end if

        ! New (coarse) side: replace the quadrant of its face
        if (rank == aOwner) then
          if (rank /= bOwner) call MPI_RECV (fr, h*h, mpi_real_kind, bOwner,  &
                                       bID, mpi_comm_world, mpistatus, ierr)
          call find (aID, localBlocks, nbMaxProc, aInd)
          call siblingCoords (bID, sx, sy, sz)
          select case (axis)
          case (1)
            oa = sy*h ; ob = sz*h
          case (2)
            oa = sx*h ; ob = sz*h
          case (3)
            oa = sx*h ; ob = sy*h
          end select
          call get_plane (aInd, axis, lc, fa)
          dmax = max( dmax, maxval(abs(fa(oa+1:oa+h, ob+1:ob+h) - fr)) )
          fa(oa+1:oa+h, ob+1:ob+h) = fr
          call put_plane (aInd, axis, lc, fa)
        end if
      end do

    end do
  end do


  ! 2. Interior faces, cell-centred B and primitives of the local new blocks
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID == -1) cycle
    if (.not.is_new(bID)) cycle
    do kc=1,ncells_z/2
      do jc=1,ncells_y/2
        do ic=1,ncells_x/2
          call reset_interior_faces (nb, 2*ic-1, 2*jc-1, 2*kc-1)
        end do
      end do
    end do
    call faces_to_cells (nb)
    call calcPrimsBlock (U, PRIM, nb, CELLS_PHYS, badcells)
  end do

#ifdef MPIP
  call mpi_reduce (dmax, gmax, 1, mpi_real_kind, mpi_max, master,             &
                   mpi_comm_world, ierr)
#else
  gmax = dmax
#endif
  if (rank == master) write(logu,'(1x,a,i8,es12.4)')                           &
       "fix_new_faces: it, max change of shared faces", it, gmax

contains

  !> True if block id was created by the refinement that is being fixed
  logical function is_new (id)
    integer, intent(in) :: id
    integer :: lev, fid
    is_new = .false.
    call meshlevel (id, lev)
    if (lev == 1) return
    call father (id, fid)
    is_new = any( fathers == fid )
  end function is_new

end subroutine fix_new_faces

!===============================================================================

!> @brief Normal-B face plane l (0..ncells) of a block, along axis
subroutine get_plane (ind, axis, l, f)

  use parameters
  use globals, only : BF
  implicit none

  integer, intent(in)  :: ind, axis, l
  real,    intent(out) :: f(ncells_x,ncells_x)

  select case (axis)
  case (1)
    f = BF(ind,1,l,1:ncells_y,1:ncells_z)
  case (2)
    f = BF(ind,2,1:ncells_x,l,1:ncells_z)
  case (3)
    f = BF(ind,3,1:ncells_x,1:ncells_y,l)
  end select

end subroutine get_plane

!===============================================================================

!> @brief Writes the normal-B face plane l of a block, along axis
subroutine put_plane (ind, axis, l, f)

  use parameters
  use globals, only : BF
  implicit none

  integer, intent(in) :: ind, axis, l
  real,    intent(in) :: f(ncells_x,ncells_x)

  select case (axis)
  case (1)
    BF(ind,1,l,1:ncells_y,1:ncells_z) = f
  case (2)
    BF(ind,2,1:ncells_x,l,1:ncells_z) = f
  case (3)
    BF(ind,3,1:ncells_x,1:ncells_y,l) = f
  end select

end subroutine put_plane

!===============================================================================

!> @brief Cell-centred B (U(6:8)) of one block = average of its two faces
subroutine faces_to_cells (ind)

  use parameters
  use globals, only : U, BF
  implicit none

  integer, intent(in) :: ind
  integer :: i, j, k

  do k=1,ncells_z
    do j=1,ncells_y
      do i=1,ncells_x
        U(ind,6,i,j,k) = 0.5*( BF(ind,1,i-1,j,k) + BF(ind,1,i,j,k) )
        U(ind,7,i,j,k) = 0.5*( BF(ind,2,i,j-1,k) + BF(ind,2,i,j,k) )
        U(ind,8,i,j,k) = 0.5*( BF(ind,3,i,j,k-1) + BF(ind,3,i,j,k) )
      end do
    end do
  end do

end subroutine faces_to_cells

!===============================================================================

!> @brief Restriction of a child's face field (Toth & Roe eq. 3, as area
!! averages) to the father's resolution, over the child's octant
!> @param cInd Local index of the child block
!> @param r1,r2,r3 Restricted Bx (faces 0:h), By, Bz in octant indices
subroutine restrict_child_faces (cInd, r1, r2, r3)

  use parameters
  use globals, only : BF
  implicit none

  integer, intent(in) :: cInd
  real, intent(out) :: r1(0:ncells_x/2,ncells_x/2,ncells_x/2)
  real, intent(out) :: r2(ncells_x/2,0:ncells_x/2,ncells_x/2)
  real, intent(out) :: r3(ncells_x/2,ncells_x/2,0:ncells_x/2)
  integer :: a, b, c, h

  h = ncells_x/2
  do c=1,h
    do b=1,h
      do a=0,h
        r1(a,b,c) = 0.25*( BF(cInd,1,2*a,2*b-1,2*c-1) + BF(cInd,1,2*a,2*b,2*c-1)  &
                         + BF(cInd,1,2*a,2*b-1,2*c  ) + BF(cInd,1,2*a,2*b,2*c  ) )
      end do
    end do
  end do
  do c=1,h
    do b=0,h
      do a=1,h
        r2(a,b,c) = 0.25*( BF(cInd,2,2*a-1,2*b,2*c-1) + BF(cInd,2,2*a,2*b,2*c-1)  &
                         + BF(cInd,2,2*a-1,2*b,2*c  ) + BF(cInd,2,2*a,2*b,2*c  ) )
      end do
    end do
  end do
  do c=0,h
    do b=1,h
      do a=1,h
        r3(a,b,c) = 0.25*( BF(cInd,3,2*a-1,2*b-1,2*c) + BF(cInd,3,2*a,2*b-1,2*c)  &
                         + BF(cInd,3,2*a-1,2*b,  2*c) + BF(cInd,3,2*a,2*b,  2*c) )
      end do
    end do
  end do

end subroutine restrict_child_faces

!===============================================================================

!> @brief Writes a restricted octant (restrict_child_faces) into the father
subroutine put_child_faces (fInd, sx, sy, sz, r1, r2, r3)

  use parameters
  use globals, only : BF
  implicit none

  integer, intent(in) :: fInd, sx, sy, sz
  real, intent(in) :: r1(0:ncells_x/2,ncells_x/2,ncells_x/2)
  real, intent(in) :: r2(ncells_x/2,0:ncells_x/2,ncells_x/2)
  real, intent(in) :: r3(ncells_x/2,ncells_x/2,0:ncells_x/2)
  integer :: h, ox, oy, oz

  h  = ncells_x/2
  ox = sx*h
  oy = sy*h
  oz = sz*h
  BF(fInd,1,ox:ox+h,  oy+1:oy+h, oz+1:oz+h) = r1
  BF(fInd,2,ox+1:ox+h, oy:oy+h,  oz+1:oz+h) = r2
  BF(fInd,3,ox+1:ox+h, oy+1:oy+h, oz:oz+h ) = r3

end subroutine put_child_faces

!===============================================================================

end module amr
