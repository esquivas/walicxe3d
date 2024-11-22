!===============================================================================
!> @file godunov.f90
!> @brief Generic wrapper for Godunov schemes (HLL1, HLL, HLLC)
!> @author Juan C. Toledo
!> @date 7/Mar/2012

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

!> @brief Godunov Module
!> @details The module contains the subroutines and utilities advance the
!> equations with Godunov type solvers, the approximate Riemann solvers are
!> included in separate source files.

module GodunovModule

  implicit none

contains

!===============================================================================

!> @brief Wrapper for Godunov-type numerical solvers
!> @details Computes intercell numerical fluxes using any Riemann solver and
!! advances the flow variables one timestep. The solver can be set to be
!! first or second order (in space and time). This is a high-level wrapper.
subroutine Godunov (order)

  use parameters
  use globals
  use tictoc
  use hydro_core, only : calcPrimsAll
  use boundaries, only : boundary
  use HLL,        only : HLLfluxes
  use HLLC,       only : HLLCfluxes
#ifdef BFIELD
  use HLLE,       only : HLLEfluxes
  use HLLD,       only : HLLDfluxes
#endif
  implicit none

  integer, intent(in) :: order

  integer :: bIndx, bID, mark, bcount
  real :: dtp

  ! -----------------------------------
  ! 1st-order timestep

  ! Halve timestep for first part of second-order Godonov type solver
  if (order.eq.2) then
    dtp = dt/2.0
  else
    dtp = dt
  end if

  ! Exchange 1-deep boundaries to fill ghost cells (all blocks)
    if (verbosity > 3) call tic(mark)
    if (verbosity > 1) then
      write(logu,*) ""
      write(logu,'(1x,a)') "> Exchanging 1-deep boundary layers ..."
    end if
  call boundary (1, U)
  if (verbosity > 3) write(logu,*) "Boundaries exchanged in", nicetoc(mark)

  ! Update primitives in ghost cells
  call calcPrimsAll (U, PRIM, CELLS_GHOST)

  ! First-order integration of all local blocks
  bcount = 0
  call tic(mark)
  if ( verbosity > 1 ) write(logu,*) ""
  if (order.eq.1) then
      if (verbosity > 1) write(logu,'(1x,a)') "> Integrating blocks ..."
  else if (order.eq.2) then
      if (verbosity > 1) write(logu,'(1x,a)') "> Integrating blocks (1st order half step) ..."
  end if

  do bIndx=1,nbMaxProc
    bID = localBlocks(bIndx)
    if (bID.ne.-1) then

      ! Compute numerical fluxes  for this block
      select case (solver_type)

        case (SOLVER_HLL1)
          call HLLfluxes (bIndx, 1)

        case (SOLVER_HLL)
          call HLLfluxes (bIndx, 1)

        case (SOLVER_HLLC)
          call HLLCfluxes (bIndx, 1)

        case (SOLVER_HLLE)
          call HLLEfluxes (bIndx, 1)

        case (SOLVER_HLLD)
          call HLLDfluxes (bIndx, 1)

      end select

      ! Apply conservative formula
      call upwindStep (bIndx, dtp)

      bcount = bcount + 1

    end if
  end do

    if (verbosity > 3) write(logu,'(1x,a,i0,a,a)') "Integrated ", bcount, " blocks in ", nicetoc(mark)

  ! -----------------------------------
  ! 2nd-order full timestep (skipped in 1st-order schemes)

  if (order.eq.2) then

    ! Exchange 2-deep boundaries to fill ghost cells (all blocks)
    if (verbosity > 3) call tic(mark)
    if (verbosity > 1) then
      write(logu,*) ""
      write(logu,'(1x,a)') "> Exchanging 2-deep boundary layers ..."
    end if

    call boundary (2, UP)

    if (verbosity > 3) write(logu,*) "Boundaries exchanged in", nicetoc(mark)

    ! Update primitives in ghost cells
    call calcPrimsAll (UP, PRIM, CELLS_ALL)

    ! Second-order integration of all local blocks
    bcount = 0
    if (verbosity > 3) call tic(mark)
    if (verbosity > 2) then
      write(logu,*) ""
      write(logu,'(1x,a)') "> Integrating blocks (2nd order full step) ..."
    endif
    do bIndx=1,nbMaxProc

      bID = localBlocks(bIndx)
      if (bID.ne.-1) then

        ! Compute numerical fluxes for this block
        select case (solver_type)

          case (SOLVER_HLL)
            call HLLfluxes (bIndx, 2)

          case (SOLVER_HLLC)
            call HLLCfluxes (bIndx, 2)

          case (SOLVER_HLLE)
            call HLLEfluxes (bIndx, 2)

          case (SOLVER_HLLD)
            call HLLDfluxes (bIndx, 2)

        end select

        ! Apply conservative formula
        call upwindStep (bIndx, dt)

        bcount = bcount + 1

      end if
    end do

      if (verbosity > 3) write(logu,'(1x,a,i0,a,a)') "Integrated ", bcount, " blocks in ", &
    nicetoc(mark)

  end if

end subroutine Godunov

!===============================================================================

!> @brief Applies conservative formula to calculate stepped U variables
!> @details This is the standard upwing Godunov step. Assumes the numerical
!! intercell fluxes FC have been calculated for this block
!> @param locIndx Block's local index
!> @param dtp Timestep to use
subroutine upwindStep (locIndx, dtp)

  use parameters
  use globals
  use clean_quit, only : clean_abort
  use sources
  use amr, only : meshlevel
  implicit none

  integer, intent(in) :: locIndx
  real,    intent(in) :: dtp

  integer :: i, j, k, bID, lev
  integer :: ieq   ! DEBUG
  logical, parameter :: debug_mode = .false.
  real :: dtdx, dtdy, dtdz
  real :: s(neqtot)

  bID = localBlocks(locIndx)
  call meshlevel (bID, lev)

  dtdx = dtp/dx(lev)
  dtdy = dtp/dy(lev)
  dtdz = dtp/dz(lev)

  ! Apply conservative formula (10.2 of Toro) to obtain UPs
  do i=1,ncells_x
    do j=1,ncells_y
      do k=1,ncells_z

! DEBUG Check FC, GC and HC values for NaNs
if (debug_mode) then
  do ieq=1,neqtot
    if (FC(ieq,i-1,j,k).ne.FC(ieq,i-1,j,k)) then
      write(logu,'(a,i2,i2,i2)') "Nan in FC at ", i-1, j, k
      call clean_abort(ERROR_GENERIC)
    end if
    if (FC(ieq,i,j,k).ne.FC(ieq,i,j,k)) then
      write(logu,'(a,i2,i2,i2)') "Nan in FC at ", i, j, k
      call clean_abort(ERROR_GENERIC)
    end if
    if (GC(ieq,i,j-1,k).ne.GC(ieq,i,j-1,k)) then
      write(logu,'(a,i2,i2,i2)') "Nan in GC at ", i, j-1, k
      call clean_abort(ERROR_GENERIC)
    end if
    if (GC(ieq,i,j,k).ne.GC(ieq,i,j,k)) then
      write(logu,'(a,i2,i2,i2)') "Nan in GC at ", i, j, k
      call clean_abort(ERROR_GENERIC)
    end if
    if (HC(ieq,i,j,k-1).ne.HC(ieq,i,j,k-1)) then
      write(logu,'(a,i2,i2,i2)') "Nan in HC at ", i, j, k-1
      call clean_abort(ERROR_GENERIC)
    end if
    if (HC(ieq,i,j,k).ne.HC(ieq,i,j,k)) then
      write(logu,'(a,i2,i2,i2)') "Nan in HC at ", i, j, k
      call clean_abort(ERROR_GENERIC)
    end if
  end do
end if

        UP(locIndx,:,i,j,k) = U(locIndx,:,i,j,k)                              &
                            + dtdx*(FC(:,i-1,j,k)-FC(:,i,j,k))                &
                            + dtdy*(GC(:,i,j-1,k)-GC(:,i,j,k))                &
                            + dtdz*(HC(:,i,j,k-1)-HC(:,i,j,k))

        !! Include source terms S
        !! (computed with the same primitives used in the fluxes calculation)
        if (eight_wave .or. user_source_terms ) then

          call source_function(locIndx,lev,i,j,k,s)
          UP(locIndx,:,i,j,k) = UP(locIndx,:,i,j,k) + dtp*s(:)

        end if

      end do
    end do
  end do

end subroutine upwindStep

!===============================================================================

end module GodunovModule
