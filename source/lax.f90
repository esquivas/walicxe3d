!===============================================================================
!> @file lax.f90
!> @brief Lax-Friedrichs first-order solver (deprecated)
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

!> @brief Lax Friedrichs module
!> @details The module contains upper level wrapper to advance the equations
!> in time using a 1st order Lax Friedrichs method, very simple and robust
!> but exteremelly diffusive
module Lax

  implicit none

contains

!===============================================================================

!> @brief First order Lax-Friefrichs scheme
subroutine LaxFriedrichs ()

  use parameters
  use globals
  use tictoc
  use hydro_core, only : prim2fluxes, calcPrimsAll
  use amr,        only : meshlevel
  use boundaries, only : boundary
  implicit none

  integer :: nb, bID, lev, i, j, k, mark
  real :: pvars(neqtot)

  ! Exchange boundaries to fill ghost cells (all blocks)
  call boundary (1, U)

  ! Update primitives in ghost cells
  call calcPrimsAll (U, PRIM, CELLS_GHOST)

  ! Integrate all local blocks
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then

      call tic (mark)

      ! Compute physical fluxes for all physical cells, plus one additional
      ! ghost cell to each side
      do i=0,ncells_x+1
        do j=0,ncells_y+1
          do k=0,ncells_z+1
            pvars(:) = PRIM(nb,:,i,j,k)
            call prim2fluxes (pvars, DIM_X, FC(:,i,j,k))
            call prim2fluxes (pvars, DIM_Y, GC(:,i,j,k))
            call prim2fluxes (pvars, DIM_Z, HC(:,i,j,k))
          end do
        end do
      end do

      ! Compute UPs using Lax-Friedrich's upwind scheme
      call meshlevel (bID,lev)
      do i=1,ncells_x
        do j=1,ncells_y
          do k=1,ncells_z
            UP(nb,:,i,j,k) = &
              ( U(nb,:,i+1,j,k) + U(nb,:,i-1,j,k) &
              + U(nb,:,i,j+1,k) + U(nb,:,i,j-1,k) &
              + U(nb,:,i,j,k+1) + U(nb,:,i,j,k-1) ) / 6.0 &
              - 0.5*dt/dx(lev)*(FC(:,i+1,j,k)-FC(:,i-1,j,k)) &
              - 0.5*dt/dy(lev)*(GC(:,i,j+1,k)-GC(:,i,j-1,k)) &
              - 0.5*dt/dz(lev)*(HC(:,i,j,k+1)-HC(:,i,j,k-1))
          end do
        end do
      end do

!      write(logu,*) "Integrated block", bID, "in", nicetoc(mark)

    end if
  end do

end subroutine LaxFriedrichs

!===============================================================================

end module Lax
