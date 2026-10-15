!===============================================================================
!> @file ghosts.f90
!> @brief Refresh of the ghost cells and ghost faces of all local blocks
!> @author A. Esquivel
!> @date 14/Oct/2026

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

!> @brief Ghost-cell refresh, used by the main loop and (through the
!! amr refresh_ghosts hook) by refineZone
module ghosts

  implicit none

contains

!===============================================================================

!> @brief Refreshes one layer of ghost cells (U), ghost faces (BF) and the
!! primitives (including ghost cells) of all local blocks
!> @details Must be called by ALL ranks.
subroutine refresh_ghosts_all ()

  use parameters
  use globals,    only : U, PRIM
  use boundaries, only : boundary
  use hydro_core, only : calcPrimsAll
#ifdef BFIELD
  use globals,    only : BF
  use flux_ct,    only : boundary_bf
#endif
  implicit none

  call boundary (1, U)
#ifdef BFIELD
  if (enable_flux_ct) call boundary_bf (1, BF, U)
#endif
  call calcPrimsAll (U, PRIM, CELLS_GHOST)

end subroutine refresh_ghosts_all

!===============================================================================

end module ghosts
