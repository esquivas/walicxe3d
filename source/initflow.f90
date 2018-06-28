!===============================================================================
!> @file initflow.f90
!> @brief Initial flow conditions module
!> @author Juan C. Toledo
!> @date 2/Feb/2012

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

!> @brief Sets the initial flow conditions
!> @details This wrapper subroutine is called to set the initial flow
!! conditions on the mesh. It will first initialize the grid with a
!! uniform IC and then call userIC(), where user-defined custom
!! initial conditions can be specified.
subroutine initflow ()

  use parameters
  use globals
  use userconds
  implicit none

  write(logu,*) ""
  write(logu,*) "============================================"
  write(logu,'(1x,a)') " Setting Initial Conditions  ..."
  write(logu,*) "============================================"
  write(logu,*) ""

  ! IC (defined in user.f90)
  call setInitialCondition (U)

  write(logu,*) ""
  write(logu,'(1x,a)') "> Done setting ICs"

end subroutine initflow
