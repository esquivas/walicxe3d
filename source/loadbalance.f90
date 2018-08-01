!===============================================================================
!> @file loadbalance.f90
!> @brief Parallel execution load balancing
!> @author Juan C. Toledo
!> @date 13/Jan/2012

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

!> @brief Load balance wrapper
!> @details This module calls subroutines and utilities needed to balance
!> the workload among different processors, following a hilbert courve.

module load_balance

  use amr, only : doBalance
  implicit none

contains

!===============================================================================


!> @brief Performs load balancing among all processes (wrapper)
!> @details Redistributes the blocks among all processes, with the purpose of
!! balancing the load. This is the high-level wrapper subroutine, intended
!! to be called by the main program.
subroutine loadBalance ()

  use parameters
  use globals
  use tictoc
  implicit none

  integer :: mark

  call tic(mark)

  write(logu,*) ""
  write(logu,'(1x,a)') "============================================"
  write(logu,'(1x,a)') " Performing Load Balance ..."
  write(logu,'(1x,a)') "============================================"
  write(logu,*) ""

  ! Main load balacing routine
  call doBalance ()

  ! Done
  write(logu,*) ""
  write(logu,'(1x,a,a)') "> Load balance completed in ", nicetoc(mark)
  write(logu,*) ""

end subroutine loadBalance


!===============================================================================

end module load_balance
