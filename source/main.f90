!===============================================================================
!> @file main.f90
!> @brief Walixce3D main program unit
!> @author Juan C. Toledo
!> @date 2/Jun/2011

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
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see http://www.gnu.org/licenses/.

!===============================================================================

!> @brief Walicxe3D Main Program

!> This is the main program unit of the Walicxe3D code. It manages all the
!! high-level operations of the code, which include initializations, initial
!! flow conditions, adaptive mesh refinements, flow variables integration,
!! boundary conditions calculation, parallel load balancing.
program Walicxe3D

  use parameters    ! Code parameters
  use globals       ! Runtime global variables
  use tictoc,        only : tic, nicetoc, stamp       ! Time-measuring library
  use clean_quit,    only  : deinit
  use init
  use hydro_core,    only : updatePrims
  use amr,           only : admesh
  use load_balance,  only : loadBalance
  use output,        only : writeOutput
  use hydro_solver,  only : hydroSolver
  use coolingModule, only : cooling
  use report,        only : main_report
  use hrate,         only : updateNeutralFraction
  implicit none    ! ALWAYS mandatory

  ! Timing mark
  call tic(start_mark)

  ! Initialize global arrays and variables
  call initmain ()

  ! Initialize base grid, and refine it (cold start)
  call basegrid ()

  ! Impose initial conditions (cold start) or load them from file (warm start)
  if (dowarm) then
    call warmstart ()
  else
    call initflow ()
  end if

  ! Distribute base load among all processors
  call loadBalance ()

  ! Update primitives (with Us array in all cells)
  call updatePrims ()

  ! Write initial condition to disk
  if (nextout.eq.0) then
    call writeOutput (0)
    nextout = 1
  end if

  ! Main Loop
  !do while(.false.)
  do while (time <= tfin/t_sc)

    it = it + 1

    if (verbosity > 2) then
      write(logu,'(a)') "================================================================================"
      write(logu,'(1x,a,i0)') "Starting Iteration " , it
      if (verbosity > 3) write(logu,'(1x,a)') stamp()
      write(logu,'(a)') "================================================================================"
    end if

    ! Hydro Solver
    call hydroSolver ()

    ! Update primitives in all blocks
    call updatePrims ()

    ! Update ionization fraction / chemistry
    if (eos_type == EOS_H_RATE) call updateNeutralFraction ()

    ! Radiative cooling
    call cooling ()

    ! Update AMR grid
    call admesh ()

    ! Load balance
    call loadBalance ()

    ! Data output (if scheduled)
    if (dumpout) then
      call writeOutput(nextout)
      nextout = nextout + 1
      dumpout = .false. ! until next time

    end if

    ! Report progress
    call main_report ()

    ! Everyone stop here
    call mpi_barrier (mpi_comm_world, ierr)

  end do

  ! Deallocate globals and terminate execution
  if (verbosity > 0) then
    write(logu,*) ""
    write(logu,'(a)') "================================================================================"
    write(logu,'(a)') STAMP()
    write(logu,'(a)') 'Execution complete!'
    write(logu,'(a,a)') 'Total elapsed time: ', nicetoc(start_mark)
  end if
  call deinit()

end program Walicxe3D
