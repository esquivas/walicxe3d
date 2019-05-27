!===============================================================================
!> @file uniformISM.f90
!> @brief Initial condition module: Uniform ISM
!> @author A. Esquivel
!> @date 28/Jun/2018

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

!> @brief Imposes a uniform ISM in the entire domain
!> @details This module contains code to set up initial uniform conditions
!!
!! NOTE: do not modify the code in this module! Import it in your user
!! module and utilize the provided code.
!!
!! How to use this module. To detonate a SNR:
!! 1) Import the 'uniformISM' module in your user module
!! 2) Create an object of type ism_params_type to hold the parameters of
!!    the uniform ISM you want to impose. The descriptions of the object's
!!    fields are given below. All fields are mandatory.
!! 3) Whenever you want to impose the uniform ISM (most likely at t=0) call the !!    impose_uniform_ism() subroutine,  passing as arguments the snrparams !!    object and the flow variables array.
module brio_wu

  implicit none

  ! ============================================================================

contains

  ! ============================================
  !> @brief Sets uniform background
  !> @details Fills the domain with uniform initial conditions.
  !> @params ism_params A ism_params_type object containing physical parameters.
  !> @param uvars Full flow variables array to be modified (U or UP)
  subroutine impose_brio_wu (uvars)

    use constants,  only : pi
    use globals,    only : localBlocks, logu, rank
    use parameters, only : nbMaxProc, nxmin, nxmax, nymin, nymax, nzmin, nzmax,&
                           d_sc, v_sc, p_sc, b_sc, neqtot, verbosity, logged, master
    use hydro_core, only : prim2flow
    use amr,        only : cellPos

    implicit none
    real, intent(inout) :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    real, parameter :: twopi = 2.*pi

    integer :: nb, bID, i, j, k
    real    :: dens, pres, vx, vy, vz, bx, by, bz, x, y, z

    real    :: primit(neqtot), flowvars(neqtot)

  if ( (verbosity > 0).and.(logged.or.(rank==master)) ) &
     write(logu,'(1x,a)') "> Imposing Brio Wu initial conditions..."

    ! Impose flow conditions
    do nb=1,nbMaxProc
      bID = localBlocks(nb)

      if (bID.ne.-1) then

        do k=nzmin,nzmax
          do j=nymin,nymax
            do i=nxmin,nxmax

              call cellPos(bID, i, j, k, x, y, z)

              if (x <= 0.5) then
                dens = 1.
                pres = 1.
                vx   = 0.
                vy   = 0.
                vz   = 0.
                Bx   = 0.75
                By   = 1.
                Bz   = 0.
              else
                dens = 0.125
                pres = 0.1
                vx   = 0.
                vy   = 0.
                vz   = 0.
                Bx   = 0.75
                By   = -1.
                Bz   = 0.
              end if

              ! set ansd scale primitives
                primit(1) = dens / d_sc
                primit(2) = vx   / v_sc
                primit(3) = vy   / v_sc
                primit(4) = vz   / v_sc
                primit(5) = pres / p_sc
                primit(6) = bx  ! / b_sc
                primit(7) = by  ! / b_sc
                primit(8) = bz  ! / b_sc
                !primit(9) = dens / d_sc

                ! Convert primitives and set flow vars for this cell
                call prim2flow (primit, flowvars)
                uvars(nb,:,i,j,k) = flowvars(:)
            end do
          end do
        end do

      end if
    end do
    if ( (verbosity > 0).and.(logged.or.(rank==master)) )  then
      write(logu,'(1x,a)') ""
      write(logu,'(1x,a)') "> Finished imposing Brio wu"
    end if

end subroutine impose_brio_wu

end module brio_wu
