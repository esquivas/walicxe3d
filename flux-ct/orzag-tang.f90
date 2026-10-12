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
!! 3) Whenever you want to impose the uniform ISM (most likely at t=0) call the
!!    impose_uniform_ism() subroutine,  passing as arguments the snrparams
!!    object and the flow variables array.
module orzag_tang

  implicit none

  !=============================================================================

contains

  ! ============================================
  !> @brief Sets uniform background
  !> @details Fills the domain with uniform initial conditions.
  !> @params ism_params A ism_params_type object containing physical parameters.
  !> @param uvars Full flow variables array to be modified (U or UP)
  subroutine impose_orzag_tang (uvars)

    use constants,  only : pi
    use globals,    only : localBlocks, logu, BF, dx, dy
    use parameters, only : nbMaxProc, nxmin, nxmax, nymin, nymax, nzmin, nzmax,&
                           d_sc, v_sc, p_sc, b_sc, neqtot, firstpas,           &
                           enable_flux_ct
    use hydro_core, only : prim2flow
    use amr,        only : cellPos, meshlevel

    implicit none
    real, intent(inout) :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    real, parameter :: twopi = 2.*pi
    integer :: nb, bID, i, j, k, lev
    real    :: dens, pres, vx, vy, vz, bx, by, bz, x, y, z
    real    :: hx, hy, bxm, bxp, bym, byp
    real    :: primit(neqtot), flowvars(neqtot)

    write(logu,'(1x,a)') "> Imposing Orzag-Tang initial conditions..."

    dens  = 25.0/(36.0*pi)
    pres =   5.0/(12.0*pi)
    vz    = 0.0
    bz    = 0.0 / b_sc

    !call refineZone (zone, zlevel)

    ! Impose flow conditions
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then

        do i=nxmin,nxmax
          do j=nymin,nymax
            do k=nzmin,nzmax

              call cellPos(bID, i, j, k, x, y, z)
              vx = -sin( y*twopi )
              vy =  sin( x*twopi )

              if (enable_flux_ct) then
                ! Face-centred B from differences of Az along the face edges
                call meshlevel (bID, lev)
                hx = 0.5*dx(lev)
                hy = 0.5*dy(lev)
                bxp =  ( az(x+hx,y+hy) - az(x+hx,y-hy) ) / dy(lev)   ! x-face i+1/2
                bxm =  ( az(x-hx,y+hy) - az(x-hx,y-hy) ) / dy(lev)   ! x-face i-1/2
                byp = -( az(x+hx,y+hy) - az(x-hx,y+hy) ) / dx(lev)   ! y-face j+1/2
                bym = -( az(x+hx,y-hy) - az(x-hx,y-hy) ) / dx(lev)   ! y-face j-1/2
                BF(nb,1,i,j,k) = bxp
                BF(nb,2,i,j,k) = byp
                BF(nb,3,i,j,k) = 0.0
                ! Cell-centred field = average of the two faces
                bx = 0.5*(bxm + bxp)
                by = 0.5*(bym + byp)
              else
                bx = -sin( y*twopi  ) / b_sc
                by =  sin(2.*x*twopi) / b_sc
              end if

              ! set ansd scale primitives
                primit(   1    ) = dens / d_sc
                primit(   2    ) = vx   / v_sc
                primit(   3    ) = vy   / v_sc
                primit(   4    ) = vz   / v_sc
                primit(   5    ) = pres / p_sc
                primit(   6    ) = bx  ! already scaled
                primit(   7    ) = by  ! already scaled
                primit(   8    ) = bz  ! already scaled
                primit(firstpas) = dens / d_sc

                ! Convert primitives and set flow vars for this cell
                call prim2flow (primit, flowvars)
                uvars(nb,:,i,j,k) = flowvars(:)

            end do
          end do
        end do

      end if
    end do

    write(logu,'(1x,a)') ""
    write(logu,'(1x,a)') "> Finished imposing the Orzag-Tang initial conditions."

  !=============================================================================
  ! Vector potential
  contains

    !> Vector potential (code units), B0 = 1/b_sc
    real function az (xx, yy)
      real, intent(in) :: xx, yy
      az = ( cos(2.*twopi*xx)/(2.*twopi) + cos(twopi*yy)/twopi ) / b_sc
    end function az

!===============================================================================

end subroutine impose_orzag_tang

end module orzag_tang
