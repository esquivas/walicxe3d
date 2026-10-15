!===============================================================================
!> @file field_loop.f90
!> @brief Initial condition module: advection of a magnetic field loop
!> @author A. Esquivel
!> @date 12/Oct/2026

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

!> @brief Field loop advection test (Gardiner & Stone 2005; Stone et al. 2008)
!> @details A weak circular field loop (beta ~ 1e6) is advected by a uniform
!! flow across a periodic box [-1,1]x[-0.5,0.5]. With v = (2,1) the loop
!! returns to its initial position at t = 1, 2, ... The field is set from
!! the vector potential Az = A0 (R - r) for r < R (0 outside), differenced
!! along the face edges (flux-CT), so div(B) = 0 to roundoff initially.
!! All values are in code units (magnetic pressure = B^2/2). With vz0 /= 0
!! (3D variant) Bz must remain zero to roundoff.
module field_loop

  implicit none

  real, parameter :: rho0 = 1.0       !< Density
  real, parameter :: p0   = 1.0       !< Pressure
  real, parameter :: vx0  = 2.0       !< Advection velocity, x
  real, parameter :: vy0  = 1.0       !< Advection velocity, y
  real, parameter :: vz0  = 0.0       !< Advection velocity, z (3D variant)
  real, parameter :: A0   = 1.0e-3    !< Vector potential amplitude (= |B|)
  real, parameter :: R0   = 0.3       !< Loop radius
  real, parameter :: xc   = 0.0       !< Loop centre, x (box centre = 0)
  real, parameter :: yc   = 0.0       !< Loop centre, y

  !> Static refinement at t=0: zone (physical units, origin at the box corner:
  !! left, right, front, back, bottom, top) refined to level smr_level
  logical, parameter :: smr          = .true.
  logical, parameter :: smr_reimpose = .true.
  integer, parameter :: smr_level = 3
  real,    parameter :: zone_smr(6) = [0.7, 1.3, 0.2, 0.8, 0.0, 0.5]



contains

  !=============================================================================

  !> @brief Imposes the field loop initial condition
  !> @param uvars Full flow variables array to be modified (U or UP)
  subroutine set_field_loop (uvars)

    use globals,    only : localBlocks, logu, BF, dx, dy
    use parameters, only : nbMaxProc, nxmin, nxmax, nymin, nymax, nzmin, nzmax,&
                           neqtot, firstpas, npassive, enable_flux_ct
    use hydro_core, only : prim2flow
    use amr,        only : cellPos, meshlevel

    implicit none
    real, intent(inout) :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    integer :: nb, bID, i, j, k, lev
    real    :: x, y, z, r, hx, hy, bxm, bxp, bym, byp, bx, by
    real    :: primit(neqtot), flowvars(neqtot)

    write(logu,'(1x,a)') "> Imposing field loop initial conditions..."

    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then

        call meshlevel (bID, lev)
        hx = 0.5*dx(lev)
        hy = 0.5*dy(lev)

        do k=nzmin,nzmax
          do j=nymin,nymax
            do i=nxmin,nxmax

              ! Cell centre with the origin at the centre of the box
              call cellPos (bID, i, j, k, x, y, z, center=.true.)
              x = x - xc
              y = y - yc
              r = sqrt(x**2 + y**2)

              if (enable_flux_ct) then
                ! Face-centred B from differences of Az along the face edges
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
                ! Point values of B = curl(Az z)
                if (r < R0 .and. r > 0.0) then
                  bx = -A0*y/r
                  by =  A0*x/r
                else
                  bx = 0.0
                  by = 0.0
                end if
              end if

              primit(1) = rho0
              primit(2) = vx0
              primit(3) = vy0
              primit(4) = vz0
              primit(5) = p0
              primit(6) = bx
              primit(7) = by
              primit(8) = 0.0
              ! Passive scalar marks the material initially inside the loop
              if (npassive >= 1) primit(firstpas:neqtot) = merge(1.0, 0.0, r < R0)

              call prim2flow (primit, flowvars)
              uvars(nb,:,i,j,k) = flowvars(:)

            end do
          end do
        end do

      end if
    end do

    write(logu,'(1x,a)') "> Finished imposing the field loop initial conditions."

  contains

    !> Vector potential Az = A0 (R0 - r) for r < R0, zero outside
    real function az (xx, yy)
      real, intent(in) :: xx, yy
      az = A0*max(R0 - sqrt(xx**2 + yy**2), 0.0)
    end function az

  end subroutine set_field_loop

  !=============================================================================

  !> @brief Imposes the field loop IC, with optional static refinement
  !> @details refineZone needs valid U (refineBlock interpolates it), so the
  !! IC is imposed on the root mesh first, the zone refined, and the IC
  !! imposed again analytically on the refined mesh.
  subroutine impose_field_loop (uvars)

    use parameters, only : nbMaxProc, nxmin, nxmax, nymin, nymax, nzmin, nzmax, neqtot
    use amr,        only : refineZone
    implicit none
    real, intent(inout) :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    call set_field_loop (uvars)
    if (smr) then
      call refineZone (zone_smr, smr_level)
      if (smr_reimpose) call set_field_loop (uvars)
    end if

  end subroutine impose_field_loop

  !=============================================================================

end module field_loop
