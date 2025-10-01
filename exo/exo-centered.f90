!===============================================================================
!> @file exo-centered.f90
!> @brief Exoplanet module for local simulations
!> @author  M. Schneiter, C. Villarreal  D'Angelo, A. Esquivel
!> @date 25/Feb/2025

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

!> @brief Exoplanet module
!> @details Problem Module for exoplanet

module exoplanet

    use parameters
    use globals, only: RANK
    implicit none

    TYPE body
      real :: x, y, z
      real :: mass, radius, amdot
      real :: rho
      real :: temp
      real :: bfield
      real :: WR,WT,WV,WD
      real :: S0H, S0He1, S0He3
    END TYPE

    type(body) star,planet,barycenter

    real :: torb             !< planet: orbital period
    real :: rorb             !<  orbital radius
    real :: omegap           !< planet: angular velocity
    real :: rho_p            !< density at rp
    real :: cs2              !< speed of sound at rp
    real :: Kconst           !< constant tripathi eq.19
    real :: rho0, R0         !< Inner region
    real :: rho_out, Rout    !< Outer region
    real :: Rbound           !< Boundary conditions
    real :: dens_c1, dens_c2 !< constants for density computation
    real :: energy_c1        !< constant for energy computation
    real :: alpha, ymol, fHII, fHeII, fHII_s, fHeII_s

  contains

  !=======================================================================

  !> @brief Module initialization
  !> @details Here the parameters of the Star are initialized, and scaled
  !! to code units

  subroutine init_exo()

    use constants,  only : RG, MSUN, YR, GGRAV, PI, MJUP, AU, DAY, RJUP, RSUN,&
                           KB, AMU, KPS
    use parameters, only: xphystot, gamma, l_sc, v_sc, t_sc, d_sc

    implicit none

    !----------------STAR PARAMETERS (cgs) ------------------
    star%mass      = 1.132*msun !1.119*msun
    star%radius    = 1.367*rsun !1.155*rsun
  !  star%amdot     = 0.0 !1.3E-16*msun/yr    ! Stellar Mass Loss rate (g s^-1)
    star%WT        = 1.0e6 !1.35e6            ! Stellar temperature (K)
    !star%WR        = 8.16*star%radius        ! wind launching radius
    star%WV        = 0. !130.0e5              ! Stellar wind velocity (cm/s)
    star%WD        = 1.3e-21  !1.3e-19        ! Stellar wind density  (g cm^-3)
    star%bfield    = 0.0                      ! Stellar magnetic field (g)
    star%S0H       = 1.2e15*10.0              ! soft_euv lyman weinn flux
    star%S0He1     = 7.3e14*10.0              ! hard_euv (1/cm2/s)
    star%S0He3     = 2.6e18*0.1               ! soft_euv **(mid-UV)**
    !----------------PLANET PARAMETERS (cgs) ------------------
    planet%mass    = 0.685*mjup
    planet%amdot   = 1.0e9                    ! Planetary Mass Loss rate (g/s)
    planet%rho     = 5e-14                    ! divide trying to improve dt
    planet%radius  = 1.98*rjup!1.359*Rjup
    planet%temp    = 1.8370E3                 ! Planets temperature
  !  planet%WT      = 0.0                     ! Planets wind temperature
    planet%WR      = 3*rjup                   ! Planetary wind radius (cm)
    planet%WV      = 15*KPS                   ! Planets wind velocity (cm/s)
    !planet%WD      = 0.0                     ! Planetary wind density
    planet%bfield  = 0.0                      ! Planetary magnetic field (g)

    !ORBITAL PARAMETERS (cgs)
    rorb = 0.03397*AU
    torb = 2.15000820*day

    !  star wind launching radius viewed from star
    star%WR = rorb-xphystot/2.0

    !  planet density at launching radius
    planet%WD = ((planet%amdot/planet%WR)/(4.0*PI*planet%WR*planet%WV) )

    ! Below this line parametes are scaled in or to  code units
    rorb = rorb/l_sc
    torb = torb/t_sc
    omegap = 2.0*pi/torb

    star%x = 0.0
    star%y = 0.0
    star%z = -rorb

    planet%x = 0.00
    planet%y = 0.00
    planet%z = 0.00

    barycenter%x = 0.0
    barycenter%y = 0.0
    barycenter%z = -rorb*star%mass/(star%mass + planet%mass)

    star%WV       = star%WV     / v_sc
    star%radius   = star%radius / l_sc
    star%WD       = star%WD     / d_sc

    planet%radius = planet%radius / l_sc
    planet%WR     = planet%WR     / l_sc
    planet%WV     = planet%WV     / v_sc
    planet%WD     = planet%WD     / d_sc
    planet%rho    = planet%rho    / d_sc
    planet%amdot  = planet%amdot  / d_sc / v_sc / l_sc / l_sc

  end subroutine init_exo

  !=======================================================================

  !> @brief Inject sources of wind
  !> @details Imposes the sources of wond from the star and planet
  !> @param real [out] u(neq,nxmin:nxmax,nymin:nymax,nzmin:nzmax) :
  !! conserver variables
  !> @param real [time] time : current integration timr
    !--------------------------------------------------------------------
   subroutine impose_initial_exo(uvars)
    use parameters, only : verbosity, xphystot, yphystot, zphystot, l_sc, p_sc,&
                           nxmin, nxmax, nymin, nymax, nzmin, nzmax, neqtot,   &
                           TauRT, inh0
    use constants, only : pi, AMU, Kb
    use globals, only : logu, maxlev, dx, dy, dz, localBlocks
    use amr,     only : refineZone, cellPos
    use hydro_core, only : prim2flow
    implicit none
    real, intent(inout) :: uvars(nbMaxProc, neqtot,  &
                                 nxmin:nxmax, nymin:nymax, nzmin:nzmax)
    real :: zone(6)
    real :: primit(neqtot)
    real :: xc, yc, zc, x, y, z
    real :: rp, xp, yp, zp, xs, ys, zs, rs
    real :: velx, vely, velz, dens, densH0, pressure

    integer :: nb, bID, i, j, k

    if (verbosity > 0) then
      write(logu,*) ""
      write(logu,'(1x,a)') "> Imposing planet (initial) ..."
    end if

    !  planet position (code units), respect to grid center
    xc = planet%x
    yc = planet%y
    zc = planet%z

    if (verbosity > 0) write(logu,*) " Refining zone around planet to level ", &
                                     maxlev
    zone(1) = xc*l_sc + xphystot/2.0 - 2.0*planet%radius*l_sc
    zone(2) = xc*l_sc + xphystot/2.0 + 2.0*planet%radius*l_sc
    zone(3) = yc*l_sc + yphystot/2.0 - 2.0*planet%radius*l_sc
    zone(4) = yc*l_sc + yphystot/2.0 + 2.0*planet%radius*l_sc
    zone(5) = zc*l_sc + zphystot/2.0 - 2.0*planet%radius*l_sc
    zone(6) = zc*l_sc + zphystot/2.0 + 2.0*planet%radius*l_sc
    call refineZone (zone, maxlev)

    write(logu,*) " About to set ICs in all the blocks"

    ! Impose flow conditions, where applicable
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then

        do i=nxmin,nxmax
          do j=nymin,nymax
            do k=nzmin,nzmax

              !  obtain cell position
              call cellPos (bID, i, j, k, x, y, z, center=.true.)
              !  distance from planet center
              rp = sqrt( (x-xc)**2 + (y-yc)**2 + (z-zc)**2 )

              if(rp <= planet%radius ) then ! inside planet

                dens = planet%rho
                velx = 0.0
                vely = 0.0
                velz = 0.0

                !  dens below is already scaled, thus is n instead of rho
                !  (d_sc = AMU)
                pressure = KB * dens * planet%temp / p_sc
                densH0   = 0.9999*dens

              elseif(rp <= planet%WR) then

                dens = planet%amdot /planet%WV /rp /rp /4.0 /PI
                velx = planet%WV * (x-xc) / rp ! + orbital speed?
                vely = planet%WV * (y-yc) / rp ! + orbital speed?
                velz = planet%WV * (z-zc) / rp ! + orbital speed?

                pressure  = KB * dens * planet%temp / p_sc
                densH0    = 0.9999*dens

              else

                !Stellar Wind
                xp = x - barycenter%x
                yp = y - barycenter%y
                zp = z - barycenter%z

                xs = x - star%x
                ys = y - star%y
                zs = z - star%z
                rs = sqrt(xs**2+ys**2+zs**2)

                velx = star%WV * xs/rs - omegap*zp
                vely = star%WV * ys/rs
                velz = star%WV * zs/rs + omegap*xp   !code units

                dens = star%WD
                pressure = KB * dens * star%WT / p_sc
                densH0   = 1e-4*dens

              end if

              primit(1) = dens
              primit(2) = velx
              primit(3) = vely
              primit(4) = velz
              primit(5) = pressure
              primit(inh0)   = densH0
              primit(TauRT) = 0.0

              ! Convert primitives and set flow vars for this cell
              call prim2flow( primit, uvars(nb,:,i,j,k) )
              !endif

            end do   ! k
          end do     ! j
        end do       ! i

      end if         ! bID
    end do           ! nb

  end subroutine impose_initial_exo

!=======================================================================

!   subroutine impose_exo(u,time)

  subroutine impose_exo(uvars, time)

    use parameters, only : verbosity, xphystot, yphystot, zphystot, l_sc, p_sc,&
                           nxmin, nxmax, nymin, nymax, nzmin, nzmax, neqtot,   &
                           TauRT, inh0
    use constants, only : pi, AMU, Kb
    use globals, only : logu, maxlev, dx, dy, dz, localBlocks
    use amr,     only : refineZone, cellPos
    use hydro_core, only : prim2flow
    implicit none
    real, intent(inout) :: uvars(nbMaxProc, neqtot,  &
                                 nxmin:nxmax, nymin:nymax, nzmin:nzmax)
    real, intent(in)    :: time
    real :: primit(neqtot)
    real :: xc, yc, zc, x, y, z
    real :: rp, velx, vely, velz, dens, densH0, pressure
    integer :: nb, bID, i, j, k
    logical :: set_cell

    !  planet position (code units), respect to grid center
    xc = planet%x
    yc = planet%y
    zc = planet%z

    ! Impose flow conditions, where applicable
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then

        do i=nxmin,nxmax
          do j=nymin,nymax
            do k=nzmin,nzmax

              set_cell = .false.

              !  obtain cell position
              call cellPos (bID, i, j, k, x, y, z, center=.true.)
              !  distance from planet center
              rp = sqrt( (x-xc)**2 + (y-yc)**2 + (z-zc)**2 )
              !  Add some Softening
              rp = max(rp, 0.5*dx(maxlev))

              if(rp <= planet%radius ) then ! inside planet

                dens = planet%rho
                velx = 0.0
                vely = 0.0
                velz = 0.0

                !  dens below is already scaled, thus is n instead of rho
                !  (d_sc = AMU)
                pressure = KB * dens * planet%temp / p_sc
                densH0   = 0.9999*dens

                set_cell = .true.

              elseif(rp <= planet%WR) then

                dens = planet%amdot /planet%WV /rp /rp /4.0 /PI
                velx = planet%WV * (x-xc) / rp ! + orbital speed?
                vely = planet%WV * (y-yc) / rp ! + orbital speed?
                velz = planet%WV * (z-zc) / rp ! + orbital speed?

                pressure  = KB * dens * planet%temp / p_sc
                densH0    = 0.9990*dens

                set_cell = .true.

              end if

              if (set_cell) then
                primit(1) = dens
                primit(2) = velx
                primit(3) = vely
                primit(4) = velz
                primit(5)     = pressure
                primit(inH0)  = densH0
                primit(TauRT) = 0.0

                ! Convert primitives and set flow vars for this cell
                call prim2flow( primit, uvars(nb,:,i,j,k) )

              endif

            end do   ! k
          end do     ! j
        end do       ! i

      end if         ! bID
    end do           ! nb

  end subroutine impose_exo

!=======================================================================

  end module exoplanet

  !=======================================================================
