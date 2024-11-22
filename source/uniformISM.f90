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
module uniformISM
  use parameters, only : mu0
  implicit none

  !===============================
  ! Please don't modify this section

  ! This custom type contains the data required by the detonateSNR to
  ! set off a supernova remnant.
  ! The variables have the following meanings:
  ! mu = mean atomic mass per particle (amu)
  ! metal = ISM metallicity (ignored if cooling type is not COOL_TABLE_METAL)
  ! y0   = ISM neutral fraction (used if EOS_H_RATE is enabled)
  ! dens = ISM mass density (g/cm^3)
  ! temp = ISM temperature (K)
  ! vx   = ISM velocity x-component (cm/s)
  ! vy   = ISM velocity y-component (cm/s)
  ! vz   = ISM velocity z-component (cm/s)
  ! bx   = ISM magnetic field x-component (G)
  ! by   = ISM magnetic field y-component (G)
  ! bz   = ISM magnetic field z-component (G)
  ! All physical values must be given in *CGS*.
  type ism_params_type

    real :: mu    = mu0
    real :: metal = 1.0
    real :: y0    = 0.9999  ! (neutral, only a small ion seed)
    real :: dens
    real :: temp
    real :: vx
    real :: vy
    real :: vz
#ifdef BFIELD
    real :: bx = 0.
    real :: by = 0.
    real :: bz = 0.
#endif

  end type ism_params_type

  ! ============================================================================

contains

  ! ============================================
  !> @brief Sets uniform background
  !> @details Fills the domain with uniform initial conditions.
  !> @params ism_params A ism_params_type object containing physical parameters.
  !> @param uvars Full flow variables array to be modified (U or UP)
  subroutine impose_uniform_ism (ism_params, uvars)

    use constants,  only : AMU, KB, COOL_TABLE_METAL, EOS_H_RATE
    use globals,    only : localBlocks, logu
    use parameters, only : nbMaxProc, neqtot, nxmin, nxmax, nymin, nymax, nzmin,&
                           nzmax, d_sc, v_sc, P_sc, B_sc, cooling_type, metalpas,     &
                           verbosity, firstpas, eos_type
    use hydro_core, only : prim2flow

    implicit none
    type(ism_params_type), intent(inout) :: ism_params
    real, intent(inout) :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    integer :: nb, bID, i, j, k
    real    :: dens, pres, temp, vx, vy, vz, mu, y0, metal
#ifdef BFIELD
    real    :: bx, by, bz
#endif
    real    :: primit(neqtot), flowvars(neqtot)

        if (verbosity > 0) write(logu,'(1x,a)') "> Imposing uniform background medium ..."

    !  unpack ISM paraeters
    mu    = ism_params%mu
    dens  = ism_params%dens
    temp  = ism_params%temp
    vx    = ism_params%vx
    vy    = ism_params%vy
    vz    = ism_params%vz
#ifdef BFIELD
    bx    = ism_params%bx
    by    = ism_params%by
    bz    = ism_params%bz
#endif
    metal = ism_params%metal
    y0    = ism_params%y0

    pres = dens/(mu*AMU)*KB*temp

    !call refineZone (zone, zlevel)

    ! Impose flow conditions
    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then

        do i=nxmin,nxmax
          do j=nymin,nymax
            do k=nzmin,nzmax

                ! set ansd scale primitives
                primit(1) = mu*dens/d_sc
                primit(2) = vx/v_sc
                primit(3) = vy/v_sc
                primit(4) = vz/v_sc
                primit(5) = pres/p_sc

                ! Magnetic field (if applicable)
#ifdef BFIELD
                primit(6) = bx/B_sc
                primit(7) = by/B_sc
                primit(8) = bz/B_sc
#endif
                ! Passive scalar for metalicity
                if (cooling_type.eq.COOL_TABLE_METAL) then
                  primit(metalpas) = metal*primit(1)
                end if

                !  Passive scalar for H_rate (neutral fraction)
                if (eos_type == EOS_H_RATE) then
                  primit(firstpas) = y0*dens/d_sc
                end if

                ! Convert primitives and set flow vars for this cell
                call prim2flow (primit, flowvars)
                uvars(nb,:,i,j,k) = flowvars(:)

            end do
          end do
        end do

      end if
    end do
    if (verbosity > 0) then
      write(logu,'(1x,a)') "Ambient medium parameters:"
      write(logu,'(1x,a,es12.5,a,es12.5,a)') "Dens= ", dens, " g cm^-3   (", dens/d_sc, " code units)"
      write(logu,'(1x,a,es12.5,a,es12.5,a)') "Pres= ", pres, " dyn cm^-2 (", pres/p_sc, ")"
      write(logu,'(1x,a,es12.5,a)') "Temp= ", temp, " K         (no code units)"
      write(logu,'(1x,a,es12.5,a,es12.5,a)') "Velx= ", vx, " cm s^-1   (", vx/v_sc, ")"
      write(logu,'(1x,a,es12.5,a,es12.5,a)') "Vely= ", vy, " cm s^-1   (", vy/v_sc, ")"
      write(logu,'(1x,a,es12.5,a,es12.5,a)') "Velz= ", vz, " cm s^-1   (", vz/v_sc, ")"
#ifdef BFIELD
      write(logu,'(1x,a,es12.5,a)') "Bx= ", bx, " G"
      write(logu,'(1x,a,es12.5,a)') "By= ", by, " G"
      write(logu,'(1x,a,es12.5,a)') "Bz= ", bz, " G"
#endif
  end if

end subroutine impose_uniform_ism

end module uniformISM
