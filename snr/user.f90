!===============================================================================
!> @file user.f90
!> @brief User-specified initial and boundary conditions
!> @author Juan C. Toledo
!> @date 20/May/2013

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

! EXAMPLE FILE: snr

! This example shows how to impose supernova remnants. We simulate a 10^3 pc
! box and impose two snr at the start (one with a "simple" model and the
! other with a more realistic Type Ia ejecta model). At t=500 yr we detonate
! a third snr in the center. We then let the simulation evole until 1 kyr.
! This is done using an equivalent max-level resolution of 128^3 and using
! four parallel processes.

! In the setInitialConditions subroutine we use imposeSNR and imposeSNRIa
! subroutines to detonate the initial remnants. These subroutines must be
! passed two arguments: a snr_params_type data object which contains the
! parameters of the snr, and the flow variables vector (uvars), which are
! received by the user IC routine from its caller.

! Note that the snr module is imported at the top of the file.

! The snr_params data object contains several fields that must be filled
! for the snr to be properly defined:
!  xc: x-coordinate of the center of the remnant (cm)
!  yc: y-coordinate of the center of the remnant (cm)
!  zc: z-coordinate of the center of the remnant (cm)
!  radius: the radius of the remnant (cm)
!  mass: the total mass of the ejecta (g)
!  energy: the total energy (kinetic+thermal) of the explosion (erg)
!  chi: the fraction of the total energy deposited as kinetic energy (0.0-1.0)
!  bx, by, bz: magnetic field components inside the remnant (G)
!  time: the time at which the snr is to be detonated. This is only for the
!   user to use if a snr is to be detonated at a later point in time.
!  armed: a logical flag that is set to .false. after the snr has been
!   detonated. This can be used by the user to prevent multiple detonations
!   of the same remnant.

! We also use the userBoundary subroutine, which is called at the end of the
! boundary conditions at each timestep, to detonate a SNR at a poin in time
! after the start of the simulation. Note how the time and armed fields of
! the snr parameters data object are used.

module userconds
! ============================================
!< @brief User-specified initial and boundary conditions
!< @details This is the file where the user sets the initial conditions of
!! the simulation, as well as any special boundary conditions or custom
!! modifications to the simulation.
!!
!! How to use this module:
!!   1) Add "use" statements to include any modules you require in the
!!      section marked [1].
!!   2) You can define aditional module-wide global variables and parameters
!!      in the section marked [2].
!!   3) Fill in the subroutine setInitialCondition(), marked [3], which
!!      is called at the beginning of the simulation.
!!   4) Optionally, fill in the subroutine userBoundary(), marked [4],
!!      which is called at the end of each boundary exchange operation.
!!   5) Optionally, fill the subroutine get_user_source_terms, marked [5],
!!      which is called at the end of the upwind step. (in each half-timestep).
!!
!!
!! All subroutines in this module automatically have access to the global
!! parameters and variables.
! ============================================

  use parameters
  use globals
  ! ============================================
  ! [1] Add HERE any aditional modules required by your simulation

  use uniformISM
  use snr

  ! ============================================
  implicit none

  ! ============================================
  ! [2] Define HERE any additional parameters or variables required
  ! by the user subroutines below if they they are not provided by
  ! an external module

  ! We can declare the snr parameters objects here. We'll fill them out
  ! during the initial conditions.
  type(ism_params_type) :: ism

  type(snr_params_type) :: snr1
  !type(snr_params_type) :: snr2
  !type(snr_params_type) :: snr3

  ! ============================================

contains

  subroutine setInitialCondition (uvars)
  ! ============================================
  ! [3] USER-DEFINED INITIAL CONDITIONS
  !
  !< @brief User-defined Initial Conditions
  !< @details This subroutine is called at the beginning of the simulation,
  !! after the base grid is built and a basic uniform initial condition
  !! is imposed. It is to be modified by the user to define the problem-
  !! specific Initial Condition.
  !!
  !! IMPORTANT: This subroutine receives the FLOW variables array to be
  !! modified as argument 'uvars'. The subroutine must modify this array,
  !! *NOT* the global arrays U, UP or PRIMS.
  !!
  !! The array has the following structure:
  !!   uvars ( block ID, equation number, cell_i, cell_j, cell_k )
  !! where the equation numbering is:
  !!   1: rho
  !!   2: rho*u
  !!   3: rho*v
  !!   4: rho*w
  !!   5: E (kinetic+thermal)
  !! If the passive magnetic field is enabled:
  !!   6: B_x
  !!   7: B_y
  !!   8: B_z
  !! If passive scalars are enabled, they begin after all other flow
  !! variables (i.e., 9+ if passive B-field enabled, 6+ if not).
  !!
  !! Note that the cell indices include ghost cells. For instance, if we
  !! had a 1-deep ghost cell layer then cells from 1 to ncells_x would be
  !! physical cells, while cells 0 and ncells_x+1  would be ghost cells.
  ! ============================================

    implicit none
    real, intent(inout) :: uvars (nbMaxProc, neqtot, &
                           nxmin:nxmax, nymin:nymax, nzmin:nzmax)

    ! ============================================

    !  uniform ISM parameters
    ism%mu   = mu0
    ism%dens = 1.0 * AMU * ism%mu
    ism%temp = 1.e3
    ism%y0   = 0.9999  ! (neutral, only a small ion seed)
    ism%vx   = 0.
    ism%vy   = 0.
    ism%vz   = 0.
#ifdef BFIELD
    ism%bx   = 0.
    ism%by   = 0.
    ism%bz   = 0.
#endif

    ! Parameters of the first SNR
    snr1%xc      = 25.0*PC
    snr1%yc      = 25.0*PC
    snr1%zc      = 25.0*PC
    snr1%radius  = 1.0*PC  ! originally 1pc
    snr1%mass    = 4*MSUN
    snr1%energy  = 1E51
    snr1%chi     = 0.5
    snr1%rho_env = ism%dens
    snr1%time    = 0.0
    snr1%y0      = 0.0   ! (fully ionized)

    ! fill the domain with a uniform ism
    call impose_uniform_ism(ism,uvars)

    ! The SN is detonated as initial conditions
    call detonateSNR(snr1, uvars)

    ! ============================================

  end subroutine setInitialCondition

  !=============================================================================

  subroutine userBoundary (uvars)
  ! ============================================
  ! [4] USER-DEFINED BOUNDARY CONDITIONS
  !
  !< @brief User-defined Boundary Conditions
  !< @details This subroutine is called once per timestep *after* standard
  !! boundary have been applied to all blocks. It allows the user to
  !! to impose an arbitrary boundary condition on the simulation.
  !!
  !! IMPORTANT: This subroutine receives the FLOW variables array to be
  !! modified as argument 'uvars'. The subroutine must modify this array,
  !! *NOT* the global arrays U and UP.
  !!
  !! The structure of this array is described in the setInitialConditions()
  !! subroutine documentation above.
  ! ============================================

    implicit none
    real, intent(inout) :: uvars (nbMaxProc, neqtot, &
                           nxmin:nxmax, nymin:nymax, nzmin:nzmax)
    ! ============================================

    ! We check both for time and armed status, and detonate when
    ! both become true. The armed property is set to false automatically
    ! by the subroutine (but can be reset by the user, if he so wishes).
    ! Also note how we de-scale the time variable to physical units.
    !if ((time*t_sc.ge.snr3%time).and.(snr3%armed)) then
    !  call detonateSNRIa(snr3,uvars)
    !end if


    ! ============================================

  end subroutine userBoundary

  !=============================================================================

  subroutine get_user_source_terms (pp, s, i, j , k)
    ! ============================================
    ! [5] USER-DEFINED SOURCE TERMS
    !
    !> @brief User-defined Source Terms
    !> @details This subroutine is called once per half timestep at the end
    !! of the upwind timestep. This allows ther user to add custom S terms,
    !! of the form:  dU/dt+dF/dx+dG/dy+dH/dz=S
    !! useful for instance to include gravity, tidal and/or inertial forces.
    !> @param real [in]   pp(neq) : vector of primitive variables
    !> @param real [inout] s(neq) : vector with source terms, has to add to
    !>  whatever is there, as other modules can add their own sources
    !> @param integer [in] i : cell index in the X direction
    !> @param integer [in] j : cell index in the Y direction
    !> @param integer [in] k : cell index in the Z direction
    !use constants,  only : Ggrav
    use parameters, only : neqtot
    !use globals,    only : dx, dy, dz, coords
    !use exoplanet
    !use radpress
    implicit none
    real,    intent(in   ) :: pp(neqtot)
    real,    intent(inout) :: s (neqtot)
    integer, intent(in   ) :: i, j, k

    ! integer, parameter  :: nb=3
    ! real    :: x(nb),y(nb),z(nb), GM(nb), rad2(nb)
    ! integer :: index
    ! real    :: xc ,yc, zc
    ! real    :: GradPhi(nb), OmegaSq
    ! real    :: rsoft
    ! real    :: v, fracv, frac_neutro

    ! rsoft = (dx*0.1)**2

    ! GM(2) = Ggrav*Star%mass/rsc/vsc2
    ! GM(1) = Ggrav*Planet%mass/rsc/vsc2

    ! !   get cell position
    ! xc = (float(i + coords(0)*nx - nxtot/2) - 0.5)*dx
    ! yc = (float(j + coords(1)*ny - nytot/2) - 0.5)*dy
    ! zc = (float(k + coords(2)*nz - nztot/2) - 0.5)*dz

    ! ! calculate distance from the sources
    ! ! planet
    ! x(1) = xc - Planet%x
    ! y(1) = yc - Planet%y
    ! z(1) = zc - Planet%z
    ! rad2(1) = x(1)**2 + y(1)**2 + z(1)**2

    ! if(rad2(1) < rsoft)then
    !   rad2(1) = rsoft
    ! endif

    ! ! star
    ! x(2) = xc - Star%x
    ! y(2) = yc - Star%y
    ! z(2) = zc - Star%z
    ! rad2(2) = x(2)**2 + y(2)**2 + z(2)**2

    ! if(rad2(2) < rsoft)then
    !   rad2(2) = rsoft
    ! endif

    ! ! barycenter
    ! x(3) = xc - Barycenter%x
    ! y(3) = 0.0 ! porque queremos la distancia en el plano orbital x,z
    ! z(3) = zc - Barycenter%z
    ! rad2(3) = x(3)**2 + y(3)**2 + z(3)**2
    ! if(rad2(3) < rsoft)then
    !   rad2(3) = rsoft
    ! endif

    ! if ( beta_pressure ) then
    !     beta(i,j,k) = 0.
    !     !  do only outside the planet
    !     if( rad2(1) >= planet%radius**2 ) then

    !       !Each cell feels a given pressure proportional to the neutrals fraction
    !       frac_neutro = pp(neqdyn+1)/pp(1)
    !       !  Radial velocity in km s^-1 at the stellar frame
    !       v = (((pp(2)+omegap*rorb)*x(2) + pp(3)*y(2) + pp(4)*z(2))/sqrt(rad2(2)))* (vsc)

    !       fracv = (v-vr(1))/(vr(Nr)-vr(1))*Nr
    !       index = int(fracv)+1

    !       if (index < 1) then
    !         index = 1
    !       else if ( index > Nr-1 ) then
    !         index = Nr-1
    !       end if
    !       !Linear interpolation for Beta
    !       Beta(i,j,k) = (Br(index) + (v-vr(index))*(Br(index+1)-Br(index))/(vr(index+1)-vr(index)))*frac_neutro
    !       !Update scale factor GM
    !       GM(2)=GM(2)*(1.-Beta(i,j,k))

    !     end if
    !   endif


    ! OmegaSq =  ( GM(2) + GM(1) )/rorb**3

    ! GradPhi(1) = GM(1)*x(1)/rad2(1)**1.5 + GM(2)*x(2)/rad2(2)**1.5 - OmegaSq*x(3)
    ! GradPhi(2) = GM(1)*y(1)/rad2(1)**1.5 + GM(2)*y(2)/rad2(2)**1.5 - OmegaSq*y(3)
    ! GradPhi(3) = GM(1)*z(1)/rad2(1)**1.5 + GM(2)*z(2)/rad2(2)**1.5 - OmegaSq*z(3)

    ! !  update source terms with gravity
    ! s(2)= s(2) - pp(1)*GradPhi(1) - 2*pp(1) * sqrt(OmegaSq) * pp(4)
    ! s(3)= s(3) - pp(1)*GradPhi(2) ! 0
    ! s(4)= s(4) - pp(1)*GradPhi(3) + 2*pp(1) * sqrt(OmegaSq) * pp(2)

    ! ! energy
    ! s(5)= s(5) - pp(1)*(pp(2)*GradPhi(1) + pp(3)*GradPhi(2) + pp(4)*GradPhi(3))

    ! ============================================

  end subroutine get_user_source_terms

end module userconds
