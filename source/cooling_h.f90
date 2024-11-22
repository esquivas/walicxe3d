!===============================================================================
!> @file cooling_H.f90
!> @brief Radiative cooling with Biro et al. 95 prescription
!> @author C. Villarreal & A. Esquivel
!> @date 7/Nov/2024

! Copyright (c) 2024 Juan C. Toledo and Alejandro Esquivel
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

!> @brief Cooling with parametrized cooling and H rate equation
!> @details Cooling with parametrized cooling and H rate equation
!> as prescribed by Biro et al. 1999

module cooling_H

!#ifdef PASSIVES

  implicit none

contains

  !======================================================================
  !> @brief betaH(T)
  !> @details @f$ \beta_H(T) @f$
  !> @param real 8[in] T : Temperature K
  function betah(T)

    implicit none

    real (kind=8) ::  betah
    real (kind=8), intent(in) ::  T
    real (kind=8)             ::  a

    a=157890./T
    betah=1.133D-24/sqrt(a)*(-0.0713+0.5*log(a)+0.640*a**(-0.33333))

  end function betah

  !======================================================================
  !> @brief Non equilibrium cooling
  !> @details   Non-equilibrium energy loss for low temperatures
  !>     considering the collisional excitation of [O I] and
  !>   [O II] lines and radiative recombination of H. This
  !>   cooling rate is multiplied by a factor of 7.033 so
  !>   that it has the same value as the "coronal equilibrium"
  !>   cooling rate at a temperature of 44770 K (at temperatures
  !>   higher than this value, the equilibrium cooling rate is
  !>   used). The collisional ionization of H and excitation
  !>   of Lyman-alpha are computed separately, and added to
  !>   the cooling rate.
  !> @param real8 [in] x1  : initial H ionization fraction
  !> @param real8 [in] x2  : final H ionization fraction
  !> @param real [in] dt  : timestep
  !> @param real8 [in] den : total density of hydrogen
  !> @param real8 [in] dh0 : density of neutral hydrogen
  !> @param real8 [in] Te0 : Temperature
  function aloss(x1,x2,den,dH0,Te0)

    implicit none

    real (kind=8)             :: aloss
!    real, intent(in)          :: dt
    real (kind=8), intent(in) :: x1,x2,den,dH0,Te0
    real, parameter :: XION=2.179e-11, XH=0.9,XO=1.e-3
    real, parameter :: C0=0.5732,C1=1.8288e-5,C2=-1.15822e-10,C3=9.4288e-16
    real, parameter :: D0=0.5856,D1=1.55083e-5,D2=-9.669e-12, D3=5.716e-19
    real, parameter :: ENK=118409.,EN=1.634E-11
    real (kind=8)   :: Te, dHp,de,dOI,dOII,omega,omegaL,omegaH,frac,qla
    real (kind=8)   :: ecoll,cion,eion,erec,Tm,Tm2,eOI,eOII,equil,fr,ex2,tanh
    real (kind=8)   :: betaf, HIIcool

    Te  = max(Te0,10.)
    dHp = (1.-X1)*den    ! hydrogen density
    de  = dHp+1.E-4*den  ! electron density
    dOI = XO*dH0         ! oxigen density
    dOII= XO*dHp         ! ionized oxigen density

    if(Te <= 1e4 ) then
      aloss = 1e-30
      return
    end if

    !   Collisionally excited Lyman alpha
    if(Te .le. 55000.) omega = C0 + Te*(C1 + Te*(C2 + Te*C3))
    if(Te .ge. 72000.) omega = D0 + Te*(D1 + Te*(D2 + Te*D3))
    if(Te .gt. 55000. .and. Te .lt. 72000.) then
      omegaL = C0 + Te*(C1 + Te*(C2 + Te*C3))
      omegaH = D0 + Te*(D1 + Te*(D2 + Te*D3))
      frac   = (Te-55000.)/17000.
      omega  = (1.-frac)*omegaL + frac*omegaH
    end if
    qla   = 8.6287E-6/(2.*sqrt(Te))*omega*exp(-ENK/Te)
    ecoll = de*dH0*qla*EN
    ecoll = max(ecoll,0.)

    !   Hydrogen recombination and collisional ionization
    cion = 5.834E-11*sqrt(Te)*exp(-1.579E5/Te)
    eion = de*dH0*cion*XION
    erec = de*dHp*(betah(Te))
    erec = max(erec,0.)

    !   [O I] and [O II] coll. excited lines
    Tm   = 1./Te
    Tm2  = Tm*Tm
    eOI  = de*dOI*10.**(1381465*Tm2-12328.69*Tm-19.82621)
    eOII = de*dOII*10.**(-2061075.*Tm2-14596.24*Tm-19.01402)
    eOI  = max(eOI,0.)
    eOII = max(eOII,0.)

    !   free-free cooling
    betaf  = 1.3*1.42E-27*Te**0.5
    HIIcool= de*dHp*betaf

    !   Equilibrium cooling (for high Te)
    equil = (1.0455E-18/Te**0.63)*(1. - exp(-(Te*1.E-5)**1.63))*de*den + HIIcool

    !   switch between non-equil. and equil. ionization cooling
    if(Te .le. 44770.) fr = 0.
    if(Te .ge. 54770.) fr = 1.
    if(Te .gt. 44770. .and. Te .lt. 54770.) then
      ex2  = exp(-2.*(Te-49770.)/500.)
      tanh = (1. - ex2)/(1. + ex2)
      fr   = 0.5*(1. + tanh)
    end if

    aloss = ecoll + eion + (erec + 7.033*(eOI + eOII))*(1.-fr) + equil*fr
    !  aloss in cgs (cm^-3 s^-1)

  end function aloss

  !=======================================================================
  !> @brief
  !> @details
  !> @param real [in] uu(neq) : primitive variablas in one cell
  !> @param real [in] uu(neq) : conserved variablas in one cell
  !> @param real [in] dt      : timestep (seconds)
  !# !> @param real [in] radphi  : photoionizing rate
  subroutine  apply_cooling_h_neq(bIndx, maxloss)

    use parameters   !#add energy per ionzation as a parameter?
    use globals,    only : U, PRIM, dt
    use hydro_core, only : calcTemp
    use constants
    implicit none

    integer, intent(in):: bIndx
    real, intent(out)  :: maxloss
    real, parameter    :: T_floor = 10.0, Tmin_cool = 1.e4

    real(kind = 8)     :: y0, y1, dh, dh0, gain, tprime, L0, ce, Temp, T1
    real               :: vel2, ETH, EK, cool_factor, dt_seconds
    real               :: frac_loss!, metal
    integer            :: i, j, k

    maxloss = 0.0

    dt_seconds = dt*t_sc

    do i=1,ncells_x
      do j=1,ncells_y
        do k=1,ncells_z

          !# H total density n_HT
          dh  = real( PRIM(bIndx,1,i,j,k) / mu0 ,  8)

          !# neutrals density
          dh0 = real( PRIM(bIndx,firstpas,i,j,k) , 8)

          !# neutral H fraction (t0)
          y0 =  real( mu0 *  U(bIndx,firstpas,i,j,k) /    U(bIndx,1,i,j,k), 8)

          !# neutral H fraction (t0+dt) fraccion actualizada
          y1  = real( mu0*PRIM(bIndx,firstpas,i,j,k) / PRIM(bIndx,1,i,j,k) , 8)

          !# update the neutral fraction for the conserved arrays Us
          U (bIndx,firstpas,i,j,k) = PRIM(bIndx,firstpas,i,j,k)

          ! Calculate temperature of this cell
          call calcTemp (PRIM(bIndx,:,i,j,k), Temp)

          ! Cooling not applied below cool_Tmin
          if (Temp > Tmin_cool) then

            !  get the energy losses L_0 ([cm^-3 s^-1])
            L0 = aloss(y0,y1,dh,dh0,real(Temp,8))    !/dh**2

            !if (dif_rad) then
            !gain   = real(radphi(bIndx,i,j,k),8)*dh0*Kb*energy_per_ionization
            !tprime = max( gain*real(Temp,8)/L0, 7000.)
            !else
            !tprime=10.
            !end if
            gain = 0.0   !# add later
            Tprime = T_floor

            ce = (2.0*L0)/(3.0*Kb*dH*real(Temp,8))

            T1=Tprime+(Temp-Tprime)*exp(-ce*dt_seconds) !# new temperature
            cool_factor = T1/Temp

            frac_loss = 1.0-cool_factor

            ! Record maximum cooling for this block before limiting
            maxloss = max(maxloss, frac_loss)

            ! Limit cool_factor directly, if needed
            if (cool_factor.lt.1.0-cooling_limit) then
              cool_factor = 1.0-cooling_limit
            end if

            ! Impose a temperature floor by adjusting cool_factor, if needed
            ! Set to 10 K by default
            T1 = Temp * cool_factor

            if (T1 < T_floor) then
              T1 = T_floor
              cool_factor = T_floor / Temp
            end if

            ! Update pressure and total energy
            PRIM(bIndx,5,i,j,k) = PRIM(bIndx,5,i,j,k) * cool_factor

            ETH = CV * PRIM(bIndx,5,i,j,k)
            vel2 = PRIM(bIndx,2,i,j,k)**2                                    &
                 + PRIM(bIndx,3,i,j,k)**2                                    &
                 + PRIM(bIndx,4,i,j,k)**2
            EK = 0.5 * PRIM(bIndx,1,i,j,k) * vel2
            U(bIndx,5,i,j,k) = EK + ETH

          end if

        end do
      end do
    end do

    end subroutine apply_cooling_h_neq

!======================================================================

!#endif

end module cooling_H
