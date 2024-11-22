!===============================================================================
!> @file cooling_schure.f90
!> @brief Radiative cooling with Schure et al cooling curves
!> @author C. Villarreal & A. Esquivel
!> @date 31/Oct/2024

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

!> @brief Module for inclussion of cooling by radiative losses

module cooling_schure

  !> The following selects the column of the ionization fraction used for low T
  !> it is set in parameters
  !> dmc_f = 1 : f_i = 1e-4
  !>         2 : f_i = 1e-3
  !>         3 : f_i = 1e-2
  !>         4 : f_i = 1e-1


contains

  !===============================================================================

  !> @brief Loads radiative cooling coefficients the Schure et al. tables
  !> @details Loads radiative cooling coefficients from a data file into the
  !! cooldata global array. The number of temperature and
  !! metallicity data points must defined in the first line of the data file,
  !! separated by blank. The second line must contain a blank-separated list
  !! of columns (for the Schure et al. cooling they are the coefficients for
  !! different ionization fractions at low T).
  !! Then, every subsequent line must contain, first, the temperature, and then the
  !! corresponding cooling coefficients for each metallicity value.

  subroutine loadcooldata_schure ()

  use parameters
  use globals, only: cooltable, rank, logu, ierr
  use clean_quit, only : clean_abort
  implicit none

  integer :: i, istat
  real (kind=8) :: data(5)  ! Table has five columns

  if (rank.eq.master) then

    if (verbosity > 1) write(logu,'(2x,a)') "Loading cooling data from file ..."
    open (unit=99, file=cooling_file, status='old', iostat=istat)
    if (istat.ne.0) then
      write(logu,'(a,a,a)') "Could not open the file ", trim(cooling_file), " !"
      write(logu,*) "***ABORTING***"
      close(99)
      call clean_abort (ERROR_COOLING_LOAD_COEFS)
    end if

    allocate( cooltable(2,180) )

    do i=1,180
      read(99,*) data(:)
      cooltable(1,i)= 10.0**data(1)
      cooltable(2,i)= 10.0**(- data(1+dmc_f) ) !data(1+dmc_f)
    end do
    close (unit=99)

    end if

      ! Broadcast cooling data to all processes
      if (rank.ne.master) then
        if (verbosity > 2) write (logu,'(2x,a)') "Receiving cooling data from master process ..."
        allocate (cooltable(2,180))
      end if
      call mpi_bcast(cooltable, 180*2, mpi_real_kind, 0, mpi_comm_world, ierr)
      if (verbosity > 1) write(logu,'(2x,a,i0,a)') "Loaded ", 180, " cooling coefficients."

  ! Set global vars for minimum and maximum temperatures
  ! Note that the temperatures are logarithmic
  ! cool_Tmin = cooltable(1,1)
  ! cool_Tmax = cooltable(1,nptsT)

    end subroutine loadcooldata_schure

  !=======================================================================
  !> @brief Returns the cooling coefficient interpolating the table
  !> @param real [in] T : Temperature K
  function get_lambda(T)
  use globals, only : cooltable
  implicit none
  real , intent(in) :: T
  integer           :: if1
  real, parameter   :: deltaTemp=0.04 !  spacing of T in tables
  real, parameter   :: logTmin = 1.0 !  log of minimum value of temp
  real (kind=8)     :: get_lambda, T0, T1, C0, C1


  if(T.gt.1e8) then
    get_lambda=0.21e-26*sqrt(T)
  else
    if1=int(( log10(T)- logTmin) /deltaTemp) +1
    if (if1 < 1) then
      get_lambda = cooltable(2,1)
      return
    end if
    T0=cooltable(1,if1)
    c0=cooltable(2,if1)
    T1=cooltable(1,if1+1)
    c1=cooltable(2,if1+1)
    get_lambda=(c1-c0)*(T-T0)/(T1-T0)+c0
  end if

  end function get_lambda

  !=======================================================================
  !> @brief High level wrapper to apply cooling with CHIANTI tables
  !> @details High level wrapper to apply cooling with CHIANTI tables
  !> @n cooling is applied in the entire domain and updates both the
  !! conserved and primitive variables
  subroutine apply_cooling_schure(bIndx,maxloss)

  use parameters !, only : nx, ny, nz, cv, Psc, tsc, dif_rad, mhd, n1_chem
  use constants,  only : Kb
  use globals!,     only : u, primit, dt_CFL
  use hydro_core!, only : u2prim
  !use difradHe
  !use network

  implicit none

  integer, intent(in)  :: bIndx
  real, intent(out)    :: maxloss
  real                 :: T , dens
  real, parameter      :: Tmin=10.
  real                 :: frac_loss
  real (kind=8)        :: Lambda0, emtauC, gain
  integer              :: i, j, k   !, iiHI, iiHeI, iiHeII
  real                 :: dt_seconds, cool_factor

  maxloss = 0.0

  dt_seconds = dt*t_sc  ! chequear este dt
  !iiHI   = n1_chem - 1 + iHI
  !iiHeI  = n1_chem - 1 + iHeI
  !iiHeII = n1_chem - 1 + iHeII

  do i=1,ncells_x
    do j=1,ncells_y
      do k=1,ncells_z


        ! Calculate temperature of this cell
        call calcTemp (PRIM(bIndx,:,i,j,k), T)

        if(T > Tmin) then

          Lambda0=get_lambda(T)

!          if (dif_rad) then
!            !  energy per photo ionization from Black 1981 (in erg)
!            gain = phHI(bIndx,i,j,k)   * U(bIndx, iiHI  ,i,j,k) * 7.75e-12     &
!                 + phHeI(bIndx,i,j,k)  * U(bIndx, iiHeI ,i,j,k) * 2.19e-11     &
!                 + phHeII(bIndx,i,j,k) * U(bIndx, iiHeII,i,j,k) * 3.10e-11
!
!          else

            gain=0.0

!          end if

          dens=PRIM(bIndx,1,i,j,k) / mu0
          !  e^{-dt/Tau}=e^{-2. L0 dt/(3 n K T)}
          emtauC = exp( -2.0*dt_seconds*dens*Lambda0/(3.0*Kb*T) )
          !  this is the Temperature factor of change
          cool_factor = (gain/(dens**2*Lambda0))*(1.0-emtauC)+emtauC

          frac_loss = 1.0-cool_factor
          ! Record maximum cooling for this block before limiting
          maxloss = max(maxloss, frac_loss)

          ! Limit cool_factor directly, if needed
          if (cool_factor.lt.1.0-cooling_limit) then
            cool_factor = 1.0-cooling_limit
          end if

          !!! DEBUG 
          !cool_factor = 1.0

          ! !  limit changes to avoid catastrophic cooling  ! uso este o el de arriba?
          ! ch_factor = min(ch_factor,0.5)
          ! ch_factor = max(ch_factor,2.0)

          ! Update pressure and total energy
          PRIM(bIndx,5,i,j,k) = PRIM(bIndx,5,i,j,k) * cool_factor

          U(bIndx,5,i,j,k) = CV *PRIM(bIndx,5,i,j,k)                          &
                           + 0.5*PRIM(bIndx,1,i,j,k)*(PRIM(bIndx,2,i,j,k)**2  &
                                                    + PRIM(bIndx,3,i,j,k)**2  &
                                                    + PRIM(bIndx,4,i,j,k)**2)
  ! #ifdef BFIELD
  !         if (mhd) then
  !           u(5,i,j,k) = u(5,i,j,k) + 0.5*(  primit(6,i,j,k)**2                  &
  !                                          + primit(7,i,j,k)**2                  &
  !                                          + primit(8,i,j,k)**2  )
  !         end if
  ! #endif

        end if
      end do
    end do
  end do


  end subroutine apply_cooling_schure

end module cooling_schure
