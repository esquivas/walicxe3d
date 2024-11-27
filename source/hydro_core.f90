!===============================================================================
!> @file prims.f90
!> @brief Calculation of primitives and conversion from/to flow variables
!> @author Alex Esquivel
!> @date 16/Aug/2018

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

!> @brief Basic hydro (and MHD) subroutines utilities
!> @details This module contains subroutines and utilities that are the
!> core of the hydro (and MHD) that are common to most implementations
!> and will be used for the different specific solvers

module hydro_core

  use clean_quit, only : clean_abort
  implicit none

contains

!===============================================================================

!> @brief Update primitives in *all* local blocks
!> @details This high-level wrapper routine updates primitives on all
!! cells of all local blocks, using the U flow vars array as source.
subroutine updatePrims ()

  use parameters, only : verbosity
  use globals, only : logu, U, PRIM
  use constants, only : CELLS_ALL
  use tictoc, only : tic, nicetoc
  implicit none

  integer :: mark

  if (verbosity > 1) write(logu,*) ""
  if (verbosity > 1) write(logu,'(1x,a)') "> Updating primitive variables ..."
  if (verbosity > 3) call tic(mark)
  call calcPrimsAll (U, PRIM, CELLS_ALL)
  if (verbosity > 3) write(logu,'(1x,a,a)') "> Updated primitives in", nicetoc(mark)

end subroutine updatePrims

!===============================================================================

!> @brief Wrapper routine to update primitives for all local blocks
!> @details The full vectors of flow variables and primitives must be given.
!! The 'cells' argument indicates the range of cells to update. See the
!! documentation of the calcPrimsBlock routine for further details.
subroutine calcPrimsAll (uvars, pvars, cells)

  use parameters, only : nxmin, nxmax, nymin, nymax, nzmin, nzmax, &
                         nbMaxProc, neqtot, verbosity
  use globals, only : logu, localBlocks
  implicit none

  real, intent(in)  :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  real, intent(out) :: pvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  integer, intent(in) :: cells

  integer :: nb, bID, badcells

  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then

      badcells = 0
      call calcPrimsBlock (uvars, pvars, nb, cells, badcells)

      if ((verbosity > 1).and.(badcells.ne.0)) then
        write(logu,'(1x,a,i0,a,i0)') "Warning: ", badcells, &
        " pressure corrections in block ", bID
      end if

    end if
  end do

end subroutine calcPrimsAll

!===============================================================================

!> @brief Calculates and update primitives in one block
!> @details The full vectors of flow variables and primitives must be given.
!! The 'cells' argument indicates the range of cells to update. Currently
!! recognized options are:
!!  CELLS_ALL: update physical and ghost cells
!!  CELLS_PHYS: update only physical cells
!!  CELLS_GHOST: update only ghost cells
!! The 'badcells' argument returns the number of cells in which
!! pressure corrections were needed.
subroutine calcPrimsBlock (uvars, pvars, locIdx, cells, badcells)

  use parameters, only : nxmin, nxmax, nymin, nymax, nzmin, nzmax, nghost,  &
                         nbMaxProc, neqtot, ncells_x, ncells_y, ncells_z
  use constants
  use globals, only : logu
  implicit none

  real, intent(in) :: uvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  real, intent(out) :: pvars(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  integer, intent(in) :: locIdx
  integer, intent(in) :: cells
  integer, intent(out) :: badcells

  integer :: i, j, k, istat

  badcells = 0

  ! Physical cells
  if ((cells.eq.CELLS_ALL).or.(cells.eq.CELLS_PHYS)) then

    do k=1,ncells_z
      do j=1,ncells_y
        do i=1,ncells_x
          call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
          if (istat.ne.0) badcells = badcells + 1
        end do
      end do
    end do

  end if

  ! Ghost cells (will skip calculation if density is zero)
  if ((cells.eq.CELLS_ALL).or.(cells.eq.CELLS_GHOST)) then

    ! Left ghost cells
    do k=nzmin,nzmax  !1,ncells_z
      do j=nymin, nymax!1,ncells_y
        do i=1-nghost,0
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
            if (istat.ne.0) badcells = badcells + 1
          end if
        end do
      end do
    end do

    ! Right ghost cells
    do k=nzmin, nzmax!1,ncells_z
      do j=nymin, nymax!1,ncells_y
        do i=ncells_x+1,ncells_x+2
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
            if (istat.ne.0) badcells = badcells + 1
          end if
        end do
      end do
    end do

    ! Front ghost cells
    do k=nzmin, nzmax!1,ncells_z
      do j=1-nghost,0
        do i=nxmin, nxmax!1,ncells_x
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
            if (istat.ne.0) badcells = badcells + 1
          end if
        end do
      end do
    end do

    ! Back ghost cells
    do k=nzmin, nzmax!1,ncells_z
      do j=ncells_y+1,ncells_y+2
        do i=nxmin, nxmax!1,ncells_x
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
            if (istat.ne.0) badcells = badcells + 1
          end if
        end do
      end do
    end do

    ! Bottom ghost cells
    do k=1-nghost,0
      do j=nymin, nymax!1,ncells_y
        do i=nxmin, nxmax!1,ncells_x
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
            if (istat.ne.0) badcells = badcells + 1
          end if
        end do
      end do
    end do

    ! Bottom ghost cells
    do k=ncells_z+1,ncells_z+2
      do j=nymin, nymax!1,ncells_y
        do i=nxmin, nxmax!1,ncells_x
          if (uvars(locIdx,1,i,j,k).ne.0.0) then
            call flow2prim (uvars(locIdx,:,i,j,k), pvars(locIdx,:,i,j,k), istat)
            if (istat.ne.0) badcells = badcells + 1
          end if
        end do
      end do
    end do

  end if

  if ((cells.ne.CELLS_ALL).and.(cells.ne.CELLS_PHYS).and.(cells.ne.CELLS_GHOST)) then
    write(logu,*) "Invalid cell range passed to updatePrimBlock!"
    write(logu,*) "***Aborting!***"
    call clean_abort (ERROR_INVALID_CELL_RANGE)
  end if

end subroutine calcPrimsBlock

!===============================================================================

!> @brief Estimates temperature from primitives, in CGS
!> @details Will use the mean atomic mass for ionized gas if the estimated
!! temperature (using the neutral mean atomic mass) exceeds ion_thres
!> @param pvars Vector of primitives
subroutine calcTemp (pvars, temp)

  use parameters
  use globals
  implicit none

  real, intent(in) :: pvars(neqtot)
  real, intent(out) :: temp

  EOS : select case(eos_type)

    case(EOS_ADIABATIC)
      !  Ideal gas
      temp = pvars(5)/pvars(1)*(mu0*AMU*p_sc/d_sc/KB)

    case(EOS_SINGLE_SPECIE)
     !  Ideal gas
      temp = pvars(5)/pvars(1)*(mu0*AMU*p_sc/d_sc/KB)

    case(EOS_TWOTEMP)
      !  Ideal gas, two mu's above/below ionization threshold
      temp = pvars(5)/pvars(1)*(mu0*AMU*p_sc/d_sc/KB)

       if (temp.gt.ion_thres) then
         temp = pvars(5)/pvars(1)*(mui*AMU*p_sc/d_sc/KB)
       end if

    case(EOS_H_RATE)
      !  uses P = (n_HI + nHII) K T ; (assumes ne=n_HII)
      !  assumes that firstpass contains the neutral H density
      if (npassive >= 1) then
        temp   = pvars(5)/( 2.0*pvars(1)/mu0 - pvars(firstpas) ) * p_sc/KB
        ! *d_sc/mH  [cm^-3]
      end if

  end select EOS

  return

end subroutine calcTemp

!===============================================================================

!> @brief Calculates the vector of conserved (flow) vars given the vector of
!! primitive variables. Passive scalars are simply copied over.
!> @param pvars(neqtot) An input vector containing (rho,u,v,w,P,s1,s2,...)
!> @param uvars(neqtot) An output vector containing (rho,rho*u,rho*v,rho*w,E,s1,s2,...)
subroutine prim2flow (pvars, uvars)

  use parameters, only : neqtot, firstpas, CV, npassive, mhd
  implicit none

  real, intent(in) :: pvars(neqtot)
  real, intent(out) :: uvars(neqtot)

  real :: v2

  uvars(1) = pvars(1)
  uvars(2) = pvars(1)*pvars(2)
  uvars(3) = pvars(1)*pvars(3)
  uvars(4) = pvars(1)*pvars(4)


  v2 = pvars(2)**2 + pvars(3)**2 + pvars(4)**2
  uvars(5) = 0.5*pvars(1)*v2 + CV*pvars(5)

#ifdef BFIELD
  if (mhd) then
    uvars(5) = uvars(5) + 0.5 * ( pvars(6)**2 + pvars(7)**2 + pvars(8)**2 )
  end if

  uvars(6) = pvars(6)
  uvars(7) = pvars(7)
  uvars(8) = pvars(8)
#endif

  if (npassive >= 1) then
    uvars(firstpas:neqtot) = pvars(firstpas:neqtot)
  end if

  return

end subroutine prim2flow

!===============================================================================

!> @brief Calculates the vector of primitive variables given the vector of
!! conserved (flow) variables. Any passive scalars are simply copied over.
!> @param uvars(neqtot) An input vector containing (rho,rho*u,rho*v,rho*w,E,...)
!> @param pvars(neqtot) An output vector containing (rho,u,v,w,P,...)
!> @param istat An output integer indicating the success status of the
!! operation (0=all good, 1=pressure floored, 2=density floored).
subroutine flow2prim (uvars, pvars, istat)

  use parameters
  use globals   ! DEBUG
  implicit none

  real, intent(in) :: uvars(neqtot)
  real, intent(out) :: pvars(neqtot)
  integer, intent(out) :: istat

  istat = 0

  if (uvars(1)==0.0) then
    write(logu,*) "Received zero density (in flow2prim)!!"
    call clean_abort(ERROR_DIVISION_BY_ZERO)
  end if

  pvars(1) = uvars(1)
  pvars(2) = uvars(2)/uvars(1)
  pvars(3) = uvars(3)/uvars(1)
  pvars(4) = uvars(4)/uvars(1)

  pvars(5) = ( uvars(5) -                                               &
             0.5* (uvars(2)**2 + uvars(3)**2 + uvars(4)**2) /uvars(1) ) / cv

#ifdef BFIELD
  if (mhd) then
    pvars(5) = pvars(5) - 0.5*(uvars(6)**2 + uvars(7)**2 + uvars(8)**2) / cv
  end if
#endif

  ! Floor on pressure
  if (pvars(5) < 1.0e-30) then
    pvars(5) = 1.0e-30
    istat = 1
  end if

  ! Floor on density
  if (pvars(1) < 1.0e-40) then
    pvars(1) = 1.0e-40
    istat = 2
  end if

#ifdef BFIELD
  pvars(6) = uvars(6)
  pvars(7) = uvars(7)
  pvars(8) = uvars(8)
#endif

  if (npassive >= 1) then
    pvars(firstpas:neqtot) = uvars(firstpas:neqtot)
  end if

  return

end subroutine flow2prim

!===============================================================================

!> @brief Calculates the vectors of physical fluxes along a specific dimension
!! given the vector of primitive variables.
!> @param pvars(neqtot) An input vector containing primitives
!> @param dimf Dimension along which fluxes are to be calculated
!> @param F(neqtot) An output vector containing fluxes along dimf
subroutine prim2fluxes (pvars, dimens, flux)

  use parameters
  implicit none

  real, intent(in) :: pvars(neqtot)
  integer, intent(in) :: dimens
  real, intent(out) :: flux(neqtot)

  real :: pvars1(neqtot)

  ! Make a copy of the primitives so the original ones are not split
  pvars1(:) = pvars(:)

  ! Swap velocity components when dimf=2 (Y) or dimf=3(Z)
  if (dimens.eq.DIM_Y) then
    call swapxy(pvars1)
  else if (dimens.eq.DIM_Z) then
    call swapxz(pvars1)
  end if

  ! Calculate Fluxes (formulae for X)
  flux(1) = pvars1(1)*pvars1(2)
  flux(2) = pvars1(1)*pvars1(2)**2 + pvars1(5)
  flux(3) = pvars1(1)*pvars1(2)*pvars1(3)
  flux(4) = pvars1(1)*pvars1(2)*pvars1(4)
  flux(5) = pvars1(2)*( 0.5*pvars1(1)*(pvars1(2)**2+pvars1(3)**2+pvars1(4)**2) &
            + (cv+1.)*pvars1(5) )
#ifdef BFIELD
  flux(6) = 0.
  flux(7) = pvars1(2)*pvars1(7) - pvars1(3)*pvars1(6)
  flux(8) = pvars1(2)*pvars1(8) - pvars1(4)*pvars1(6)

  if (mhd) then
    flux(2) = flux(2) + 0.5 *( pvars1(7)**2 + pvars1(8)**2 - pvars1(6)**2 )
    flux(3) = flux(3) - pvars1(6)*pvars1(7)
    flux(4) = flux(4) - pvars1(6)*pvars1(8)
    flux(5) = flux(5) + pvars1(2)*(pvars1(6)**2+pvars1(7)**2+pvars1(8)**2)     &
    - pvars1(6)*( pvars1(2)*pvars1(6)+pvars1(3)*pvars1(7)+pvars1(4)*pvars1(8) )
  end if
#endif

  if (npassive >= 1) then
    flux(firstpas:neqtot) = pvars1(2)*pvars1(firstpas:neqtot)
  end if

  ! Swap flux components to get correct fluxes for Y or Z
  if (dimens.eq.DIM_Y) then
    call swapxy(flux)
  else if (dimens.eq.DIM_Z) then
    call swapxz(flux)
  end if

  return

end subroutine prim2fluxes

!===============================================================================

!> @brief Swaps the x and y components of a vector of hydro variables
!> @param vec The vector to be swapped
subroutine swapxy (vec)

  use parameters
  implicit none

  real, intent(inout) :: vec(neqtot)

  real :: temp

  temp = vec(2)
  vec(2) = vec(3)
  vec(3) = temp

#ifdef BFIELD
  temp = vec(6)
  vec(6) = vec(7)
  vec(7) = temp
#endif

  return

end subroutine swapxy

!===============================================================================

!> @brief Swaps the x and z components of a vector of hydro variables
!> @param vec The vector to be swapped
subroutine swapxz (vec)

  use parameters
  implicit none

  real, intent(inout) :: vec(neqtot)

  real :: temp

  temp = vec(2)
  vec(2) = vec(4)
  vec(4) = temp

#ifdef BFIELD
  temp = vec(6)
  vec(6) = vec(8)
  vec(8) = temp
#endif

  return

end subroutine swapxz

!===============================================================================

!> @brief Calculates the soundspeed given a vector of primitives
!> @param pvars Vector of primitive variables
!> @param csound The hydrodynamical speed of sound
subroutine sound (pvars, csound)

  use parameters
  implicit none

  real, intent(in) :: pvars(neqtot)
  real, intent(out) :: csound

  csound = sqrt(gamma*pvars(5)/pvars(1))

  return

end subroutine sound

!=======================================================================

!> @brief Computes the fast magnetosonic speeds  in the 3 coordinates
!> @details Computes the fast magnetosonic speeds  in the 3 coordinates
!> @param real [in] p  : value of pressure
!> @param real [in] d  : value of density
!> @param real [in] Bx : value of the x component of the magnetic field
!> @param real [in] By : value of the y component of the magnetic field
!> @param real [in] Bz : value of the z component of the magnetic field
!> @param real [out] csx : fast magnetosonic speed in x
!> @param real [out] csy : fast magnetosonic speed in y
!> @param real [out] csz : fast magnetosonic speed in z

subroutine cfast(p,d,bx,by,bz,cfx,cfy,cfz)

  use parameters, only : gamma
  implicit none
  real, intent(in) :: p, d, bx, by, bz
  real, intent(out) ::cfx,cfy,cfz
  real :: b2

  b2=bx*bx+by*by+bz*bz
  cfx=sqrt(0.5*((gamma*p+b2)+sqrt((gamma*p+b2)**2-4.*gamma*p*bx*bx))/d)
  cfy=sqrt(0.5*((gamma*p+b2)+sqrt((gamma*p+b2)**2-4.*gamma*p*by*by))/d)
  cfz=sqrt(0.5*((gamma*p+b2)+sqrt((gamma*p+b2)**2-4.*gamma*p*bz*bz))/d)

end subroutine cfast

!=======================================================================

!> @brief Computes the fast magnetosonic speed in the x direction
!> @details Computes the fast magnetosonic speed in the x direction
!> @param real [in] prim(neq) : vector with the primitives in one cell

#ifdef BFIELD

subroutine cfastX(prim,cfX)

  use parameters, only : neqtot, gamma
  implicit none
  real, intent(in) :: prim(neqtot)
  real, intent(out) ::cfX
  real :: b2, cs2va2

  b2=prim(6)**2+prim(7)**2+prim(8)**2
  cs2va2 = (gamma*prim(5)+b2)/prim(1)   ! cs^2 + ca^2

  cfx=sqrt(0.5*(cs2va2+sqrt(cs2va2**2-4.*gamma*prim(5)*prim(6)**2/prim(1)/prim(1) ) ) )

end subroutine cfastX

#endif

!===============================================================================

!> @brief Applies artificial viscosity to a block's flow variables
!> @details Viscosity: up(n)=up(n)+eta*\nabla^2(u(n))
!> @param bIndx Block's local index
!> @param U Flow variables at t^n
!> @param UP Flow variables at t^(n+1)
!!subroutine viscosity (bIndx, U, UP)
!!
!!  use parameters
!!  implicit none
!!
!!  integer, intent(in) :: bIndx
!!  real, intent(in) :: U (nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
!!  real, intent(inout) :: UP (nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
!!
!!  integer :: i, j, k
!!
!!  do i=1,ncells_x
!!    do j=1,ncells_y
!!      do k=1,ncells_z
!!        UP(bIndx,:,i,j,k)   &
!!          = UP(bIndx,:,i,j,k) + visc_eta*(   &
!!          + U(bIndx,:,i+1,j,k) + U(bIndx,:,i-1,j,k)   &
!!          + U(bIndx,:,i,j+1,k) + U(bIndx,:,i,j-1,k)   &
!!          + U(bIndx,:,i,j,k+1) + U(bIndx,:,i,j,k-1)   &
!!          - 6.0*U(bIndx,:,i,j,k) )
!!      end do
!!    end do
!!  end do
!!
!!end subroutine viscosity

!===============================================================================

!> @brief Applies the flux limiter to average left/right states
!> @details This routine uses the X-component of speed.
!> @param pll Pressure two cells "left" of interface
!> @param pl Pressure one cell "left" of interface
!> @param prr Pressure two cells "right" of interface
!> @param pr Pressure once cell "right" of interface
!> @param lim Limiter to use. See parameters.f90 for supported options.
!> @param neqs Number of equations
subroutine limiter (pll,pl,pr,prr,lim,neqs)

  implicit none

  integer, intent(in) :: neqs
  integer, intent(in) :: lim
  real, intent(in) :: pll(neqs)
  real, intent(in) :: prr(neqs)
  real, intent(inout) :: pl(neqs)
  real, intent(inout) :: pr(neqs)


  real :: dl, dm, dr, al, ar
  integer :: ieq

  do ieq=1,neqs
    dl = pl(ieq) - pll(ieq)
    dm = pr(ieq) - pl(ieq)
    dr = prr(ieq) - pr(ieq)
    al = average(dl, dm, lim)
    ar = average(dm, dr, lim)
    pl(ieq) = pl(ieq) + 0.5*al
    pr(ieq) = pr(ieq) - 0.5*ar
  end do

contains

  real function average (a,b,opt)

    use constants
    implicit none

    real, intent(in) :: a, b
    integer, intent(in) :: opt

    real :: s, c, d, eps

    select case (opt)

    case (LIMITER_NO_AVERAGE)
      average = 0.

    case (LIMITER_NONE)
      average = 0.5*(a+b)

    case (LIMITER_VANLEER)
      if (a*b.le.0.0) then
        average = 0.0
      else
        average = a*b*(a+b)/(a*a+b*b)
      end if

    case (LIMITER_MINMOD)
      s = sign(1.0,a)
      average = s*max(0.0, min(abs(a), s*b))

    case (LIMITER_ALBADA)   ! NOT WORKING
      eps = 1.0e-7
      average = (a*(b*b+eps)+b*(a*a+eps))/(a*a+b*b*eps)

    case (LIMITER_UMIST)
      s = sign(1.0,a)
      c = 0.25*a + 0.75*b
      d = 0.75*a + 0.25*b
      average = min(2.0*abs(a), 2.0*s*b, s*c, s*d)
      average = s*max(0.0, average)

    case (LIMITER_WOODWARD)
      s = sign(1.0,a)
      c = 0.5*(a+b)
      average = min(2.0*abs(a), 2*s*b, s*c)
      average = s*max(0.0, average)

    case (LIMITER_SUPERBEE)
      s = sign(1.0,b)
      c = min(2.0*abs(b), s*a)
      d = min(abs(b),2.0*s*a)
      average = s*max(0.0,c,d)

    case default
      average = 0.0
      write(*,'(a)') "WARNING: no averaging in limiter!"
      write(*,'(a,i2)') "Passed limiter value: ", opt

    end select

  end function average

end subroutine limiter

!===============================================================================

end module hydro_core
