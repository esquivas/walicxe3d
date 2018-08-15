!===============================================================================
!> @file hlle.f90
!> @brief HLLE Riemann solver
!> @author Juan C. Toledo & Alex Esquivel
!> @date 14/Ago/2018

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

!> @brief Harten, Lax, van Leer (HLLE) Approximate Riemann Solver Module
!> @details Computes the  intercell numerical fluxes for every cell interface
!! in a block using the HLLE solver.

module HLLE

#ifdef BFIELD

  implicit none

contains

!===============================================================================

!> @brief Harten, Lax, van Leer (HLLE) Approximate Riemann Solver
!> @details Computes the intercell numerical fluxes for every cell interface
!! in a block using the HLLE solver. This routine assumes the global vector of
!! primitives is up-to-date (including enough ghost cells) and exits with
!! computed values in the intercell numerical fluxes (FC, GC, HC). Interfaces
!! are defined to the right, that is:
!! FC(:,i,j,k): intercell flux between cell (i,j,k) and cell (i+1,j,k)
!! GC(:,i,j,k): intercell flux between cell (i,j,k) and cell (i,j+1,k)
!! HC(:,i,j,k): intercell flux between cell (i,j,k) and cell (i,j,k+1)
!! The order of the solver must be specified. First-order employs piece-wise
!! constant primitives in each cell, while second-ordr uses linear
!! reconstruction to interpolate the primitives at cell interfaces.
!> @param locIndx Local index of the block
!> @param order Order of solver (interpolation for primitives)
subroutine HLLEfluxes (locIndx, order)

  use parameters
  use globals
  use hydro_core, only : swapxy, swapxz, limiter
  use amr, only : validCell
  implicit none

  integer, intent(in) :: locIndx
  integer, intent(in) :: order

  integer :: i, j, k
  real :: pl(neqtot), pr(neqtot), pll(neqtot), prr(neqtot), ff(neqtot)
  logical :: valid

  ! First-order method
  select case (order)

  ! -------------------------------------------

  case (1)   ! First-order

    do i=0,ncells_x
      do j=0,ncells_y
        do k=0,ncells_z

          call validCell (i,j,k,valid)
          if (valid) then

            ! X dimension
            pl(:) = PRIM(locIndx,:,i,j,k)
            pr(:) = PRIM(locIndx,:,i+1,j,k)
            call primfhlle (pL, pR, ff)
            FC(:,i,j,k) = ff(:)

            ! Y dimension
            pL(:) = PRIM(locIndx,:,i,j,k)
            pR(:) = PRIM(locIndx,:,i,j+1,k)
            call swapxy (pL)
            call swapxy (pR)
            call primfhlle (pL, pR, ff)
            call swapxy (ff)
            GC(:,i,j,k) = ff(:)

            ! Z dimension
            pL(:) = PRIM(locIndx,:,i,j,k)
            pR(:) = PRIM(locIndx,:,i,j,k+1)
            call swapxz (pL)
            call swapxz (pR)
            call primfhlle (pL, pR, ff)
            call swapxz (ff)
            HC(:,i,j,k) = ff(:)

          end if

        end do
      end do
    end do

  ! -------------------------------------------

  case (2)  ! Second-order - requires limiter

    do i=0,ncells_x
      do j=0,ncells_y
        do k=0,ncells_z

          call validCell (i,j,k,valid)
          if (valid) then

            ! X dimension
            pll(:) = PRIM(locIndx,:,i-1,j,k)
            pl(:)  = PRIM(locIndx,:,i,  j,k)
            pr(:)  = PRIM(locIndx,:,i+1,j,k)
            prr(:) = PRIM(locIndx,:,i+2,j,k)
            call limiter (pll,pl,pr,prr,limiter_type,neqtot)
            call primfhlle (pl, pr, ff)
            FC(:,i,j,k) = ff(:)

            ! Y dimension
            pll(:) = PRIM(locIndx,:,i,j-1,k)
            pl(:)  = PRIM(locIndx,:,i,j  ,k)
            pr(:)  = PRIM(locIndx,:,i,j+1,k)
            prr(:) = PRIM(locIndx,:,i,j+2,k)
            call swapxy (pll)
            call swapxy (pl)
            call swapxy (pr)
            call swapxy (prr)
            call limiter (pll,pl,pr,prr,limiter_type,neqtot)
            call primfhlle (pl, pr, ff)
            call swapxy (ff)
            GC(:,i,j,k) = ff(:)

            ! Z dimension
            pll(:) = PRIM(locIndx,:,i,j,k-1)
            pl(:)  = PRIM(locIndx,:,i,j,k  )
            pr(:)  = PRIM(locIndx,:,i,j,k+1)
            prr(:) = PRIM(locIndx,:,i,j,k+2)
            call swapxz (pll)
            call swapxz (pl)
            call swapxz (pr)
            call swapxz (prr)
            call limiter (pll,pl,pr,prr,limiter_type,neqtot)
            call primfhlle (pl, pr, ff)
            call swapxz (ff)
            HC(:,i,j,k) = ff(:)

          end if

        end do
      end do
    end do

  end select

end subroutine HLLEfluxes

!===============================================================================

!> @brief HLLE intercell fluxes along X using primitives at an interface
!> @param pL Vector of primitives at LEFT of interface
!> @param pR Vector of primitives at RIGHT of interface
!> @param ff Returned vector of HLLE intercell fluxes
subroutine primfhlle (pL, pR, ff)

  use parameters
  use hydro_core, only : prim2fluxes, prim2flow
  implicit none

  real, intent(in) :: pL(neqtot)
  real, intent(in) :: pR(neqtot)
  real, intent(out) :: ff(neqtot)

  real :: sl, sr
  real :: uL(neqtot), uR(neqtot)
  real :: fL(neqtot), fR(neqtot)

  ! Calculate wavespeeds
  call wavespeedHLLE (pL, pR, sl, sr)

!  write(*,*) sl, sr

  ! Obtain intercell fluxes as given by 10.21 of Toro
  if (sl.ge.0) then
    call prim2fluxes (pL, DIM_X, fL)
    ff(:) = fL(:)
    return
  else if (sr.le.0) then
    call prim2fluxes (pR, DIM_X, fR)
    ff(:) = fR(:)
    return
  else
    call prim2fluxes (pL, DIM_X, fL)
    call prim2fluxes (pR, DIM_X, fR)
    call prim2flow (pL, uL)
    call prim2flow (pR, uR)
    ff(:) = (sr*fl(:)-sl*fr(:)+sl*sr*(uR(:)-uL(:)))/(sr-sl)
    return
  end if

end subroutine primfhlle

!===============================================================================

!> @brief Obtains wavespeed estimates given primitives at interface
!> @details Computes wavespeed estimates for the two outer shock waves in the
!! Riemann fan. This routine employs the x-component of speed.
!> @param primL(neqtot) Primitives "left" of interface
!> @param primR(neqtot) Primitives "right" of interface
!> @param sl Wavespeed estimate for "left" moving wave
!> @param sr Wavespeed estimate for "right" moving wave
subroutine wavespeedHLLE (primL, primR, sl, sr)

  use parameters, only : neqtot
  use hydro_core, only : cfastX
  implicit none

  real, intent(in) :: primL(neqtot)
  real, intent(in) :: primR(neqtot)
  real, intent(out) :: sl
  real, intent(out) :: sr
  real              :: cfL, cfR

  ! Sound speeds
  call cfastX(primL,cfL)
  call cfastX(primR,cfR)

  ! Compute SL and SR
  ! Davis direct bounded, 10.38 of Toro
  sl = min( primL(2)-cfL, primR(2)-cfR )
  sr = max( primL(2)+cfL, primR(2)+cfR )

  return

end subroutine wavespeedHLLE
!===============================================================================

#endif

end module HLLE
