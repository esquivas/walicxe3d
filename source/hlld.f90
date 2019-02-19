!===============================================================================
!> @file hlld.f90
!> @brief HLLD Riemann solver
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

!> @brief Harten, Lax, van Leer (HLLD) Approximate Riemann Solver Module
!> @details Computes the  intercell numerical fluxes for every cell interface
!! in a block using the HLLD solver.

module HLLD

#ifdef BFIELD

  implicit none

contains

!===============================================================================

!> @brief Harten, Lax, van Leer (HLLD) Approximate Riemann Solver
!> @details Computes the intercell numerical fluxes for every cell interface
!! in a block using the HLLD solver. This routine assumes the global vector of
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
subroutine HLLDfluxes (locIndx, order)

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
            call primfhlld (pL, pR, ff)
            FC(:,i,j,k) = ff(:)

            ! Y dimension
            pL(:) = PRIM(locIndx,:,i,j,k)
            pR(:) = PRIM(locIndx,:,i,j+1,k)
            call swapxy (pL)
            call swapxy (pR)
            call primfhlld (pL, pR, ff)
            call swapxy (ff)
            GC(:,i,j,k) = ff(:)

            ! Z dimension
            pL(:) = PRIM(locIndx,:,i,j,k)
            pR(:) = PRIM(locIndx,:,i,j,k+1)
            call swapxz (pL)
            call swapxz (pR)
            call primfhlld (pL, pR, ff)
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
            call primfhlld (pl, pr, ff)
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
            call primfhlld (pl, pr, ff)
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
            call primfhlld (pl, pr, ff)
            call swapxz (ff)
            HC(:,i,j,k) = ff(:)

          end if

        end do
      end do
    end do

  end select

end subroutine HLLDfluxes

!===============================================================================

!> @brief HLLD intercell fluxes along X using primitives at an interface
!> @param primL Vector of primitives at LEFT of interface
!> @param primR Vector of primitives at RIGHT of interface
!> @param ff Returned vector of HLLD intercell fluxes
subroutine primfhlld(priml,primr,ff)

  use parameters, only : neqtot, cv, npassive, firstpas
  use globals,    only : logu
  use hydro_core, only : cfastX, prim2fluxes
  use clean_quit, only : clean_abort
  use constants,  only  : DIM_X, ERROR_NAN_IN_FLUXES
  implicit none
  real, dimension(neqtot),intent(in   ) :: priml, primr
  real, dimension(neqtot),intent(inout) :: ff

  real :: csl, csr, sl, sr, slmul, srmur, rholul, rhorur, sM
  real :: pTL, pTR, Bx, signBx
  real :: slmsM, srmsM, rhostl, rhostr, sstl, sstr
  real :: pst, el, er, denl, denr, sMmul, sMmur
  real :: vstl, wstl, bystl, bzstl, estl, vdotbl, vstdotbstl
  real :: vstr, wstr, bystr, bzstr, estr, vdotbr, vstdotbstr
  real :: dd, vstst, wstst, bystst, bzstst
  real ::  vststdotbstst, eststl, eststr

  call cfastX(priml,csl)
  call cfastX(primr,csr)

  sr=max(priml(2)+csl,primr(2)+csr)
  sl=min(priml(2)-csl,primr(2)-csr)

  ! UL region -----------------------------------
  if (sl > 0) then
    call prim2fluxes(priml,DIM_X,ff)
    return
  endif

  ! UR region -----------------------------------
  if (sr < 0) then
    call prim2fluxes(primr,DIM_X,ff)
    return
  endif

  Bx= 0.5* (primL(6)+primR(6) )
  signBx= sign(1.,Bx)

  !  Total pressure
  pTL=primL(5) + 0.5*( bx**2+primL(7)**2+primL(8)**2 )
  pTR=primR(5) + 0.5*( bx**2+primR(7)**2+primR(8)**2 )

  slmul=sl-priml(2)
  srmur=sr-primr(2)

  rholul=priml(1)*priml(2)
  rhorur=primr(1)*primr(2)

  sM = (srmur*rhorur-slmul*rholul-pTR+pTL)/( srmur*primr(1)-slmul*priml(1) )

  srmsM=sr-sM
  slmsM=sl-sM

  rhostl=priml(1)*slmul/slmsM      !rhoL*
  rhostr=primr(1)*srmur/srmsM      !rhoR*

  sstl=sM - abs(bx)/sqrt(rhostl)  !SL*
  sstr=sM + abs(bx)/sqrt(rhostr)  !SR*

  pst= (srmur*primr(1)*pTL - slmul*priml(1)*pTR               &
  + priml(1)*primr(1)*srmur*slmul*(primr(2)-priml(2)) ) &
  /( srmur*primr(1)-slmul*priml(1) )

  ! UL* region -----------------------------------
  if(sstl >= 0) then
    !
    el=0.5*priml(1)*(priml(2)**2+priml(3)**2+priml(4)**2)+cv*priml(5) &  !eL
    +0.5*(bx**2+priml(7)**2+priml(8)**2)

    sMmuL=sM - priml(2)
    denl=priml(1)*slmul*slmsM-bx**2

    if(denl == 0) then
      vstl = primL(3)
      wstl = primL(4)
      bystl= 0.!primL(7)
      bzstl= 0.!primL(8)
      !print*,'stopped @ HLLD'
      !stop
    else

      vstl = priml(3) - bx*priml(7)*sMmul/denl                      !vL*
      wstl = priml(4) - bx*priml(8)*sMmul/denl                      !wL*
      bystl= priml(7)*( priml(1)*slmul**2 - bx**2 )/denl            !byL*
      bzstl= priml(8)*( priml(1)*slmul**2 - bx**2 )/denl            !bzL*

    end if

    vdotbl    = priml(2)*bx + priml(3)*priml(7) + priml(4)*priml(8) !vL dot BL
    vstdotbstl= sM*bx + vstl*bystl + wstl*bzstl                     !vL* dot BL*

    estl= ( slmul*el -pTL*priml(2) +pst*sM +bx*(vdotbl-vstdotbstl) )/slmsM !eL*

    ff(1) = rhostl*sM
    ff(2) = rhostl*SM**2+pst-bx**2
    ff(3) = rhostl*sM*vstl-bx*bystl
    ff(4) = rhostl*sM*wstl-bx*bzstl
    ff(5) = sM*(estl+pst)-bx*(vstdotbstl)
    ff(6) = 0.
    ff(7) = bystl*sM-bx*vstl
    ff(8) = bzstl*sM-bx*wstl

#ifdef PASSIVES
    if (npassive >= 1) then
      ff(firstpas:neqtot)=sM*priml(firstpas:neqtot)*slmul/slmsM
    end if
#endif

    return
  endif

  ! UR* region -----------------------------------
  if(sstr <= 0) then

    er=0.5*primr(1)*(primr(2)**2+primr(3)**2+primr(4)**2)+cv*primr(5) &  !eR
    +0.5*(bx**2+primr(7)**2+primr(8)**2)

    sMmuR=sM - primr(2)
    denr=primr(1)*srmur*sRmsM-bx**2

    if(denr == 0) then
      vstl = primL(3)
      wstl = primL(4)
      bystl= 0.!primL(7)
      bzstl= 0.!primL(8)
      !print*,'stopped @ HLLD'
      !stop
    else

      vstr = primr(3) - bx*primr(7)*sMmuR/denr                      !vR*
      wstr = primr(4) - bx*primr(8)*sMmuR/denr                      !wR*
      bystr= primr(7)*( primr(1)*srmur**2 - bx**2 )/denr            !byR*
      bzstr= primr(8)*( primr(1)*srmur**2 - bx**2 )/denr            !bzR*

    end if

    vdotbr    = primr(2)*bx + primr(3)*primr(7) + primr(4)*primr(8) !vR dot BR
    vstdotbstr= sM*bx + vstr*bystr + wstr*bzstr                     !vR* dot BR*

    estr= ( srmur*er -pTR*primr(2) +pst*sM +bx*(vdotbr-vstdotbstr) )/srmsM !eR*

    ff(1) = rhostr*sM
    ff(2) = rhostr*SM**2+pst-bx**2
    ff(3) = rhostr*sM*vstr-bx*bystr
    ff(4) = rhostr*sM*wstr-bx*bzstr
    ff(5) = sM*(estr+pst)-bx*(vstdotbstr)
    ff(6) = 0.
    ff(7) = bystr*sM-bx*vstr
    ff(8) = bzstr*sM-bx*wstr

#ifdef PASSIVES
    if (npassive >= 1) then
      ff(firstpas:neqtot)=sM*primr(firstpas:neqtot)*srmur/srmsM
    end if
#endif

    return
  endif

  !   All this are needed on both the UL** and UR** regions
  sMmul= sM - priml(2)
  sMmur= sM - primr(2)

  denl=priml(1)*slmul*slmsM-bx**2
  denr=primr(1)*srmur*srmsM-bx**2

  if(denl == 0) then
    vstl =priml(3)
    wstl =priml(4)
    bystl=0.!priml(7)
    bzstl=0.!priml(8)
    !print*,'stopped @ HLLD'
    !stop
  else
    vstl = priml(3) - bx*priml(7)*sMmul/denl                      !vL*
    wstl = priml(4) - bx*priml(8)*sMmul/denl                      !wL*
    bystl= priml(7)*( priml(1)*slmul**2 - bx**2 )/denl            !byL*
    bzstl= priml(8)*( priml(1)*slmul**2 - bx**2 )/denl            !bzL*
  endif

  if(denr == 0) then
    vstr =primr(3)
    wstr =primr(4)
    bystr=0.!primr(7)
    bzstr=0.!primr(8)
    !print*,'stopped @ HLLD'
    !stop
  else
    vstr = primr(3) - bx*primr(7)*sMmuR/denr                      !vR*
    wstr = primr(4) - bx*primr(8)*sMmuR/denr                      !wR*
    bystr= primr(7)*( primr(1)*srmur**2 - bx**2 )/denr            !byR*
    bzstr= primr(8)*( primr(1)*srmur**2 - bx**2 )/denr            !bzR*
  endif

  dd=sqrt(rhostl)+sqrt(rhostr)

  vstst =(sqrt(rhostl)*vstl + sqrt(rhostr)*vstr + (bystr-bystl)*signBx )/dd  !v**
  wstst =(sqrt(rhostl)*wstl + sqrt(rhostr)*wstr + (bzstr-bzstl)*signBx )/dd  !w**

  bystst=(sqrt(rhostl)*bystr + sqrt(rhostr)*bystl +  &      !by**
  sqrt(rhostl*rhostr)*(vstr-vstl)*signBx )/dd

  bzstst=(sqrt(rhostl)*bzstr + sqrt(rhostr)*bzstl +  &      !bz**
  sqrt(rhostl*rhostr)*(wstr-wstl)*signBx )/dd

  vststdotbstst= sM*bx + vstst*bystst + wstst*bzstst        !v** dot B**

  ! UL** region -----------------------------------
  if(sM >= 0 ) then
    !
    el=0.5*priml(1)*(priml(2)**2+priml(3)**2+priml(4)**2)+cv*priml(5) &  !eL
    +0.5*(bx**2+priml(7)**2+priml(8)**2)

    vdotbl    = priml(2)*bx + priml(3)*priml(7) + priml(4)*priml(8) !vL dot BL
    vstdotbstl= sM*bx+ vstl*bystl + wstl*bzstl                      !vL* dot BL*

    estl= ( slmul*el -pTL*priml(2) +pst*sM +bx*(vdotbl-vstdotbstl) )/slmsM !eL*

    eststl= estl - sqrt(rhostl)*(vstdotbstl-vststdotbstst)*signBx !eL**

    ff(1) = rhostl*sM
    ff(2) = rhostl*SM**2+pst-bx**2
    ff(3) = rhostl*sM*vstst-bx*bystst
    ff(4) = rhostl*sM*wstst-bx*bzstst
    ff(5) = sm*(eststl+pst)-bx*(vststdotbstst)
    ff(6) = 0.
    ff(7) = bystst*sM-bx*vstst
    ff(8) = bzstst*sM-bx*wstst

#ifdef PASSIVES
    if (npassive >= 1) then
      ff(firstpas:neqtot) = sM*priml(nfirstpas:neqtot)*slmul/slmsM
    end if
#endif

    return
  endif

  ! UR** region -----------------------------------
  if(sM <= 0 ) then

    er=0.5*primr(1)*(primr(2)**2+primr(3)**2+primr(4)**2)+cv*primr(5) &  !eR
    +0.5*(bx**2+primr(7)**2+primr(8)**2)

    vdotbr    = primr(2)*bx + primr(3)*primr(7) + primr(4)*primr(8) !vR dot BR
    vstdotbstr= sM*bx+ vstr*bystr + wstr*bzstr                      !vR* dot BR*

    estr= ( srmur*er -pTR*primr(2) +pst*sM +bx*(vdotbr-vstdotbstr) )/srmsM !eR*

    eststr= estr + sqrt(rhostr)*(vstdotbstr-vststdotbstst)*signBx !eR**

    ff(1) = rhostr*sM
    ff(2) = rhostr*SM**2+pst-bx**2
    ff(3) = rhostr*sM*vstst-bx*bystst
    ff(4) = rhostr*sM*wstst-bx*bzstst
    ff(5) = sm*(eststr+pst)-bx*(vststdotbstst)
    ff(6) = 0.
    ff(7) = bystst*sM-bx*vstst
    ff(8) = bzstst*sM-bx*wstst

#ifdef PASSIVES
    if (npassive >= 1) then
      ff(firstpas:neqtot) = sM*primr(firstpas:neqtot)*srmur/srmsM
    end if
#endif

    return
  endif

  write(logu,'(a,5es12.3)') 'Error in HLLD routine', sM,sl,sr, csl, csr
  call clean_abort(ERROR_NAN_IN_FLUXES)
  !stop

end subroutine primfhlld

!===============================================================================

#endif

end module HLLD
