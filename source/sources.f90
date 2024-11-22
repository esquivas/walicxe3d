!===============================================================================
!> @file sources.f90
!> @brief Source terms for the hydrodynamical solver
!> @author Alex Esquivel
!> @date 14/Aug/2018

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

!> @brief Source terms module
!> @details The module applies most source terms after the hydro (or MHD)
!! timestep
module sources

  implicit none

contains

  !=======================================================================
#ifdef BFIELD

  subroutine divergence_B(locIndx,lev,i,j,k,div)

    !> @brief Computes div(B)
    !> @details Computes div(B)
    !> @param integer [in] i : cell index in the X direction
    !> @param integer [in] j : cell index in the Y direction
    !> @param integer [in] k : cell index in the Z direction
    !> @param real [out] d :: div(B)
    use globals
    implicit none
    integer, intent(in) :: locIndx,lev,i,j,k
    real, intent(out)   :: div

    div=  (prim(locIndx,6,i+1,j,k)-prim(locIndx,6,i-1,j,k))/(2.*dx(lev) )  &
        + (prim(locIndx,7,i,j+1,k)-prim(locIndx,7,i,j-1,k))/(2.*dy(lev) )  &
        + (prim(locIndx,8,i,j,k+1)-prim(locIndx,8,i,j,k-1))/(2.*dz(lev) )

  end subroutine divergence_B

#endif

  !=======================================================================
#ifdef BFIELD

  subroutine divbcorr_8w_source(locIndx,lev,i,j,k,s)

    !> @brief 8 Wave source terms for div(B) correction
    !> @details  Adds terms proportional to div B in Faraday's Law,
    !! momentum equation and energy equation as propoes in Powell et al. 1999
    !> @param integer [in] i : cell index in the X direction
    !> @param integer [in] j : cell index in the Y direction
    !> @param integer [in] k : cell index in the Z direction
    !> @param real [out] s(neq) : vector with source terms
    use parameters, only : neqtot
    use globals,    only : prim
    implicit none
    integer, intent(in) :: locIndx, lev, i, j, k
    real, intent(inout) :: s(neqtot)
    real                :: divB, pp(neqtot)

    call divergence_B(locIndx,lev,i,j,k,divB)

    pp(:) = prim(locIndx,:,i,j,k)
    ! update source terms
    ! momenta
    s(2)= s(2)-divB*pp(6)
    s(3)= s(3)-divB*pp(7)
    s(4)= s(4)-divB*pp(8)

    ! energy
    s(5)= s(5)-divB*(pp(2)*pp(6)+pp(3)*pp(7)+pp(4)*pp(8))

    ! Faraday law
    s(6)=s(6)-divB*pp(2)
    s(7)=s(7)-divB*pp(3)
    s(8)=s(8)-divB*pp(4)

  end subroutine divbcorr_8w_source

#endif

  !=======================================================================
  subroutine source_function(locIndx,lev,i,j,k,s)

    !> @brief Upper level wrapper for sources
    !> @details Upper level wrapper for sources
    !! @n Main driver, this is called from the upwind stepping
    !> @param integer [in] locIndx : local index of current block
    !> @param integer [in] lev     : level of refinement of currrent block
    !> @param integer [in] i       : cell index in the X direction
    !> @param integer [in] j       : cell index in the Y direction
    !> @param integer [in] k       : cell index in the Z direction
    !> @param real [in] prim(neq) : vector of primitive variables
    !> @param real [out] s(neq) : vector with source terms
    use parameters, only : neqtot, user_source_terms, eight_wave
    use userconds, only : get_user_source_terms
    use globals,   only : PRIM
    implicit none
    integer, intent(in)  :: locIndx, lev, i, j, k
    real, intent(out)    :: s(neqtot)
    !real :: x, y, z, r

    ! resets the source terms
    s(:) = 0.

    !  user source terms (such as gravity, tidal or inertial forces)
    if (user_source_terms) call get_user_source_terms(PRIM(locIndx,:,i,j,k),s,i,j,k)

#ifdef BFIELD
    !  divergence correction Powell et al. 1999
    if (eight_wave) call divbcorr_8w_source(locIndx,lev,i,j,k,s)
#endif

    return

  end subroutine source_function

!===============================================================================

end module sources
