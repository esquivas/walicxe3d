!===============================================================================
!> @file flux_cd.f90
!> @brief Field CD Module
!> @author Alex Esquivel
!> @date 14/May/2019

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

!> @brief Field CD Module
!> @details Implements the Field-Interpolated centrered
!! diferences method to update the B field to enforce
!! the div(B)=0 constrain (Toth 2000)

module flux_cd

#ifdef BFIELD

use parameters, only : ncells_x, ncells_y, ncells_z
implicit none
real :: efield(3,0:ncells_x+1,0:ncells_y+1,0:ncells_z+1)


contains

!=======================================================================

!> @brief Upper level wrapper for field-CD update
!> @details Upper level wrapper for field-CD, updates the
!> hydro variables with upwind scheme and the field as field-CD
!> @param integer [in] i : cell index in the X direction
!> @param integer [in] j : cell index in the Y direction
!> @param integer [in] k : cell index in the Z direction
!> @param real [in] dt : timestep

subroutine flux_cd_update(locIndx,i,j,k,dtdx,dtdy,dtdz)

  use parameters, only : neqhydro,firstpas, npassive
  use globals, only : up, u, fc, gc, hc
  implicit none
  integer, intent(in)  :: locIndx, i, j, k
  real, intent (in)    :: dtdx, dtdy, dtdz

  !   hydro variables (and passive eqs.)
  up(locIndx,:neqhydro,i,j,k)= u(locIndx,:neqhydro,i,j,k)             &
                   - dtdx*(fc(:neqhydro,i,j,k)-fc(:neqhydro,i-1,j,k)) &
                   - dtdy*(gc(:neqhydro,i,j,k)-gc(:neqhydro,i,j-1,k)) &
                   - dtdz*(hc(:neqhydro,i,j,k)-hc(:neqhydro,i,j,k-1))

#ifdef PASSIVES
  if (npassive >= 1) &
  up(locIndx,firstpas:,i,j,k)=u(locIndx,firstpas:,i,j,k)      &
           - dtdx*(fc(firstpas:,i,j,k)-fc(firstpas:,i-1,j,k)) &
           - dtdy*(gc(firstpas:,i,j,k)-gc(firstpas:,i,j-1,k)) &
           - dtdz*(hc(firstpas:,i,j,k)-hc(firstpas:,i,j,k-1))
#endif
  ! evolution of B with field-CD
  up(locIndx,6,i,j,k)=u(locIndx,6,i,j,k)                    &
         - 0.5*dtdy* ( efield(3,i,j+1,k)-efield(3,i,j-1,k)) &
         - 0.5*dtdz* (-efield(2,i,j,k+1)+efield(2,i,j,k-1))

   up(locIndx,7,i,j,k)=u(locIndx,7,i,j,k)                   &
         - 0.5*dtdx* (-efield(3,i+1,j,k)+efield(3,i-1,j,k)) &
         - 0.5*dtdz* ( efield(1,i,j,k+1)-efield(1,i,j,k-1))

   up(locIndx,8,i,j,k)=u(locIndx,8,i,j,k)                   &
         - 0.5*dtdx* ( efield(2,i+1,j,k)-efield(2,i-1,j,k)) &
         - 0.5*dtdy* (-efield(1,i,j+1,k)+efield(1,i,j-1,k))

end subroutine flux_cd_update

!===============================================================================

#endif

end module flux_cd
