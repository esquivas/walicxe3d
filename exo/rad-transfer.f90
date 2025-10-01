!===============================================================================
!> @file radTransfer.f90
!> @brief Photoionization radiation transfer module
!> @author A. Esquivel
!> @date May/26/2015

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

!> @brief Computes the radiation transfer of ionizing photons
module radTransfer

  implicit none

  real, parameter   :: a0 = 6.3e-18       !< Fotoionization cross section
  real, parameter   :: F0 = 2.e13         !< ionizing photon Flux [ cm^-2 s^-1 ]
  integer, allocatable :: SortedBlocks(:) !< array of local blocks in Z order

contains

!===============================================================================
!> @brief Initializes radiation transfer module
subroutine initRadTransfer()

  use parameters, only : nbMaxProc
  use globals,    only : rank
  implicit none
  integer :: rand_size
  integer, allocatable, dimension(:) :: rand_seed
  character (len=10) :: system_time
  real :: rtime

  !>  Allocate memory for Sorted Block list
  allocate (SortedBlocks(nbMaxProc))

  !>  Initialize random number generator
  call random_seed(size=rand_size)
  allocate(rand_seed(1:rand_size))
  call date_and_time(time=system_time)
  read(system_time,*) rtime
  rand_seed=int(rtime*1000.0)
#ifdef MPIP
  rand_seed=rand_seed*rank
#endif
  call random_seed(put=rand_seed)
  deallocate(rand_seed)

end subroutine initRadTransfer

!===============================================================================
!> @brief Retuirns the z position of each block reference corner
subroutine getZ(bID, z)

  use parameters, only : ncells_z
  use globals,    only : maxlev
  use amr,        only : bcoords, meshlevel
  implicit none
  integer, intent(in)  :: bID
  integer, intent(out) :: z
  integer              :: ilev, bx, by, bz

   ! Obtain block coords
  call meshlevel (bID, ilev)
  call bcoords(bID, bx, by, bz)

  z = (bz-1)*ncells_z*2**(maxlev-ilev)

  return

end subroutine

!===============================================================================
!> @brief Numerical Recipes heapsort index
subroutine indexx(N,ARRIN,INDX)
  implicit none
  !
  ! Array dimensions
  integer N
  !
  ! Integers
  integer I,INDX(N),INDXT,IR,J,L
  integer ARRIN(N),Q
  !
  ! Code
  forall (J=1:N)
     INDX(J)=J
  end forall
  if (n.le.1) return
  L=N/2+1
  IR=N
  do while (.true.)
     if (L.gt.1) then
        L=L-1
        INDXT=INDX(L)
        Q=ARRIN(INDXT)
     else
        INDXT=INDX(IR)
        Q=ARRIN(INDXT)
        INDX(IR)=INDX(1)
        IR=IR-1
        if (IR.eq.1) then
           INDX(1)=INDXT
           RETURN
        end if
     end if
     I=L
     J=L+L
     do while (J.le.IR)
        if (J.lt.IR) then
           if (ARRIN(INDX(J)).lt.ARRIN(INDX(J+1))) J=J+1
        endif
        if (Q.lt.ARRIN(INDX(J))) then
           INDX(I)=INDX(J)
           I=J
           J=J+J
        else
           J=IR+1
        end if
     end do
     INDX(I)=INDXT
  end do
end subroutine indexx

!===============================================================================
!> @brief Create an index table sorted by the z position
subroutine ZIndex(nbMaxProc, localBlocks, SortedBlocks, Zcoord)

  use amr,        only : bcoords, meshlevel
  implicit none
  integer, intent(in)  :: nbMaxProc
  integer, intent(in)  :: localBlocks(nbMaxProc)
  integer, intent(out) :: SortedBlocks(nbMaxProc)
  integer, intent(out) :: Zcoord(nbMaxProc)
  integer nb, bID
  Zcoord(:) = -1

  !  Store the Z position for all local blocks
  do nb=1, nbMaxProc
    bID = localBlocks(nb)
    if (bID.ne.-1) then
      call getZ(bID, Zcoord(nb) )
      cycle
    end if
  end do

  !  create index table
  call indexx(nbMaxProc, Zcoord, SortedBlocks)

  return

end subroutine Zindex

!===============================================================================
!> @brief Driver module for the radiation transfer
subroutine DoRadTransfer()

  use globals,    only : localBlocks, logu, rank, U, PRIM, dz, maxlev, ierr
  use parameters, only : nbMaxProc, TauRT, ncells_x, ncells_y, ncells_z, nzmax,&
                         l_sc, mpi_real_kind, inH0
  use constants,  only : NEIGH_SAME,NEIGH_COARSER,NEIGH_FINER,NEIGH_BOUNDARY,  &
                         TOP, BOTTOM
  use amr,        only : meshlevel, neighbors, getOwner, absCoords,            &
                         siblingCoords
  use utils,      only : find
  use mpi
  implicit none
  integer :: nb, bID, ilev, iID, neighI, sx, sy, sz
  integer :: i, j, k, ip, jp, i1, j1, i2, j2
  integer :: ntypeB, ntypeT, neighT(4), neighB(4), rankB, rankT
  integer :: Zcoord(nbMaxProc)
  real    :: sendBuf(ncells_x,ncells_y), recvBuff(ncells_x,ncells_y)
  integer, parameter :: nx2 = ncells_x/2, ny2=ncells_y/2
  real    :: sendHalf(nx2,ny2), recvHalf(nx2,ny2)
  integer :: mpistatus(MPI_STATUS_SIZE), mpirequest
  real    :: dtau, dl
  ! DEBUG
  !integer ::  bx, by, bz

  !write(logu,'(a,i0,a,i0)') 'rank: ', rank, ' commencing radiation transfer ...'
  !write(logu,'(a)') ' '
  !write(logu,'(a)') &
  !  '~Active   nb   bID lev    Z |        neighB      :      neighT           |  Coords:'

  !  get index table based on each block position in Z
  call Zindex ( nbMaxProc, localBlocks, sortedBlocks, Zcoord)

  !do nb = 1, nbMaxProc
  !  iID = SortedBlocks(nb)
  !  bID = localBlocks(iID)
  !  if (bID.ne.-1) then
  !    call meshlevel(bID, ilev)
  !    call absCoords( bID, 1, 1, 1, bx, by, bz)
  !    call neighbors(bID, BOTTOM, ntypeB, neighB )
  !    call neighbors(bID, TOP   , ntypeT, neighT )
  !    write(logu,'(a, i5,i6,i3, i6,a,4i5,a,4i5,a,3i4)') &
  !        '~Active', nb, bID, ilev, Zcoord(iID) , ' |',&
  !                                     neighB, ':', neighT,'   |', bx, by, bz
  !  end if
  !end do

  !  main loop
  do nb = 1, nbMaxProc
    iID = SortedBlocks(nb)
    bID = localBlocks(iID)
    if (bID.ne.-1) then

      !  Determine neighbors in the Z direction
      call neighbors(bID, BOTTOM, ntypeB, neighB )
      call neighbors(bID, TOP   , ntypeT, neighT )

      call meshlevel(bID, ilev)
      dl = dz(ilev)*l_sc    !  [ cgs ]

      !-------------------------------------------------------------------------
      !  FIRST RECEIVE, looking at bottom neighbors
      !-------------------------------------------------------------------------
      !-------- DOMAIN BOUNDARY-------------------------------------------------
      if (ntypeB == NEIGH_BOUNDARY) then

        !  Only do in blocks @ z=0
        if (Zcoord(iID).eq.0) then
          U( iID, TauRT, :, :, 0 ) = 0.0 !U( iID, inH0, :,:, 0 ) * a0 * dl! dTau(0)
        end if

      !-------- SAME LEVEL NEIGHBORS -------------------------------------------
      else if (ntypeB == NEIGH_SAME) then

        call getOwner(neighB(1), rankB )

        !  Neighbor block is local & @same level
        if ( rank == rankB )  then

          call find( neighB(1), localBlocks, nbMaxProc, neighI )
          U( iID, TauRT, :, :, 0 ) = U( neighI, TauRT, :,:, ncells_z )

        !  Neighbor @same level is on other procesor
        else

          !write(logu,'(a,i5,a,i2,a,i4)') &
          !'** 1 Ready to recv ', neighB(1), ' from: ',rankB, ' to pair w: ',bID
          call MPI_RECV(recvBuff, ncells_x*ncells_y, mpi_real_kind, rankB,     &
                        neighB(1), mpi_comm_world, mpistatus, ierr )
          U ( iID, TauRT, 1:ncells_x, 1:ncells_y, 0 ) =                        &
                                                 recvBuff(1:ncells_x,1:ncells_y)

        end if

      !-------- HIGHER LEVEL NEIGHBORS -----------------------------------------
      else if (ntypeB == NEIGH_FINER) then

        !  Loop over @higher level neighbors
        do k=1,4

          !  Determine ownership, sibling location and its own block coordinates
          call getOwner(neighB(k), rankB )
          call siblingCoords(neighB(k), sx, sy, sz)

          !  Neighbor is local (@higher-level)
          if(rank == rankB) then

            call find( neighB(k), localBlocks, nbMaxProc, neighI )
            do j=1,ny2
              do i=1,nx2
                i1 = i*2-1       ;   j1 = j*2-1
                ip = i + sx*nx2  ;  jp = j + sy*ny2
                U(iID, TauRT, ip, jp, 0) =                                     &
                !sum( U(neighI, TauRT, i1:i1+1, j1:j1+1, ncells_z) )
                sum( U(neighI, TauRT, i1:i1+1, j1:j1+1, ncells_z) ) / 4.0

              end do
            end do

          !  Neighbor is in other processor (@higher-level)
          else

            i1 = 1 + nx2*sx   ; i2 = nx2*(1+sx)
            j1 = 1 + ny2*sy   ; j2 = ny2*(1+sy)

            !write(logu,'(a,i5,a,i2,a,i4)') &
            !'** 2 Ready to recv ', neighB(k),' from: ',rankB, ' to pair w: ',bID
            call MPI_RECV(recvHalf, nx2*ny2, mpi_real_kind, rankB, neighB(k),&
                          mpi_comm_world, mpistatus, ierr )
            U   (iID, TauRT, i1:i2, j1:j2, 0) = recvHalf(1:nx2, 1:ny2)

          end if
        end do

      !-------- COARSER LEVEL NEIGHBORS ----------------------------------------
      else if (ntypeB == NEIGH_COARSER) then

        !  Determine ownership, sibling location and its own block coordinates
        call getOwner(neighB(1), rankB )
        call siblingCoords(bID, sx, sy, sz)

        !  Neighbor is local (@coarser-level)
        if (rankB == rank) then

          call find( neighB(1), localBlocks, nbMaxProc, neighI )
          do j=1,ncells_y
            do i=1,ncells_x

                ip = (ncells_x/2)*sx + (i+1)/2   ! INT division
                jp = (ncells_y/2)*sy + (j+1)/2   ! INT division
                U(iID, TauRT, i, j, 0) = U(neighI, TauRT, ip, jp, ncells_z )

            end do
          end do

        !  Neighbor in other processor (@coarser-level)
        else

          !write(logu,'(a,i5,a,i2,a,i4)') &
          !'** 3 Ready to recv ', neighB(1), ' from: ',rankB, ' to pair w: ',bID
          call MPI_RECV(recvHalf, nx2*ny2 , mpi_real_kind, rankB, neighB(1),  &
                        mpi_comm_world, mpistatus, ierr )
          do j=1,ncells_y
            do i=1,ncells_x

              ip = (i+1)/2   ! INT division + (ncells_x/2)*sx
              jp = (j+1)/2   ! INT division + (ncells_y/2)*sy
              U(iID, TauRT, i, j, 0) = recvHalf(ip, jp)

            end do
          end do

        end if

      end if

      !-------------------------------------------------------------------------
      !  compute tau on current block
      !-------------------------------------------------------------------------
      do k = 1, ncells_x
        do j = 1, ncells_y
          do i = 1, ncells_z

            dtau = a0 * dl * U(iID, inH0, i, j, k)
            !  U contains Tau, PRIM phi
            U(iID, TauRT, i, j, k) = U(iID, TauRT, i, j, k-1 ) + dtau

            !PRIM(iID, TauRT, i, j, k) = S0*exp(-U(iID, TauRT, i, j, k-1 )  )*  &
            !               (1.0-exp(-dtau))/U(iID, inH0, i, j, k)/dl**3
          end do
        end do
      end do

      !-------------------------------------------------------------------------
      !  After got Tau, look for TOP neighbors to SEND data to other
      !  processors if needed
      !-------------------------------------------------------------------------
      !-------- SAME LEVEL NEIGHBORS -------------------------------------------
      if (ntypeT == NEIGH_SAME) then

        call getOwner(neighT(1), rankT )
        !  Neighbor @same level is on other procesor
        if (rank /= rankT) then

          sendBuf(1:ncells_x,1:ncells_y) = &
                          U ( iID, TauRT, 1:ncells_x, 1:ncells_y, ncells_z )

          call MPI_ISEND( sendBuf, ncells_x*ncells_y, mpi_real_kind, rankT,  &
                          bID, mpi_comm_world, mpirequest, ierr )
          !write(logu,'(a,i5,a,i2)') '1 Just sent ', bID, ' to: ',rankT

        end if

      !-------- COARSER LEVEL NEIGHBOR  ----------------------------------------
      else if (ntypeT == NEIGH_COARSER) then

        call getOwner(neighT(1), rankT )
        if (rank /= rankT) then

          call siblingCoords(neighT(1), sx, sy, sz)
          do j=1,ny2
            do i=1,nx2
              i1 = i*2-1
              j1 = j*2-1
              sendHalf(i,j) =  &
              sum( U(iID, TauRT, i1:i1+1, j1:j1+1, ncells_z) )/ 4.0
            end do
          end do

          call MPI_ISEND( sendHalf, nx2*ny2, mpi_real_kind, rankT, bID,      &
                          mpi_comm_world, mpirequest, ierr )
          !write(logu,'(a,i5,a,i2)') '3 Just sent ', bID, ' to: ',rankT

        end if

      !-------- HIGHER LEVEL NEIGHBORS -----------------------------------------
      else if (ntypeT == NEIGH_FINER) then

        do k =1, 4
          call getOwner(neighT(k), rankT )
          call siblingCoords(neighT(k), sx, sy, sz)
          if (rank /= rankT) then

            i1 = 1 + nx2*sx   ; i2 = nx2*(1+sx)
            j1 = 1 + ny2*sy   ; j2 = ny2*(1+sy)
            sendHalf(1:nx2,1:ny2) = U(iID, TauRT, i1:i2, j1:j2, ncells_z)

            call MPI_ISEND( sendHalf, nx2*ny2, mpi_real_kind, rankT, bID,    &
                            mpi_comm_world, mpirequest, ierr )
            !write(logu,'(a,i5,a,i2)') '2 Just sent ', bID, ' to: ',rankT

          end if
        end do
      end if

      !  synch U / prim
      PRIM(iID,TauRT,:,:,:) = U(iID,TauRT,:,:,:)
    end if
  end do

write(logu,'(a,i0,a)') 'rank: ', rank, ' radiation transfer finished...'

end subroutine DoRadTransfer

!===============================================================================
end module radTransfer