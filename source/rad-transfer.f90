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

  use mpi
  implicit none

  real, parameter   :: a0 = 6.3e-18       !< Fotoionization cross section
  real, parameter   :: F0 = 2.e13         !< ionizing photon Flux [ cm^-2 s^-1 ]
  integer, allocatable :: SortedBlocks(:) !< array of local blocks in Z order

  !> MPI send and receive buffers
  character(len=:), allocatable :: sendBuf, sendBuf2, recvBuf, recvBuf2
  integer                       :: bufSize, bufSize2

contains

!===============================================================================
!> @brief Initializes radiation transfer module
subroutine initRadTransfer()

  use parameters, only : nbMaxProc, ncells_x, ncells_y
  use globals,    only : rank
  implicit none
  integer, parameter :: nx2 = ncells_x/2, ny2=ncells_y/2
  integer :: ps, ierr
  integer :: rand_size
  integer, allocatable, dimension(:) :: rand_seed
  character (len=10) :: system_time
  real :: rtime

  !>  Allocate memory for Sorted Block list
  allocate (SortedBlocks(nbMaxProc))

  ! ---- compute packed buffer size ----
  ! entire layer buffered
  bufsize = 0
  call MPI_Pack_size(1, MPI_INTEGER, mpi_comm_world, ps, ierr)
  bufsize = bufsize + ps
  call MPI_Pack_size(1, MPI_INTEGER, mpi_comm_world, ps, ierr)
  bufsize = bufsize + ps
    call MPI_Pack_size(1, MPI_INTEGER, mpi_comm_world, ps, ierr)
  bufsize = bufsize + ps
  call MPI_Pack_size(ncells_x * ncells_y, MPI_DOUBLE_PRECISION, mpi_comm_world,&
                     ps, ierr)
  bufsize = bufsize + ps
! quart of layer buffered
  bufsize2 = 0
  call MPI_Pack_size(1, MPI_INTEGER, mpi_comm_world, ps, ierr)
  bufsize2 = bufsize2 + ps
  call MPI_Pack_size(1, MPI_INTEGER, mpi_comm_world, ps, ierr)
  bufsize2 = bufsize2 + ps
    call MPI_Pack_size(1, MPI_INTEGER, mpi_comm_world, ps, ierr)
  bufsize2 = bufsize2 + ps
  call MPI_Pack_size(    nx2  *  ny2,     MPI_DOUBLE_PRECISION, mpi_comm_world,&
                     ps, ierr)
  bufsize2 = bufsize2 + ps

  !> Allocate buffers for communication
  allocate(character( len=bufsize  ) :: sendBuf  )
  allocate(character( len=bufsize  ) :: recvBuf  )
  allocate(character( len=bufsize2 ) :: sendBuf2 )
  allocate(character( len=bufsize2 ) :: recvBuf2 )

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
!> @brief receive non-local ray package
subroutine get_rays(rankB)

  use globals,    only : localBlocks, logu, U, dz, logu
  use parameters, only : nbMaxProc, TauRT, ncells_x, ncells_y, ncells_z,       &
                         mpi_real_kind, l_sc, inH0
  use constants,  only : NEIGH_SAME,NEIGH_COARSER,NEIGH_FINER,NEIGH_BOUNDARY,  &
                         TOP, BOTTOM
  use amr,        only : meshlevel, neighbors, siblingCoords
  use utils,      only : find
  implicit none
  integer, intent(in ) :: rankB
  integer, parameter :: nx2 = ncells_x/2, ny2=ncells_y/2
  integer :: i, j, k, ip, jp, i1, j1
  integer :: mpistatus(MPI_STATUS_SIZE), ierr, count, ilev
  integer :: bIDs, srcID, dstID, pos, ntypeT, neighI, sx, sy, sz
  real    :: recvData(ncells_x,ncells_y), recvData2(nx2,ny2), dtau, dl

! ---- find size of packed message ----
call MPI_Probe(rankB, MPI_ANY_TAG, mpi_comm_world, mpistatus, ierr)
bIDs  = mpistatus(MPI_TAG)  ! bID of source block
call MPI_Get_count(mpistatus, MPI_PACKED, count, ierr)

!-----------------------------------------------------------------------------
if (count == bufSize) then  ! receiving full layer

  call MPI_Recv(recvBuf, count, MPI_PACKED, rankB, bIDs, mpi_comm_world,       &
                mpistatus, ierr)
  ! ---- unpack ----
  pos = 0
  call MPI_Unpack(recvBuf, count, pos, srcID, 1, MPI_INTEGER,                  &
                  mpi_comm_world, ierr)
  call MPI_Unpack(recvBuf, count, pos, dstID, 1, MPI_INTEGER,                  &
                  mpi_comm_world, ierr)
  call MPI_Unpack(recvBuf, count, pos, ntypeT, 1, MPI_INTEGER,                 &
                  mpi_comm_world, ierr)

  ! this occurs only when neighbors are @same level
  call MPI_Unpack(recvBuf, count, pos, recvData, ncells_x*ncells_y,            &
                  MPI_DOUBLE_PRECISION, mpi_comm_world, ierr)

  call find( dstID, localBlocks, nbMaxProc, neighI )

  U(neighI, TauRT, 1:ncells_x, 1:ncells_y, 0 ) = recvData(1:ncells_x,1:ncells_y)

  !write(logu,'(a,i5,a,i5,a,i5,a,i5,a,i5)') &
  ! '**1** received', count, ' from rank  ', rankB, ' bID:', bIDs , ' from ', srcID, ' to ', dstID

   call meshlevel(dstID, ilev)
   dl = dz(ilev)*l_sc
  ! compute tau
  do k = 1, ncells_z
    do j = 1, ncells_y
      do i = 1, ncells_x
        dtau = a0 * dl * U(neighI, inH0, i, j, k)
        U(neighI, TauRT, i, j, k) = U(neighI, TauRT, i, j, k-1 ) + dtau
      end do
    end do
  end do

  return
  !-----------------------------------------------------------------------------
else if (count == bufSize2) then ! receiving 1/4 of layer

  call MPI_Recv(recvBuf2, count, MPI_PACKED, rankB, bIDs, mpi_comm_world,      &
                mpistatus, ierr)
  ! ---- unpack ----
  pos = 0
  call MPI_Unpack(recvBuf2, count, pos, srcID, 1, MPI_INTEGER,                 &
                  mpi_comm_world, ierr)
  call MPI_Unpack(recvBuf2, count, pos, dstID, 1, MPI_INTEGER,                 &
                  mpi_comm_world, ierr)
  call MPI_Unpack(recvBuf2, count, pos, ntypeT, 1, MPI_INTEGER,                 &
                  mpi_comm_world, ierr)

  ! this occurs only when neighbors are @same level
  call MPI_Unpack(recvBuf2, count, pos, recvData2, nx2*ny2,                    &
                  MPI_DOUBLE_PRECISION, mpi_comm_world, ierr)

  call find( dstID, localBlocks, nbMaxProc, neighI )

  if (ntypeT == NEIGH_COARSER) then

    call siblingCoords(bIDs, sx, sy, sz)
    i1 = sx*nx2 + 1      ;   j1 = sy*ny2 + 1
    ip = (sx + 1)*nx2    ;   jp = (sy + 1)*ny2
    U(neighI, TauRT, i1:ip , j1:jp, 0 ) = recvData2( 1:nx2, 1:ny2)

    call meshlevel(dstID, ilev)
    dl = dz(ilev)*l_sc
    ! compute tau
    do k = 1, ncells_z
      do j = 1, ny2
        do i = 1, nx2
          ip = i + sx*nx2  ;  jp =  j  + sy*ny2
          dtau = a0 * dl * U(neighI, inH0, ip, jp, k)
          U(neighI, TauRT, ip, jp, k) = U(neighI, TauRT, ip, jp, k-1 ) + dtau
        end do
      end do
    end do

    !write(logu,'(a,i5,a,i5,a,i5,a,i5,a,i5,a,3i3)') &
    !  '*3* received', count, ' from rank  ', rankB, ' bID:', bIDs , ' from ',  &
    !  srcID, ' to ', dstID, ' | ', sx, sy, sz

  else if (ntypeT == NEIGH_FINER) then

    do j=1,ncells_y
      do i=1,ncells_x
        ip = (i+1)/2
        jp = (j+1)/2
        U(neighI, TauRT, i, j, 0) = recvData2(ip, jp)
      end do
    end do

    call meshlevel(dstID, ilev)
    dl = dz(ilev)*l_sc
    ! compute tau
    do k = 1, ncells_z
      do j = 1, ncells_y
        do i = 1, ncells_x
        dtau = a0 * dl * U(neighI, inH0, i, j, k)
        U(neighI, TauRT, i, j, k) = U(neighI, TauRT, i, j, k-1 ) + dtau
        end do
      end do
    end do

    !write(logu,'(a,i5,a,i5,a,i5,a,i5,a,i5,a,3i3)') &
    !  '*2* received', count, ' from rank  ', rankB, ' bID:', bIDs , ' from ',  &
    !  srcID, ' to ', dstID, ' | ', sx, sy, sz, ' | ', sx, sy, sz

  end if

  return
end if

write(logu,*) ' ERROR HERE', rankB, bIDS, dstID

end subroutine get_rays

!===============================================================================
!> @brief Driver module for the radiation transfer
subroutine DoRadTransfer()

  use globals,    only : localBlocks, logu, rank, U, PRIM, dz, maxlev, ierr
  use parameters, only : nbMaxProc, TauRT, ncells_x, ncells_y, ncells_z, nzmax,&
                         l_sc, mpi_real_kind, inH0
  use constants,  only : NEIGH_SAME,NEIGH_COARSER,NEIGH_FINER,NEIGH_BOUNDARY,  &
                         TOP, BOTTOM
  use amr,        only : meshlevel, neighbors, getOwner, absCoords,            &
                         siblingCoords, getRefCorner
  use utils,      only : find
  implicit none
  integer :: nb, bID, ilev, iID, neighI, sx, sy, sz, mpirequest
  integer :: i, j, k, l, ip, jp, i1, j1, i2, j2, pos
  integer :: ntypeB, ntypeT, neighT(4), neighB(4), rankB, rankT
  integer :: Zcoord(nbMaxProc)
  integer, parameter :: nx2 = ncells_x/2, ny2=ncells_y/2
  real    :: sendData(ncells_x,ncells_y), sendData2(nx2,ny2)
  real    :: dtau, dl

  !  get index table based on each block position in Z
  call Zindex ( nbMaxProc, localBlocks, sortedBlocks, Zcoord)

  ! !  Report active block info in rank
  ! do nb = 1, nbMaxProc
  !   iID = SortedBlocks(nb)
  !   bID = localBlocks(iID)
  !
  !   !rayno(iID)=0
  !   if (bID.ne.-1) then
  !   call getrefcorner(bID,xgb,ygb,zgb)
  !   call meshlevel(bID, ilev)
  !
  !     write( logu, '(a,i5, a, i5,a,i5, a, i5, a, 3f7.2,a,i0)')                 &
  !     'rank: ', rank, ' nb ', nb, ' iID: ', iiD, ' bID: ', bID, ' | ',         &
  !      xgb, ygb, zgb, ' | level: ', ilev
  !
  !   end if
  ! end do
  ! write(logu,'(a)') '--------'

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

          ! compute tau
          do k = 1, ncells_z
            do j = 1, ncells_y
              do i = 1, ncells_x
                dtau = a0 * dl * U(iID, inH0, i, j, k)
                U(iID, TauRT, i, j, k) = U(iID, TauRT, i, j, k-1 ) + dtau
              end do
            end do
          end do

        end if

      !-------- SAME LEVEL NEIGHBORS -------------------------------------------
      else if (ntypeB == NEIGH_SAME) then

        call getOwner(neighB(1), rankB )

        !  Neighbor block is local & @same level
        if ( rank == rankB )  then

          call find( neighB(1), localBlocks, nbMaxProc, neighI )
          U( iID, TauRT, :, :, 0 ) = U( neighI, TauRT, :,:, ncells_z )

          ! compute tau
          do k = 1, ncells_z
            do j = 1, ncells_y
              do i = 1, ncells_x
                dtau = a0 * dl * U(iID, inH0, i, j, k)
                U(iID, TauRT, i, j, k) = U(iID, TauRT, i, j, k-1 ) + dtau
              end do
            end do
          end do

        !  Neighbor @same level is on other procesor
        else

          call get_rays(rankB)

        end if

      !-------- HIGHER LEVEL NEIGHBORS -----------------------------------------
      else if (ntypeB == NEIGH_FINER) then

        !  Loop over 4 higher level neighbors
        do l=1,4

          !  Determine ownership, sibling location and its own block coordinates
          call getOwner(neighB(l), rankB )

          !  Neighbor is local (@higher-level)
          if(rank == rankB) then

            call find( neighB(l), localBlocks, nbMaxProc, neighI )
            call siblingCoords(neighB(l), sx, sy, sz)

            do j=1,ny2
              do i=1,nx2
                i1 = i*2-1       ;  j1 = j*2-1
                ip = i + sx*nx2  ;  jp =  j  + sy*ny2
                U(iID, TauRT, ip, jp, 0) =                                     &
                      sum( U(neighI, TauRT, i1:i1+1, j1:j1+1, ncells_z) ) / 4.0
              end do
            end do

            ! compute tau
            do k = 1, ncells_z
              do j = 1, ny2
                do i = 1, nx2
                  ip = i + sx*nx2  ;  jp =  j  + sy*ny2
                  dtau = a0 * dl * U(iID, inH0, ip, jp, k)
                  U(iID, TauRT, ip, jp, k) = U(iID, TauRT, ip, jp, k-1 ) + dtau
                end do
              end do
            end do

          !  Neighbor is in other processor (@higher-level)
          else

            call get_rays(rankB)

          end if

        end do

      !-------- COARSER LEVEL NEIGHBORS ----------------------------------------
      else if (ntypeB == NEIGH_COARSER) then

        !  Determine ownership, sibling location and its own block coordinates
        call getOwner(neighB(1), rankB )

        !  Neighbor is local (@coarser-level)
        if (rankB == rank) then

          call siblingCoords(bID, sx, sy, sz)
          call find( neighB(1), localBlocks, nbMaxProc, neighI )

          do j=1,ncells_y
            do i=1,ncells_x
                ip = nx2*sx + (i+1)/2
                jp = ny2*sy + (j+1)/2
                U(iID, TauRT, i, j, 0) = U(neighI, TauRT, ip, jp, ncells_z )
            end do
          end do

          ! compute tau
          do k = 1, ncells_z
            do j = 1, ncells_y
              do i = 1, ncells_x
                dtau = a0 * dl * U(iID, inH0, i, j, k)
                U(iID, TauRT, i, j, k) = U(iID, TauRT, i, j, k-1 ) + dtau
              end do
            end do
          end do

        !  Neighbor in other processor (@coarser-level)
        else

          call get_rays(rankB)

        end if

      end if

      !-------------------------------------------------------------------------
      !  After got Tau, look for TOP neighbors to SEND data to other
      !  processors if needed
      !-------------------------------------------------------------------------

      !-------- SAME LEVEL NEIGHBORS -------------------------------------------
      if (ntypeT == NEIGH_SAME) then

        call getOwner(neighT(1), rankT )
        !  Neighbor @same level is on other procesor
        if (rank /= rankT) then

          sendData(1:ncells_x,1:ncells_y) = &
                          U ( iID, TauRT, 1:ncells_x, 1:ncells_y, ncells_z )

          ! --- pack data ----
          pos = 0
          call MPI_Pack(bID,       1, MPI_INTEGER, sendBuf, bufsize, pos,      &
                        mpi_comm_world, ierr)
          call MPI_Pack(neighT(1), 1, MPI_INTEGER, sendBuf, bufsize, pos,      &
                        mpi_comm_world, ierr)
          call MPI_Pack(ntypeT,    1, MPI_INTEGER, sendBuf, bufsize, pos,      &
                        mpi_comm_world, ierr)

          call MPI_Pack(sendData, ncells_x*ncells_y , MPI_DOUBLE_PRECISION,    &
                        sendBuf, bufsize, pos, mpi_comm_world, ierr)

          call MPI_ISEND( sendBuf, pos, MPI_PACKED, rankT, bID,                &
                          mpi_comm_world, mpirequest, ierr )

          !write(logu,'(a,i5,a,i2,a,i6)')                                       &
          !'*1*  SAME  Just sent ', bID, ' to: ', rankT, ' to pair with:', neighT(1)

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
              sendData2(i,j) =  &
              sum( U(iID, TauRT, i1:i1+1, j1:j1+1, ncells_z) )/ 4.0
            end do
          end do

          ! --- pack data ----
          pos = 0
          call MPI_Pack(bID,       1, MPI_INTEGER, sendBuf2, bufsize, pos,     &
                        mpi_comm_world, ierr)
          call MPI_Pack(neighT(1), 1, MPI_INTEGER, sendBuf2, bufsize, pos,     &
                        mpi_comm_world, ierr)
          call MPI_Pack(ntypeT,    1, MPI_INTEGER, sendBuf2, bufsize, pos,     &
                        mpi_comm_world, ierr)

          call MPI_Pack(sendData2, nx2*ny2 , MPI_DOUBLE_PRECISION,  sendBuf2,  &
                        bufsize2, pos, mpi_comm_world, ierr)

          call MPI_ISEND( sendBuf2, pos, MPI_PACKED, rankT, bID,               &
                          mpi_comm_world, mpirequest, ierr )

          !write(logu,'(a,i5,a,i2,a,i6)')                                       &
          !'*3* COARSE Just sent ', bID, ' to: ', rankT, ' to pair with:', neighT(1)

        end if

      !-------- HIGHER LEVEL NEIGHBORS -----------------------------------------
      else if (ntypeT == NEIGH_FINER) then

        do k =1, 4
          call getOwner(neighT(k), rankT )
          call siblingCoords(neighT(k), sx, sy, sz)
          if (rank /= rankT) then

            i1 = 1 + nx2*sx   ; i2 = nx2*(1+sx)
            j1 = 1 + ny2*sy   ; j2 = ny2*(1+sy)
            sendData2(1:nx2,1:ny2) = U(iID, TauRT, i1:i2, j1:j2, ncells_z)

            ! --- pack data ----
            pos = 0
            call MPI_Pack(bID,       1, MPI_INTEGER, sendBuf2, bufsize, pos,   &
                          mpi_comm_world, ierr)
            call MPI_Pack(neighT(k), 1, MPI_INTEGER, sendBuf2, bufsize, pos,   &
                          mpi_comm_world, ierr)
            call MPI_Pack(ntypeT,    1, MPI_INTEGER, sendBuf2, bufsize, pos,    &
                        mpi_comm_world, ierr)

            call MPI_Pack(sendData2, nx2*ny2 , MPI_DOUBLE_PRECISION,  sendBuf2,&
                          bufsize2, pos, mpi_comm_world, ierr)

            call MPI_ISEND( sendBuf2, pos, MPI_PACKED, rankT, bID,             &
                          mpi_comm_world, mpirequest, ierr )

            !write(logu,'(a,i5,a,i2,a,i6,a,f)')                                 &
            !'*2*  FINE  Just sent ', bID, ' to: ', rankT, ' to pair with',     &
            ! neighT(K), ' | ', sendData2(2,2)

          end if
        end do
      end if

      !  synch U / prim
      PRIM(iID,TauRT,:,:,:) = U(iID,TauRT,:,:,:)
    end if
  end do

!write(logu,'(a,i0,a)') 'rank: ', rank, ' radiation transfer finished...'

end subroutine DoRadTransfer

!===============================================================================
end module radTransfer
