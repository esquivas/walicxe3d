!===============================================================================
!> @file flux_correction.f90
!> @brief Flux correction at fine-coarse block interfaces
!> @author A. Esquivel
!> @date 06/Oct/2026

!> @brief Flux correction (refluxing) module
!> @details All levels advance with the same global dt, so conservation across
!! fine-coarse interfaces is restored by replacing, in the coarse cells next to
!! a finer neighbor, the coarse face flux with the area average of the four
!! fine face fluxes. Face fluxes are stored during the final (2nd order) stage
!! with store_face_fluxes; apply_flux_correction then corrects UP once all
!! blocks have been updated (before UP is copied into U).
module flux_correction

  implicit none

  !> Face fluxes of each local block: (neqtot, n, n, face, local index)
  !!  face with a COARSER neighbor: 2x2-averaged fine fluxes in (1:n/2,1:n/2)
  !!  face with a FINER neighbor  : own (coarse) fluxes in (1:n,1:n)
  real, allocatable, save :: faceflux(:,:,:,:,:)

contains

!===============================================================================

!> @brief Stores the face fluxes of one block at fine-coarse interfaces
!> @details Must be called right after the block's fluxes (FC, GC, HC) are
!! computed and before the next block overwrites them.
!> @param locIndx Local index of the block
subroutine store_face_fluxes (locIndx)

  use parameters
  use globals, only : localBlocks, FC, GC, HC
  use amr,     only : neighbors
  implicit none

  integer, intent(in) :: locIndx

  integer :: bID, dir, ntype, neighs(4), a, b, n, nh
  real    :: face(neqtot, ncells_x, ncells_x)

  n  = ncells_x
  nh = n/2

  if (.not.allocated(faceflux)) then
    allocate( faceflux(neqtot, n, n, 6, nbMaxProc) )
    faceflux(:,:,:,:,:) = 0.0
  end if

  bID = localBlocks(locIndx)

  do dir=1,6

    call neighbors (bID, dir, ntype, neighs)
    if ((ntype /= NEIGH_COARSER).and.(ntype /= NEIGH_FINER)) cycle

    ! Fluxes on this face, transverse indices in ascending axis order
    select case (dir)
    case (LEFT)
      face(:,:,:) = FC(:,0,1:n,1:n)
    case (RIGHT)
      face(:,:,:) = FC(:,n,1:n,1:n)
    case (FRONT)
      face(:,:,:) = GC(:,1:n,0,1:n)
    case (BACK)
      face(:,:,:) = GC(:,1:n,n,1:n)
    case (BOTTOM)
      face(:,:,:) = HC(:,1:n,1:n,0)
    case (TOP)
      face(:,:,:) = HC(:,1:n,1:n,n)
    end select

    if (ntype == NEIGH_COARSER) then
      ! Fine side: average each 2x2 group of fine faces (= one coarse face)
      do b=1,nh
        do a=1,nh
          faceflux(:,a,b,dir,locIndx) = 0.25*( face(:,2*a-1,2*b-1)            &
                                             + face(:,2*a  ,2*b-1)            &
                                             + face(:,2*a-1,2*b  )            &
                                             + face(:,2*a  ,2*b  ) )
        end do
      end do
    else
      ! Coarse side: keep own fluxes
      faceflux(:,:,:,dir,locIndx) = face(:,:,:)
    end if

  end do

end subroutine store_face_fluxes

!===============================================================================

!> @brief Corrects UP in coarse cells adjacent to finer neighbors
!> @details Must be called by ALL ranks. Every rank sweeps all blocks in the
!! same order (as in boundary.f90), so blocking send/recv pairs match.
subroutine apply_flux_correction ()

  use parameters
  use globals,    only : globalBlocks, localBlocks, UP, dt, dx, dy, dz,       &
                         rank, ierr, logu
  use amr,        only : neighbors, getOwner, siblingCoords, meshlevel
  use boundaries, only : opposite
  use utils,      only : find
  use clean_quit, only : clean_abort
  implicit none

  integer :: dir, oppdir, nb, b, destID, destOwner, destInd
  integer :: srcID, srcOwner, srcInd, ntype, neighs(4)
  integer :: sx, sy, sz, oa, ob, a, c, ia, ib, lev, n, nh, nData
  integer :: mpistatus(MPI_STATUS_SIZE)
  real    :: buf(neqtot, ncells_x/2, ncells_x/2), dF(neqtot), dtdl

  n     = ncells_x
  nh    = n/2
  nData = neqtot*nh*nh
  oa    = 0
  ob    = 0
  dtdl  = 0.0

  ! dir is the face of the (coarse) destination block
  do dir=1,6
    call opposite (dir, oppdir)

    do nb=1,nbMaxGlobal
      destID = globalBlocks(nb)
      if (destID == -1) cycle

      call neighbors (destID, dir, ntype, neighs)
      if (ntype /= NEIGH_FINER) cycle

      call getOwner (destID, destOwner)

      do b=1,4
        srcID = neighs(b)
        call getOwner (srcID, srcOwner)
        if (srcOwner == MPI_PROC_NULL) then
          write(logu,'(a,i0,a,i0)') "Flux correction: fine block ", srcID,    &
                                    " not found for block ", destID
          call clean_abort (ERROR_GENERIC)
        end if

        ! Sender: owns the fine block but not the coarse one
        if ((rank == srcOwner).and.(rank /= destOwner)) then
          call find (srcID, localBlocks, nbMaxProc, srcInd)
          buf(:,:,:) = faceflux(:,1:nh,1:nh,oppdir,srcInd)
          call MPI_SEND (buf, nData, mpi_real_kind, destOwner, srcID,         &
                         mpi_comm_world, ierr)
        end if

        ! Receiver: owns the coarse block
        if (rank == destOwner) then

          if (srcOwner == rank) then
            call find (srcID, localBlocks, nbMaxProc, srcInd)
            buf(:,:,:) = faceflux(:,1:nh,1:nh,oppdir,srcInd)
          else
            call MPI_RECV (buf, nData, mpi_real_kind, srcOwner, srcID,        &
                           mpi_comm_world, mpistatus, ierr)
          end if

          call find (destID, localBlocks, nbMaxProc, destInd)
          call meshlevel (destID, lev)
          call siblingCoords (srcID, sx, sy, sz)

          ! Offset of the fine block's quadrant on the coarse face
          select case (dir)
          case (LEFT, RIGHT)
            oa = sy*nh ; ob = sz*nh ; dtdl = dt/dx(lev)
          case (FRONT, BACK)
            oa = sx*nh ; ob = sz*nh ; dtdl = dt/dy(lev)
          case (BOTTOM, TOP)
            oa = sx*nh ; ob = sy*nh ; dtdl = dt/dz(lev)
          end select

          do c=1,nh
            do a=1,nh
              ia = oa + a
              ib = ob + c
              ! (coarse flux) - (averaged fine flux)
              dF(:) = faceflux(:,ia,ib,dir,destInd) - buf(:,a,c)
              select case (dir)
              case (LEFT)
                UP(destInd,:,1,ia,ib) = UP(destInd,:,1,ia,ib) - dtdl*dF(:)
              case (RIGHT)
                UP(destInd,:,n,ia,ib) = UP(destInd,:,n,ia,ib) + dtdl*dF(:)
              case (FRONT)
                UP(destInd,:,ia,1,ib) = UP(destInd,:,ia,1,ib) - dtdl*dF(:)
              case (BACK)
                UP(destInd,:,ia,n,ib) = UP(destInd,:,ia,n,ib) + dtdl*dF(:)
              case (BOTTOM)
                UP(destInd,:,ia,ib,1) = UP(destInd,:,ia,ib,1) - dtdl*dF(:)
              case (TOP)
                UP(destInd,:,ia,ib,n) = UP(destInd,:,ia,ib,n) + dtdl*dF(:)
              end select
            end do
          end do

        end if
      end do
    end do
  end do

   !write(logu,*) 'applied flux correction', dtdl, df(1)

end subroutine apply_flux_correction

!===============================================================================

end module flux_correction
