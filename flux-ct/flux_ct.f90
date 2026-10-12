!===============================================================================
!> @file flux_ct.f90
!> @brief Constrained transport (flux-CT) of the face-centred magnetic field
!> @author A. Esquivel
!> @date 12/Oct/2026

!> @brief Flux-CT module
!> @details Evolves the face-centred field BF with constrained transport
!! (Evans & Hawley 1988), using edge EMFs built from the Godunov fluxes with
!! the upwind averaging of Gardiner & Stone (2005, JCP 205, 509), as in
!! Athena (Stone et al. 2008, sect. 5.3). The face-based div(B) stays at
!! roundoff, and U/UP(6:8) is always the average of the two faces.
!! Per block and stage, right after its fluxes and upwindStep:
!!   flux_ct_update(locIndx, dtp) : EMFs, BFP = BF + CT update, UP(6:8)
!! After the final stage, for all local blocks:
!!   flux_ct_advance()            : BF = BFP
!! Needs edge/corner ghost cells and the fluxes on one extra transverse
!! layer (see HLLDfluxes).
module flux_ct

#ifdef BFIELD

  implicit none

  !> Edge EMFs of the current block, "right" indexing like the fluxes:
  !! EX(i,j,k) at (i,j+1/2,k+1/2), EY(i,j,k) at (i+1/2,j,k+1/2),
  !! EZ(i,j,k) at (i+1/2,j+1/2,k)
  real, allocatable, save :: EX(:,:,:), EY(:,:,:), EZ(:,:,:)
  !> Cell-centred reference field E = -v x B of the current block
  real, allocatable, save :: ER(:,:,:,:)
  !> Largest difference found between the two copies of a shared block face
  !! (set by boundary_bf; should be exactly 0)
  real, save :: bf_shared_err = 0.0

contains

!===============================================================================

!> @brief EMFs, CT update of BFP and cell-centred B in UP for one block
!> @param locIndx Local index of the block
!> @param dtp Timestep of this stage (as passed to upwindStep)
subroutine flux_ct_update (locIndx, dtp)

  use parameters
  use globals, only : localBlocks, PRIM, FC, GC, HC, BF, BFP, UP, dx, dy, dz
  use amr,     only : meshlevel
  implicit none

  integer, intent(in) :: locIndx
  real,    intent(in) :: dtp

  integer :: i, j, k, n, lev
  real    :: dtdx, dtdy, dtdz

  n = ncells_x
  if (.not.allocated(EX)) then
    allocate( EX(1:n,0:n,0:n), EY(0:n,1:n,0:n), EZ(0:n,0:n,1:n) )
    allocate( ER(3,0:n+1,0:n+1,0:n+1) )
  end if

  call meshlevel (localBlocks(locIndx), lev)
  dtdx = dtp/dx(lev)
  dtdy = dtp/dy(lev)
  dtdz = dtp/dz(lev)

  ! 1) Reference E = -v x B at cell centres, including one ghost layer
  do k=0,n+1
    do j=0,n+1
      do i=0,n+1
        ER(1,i,j,k) = PRIM(locIndx,4,i,j,k)*PRIM(locIndx,7,i,j,k)              &
                    - PRIM(locIndx,3,i,j,k)*PRIM(locIndx,8,i,j,k)
        ER(2,i,j,k) = PRIM(locIndx,2,i,j,k)*PRIM(locIndx,8,i,j,k)              &
                    - PRIM(locIndx,4,i,j,k)*PRIM(locIndx,6,i,j,k)
        ER(3,i,j,k) = PRIM(locIndx,3,i,j,k)*PRIM(locIndx,6,i,j,k)              &
                    - PRIM(locIndx,2,i,j,k)*PRIM(locIndx,7,i,j,k)
      end do
    end do
  end do

  ! 2) Edge EMFs. Face values of E from the induction fluxes:
  !    x-faces: Ez = -F(By), Ey = +F(Bz)
  !    y-faces: Ez = +G(Bx), Ex = -G(Bz)
  !    z-faces: Ey = -H(Bx), Ex = +H(By)

  ! Ez at (i+1/2,j+1/2,k): a = x, b = y
  do k=1,n
    do j=0,n
      do i=0,n
        EZ(i,j,k) = emf_edge(                                                  &
          -FC(7,i,j,k), -FC(7,i,j+1,k),  FC(1,i,j,k),  FC(1,i,j+1,k),          &
           GC(6,i,j,k),  GC(6,i+1,j,k),  GC(1,i,j,k),  GC(1,i+1,j,k),          &
           ER(3,i,j,k),  ER(3,i+1,j,k),  ER(3,i,j+1,k), ER(3,i+1,j+1,k) )
      end do
    end do
  end do

  ! Ex at (i,j+1/2,k+1/2): a = y, b = z
  do k=0,n
    do j=0,n
      do i=1,n
        EX(i,j,k) = emf_edge(                                                  &
          -GC(8,i,j,k), -GC(8,i,j,k+1),  GC(1,i,j,k),  GC(1,i,j,k+1),          &
           HC(7,i,j,k),  HC(7,i,j+1,k),  HC(1,i,j,k),  HC(1,i,j+1,k),          &
           ER(1,i,j,k),  ER(1,i,j+1,k),  ER(1,i,j,k+1), ER(1,i,j+1,k+1) )
      end do
    end do
  end do

  ! Ey at (i+1/2,j,k+1/2): a = z, b = x
  do k=0,n
    do j=1,n
      do i=0,n
        EY(i,j,k) = emf_edge(                                                  &
          -HC(6,i,j,k), -HC(6,i+1,j,k),  HC(1,i,j,k),  HC(1,i+1,j,k),          &
           FC(8,i,j,k),  FC(8,i,j,k+1),  FC(1,i,j,k),  FC(1,i,j,k+1),          &
           ER(2,i,j,k),  ER(2,i,j,k+1),  ER(2,i+1,j,k), ER(2,i+1,j,k+1) )
      end do
    end do
  end do

  ! 3) CT update of the faces, always from BF at time n (Athena eqs. 16-18)
  ! Bx at x-faces i+1/2: dBx/dt = -(dEz/dy - dEy/dz)
  do k=1,n
    do j=1,n
      do i=0,n
        BFP(locIndx,1,i,j,k) = BF(locIndx,1,i,j,k)                             &
                             - dtdy*( EZ(i,j,k) - EZ(i,j-1,k) )                &
                             + dtdz*( EY(i,j,k) - EY(i,j,k-1) )
      end do
    end do
  end do
  ! By at y-faces j+1/2: dBy/dt = -(dEx/dz - dEz/dx)
  do k=1,n
    do j=0,n
      do i=1,n
        BFP(locIndx,2,i,j,k) = BF(locIndx,2,i,j,k)                             &
                             + dtdx*( EZ(i,j,k) - EZ(i-1,j,k) )                &
                             - dtdz*( EX(i,j,k) - EX(i,j,k-1) )
      end do
    end do
  end do
  ! Bz at z-faces k+1/2: dBz/dt = -(dEy/dx - dEx/dy)
  do k=0,n
    do j=1,n
      do i=1,n
        BFP(locIndx,3,i,j,k) = BF(locIndx,3,i,j,k)                             &
                             - dtdx*( EY(i,j,k) - EY(i-1,j,k) )                &
                             + dtdy*( EX(i,j,k) - EX(i,j-1,k) )
      end do
    end do
  end do

  ! 4) Cell-centred field = average of the two faces (Athena eqs. 26-28)
  do k=1,n
    do j=1,n
      do i=1,n
        UP(locIndx,6,i,j,k) = 0.5*( BFP(locIndx,1,i-1,j,k) + BFP(locIndx,1,i,j,k) )
        UP(locIndx,7,i,j,k) = 0.5*( BFP(locIndx,2,i,j-1,k) + BFP(locIndx,2,i,j,k) )
        UP(locIndx,8,i,j,k) = 0.5*( BFP(locIndx,3,i,j,k-1) + BFP(locIndx,3,i,j,k) )
      end do
    end do
  end do

end subroutine flux_ct_update

!===============================================================================

!> @brief Accepts the stepped face field: BF = BFP on all local blocks
subroutine flux_ct_advance ()

  use parameters, only : nbMaxProc
  use globals,    only : localBlocks, BF, BFP
  implicit none

  integer :: nb

  do nb=1,nbMaxProc
    if (localBlocks(nb) /= -1) BF(nb,:,:,:,:) = BFP(nb,:,:,:,:)
  end do

end subroutine flux_ct_advance

!===============================================================================

!> @brief GS05 edge EMF in the (a,b) plane of the edge (a+1/2, b+1/2)
!> @details ea0, ea1: E on the a-faces (a+1/2,b) and (a+1/2,b+1), with mass
!! fluxes ma0, ma1; eb0, eb1: E on the b-faces (a,b+1/2) and (a+1,b+1/2),
!! with mass fluxes mb0, mb1; rXY: reference E at cell (a+X, b+Y).
!! Average of the four face values plus four upwinded (face - cell)
!! corrections (Athena eqs. 79-81).
pure real function emf_edge (ea0, ea1, ma0, ma1, eb0, eb1, mb0, mb1,         &
                             r00, r10, r01, r11)

  implicit none
  real, intent(in) :: ea0, ea1, ma0, ma1, eb0, eb1, mb0, mb1
  real, intent(in) :: r00, r10, r01, r11

  emf_edge = 0.25*( ea0 + ea1 + eb0 + eb1                                      &
                  + upwind(ma0, eb0 - r00, eb1 - r10)                          &
                  + upwind(ma1, eb0 - r01, eb1 - r11)                          &
                  + upwind(mb0, ea0 - r00, ea1 - r01)                          &
                  + upwind(mb1, ea0 - r10, ea1 - r11) )

end function emf_edge

!===============================================================================

!> @brief Upwind selection by the sign of the mass flux m (average if m=0)
pure real function upwind (m, lo, hi)

  implicit none
  real, intent(in) :: m, lo, hi

  if (m > 0.0) then
    upwind = lo
  else if (m < 0.0) then
    upwind = hi
  else
    upwind = 0.5*(lo + hi)
  end if

end function upwind

!===============================================================================

!> @brief Fills the ghost faces of a face-centred field (BF or BFP)
!> @details Same sweep as normalBoundary (x, then y, then z, with
!! edges=.true.), so edge and corner ghost faces are filled as well. With
!! the face index convention (index i is face i+1/2) a ghost layer is the
!! source layer shifted by n cells, as for cell data. The low ghost layer
!! includes index 0, the block's own low face, which is overwritten by the
!! neighbour's copy; the difference is recorded in bf_shared_err.
!! Only same-level neighbours (incl. periodic) are supported for now.
!> @param depth Number of ghost layers
!> @param bfield Face-centred field (BF or BFP)
subroutine boundary_bf (depth, bfield)

  use parameters
  use globals,    only : globalBlocks, localBlocks, rank, ierr, logu
  use amr,        only : neighbors, getOwner
  use boundaries, only : layerLimits, opposite
  use utils,      only : find
  use clean_quit, only : clean_abort
  implicit none

  integer, intent(in) :: depth
  real, intent(inout) :: bfield(nbMaxProc, 3, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

  integer :: direction, src_face, axis, nb, destID, destOwner, destInd
  integer :: srcID, srcOwner, srcInd, ntype, neighs(4), nData
  integer :: i1, i2, j1, j2, k1, k2, i3, i4, j3, j4, k3, k4
  integer :: mpistatus(MPI_STATUS_SIZE)
  real    :: saved(ncells_x, ncells_x), now(ncells_x, ncells_x)

  ! LEFT=1, RIGHT=2 (x), FRONT=3, BACK=4 (y), BOTTOM=5, TOP=6 (z)
  do direction=1,6

    call opposite (direction, src_face)
    axis = (direction+1)/2

    do nb=1,nbMaxGlobal
      destID = globalBlocks(nb)
      if (destID == -1) cycle

      call getOwner (destID, destOwner)
      call neighbors (destID, direction, ntype, neighs)

      select case (ntype)

      case (NEIGH_SAME)

        srcID = neighs(1)
        call getOwner (srcID, srcOwner)
        call layerLimits (src_face,  depth, .false., i3, i4, j3, j4, k3, k4, edges=.true.)
        call layerLimits (direction, depth, .true.,  i1, i2, j1, j2, k1, k2, edges=.true.)
        nData = 3*(i2-i1+1)*(j2-j1+1)*(k2-k1+1)

        ! Source side
        if ((rank == srcOwner).and.(rank /= destOwner)) then
          call find (srcID, localBlocks, nbMaxProc, srcInd)
          call MPI_SEND (bfield(srcInd,1:3,i3:i4,j3:j4,k3:k4), nData,         &
                         mpi_real_kind, destOwner, srcID, mpi_comm_world, ierr)
        end if

        ! Destination side
        if (rank == destOwner) then
          call find (destID, localBlocks, nbMaxProc, destInd)
          if (mod(direction,2) == 1) call own_low_face (bfield, destInd, axis, saved)
          if (rank == srcOwner) then
            call find (srcID, localBlocks, nbMaxProc, srcInd)
            bfield(destInd,1:3,i1:i2,j1:j2,k1:k2) = bfield(srcInd,1:3,i3:i4,j3:j4,k3:k4)
          else
            call MPI_RECV (bfield(destInd,1:3,i1:i2,j1:j2,k1:k2), nData,      &
                           mpi_real_kind, srcOwner, srcID, mpi_comm_world,     &
                           mpistatus, ierr)
          end if
          if (mod(direction,2) == 1) then
            call own_low_face (bfield, destInd, axis, now)
            bf_shared_err = max(bf_shared_err, maxval(abs(now - saved)))
          end if
        end if

      case default

        ! Level jumps and non-periodic domain boundaries: not implemented yet
        if (rank == destOwner) then
          write(logu,'(a,i0,a,i0,a,i0)') "boundary_bf: neighbour type ", ntype, &
                 " (dir ", direction, ") not supported yet, block ", destID
          call clean_abort (ERROR_GENERIC)
        end if

      end select

    end do
  end do

end subroutine boundary_bf

!===============================================================================

!> @brief Copies the block's own low face (index 0) of component axis
subroutine own_low_face (bfield, locIndx, axis, face)

  use parameters
  implicit none

  real, intent(in)    :: bfield(nbMaxProc, 3, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  integer, intent(in) :: locIndx, axis
  real, intent(out)   :: face(ncells_x, ncells_x)

  integer :: n

  n = ncells_x
  select case (axis)
  case (1)
    face(:,:) = bfield(locIndx,1,0,1:n,1:n)
  case (2)
    face(:,:) = bfield(locIndx,2,1:n,0,1:n)
  case (3)
    face(:,:) = bfield(locIndx,3,1:n,1:n,0)
  end select

end subroutine own_low_face

!===============================================================================

#endif

end module flux_ct
