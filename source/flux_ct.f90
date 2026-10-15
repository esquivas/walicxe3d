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

  !> Edge EMFs of every local block, "right" indexing like the fluxes:
  !! EMF(1,i,j,k,nb) = Ex at (i,j+1/2,k+1/2)   (i=1..n,  j,k=0..n)
  !! EMF(2,i,j,k,nb) = Ey at (i+1/2,j,k+1/2)   (j=1..n,  i,k=0..n)
  !! EMF(3,i,j,k,nb) = Ez at (i+1/2,j+1/2,k)   (k=1..n,  i,j=0..n)
  real, allocatable, save :: EMF(:,:,:,:,:)
  !> Cell-centred reference field E = -v x B of the current block
  real, allocatable, save :: ER(:,:,:,:)
  !> Largest |coarse face - average of its 4 fine faces| (check_cf_faces)
  real, save :: cf_face_err = 0.0
  !> Largest difference found between the two copies of a shared block face
  !! (set by boundary_bf; should be exactly 0)
  real, save :: bf_shared_err = 0.0
  !> Weights of the edge EMFs during sync_emf: number of block copies averaged
  !! into the value, or WFINE for values restricted from finer blocks
  real, allocatable, save :: EW(:,:,:,:,:)
  real,    parameter :: WFINE = 1000.0
  integer, parameter :: SYNC_AVERAGE = 1, SYNC_MAXW = 2

contains

!===============================================================================

!> @brief Edge EMFs of one block (GS05), stored in EMF(:,:,:,:,locIndx)
!> @details Call right after the block's fluxes are computed (pass A). The
!! EMFs are synchronised across blocks (sync_emf) before the CT update.
!> @param locIndx Local index of the block
subroutine flux_ct_emf (locIndx)

  use parameters
  use globals, only : PRIM, FC, GC, HC
  implicit none

  integer, intent(in) :: locIndx

  integer :: i, j, k, n

  n = ncells_x
  if (.not.allocated(EMF)) then
    allocate( EMF(3,0:n,0:n,0:n,nbMaxProc) )
    EMF(:,:,:,:,:) = 0.0
    allocate( EW(3,0:n,0:n,0:n,nbMaxProc) )
    EW(:,:,:,:,:) = 1.0
    allocate( ER(3,0:n+1,0:n+1,0:n+1) )
  end if

  ! 1) Reference E = -v x B
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

  ! 2) Edge EMFs
  ! Ez at (i+1/2,j+1/2,k): a = x, b = y
  do k=1,n
    do j=0,n
      do i=0,n
        EMF(3,i,j,k, locIndx) = emf_edge(                                      &
          -FC(7,i,j,k), -FC(7,i,j+1,k),  FC(1,i,j,k),   FC(1,i,j+1,k),         &
           GC(6,i,j,k),  GC(6,i+1,j,k),  GC(1,i,j,k),   GC(1,i+1,j,k),         &
           ER(3,i,j,k),  ER(3,i+1,j,k),  ER(3,i,j+1,k), ER(3,i+1,j+1,k) )
      end do
    end do
  end do

  ! Ex at (i,j+1/2,k+1/2): a = y, b = z
  do k=0,n
    do j=0,n
      do i=1,n
        EMF(1,i,j,k,locIndx) = emf_edge(                                       &
          -GC(8,i,j,k), -GC(8,i,j,k+1),  GC(1,i,j,k),   GC(1,i,j,k+1),         &
           HC(7,i,j,k),  HC(7,i,j+1,k),  HC(1,i,j,k),   HC(1,i,j+1,k),         &
           ER(1,i,j,k),  ER(1,i,j+1,k),  ER(1,i,j,k+1), ER(1,i,j+1,k+1) )
      end do
    end do
  end do

  ! Ey at (i+1/2,j,k+1/2): a = z, b = x
  do k=0,n
    do j=1,n
      do i=0,n
        EMF(2,i,j,k,locIndx) = emf_edge(                                       &
          -HC(6,i,j,k), -HC(6,i+1,j,k),  HC(1,i,j,k),   HC(1,i+1,j,k),         &
           FC(8,i,j,k),  FC(8,i,j,k+1),  FC(1,i,j,k),   FC(1,i,j,k+1),         &
           ER(2,i,j,k),  ER(2,i,j,k+1),  ER(2,i+1,j,k), ER(2,i+1,j,k+1) )
      end do
    end do
  end do

end subroutine flux_ct_emf

!===============================================================================

!> @brief EMFs, CT update of BFP and cell-centred B in UP for one block
!> @param locIndx Local index of the block
!> @param dtp Timestep of this stage (as passed to upwindStep)
subroutine flux_ct_update (locIndx, dtp)

  use parameters
  use globals, only : localBlocks, BF, BFP, UP, dx, dy, dz
  use amr,     only : meshlevel
  implicit none

  integer, intent(in) :: locIndx
  real,    intent(in) :: dtp

  integer :: i, j, k, n, lev
  real    :: dtdx, dtdy, dtdz

  n = ncells_x

  call meshlevel (localBlocks(locIndx), lev)
  dtdx = dtp/dx(lev)
  dtdy = dtp/dy(lev)
  dtdz = dtp/dz(lev)

  ! 3) CT update of the faces, always from BF at time n (Athena eqs. 16-18)
  ! Bx at x-faces i+1/2: dBx/dt = -(dEz/dy - dEy/dz)
  do k=1,n
    do j=1,n
      do i=0,n
        BFP(locIndx,1,i,j,k) = BF(locIndx,1,i,j,k)                             &
                    - dtdy*( EMF(3,i,j,k,locIndx) - EMF(3,i,j-1,k,locIndx) )   &
                    + dtdz*( EMF(2,i,j,k,locIndx) - EMF(2,i,j,k-1,locIndx) )
      end do
    end do
  end do
  ! By at y-faces j+1/2: dBy/dt = -(dEx/dz - dEz/dx)
  do k=1,n
    do j=0,n
      do i=1,n
        BFP(locIndx,2,i,j,k) = BF(locIndx,2,i,j,k)                             &
                    + dtdx*( EMF(3,i,j,k,locIndx) - EMF(3,i-1,j,k,locIndx) )   &
                    - dtdz*( EMF(1,i,j,k,locIndx) - EMF(1,i,j,k-1,locIndx) )
      end do
    end do
  end do
  ! Bz at z-faces k+1/2: dBz/dt = -(dEy/dx - dEx/dy)
  do k=0,n
    do j=1,n
      do i=1,n
        BFP(locIndx,3,i,j,k) = BF(locIndx,3,i,j,k)                             &
                    - dtdx*( EMF(2,i,j,k,locIndx) - EMF(2,i-1,j,k,locIndx) )   &
                    + dtdy*( EMF(1,i,j,k,locIndx) - EMF(1,i,j-1,k,locIndx) )
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

!> @brief Combines the copies of shared edge EMFs between same-level blocks
!> @details Every same-level pair of face neighbours combines the EMFs on
!! their common face (the components tangential to it), with weights EW.
!! Pairs are taken as (high face of a, low face of b) in an x, then y, then
!! z sweep. Each edge is exchanged in the two sweeps transverse to it, so
!! the up to four blocks around it are combined in two steps. Used by
!! sync_emf on one level at a time:
!!  - SYNC_AVERAGE: weighted average, weights added; values restricted from
!!    the next finer level (weight WFINE) win. Pairs always merge disjoint
!!    sets of copies, so with 4 blocks around an edge all copies end with
!!    the full average (weight 4). With 3 (a coarser block at an inner
!!    corner of a finer region) one copy is left with weight 2.
!!  - SYNC_MAXW: the copy with the larger weight wins, which completes the
!!    3-block edges and spreads any remaining WFINE values.
!! With equal weights the result is bit-identical to a plain average, so
!! copies that are already identical are left unchanged. Must be called by
!! ALL ranks, after flux_ct_emf on every local block.
!> @param mode SYNC_AVERAGE or SYNC_MAXW
!> @param lev refinement level to synchronize
subroutine sync_same_level (mode, lev)

  use parameters
  use globals, only : globalBlocks, localBlocks, rank, ierr
  use amr,     only : neighbors, getOwner, meshlevel
  use utils,   only : find

  implicit none
  integer, intent(in) :: mode, lev

  integer :: axis, nb, aID, aOwner, aInd, bID, bOwner, bInd, ntype, neighs(4)
  integer :: nData, alev
  integer :: mpistatus(MPI_STATUS_SIZE)
  real    :: fa(4*ncells_x*(ncells_x+1)), fb(4*ncells_x*(ncells_x+1))

  nData = 4*ncells_x*(ncells_x+1)
  aInd  = -1
  bInd  = -1

  do axis=1,3              ! RIGHT/BACK/TOP face of a = LEFT/FRONT/BOTTOM of b

    do nb=1,nbMaxGlobal
      aID = globalBlocks(nb)
      if (aID == -1) cycle
      call meshlevel (aID, alev)
      if (alev /= lev) cycle

      call neighbors (aID, 2*axis, ntype, neighs)
      if (ntype /= NEIGH_SAME) cycle
      bID = neighs(1)
      call getOwner (aID, aOwner)
      call getOwner (bID, bOwner)

      if (rank == aOwner) then
        call find (aID, localBlocks, nbMaxProc, aInd)
        call pack_face_emf (aInd, axis, ncells_x, fa)       ! a's high face
      end if
      if (rank == bOwner) then
        call find (bID, localBlocks, nbMaxProc, bInd)
        call pack_face_emf (bInd, axis, 0, fb)              ! b's low face
      end if

      if ((rank == aOwner).and.(rank == bOwner)) then
        call combine_face_emf( fa, fb, mode )
        call unpack_face_emf (aInd, axis, ncells_x, fa)
        call unpack_face_emf (bInd, axis, 0, fa)
      else if (rank == aOwner) then
        call MPI_SENDRECV (fa, nData, mpi_real_kind, bOwner, aID,             &
                           fb, nData, mpi_real_kind, bOwner, bID,             &
                           mpi_comm_world, mpistatus, ierr)
        call combine_face_emf( fa, fb, mode )
        call unpack_face_emf (aInd, axis, ncells_x, fa)
      else if (rank == bOwner) then
        call MPI_SENDRECV (fb, nData, mpi_real_kind, aOwner, bID,             &
                           fa, nData, mpi_real_kind, aOwner, aID,             &
                           mpi_comm_world, mpistatus, ierr)
        call combine_face_emf( fa, fb, mode )
        call unpack_face_emf (bInd, axis, 0, fa)
      end if

    end do
  end do

end subroutine sync_same_level

!===============================================================================

!> @brief Makes all copies of every shared edge EMF identical
!> @details Each copy carries a weight EW (number of block copies averaged
!! into it). Level by level, from the finest: 1) same-level weighted
!! average, then 2) larger weight wins, so the blocks of that level agree
!! around every edge, including edges shared by only three of them (inner
!! corners of a finer region); values restricted from the next finer level
!! (weight WFINE) win in both passes; 3) the coarse edges on coarse-fine
!! faces of the next coarser level are replaced by the restricted values
!! (weight WFINE). Each level is final before it is restricted. Must be
!! called by ALL ranks after flux_ct_emf on every local block.
subroutine sync_emf ()

  use parameters, only : nbMaxProc
  use globals,    only : localBlocks, maxlev
  implicit none

  integer :: nb, lev

  do nb=1,nbMaxProc
    if (localBlocks(nb) /= -1) EW(:,:,:,:,nb) = 1.0
  end do

  do lev=maxlev,1,-1
    call sync_same_level (SYNC_AVERAGE, lev)  ! weighted average (WFINE wins)
    call sync_same_level (SYNC_MAXW, lev)     ! completes edges shared by 3 blocks
    if (lev > 1) call restrict_fine_emf (lev-1)
  end do

end subroutine sync_emf

!===============================================================================

!> @brief Replaces the coarse edge EMFs on every coarse-fine face by the
!! average of the two collinear fine edge EMFs, with weight WFINE
!> @details For each coarse face whose neighbours are four finer blocks,
!! each fine block restricts its face EMFs to the coarse edges of its
!! quadrant and sends them to the coarse owner. This keeps every coarse face
!! equal to the average of its four fine faces. The fine copies must
!! already agree (same-level passes of sync_emf done before).
!> @param clev Level of the coarse blocks (their finer neighbours are at clev+1)
subroutine restrict_fine_emf (clev)

  use parameters
  use globals,    only : globalBlocks, localBlocks, rank, ierr, logu
  use amr,        only : neighbors, getOwner, siblingCoords, meshlevel
  use utils,      only : find
  use clean_quit, only : clean_abort
  implicit none

  integer, intent(in) :: clev

  integer :: dir, axis, nb, ib, cID, cOwner, cIdx, fID, fOwner, fIdx, cl
  integer :: ntype, neighs(4), n, h, lc, lf, sx, sy, sz, oa, ob, a, b, nData
  integer :: mpistatus(MPI_STATUS_SIZE)
  real    :: t1(ncells_x,0:ncells_x), t2(0:ncells_x,ncells_x)
  real    :: w1(ncells_x,0:ncells_x), w2(0:ncells_x,ncells_x)
  real    :: r1(ncells_x/2,0:ncells_x/2), r2(0:ncells_x/2,ncells_x/2)
  real    :: buf(2*(ncells_x/2)*(ncells_x/2+1))

  n     = ncells_x
  h     = n/2
  nData = 2*h*(h+1)
  oa    = 0
  ob    = 0

  do dir=1,6
    axis = (dir+1)/2
    lc   = merge(0, n, mod(dir,2) == 1)   ! coarse face (low face for LEFT/FRONT/BOTTOM)
    lf   = n - lc                         ! the fine blocks' face, on the other side

    do nb=1,nbMaxGlobal
      cID = globalBlocks(nb)
      if (cID == -1) cycle
      call meshlevel (cID, cl)
      if (cl /= clev) cycle
      call neighbors (cID, dir, ntype, neighs)
      if (ntype /= NEIGH_FINER) cycle
      call getOwner (cID, cOwner)

      do ib=1,4
        fID = neighs(ib)
        call getOwner (fID, fOwner)
        if (fOwner == MPI_PROC_NULL) then
          write(logu,'(a,i0,a,i0)') "restrict_fine_emf: fine block ", fID,    &
                                    " not found for block ", cID
          call clean_abort (ERROR_GENERIC)
        end if

        ! Fine side: restrict to the coarse edges of this quadrant
        if (rank == fOwner) then
          call find (fID, localBlocks, nbMaxProc, fIdx)
          call get_face (fIdx, axis, lf, t1, t2, w1, w2)
          do b=0,h
            do a=1,h
              r1(a,b) = 0.5*( t1(2*a-1,2*b) + t1(2*a,2*b) )
            end do
          end do
          do b=1,h
            do a=0,h
              r2(a,b) = 0.5*( t2(2*a,2*b-1) + t2(2*a,2*b) )
            end do
          end do
          if (rank /= cOwner) then
            buf(1:h*(h+1))         = reshape( r1, [h*(h+1)] )
            buf(h*(h+1)+1:2*h*(h+1)) = reshape( r2, [h*(h+1)] )
            call MPI_SEND (buf, nData, mpi_real_kind, cOwner, fID,            &
                           mpi_comm_world, ierr)
          end if
        end if

        ! Coarse side: replace and flag the edges of the quadrant
        if (rank == cOwner) then
          if (rank /= fOwner) then
            call MPI_RECV (buf, nData, mpi_real_kind, fOwner, fID,            &
                           mpi_comm_world, mpistatus, ierr)
            r1 = reshape( buf(1:h*(h+1)),            [h, h+1] )
            r2 = reshape( buf(h*(h+1)+1:2*h*(h+1)),  [h+1, h] )
          end if
          call find (cID, localBlocks, nbMaxProc, cIdx)
          call siblingCoords (fID, sx, sy, sz)
          select case (axis)
          case (1)
            oa = sy*h ; ob = sz*h
          case (2)
            oa = sx*h ; ob = sz*h
          case (3)
            oa = sx*h ; ob = sy*h
          end select
          call get_face (cIdx, axis, lc, t1, t2, w1, w2)
          t1(oa+1:oa+h, ob:ob+h) = r1 ; w1(oa+1:oa+h, ob:ob+h) = WFINE
          t2(oa:oa+h, ob+1:ob+h) = r2 ; w2(oa:oa+h, ob+1:ob+h) = WFINE
          call put_face (cIdx, axis, lc, t1, t2, w1, w2)
        end if

      end do
    end do
  end do

end subroutine restrict_fine_emf

!===============================================================================

!> @brief Check: every coarse face on a coarse-fine interface must equal the
!! average of its four fine faces. Result (max abs difference over all such
!! faces) in cf_face_err. Must be called by ALL ranks.
subroutine check_cf_faces ()

  use parameters
  use globals,    only : globalBlocks, localBlocks, rank, ierr, BF
  use amr,        only : neighbors, getOwner, siblingCoords
  use utils,      only : find
  implicit none

  integer :: dir, axis, nb, ib, cID, cOwner, cIdx, fID, fOwner, fIdx
  integer :: ntype, neighs(4), n, h, lc, lf, sx, sy, sz, oa, ob, a, b
  integer :: mpistatus(MPI_STATUS_SIZE)
  real    :: fc(ncells_x,ncells_x), fr(ncells_x/2,ncells_x/2)

  n  = ncells_x
  h  = n/2
  oa = 0
  ob = 0
  cf_face_err = 0.0

  do dir=1,6
    axis = (dir+1)/2
    lc   = merge(0, n, mod(dir,2) == 1)
    lf   = n - lc

    do nb=1,nbMaxGlobal
      cID = globalBlocks(nb)
      if (cID == -1) cycle
      call neighbors (cID, dir, ntype, neighs)
      if (ntype /= NEIGH_FINER) cycle
      call getOwner (cID, cOwner)

      do ib=1,4
        fID = neighs(ib)
        call getOwner (fID, fOwner)

        if (rank == fOwner) then
          call find (fID, localBlocks, nbMaxProc, fIdx)
          call normal_face (fIdx, axis, lf, fc)
          do b=1,h
            do a=1,h
              fr(a,b) = 0.25*( fc(2*a-1,2*b-1) + fc(2*a,2*b-1)                &
                             + fc(2*a-1,2*b  ) + fc(2*a,2*b  ) )
            end do
          end do
          if (rank /= cOwner) call MPI_SEND (fr, h*h, mpi_real_kind, cOwner,  &
                                             fID, mpi_comm_world, ierr)
        end if

        if (rank == cOwner) then
          if (rank /= fOwner) call MPI_RECV (fr, h*h, mpi_real_kind, fOwner,  &
                                       fID, mpi_comm_world, mpistatus, ierr)
          call find (cID, localBlocks, nbMaxProc, cIdx)
          call siblingCoords (fID, sx, sy, sz)
          select case (axis)
          case (1)
            oa = sy*h ; ob = sz*h
          case (2)
            oa = sx*h ; ob = sz*h
          case (3)
            oa = sx*h ; ob = sy*h
          end select
          call normal_face (cIdx, axis, lc, fc)
          cf_face_err = max( cf_face_err,                                     &
                             maxval(abs( fc(oa+1:oa+h, ob+1:ob+h) - fr )) )
        end if

      end do
    end do
  end do

end subroutine check_cf_faces

!===============================================================================

!> @brief Normal face field of one face of a block, f(a,b) (transverse axes)
subroutine normal_face (ind, axis, l, f)

  use parameters, only : ncells_x
  use globals,    only : BF
  implicit none

  integer, intent(in)  :: ind, axis, l
  real,    intent(out) :: f(ncells_x,ncells_x)

  integer :: n

  n = ncells_x
  select case (axis)
  case (1)
    f = BF(ind,1,l,1:n,1:n)
  case (2)
    f = BF(ind,2,1:n,l,1:n)
  case (3)
    f = BF(ind,3,1:n,1:n,l)
  end select

end subroutine normal_face

!===============================================================================

!> @brief Packs the tangential EMFs of one face and their weights into f(4m),
!! m = n(n+1): t1, t2, weights(t1), weights(t2)
subroutine pack_face_emf (ind, axis, l, f)

  use parameters, only : ncells_x
  implicit none

  integer, intent(in)  :: ind, axis, l
  real,    intent(out) :: f(4*ncells_x*(ncells_x+1))

  real    :: t1(ncells_x,0:ncells_x), t2(0:ncells_x,ncells_x)
  real    :: w1(ncells_x,0:ncells_x), w2(0:ncells_x,ncells_x)
  integer :: m

  m = ncells_x*(ncells_x+1)
  call get_face (ind, axis, l, t1, t2, w1, w2)
  f(1:m)       = reshape( t1, [m] )
  f(m+1:2*m)   = reshape( t2, [m] )
  f(2*m+1:3*m) = reshape( w1, [m] )
  f(3*m+1:4*m) = reshape( w2, [m] )

end subroutine pack_face_emf

!===============================================================================

!> @brief Inverse of pack_face_emf
subroutine unpack_face_emf (ind, axis, l, f)

  use parameters, only : ncells_x
  implicit none

  integer, intent(in) :: ind, axis, l
  real,    intent(in) :: f(4*ncells_x*(ncells_x+1))

  real    :: t1(ncells_x,0:ncells_x), t2(0:ncells_x,ncells_x)
  real    :: w1(ncells_x,0:ncells_x), w2(0:ncells_x,ncells_x)
  integer :: n, m

  n = ncells_x
  m = n*(n+1)
  t1 = reshape( f(1:m),       [n, n+1] )
  t2 = reshape( f(m+1:2*m),   [n+1, n] )
  w1 = reshape( f(2*m+1:3*m), [n, n+1] )
  w2 = reshape( f(3*m+1:4*m), [n+1, n] )
  call put_face (ind, axis, l, t1, t2, w1, w2)

end subroutine unpack_face_emf


!===============================================================================
!> @brief Combines two packed copies of a face (a first, b second)
!> @details mode SYNC_AVERAGE: weighted average of the two copies, weights
!! added (values restricted from fine blocks, weight WFINE, still win).
!! mode SYNC_MAXW: the copy with the larger weight wins (equal weights:
!! average, a no-op for identical copies). Symmetric in the sense that both
!! owners call it with the same (a, b) order, so they get identical results.
!! Result in fa.
subroutine combine_face_emf (fa, fb, mode)

  use parameters, only : ncells_x
  implicit none

  real,    intent(inout) :: fa(4*ncells_x*(ncells_x+1))
  real,    intent(in)    :: fb(4*ncells_x*(ncells_x+1))
  integer, intent(in)    :: mode

  integer :: l, nv
  real    :: wa, wb

  nv = 2*ncells_x*(ncells_x+1)
  do l=1,nv
    wa = fa(nv+l)
    wb = fb(nv+l)
    if ((mode == SYNC_AVERAGE).and.(max(wa,wb) < WFINE)) then
      fa(l)    = (wa*fa(l) + wb*fb(l))/(wa + wb)
      fa(nv+l) = wa + wb
    else if (wb > wa) then
      fa(l)    = fb(l)
      fa(nv+l) = wb
    else if (wa == wb) then
      fa(l) = 0.5*(fa(l) + fb(l))
    end if
  end do

end subroutine combine_face_emf

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
!! Level jumps: tangential ghost faces from the cell-centred field
!! (faces_from_cells). Non-periodic domain boundaries: bc_faces.
!> @param depth  : Number of ghost layers
!> @param bfield : Face-centred field (BF or BFP)
!> @param uvars  : Conserved flow variables array
subroutine boundary_bf (depth, bfield, uvars)

  use parameters
  use globals,    only : globalBlocks, localBlocks, rank, ierr, logu
  use amr,        only : neighbors, getOwner
  use boundaries, only : layerLimits, opposite
  use utils,      only : find
  use clean_quit, only : clean_abort
  implicit none

  integer, intent(in) :: depth
  real, intent(inout) :: bfield(nbMaxProc, 3, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  real, intent(in)    :: uvars (nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)

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
        ! Also copy face index 0 of the transverse directions that layerLimits
        ! does not extend over ghost cells (x pass: y and z; y pass: z). In the
        ! ghost columns these are the neighbouring blocks' own low faces (normal
        ! B), needed by the fluxes on the extra transverse layer; they are not
        ! filled later when the transverse neighbour is at a level jump. The
        ! tangential components copied there are overwritten by later passes.
        if (axis == 1) then
          j1 = j1 - 1 ; j3 = j3 - 1
          k1 = k1 - 1 ; k3 = k3 - 1
        else if (axis == 2) then
          k1 = k1 - 1 ; k3 = k3 - 1
        end if
        nData = 3*(i2-i1+1)*(j2-j1+1)*(k2-k1+1)

        ! Source side
        if ((rank == srcOwner).and.(rank /= destOwner)) then
          call find (srcID, localBlocks, nbMaxProc, srcInd)
          call MPI_SEND (bfield(srcInd,1:3,i3:i4,j3:j4,k3:k4), nData,          &
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
            call MPI_RECV (bfield(destInd,1:3,i1:i2,j1:j2,k1:k2), nData,       &
                           mpi_real_kind, srcOwner, srcID, mpi_comm_world,     &
                           mpistatus, ierr)
          end if
          if (mod(direction,2) == 1) then
            call own_low_face (bfield, destInd, axis, now)
            bf_shared_err = max(bf_shared_err, maxval(abs(now - saved)))
          end if
        end if

      case (NEIGH_COARSER, NEIGH_FINER)

        ! Level jump: tangential ghost faces from the cell-centred field
        if (rank == destOwner) then
          call find (destID, localBlocks, nbMaxProc, destInd)
          call faces_from_cells (bfield, uvars, destInd, direction)
        end if

      case (NEIGH_BOUNDARY)

        ! Non-periodic domain boundary (periodic neighbours are NEIGH_SAME)
        if (rank == destOwner) then
          call find (destID, localBlocks, nbMaxProc, destInd)
          call bc_faces (bfield, destInd, direction, depth)
        end if

      case default
      ! Unknown neighbour type
        if (rank == destOwner) then
          write(logu,'(a,i0,a,i0,a,i0)') "boundary_bf: neighbour type ", ntype,&
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

!> @brief Tangential ghost faces (one layer) from the cell-centred field
!> @details Used next to level jumps: each tangential face in the ghost
!! layer gets the average of its two adjacent cells of uvars(6:8) (ghost
!! cells, including edges, already filled by boundary). The normal component
!! is not touched (index 0 is the block's own face).
subroutine faces_from_cells (bfield, uvars, ind, direction)

  use parameters
  implicit none

  real, intent(inout) :: bfield(nbMaxProc,  3,     nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  real, intent(in)    :: uvars (nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  integer, intent(in) :: ind, direction

  integer :: i, j, k, n

  n = ncells_x
  select case (direction)
  case (LEFT, RIGHT)
    i = merge(0, n+1, direction == LEFT)
    do k=0,n+1
      do j=0,n
        bfield(ind,2,i,j,k) = 0.5*( uvars(ind,7,i,j,k) + uvars(ind,7,i,j+1,k) )
      end do
    end do
    do k=0,n
      do j=0,n+1
        bfield(ind,3,i,j,k) = 0.5*( uvars(ind,8,i,j,k) + uvars(ind,8,i,j,k+1) )
      end do
    end do
  case (FRONT, BACK)
    j = merge(0, n+1, direction == FRONT)
    do k=0,n+1
      do i=0,n
        bfield(ind,1,i,j,k) = 0.5*( uvars(ind,6,i,j,k) + uvars(ind,6,i+1,j,k) )
      end do
    end do
    do k=0,n
      do i=0,n+1
        bfield(ind,3,i,j,k) = 0.5*( uvars(ind,8,i,j,k) + uvars(ind,8,i,j,k+1) )
      end do
    end do
  case (BOTTOM, TOP)
    k = merge(0, n+1, direction == BOTTOM)
    do j=0,n+1
      do i=0,n
        bfield(ind,1,i,j,k) = 0.5*( uvars(ind,6,i,j,k) + uvars(ind,6,i+1,j,k) )
      end do
    end do
    do j=0,n
      do i=0,n+1
        bfield(ind,2,i,j,k) = 0.5*( uvars(ind,7,i,j,k) + uvars(ind,7,i,j+1,k) )
      end do
    end do
  end select

end subroutine faces_from_cells

!===============================================================================

!> @brief Ghost faces at a non-periodic domain boundary
!> @details Same index set as the same-level copy of boundary_bf (layerLimits
!! with edges). Tangential components: copied from the boundary cell
!! (BC_FREE_FLOW) or its mirror cell (BC_REFLECTIVE). Normal component: zero
!! gradient (free flow) or antisymmetric about the wall (reflective, i.e. a
!! conducting wall, as in Athena). The wall face itself (index 0 on the low
!! side, ncells on the high side) is the block's own face, evolved by CT,
!! and is not modified.
!> @param bfield Face-centred field (BF or BFP)
!> @param ind Local index of the block
!> @param direction Domain boundary the block touches (LEFT ... TOP)
!> @param depth Number of ghost layers
subroutine bc_faces (bfield, ind, direction, depth)

  use parameters
  use globals,    only : logu
  use boundaries, only : layerLimits
  use clean_quit, only : clean_abort
  implicit none

  real, intent(inout) :: bfield(nbMaxProc, 3, nxmin:nxmax, nymin:nymax, nzmin:nzmax)
  integer, intent(in) :: ind, direction, depth

  integer :: bc_type, axis, c, i, j, k, ip, jp, kp, l, lp, lw, n
  integer :: i1, i2, j1, j2, k1, k2
  logical :: low
  real    :: sgn

  select case (direction)
  case (LEFT)
    bc_type = bc_left
  case (RIGHT)
    bc_type = bc_right
  case (FRONT)
    bc_type = bc_front
  case (BACK)
    bc_type = bc_back
  case (BOTTOM)
    bc_type = bc_bottom
  case default
    bc_type = bc_top
  end select
  if ((bc_type /= BC_FREE_FLOW).and.(bc_type /= BC_REFLECTIVE)) then
    write(logu,'(a,i0,a,i0)') "bc_faces: BC type ", bc_type,                  &
                              " not supported, direction ", direction
    call clean_abort (ERROR_GENERIC)
  end if

  n    = ncells_x
  axis = (direction+1)/2
  low  = (mod(direction,2) == 1)
  lw   = merge(0, n, low)                       ! wall face index
  sgn  = merge(-1.0, 1.0, bc_type == BC_REFLECTIVE)

  call layerLimits (direction, depth, .true., i1, i2, j1, j2, k1, k2, edges=.true.)

  do k=k1,k2
    do j=j1,j2
      do i=i1,i2

        ! Ghost index along the normal
        l = merge(i, merge(j, k, axis == 2), axis == 1)

        ! Tangential components: boundary cell (free flow) or mirror cell
        if (bc_type == BC_FREE_FLOW) then
          lp = merge(1, n, low)
        else
          lp = merge(1-l, 2*n+1-l, low)
        end if
        ip = merge(lp, i, axis == 1)
        jp = merge(lp, j, axis == 2)
        kp = merge(lp, k, axis == 3)
        do c=1,3
          if (c /= axis) bfield(ind,c,i,j,k) = bfield(ind,c,ip,jp,kp)
        end do

        ! Normal component (face index l = face l+1/2): zero gradient or
        ! antisymmetric about the wall face; the wall face is not touched
        if (l /= lw) then
          lp = merge(lw, 2*lw-l, bc_type == BC_FREE_FLOW)
          ip = merge(lp, i, axis == 1)
          jp = merge(lp, j, axis == 2)
          kp = merge(lp, k, axis == 3)
          bfield(ind,axis,i,j,k) = sgn*bfield(ind,axis,ip,jp,kp)
        end if

      end do
    end do
  end do

end subroutine bc_faces

!===============================================================================

!> @brief Gets the tangential edge EMFs (and weights) on one face of a block
!> @details t1(a,b): cells along a, edges along b; t2(a,b): edges along a,
!! cells along b; (a,b) are the transverse axes of the face in ascending
!! order: x-face (y,z) -> t1=Ey, t2=Ez; y-face (x,z) -> t1=Ex, t2=Ez;
!! z-face (x,y) -> t1=Ex, t2=Ey. w1, w2: their weights EW (sync_emf).
subroutine get_face (ind, axis, l, t1, t2, w1, w2)

  use parameters, only : ncells_x
  implicit none

  integer, intent(in)  :: ind, axis, l
  real,    intent(out) :: t1(ncells_x,0:ncells_x), t2(0:ncells_x,ncells_x)
  real,    intent(out) :: w1(ncells_x,0:ncells_x), w2(0:ncells_x,ncells_x)

  integer :: n

  n = ncells_x
  select case (axis)
  case (1)
    t1 = EMF(2,l,1:n,0:n,ind) ; w1 = EW(2,l,1:n,0:n,ind)
    t2 = EMF(3,l,0:n,1:n,ind) ; w2 = EW(3,l,0:n,1:n,ind)
  case (2)
    t1 = EMF(1,1:n,l,0:n,ind) ; w1 = EW(1,1:n,l,0:n,ind)
    t2 = EMF(3,0:n,l,1:n,ind) ; w2 = EW(3,0:n,l,1:n,ind)
  case (3)
    t1 = EMF(1,1:n,0:n,l,ind) ; w1 = EW(1,1:n,0:n,l,ind)
    t2 = EMF(2,0:n,1:n,l,ind) ; w2 = EW(2,0:n,1:n,l,ind)
  end select

end subroutine get_face

!===============================================================================

!> @brief Inverse of get_face
subroutine put_face (ind, axis, l, t1, t2, w1, w2)

  use parameters, only : ncells_x
  implicit none

  integer, intent(in) :: ind, axis, l
  real,    intent(in) :: t1(ncells_x,0:ncells_x), t2(0:ncells_x,ncells_x)
  real,    intent(in) :: w1(ncells_x,0:ncells_x), w2(0:ncells_x,ncells_x)

  integer :: n

  n = ncells_x
  select case (axis)
  case (1)
    EMF(2,l,1:n,0:n,ind) = t1 ; EW(2,l,1:n,0:n,ind) = w1
    EMF(3,l,0:n,1:n,ind) = t2 ; EW(3,l,0:n,1:n,ind) = w2
  case (2)
    EMF(1,1:n,l,0:n,ind) = t1 ; EW(1,1:n,l,0:n,ind) = w1
    EMF(3,0:n,l,1:n,ind) = t2 ; EW(3,0:n,l,1:n,ind) = w2
  case (3)
    EMF(1,1:n,0:n,l,ind) = t1 ; EW(1,1:n,0:n,l,ind) = w1
    EMF(2,0:n,1:n,l,ind) = t2 ; EW(2,0:n,1:n,l,ind) = w2
  end select

end subroutine put_face

!===============================================================================

#endif

end module flux_ct
