!===============================================================================
!> @file diagnostics.f90
!> @brief Global conservation and div(B) diagnostics
!> @author A. Esquivel
!> @date 06/Oct/2026

!> @brief Global diagnostics module
!> @details Computes volume integrals of the conserved variables and
!! div(B) norms over all active blocks and writes them (master only)
!! to datadir/diagnostics.dat
module diagnostics

  implicit none

  !> Compute diagnostics every diag_every iterations
  integer, parameter :: diag_every = 10

contains

!===============================================================================

!> @brief Computes and writes global conserved quantities and div(B) norms
!> @details Must be called by ALL ranks: it exchanges boundaries and does
!! MPI reductions. Integrals use U in physical cells; div(B) uses PRIM
!! with refreshed ghost cells.
subroutine conservation_diagnostics ()

  use parameters
  use globals
  use amr,        only : meshlevel
#ifdef BFIELD
  use sources,    only : divergence_B
  use flux_ct,    only : bf_shared_err, check_cf_faces, cf_face_err
#endif
  implicit none

  ! 1:mass 2:passive 3-5:momentum 6:E 7-9:B 10:Emag 11:|divB|dx dV 12:|B| dV
  integer, parameter :: nsum = 12
  real    :: loc_sum(nsum), glob_sum(nsum), loc_max, glob_max,                 &
             loc_mis, glob_mis, glob_shared, glob_cf
  integer :: loc_nlev(maxlev), glob_nlev(maxlev)
  integer :: nb, bID, lev, i, j, k, unitd
  real    :: dV, divB, bmag
  logical :: exists
  character(len=256) :: fname

  if (mod(it, diag_every) /= 0) return

  loc_sum(:)  = 0.0
  loc_max     = 0.0
  loc_mis     = 0.0
  loc_nlev(:) = 0

  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID /= -1) then

      call meshlevel (bID, lev)
      loc_nlev(lev) = loc_nlev(lev) + 1
      dV = dx(lev)*dy(lev)*dz(lev)

      do k=1,ncells_z
        do j=1,ncells_y
          do i=1,ncells_x
            loc_sum(1)   = loc_sum(1)   + U(nb,1,i,j,k)*dV
            if (npassive >= 1) &
            loc_sum(2)   = loc_sum(2)   + U(nb,firstpas,i,j,k)*dV
            loc_sum(3:5) = loc_sum(3:5) + U(nb,2:4,i,j,k)*dV
            loc_sum(6)   = loc_sum(6)   + U(nb,5,i,j,k)*dV
#ifdef BFIELD
            bmag = sqrt(U(nb,6,i,j,k)**2 + U(nb,7,i,j,k)**2 + U(nb,8,i,j,k)**2)
            loc_sum(7:9) = loc_sum(7:9) + U(nb,6:8,i,j,k)*dV
            loc_sum(10)  = loc_sum(10)  + 0.5*bmag**2*dV
            if (enable_flux_ct) then
              divB = ( BF(nb,1,i,j,k) - BF(nb,1,i-1,j,k) ) / dx(lev)           &
                   + ( BF(nb,2,i,j,k) - BF(nb,2,i,j-1,k) ) / dy(lev)           &
                   + ( BF(nb,3,i,j,k) - BF(nb,3,i,j,k-1) ) / dz(lev)
            else
              call divergence_B (nb, lev, i, j, k, divB)
            end if
            loc_sum(11)  = loc_sum(11)  + abs(divB)*dx(lev)*dV
            loc_sum(12)  = loc_sum(12)  + bmag*dV
            loc_max      = max( loc_max, abs(divB)*dx(lev) )
            if (enable_flux_ct) then
              loc_mis = max( loc_mis,                                          &
                abs(U(nb,6,i,j,k) - 0.5*(BF(nb,1,i-1,j,k) + BF(nb,1,i,j,k))),  &
                abs(U(nb,7,i,j,k) - 0.5*(BF(nb,2,i,j-1,k) + BF(nb,2,i,j,k))),  &
                abs(U(nb,8,i,j,k) - 0.5*(BF(nb,3,i,j,k-1) + BF(nb,3,i,j,k))) )
            end if
#endif
          end do
        end do
      end do

    end if
  end do

  if (enable_flux_ct) call check_cf_faces ()

#ifdef MPIP
  call mpi_reduce (loc_sum, glob_sum, nsum, mpi_real_kind, mpi_sum, master,    &
                   mpi_comm_world, ierr)
  call mpi_reduce (loc_max, glob_max, 1, mpi_real_kind, mpi_max, master,       &
                   mpi_comm_world, ierr)
  call mpi_reduce (loc_nlev, glob_nlev, maxlev, mpi_integer, mpi_sum, master,  &
                   mpi_comm_world, ierr)
  call mpi_reduce (loc_mis, glob_mis, 1, mpi_real_kind, mpi_max, master,       &
                   mpi_comm_world, ierr)
  call mpi_reduce (cf_face_err , glob_cf, 1, mpi_real_kind, mpi_max, master,   &
                   mpi_comm_world, ierr)
#else
  glob_sum  = loc_sum
  glob_max  = loc_max
  glob_nlev = loc_nlev
  glob_mis  = loc_mis
  glob_cf   = cf_face_err
#endif

#ifdef BFIELD
#ifdef MPIP
  call mpi_reduce (bf_shared_err, glob_shared, 1, mpi_real_kind, mpi_max,      &
                   master, mpi_comm_world, ierr)
#else
  glob_shared = bf_shared_err
#endif
  if (enable_flux_ct .and. rank == master) write(logu,'(a,i8,3es12.4)')        &
       " flux-CT: |B_cell-<B_face>|, shared-face diff at it ", it,             &
       glob_mis, glob_shared, glob_cf
  bf_shared_err = 0.0
#else
  glob_shared = 0.0
#endif

  if (rank == master) then
    fname = trim(diagdir)//'diagnostics.dat'
    inquire (file=trim(fname), exist=exists)
    open (newunit=unitd, file=trim(fname), status='unknown', &
          position='append', action='write')
    if (.not.exists) write(unitd,'(a)') '# it time mass passive px py pz E '// &
      'Bx By Bz Emag <|divB|dx>/<|B|> max|divB|dx nblocks(lev=1..maxlev)'
    write(unitd,'(i8,13es23.15,*(1x,i6))') it, time, glob_sum(1:10),       &
          glob_sum(11)/max(glob_sum(12),1.0e-30), glob_max, glob_nlev
    close(unitd)
  end if

end subroutine conservation_diagnostics

!===============================================================================

!> @brief Temporary test: are face and edge/corner ghost cells filled correctly?
!> @details Zeroes all ghost cells of the local blocks, calls boundary(2,U),
!! and checks the velocity against the t=0 Orszag-Tang solution
!! (vx = -sin 2pi y, vy = sin 2pi x) at every ghost cell. Reports, for face
!! ghosts and for edge/corner ghosts: number of unfilled cells (rho = 0)
!! and the maximum velocity error.
subroutine test_edge_ghosts ()

  use parameters
  use globals
  use boundaries, only : boundary
  use amr,        only : cellPos
  implicit none

  integer :: nb, bID, i, j, k, nout
  real    :: loc(4), glob(4), x, y, z, err
  real, parameter :: twopi = 2.0*pi

  ! Wipe all ghost cells
  do nb=1,nbMaxProc
    if (localBlocks(nb) == -1) cycle
    do k=nzmin,nzmax
      do j=nymin,nymax
        do i=nxmin,nxmax
          nout = count([ (i<1).or.(i>ncells_x), (j<1).or.(j>ncells_y),        &
                         (k<1).or.(k>ncells_z) ])
          if (nout > 0) U(nb,:,i,j,k) = 0.0
        end do
      end do
    end do
  end do

  call boundary (2, U)

  ! loc(1:2): unfilled count face / edge+corner; loc(3:4): max error
  loc(:) = 0.0
  do nb=1,nbMaxProc
    bID = localBlocks(nb)
    if (bID == -1) cycle
    do k=nzmin,nzmax
      do j=nymin,nymax
        do i=nxmin,nxmax
          nout = count([ (i<1).or.(i>ncells_x), (j<1).or.(j>ncells_y),        &
                         (k<1).or.(k>ncells_z) ])
          if (nout == 0) cycle
          if (U(nb,1,i,j,k) == 0.0) then
            if (nout == 1) then
              loc(1) = loc(1) + 1.0
            else
              loc(2) = loc(2) + 1.0
            end if
          else
            call cellPos (bID, i, j, k, x, y, z)
            err = max( abs(U(nb,2,i,j,k)/U(nb,1,i,j,k) + sin(twopi*y)),      &
                       abs(U(nb,3,i,j,k)/U(nb,1,i,j,k) - sin(twopi*x)) )
            if (nout == 1) then
              loc(3) = max(loc(3), err)
            else
              loc(4) = max(loc(4), err)
            end if
          end if
        end do
      end do
    end do
  end do

#ifdef MPIP
  call mpi_reduce (loc(1:2), glob(1:2), 2, mpi_real_kind, mpi_sum, master,    &
                   mpi_comm_world, ierr)
  call mpi_reduce (loc(3:4), glob(3:4), 2, mpi_real_kind, mpi_max, master,    &
                   mpi_comm_world, ierr)
#else
  glob = loc
#endif
  if (rank == master) then
    write(logu,'(a,2f10.0)') " test_edge_ghosts: unfilled face / edge+corner:", glob(1:2)
    write(logu,'(a,2es12.4)') " test_edge_ghosts: max v error face / edge+corner:", glob(3:4)
  end if

end subroutine test_edge_ghosts

!===============================================================================

end module diagnostics
