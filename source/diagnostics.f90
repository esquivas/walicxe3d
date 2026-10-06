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
  use boundaries, only : boundary
  use hydro_core, only : calcPrimsAll
  use sources,    only : divergence_B
#endif
  implicit none

  ! 1:mass 2:passive 3-5:momentum 6:E 7-9:B 10:Emag 11:|divB|dx dV 12:|B| dV
  integer, parameter :: nsum = 12
  real    :: loc_sum(nsum), glob_sum(nsum), loc_max, glob_max
  integer :: loc_nlev(maxlev), glob_nlev(maxlev)
  integer :: nb, bID, lev, i, j, k, unitd
  real    :: dV, divB, bmag
  logical :: exists
  character(len=256) :: fname

  if (mod(it, diag_every) /= 0) return

#ifdef BFIELD
  ! Refresh 1-deep ghost cells (U ghosts are stale after the step)
  call boundary (1, U)
  call calcPrimsAll (U, PRIM, CELLS_GHOST)
#endif

  loc_sum(:)  = 0.0
  loc_max     = 0.0
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
            bmag = sqrt( U(nb,6,i,j,k)**2 + U(nb,7,i,j,k)**2 + U(nb,8,i,j,k)**2 )
            loc_sum(7:9) = loc_sum(7:9) + U(nb,6:8,i,j,k)*dV
            loc_sum(10)  = loc_sum(10)  + 0.5*bmag**2*dV
            call divergence_B (nb, lev, i, j, k, divB)
            loc_sum(11)  = loc_sum(11)  + abs(divB)*dx(lev)*dV
            loc_sum(12)  = loc_sum(12)  + bmag*dV
            loc_max      = max( loc_max, abs(divB)*dx(lev) )
#endif
          end do
        end do
      end do

    end if
  end do

#ifdef MPIP
  call mpi_reduce (loc_sum, glob_sum, nsum, mpi_real_kind, mpi_sum, master, &
                   mpi_comm_world, ierr)
  call mpi_reduce (loc_max, glob_max, 1, mpi_real_kind, mpi_max, master,    &
                   mpi_comm_world, ierr)
  call mpi_reduce (loc_nlev, glob_nlev, maxlev, mpi_integer, mpi_sum, master,&
                   mpi_comm_world, ierr)
#else
  glob_sum  = loc_sum
  glob_max  = loc_max
  glob_nlev = loc_nlev
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

end module diagnostics
