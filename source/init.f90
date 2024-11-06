!===============================================================================
!> @file initmain.f90
!> @brief Basic allocations and initializations
!> @author Juan C. Toledo
!> @date 3/Jun/2011

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

!> @brief Code Initialization
!> @details This module contains subroutines and utilities to initialize the
!> code.

module init

  implicit none

contains

!===============================================================================

!> @brief Code-wide initializations
!> @details Initializes MPI, allocates big arrays, loads cooling data,
!! reports compilation parameters
subroutine initmain ()

  use parameters
  use globals
  !use userconds, only : initializeUserModule
  use tictoc
  use clean_quit, only : clean_abort
  use cooling_schure, only : loadcooldata_schure
  use coolingModule,  only : loadcooldata, loadcooldata_metal
  implicit none

  integer :: nps, istat, inext
  integer :: mark, l
  logical :: existing, success
  real :: totalsize
  character(3) :: rankstr
  character(1) :: slash

  ! =================================

  ! MPI Initialization
#ifdef MPIP
  call mpi_init (ierr)
  call mpi_comm_size (mpi_comm_world, nps, ierr)
  call mpi_comm_rank (mpi_comm_world, rank, ierr)
  if (nps.ne.nProcs) then
    if (rank.eq.master) then
      write(*,'(a,i0,a)') "Actual number of processors ( ", nps, " ) does not match value specified "
      write(*,'(a,i0,a)') "in parameters.f90 ( ", nProcs, " ) !"
      write(*,'(a)') "***ABORTING***"
      write(*,*) ""
    end if
    call clean_abort (ERROR_WRONG_NPROCS)
  end if
#else
  rank = master
#endif

  ! =================================

  ! Start logging
  if (logged) then

    ! Generate logfile name
    l = len_trim(logdir)
    if (logdir(l:l).ne.'/') then
      slash = '/'
    else
      slash = ''
    end if

    write(rankstr,'(i3.3)') rank
    logfile = trim(logdir) // trim(slash) // "rank" // rankstr // ".log"

    ! create log directory in case it's missing
#ifdef gfortran
    inquire(     file=trim(logdir)//'/.', exist=existing )
#endif
#ifdef ifort
    inquire(directory=trim(logdir),exist=existing)
#endif
    if (.not.existing .and. rank ==master ) then
      write(*,'(a)') "Could not find logdir, creating it anew"
      call system('mkdir -p ' // trim(logdir) )
    end if
    call mpi_barrier(mpi_comm_world, ierr)

    ! Check if file already exists
    success = .false.
    inext = 1
    do while(.not.success)
      inquire (file=logfile, exist=existing)
      if (existing) then
        ! file exists - append a number, try again
        write(logfile,'(a,a,a,a,i0,a)') trim(logdir) // trim(slash), &
          "rank", rankstr, "-", inext, ".log"
        inext = inext + 1
      else
        ! Open logfile
        logu = 10 + nProcs + rank
        open (unit=logu, file=logfile, status='unknown', position="append", &
          iostat=istat)
        success = .true.
      end if
    end do

    if (istat.ne.0) then
      if (rank.eq.master) then
        write(*,'(a)') "Could not open the log file!"
        write(*,'(a,a)') "Tried to open: ", logfile
        write(*,'(a,a,a)') "Does the logdir '", trim(logdir), "' exist ?"
        write(*,'(1x,a)') "***ABORTING***"
      end if
      close(logu)
      call clean_abort (ERROR_NO_LOGFILE)
    end if

  else

    logu = 6

  end if

  ! =================================
  if (verbosity > 3) call tic(mark)
  if ( (verbosity > 0).and.(logged.or.(rank==master)) ) then
    write(logu,*) ""
    write(logu,'(a)') "================================================================================"
    if (logged) then
      write(logu,*) ""
      write(logu,'(1x,a,a)') "> Started logfile ", trim(logfile)
    end if
    write(logu,*) ""
    write(logu,'(1x,a)') "> Initializing ... "

    if (dowarm) then
      write(logu,'(1x,a)') "A WARM START has been scheduled"
    end if
  endif
  ! Get hostname
  call HostNm(host)

  if ( (verbosity > 0).and.(logged.or.(rank==master)) ) then
    if(logged.or.rank.eq.master) then
      write(logu,*)
      write(logu,'(1x,a)') "***************************************************"
      write(logu,'(1x,a)') "*   __    __      _ _               _____ ____    *"
      write(logu,'(1x,a)') "*  / / /\ \ \__ _| (_) _____  _____|___ /|  _ \   *"
      write(logu,'(1x,a)') "*  \ \/  \/ / _` | | |/ __\ \/ / _ \ |_ \| | | |  *"
      write(logu,'(1x,a)') "*   \  /\  / (_| | | | (__ |  |  __/___) | |_| |  *"
      write(logu,'(1x,a)') "*    \/  \/ \__,_|_|_|\___/_/\_\___|____/|____/   *"
      write(logu,'(1x,a)') "*                                                 *"
      write(logu,'(1x,a)') "*         Version 1.2 ($Revision:: 76  $)         *"
      write(logu,'(1x,a)') "*                                                 *"
#ifdef MPIP
      write(logu,'(1x,a,i3,a)') "*        Running with MPI on ", nProcs , " processors       *"
#else
      write(logu,'(1x,a)') "*               Running serially                  *"
#endif
      write(logu,'(1x,a,a,a)') "*        Hostname: ", host, "                *"
      write(logu,'(1x,a,a,a)') "*        Start: ", STAMP(), "            *"
      write(logu,'(1x,a)') "*                                                 *"
      write(logu,'(1x,a)') "***************************************************"
    end if
  end if
  if (verbosity > 3) call tic (start_mark)
  if (verbosity > 0) then
    write(logu,*) ""
    write(logu,'(1x,a,i0,a)') "Processor ", rank, " ready."
  end if
  if ( (verbosity > 0).and.(logged.or.(rank==master)) ) then
    ! =================================

    ! Report parameters. For warm starts, critical parameters are also verified.

    write(logu,*) ""
    write(logu,'(1x,a)') "============================================"
    write(logu,'(1x,a)') " Doing basic initializations ..."
    write(logu,'(1x,a)') "============================================"

    write(logu,*) ""
    write(logu,'(1x,a,es12.5,a,f7.3,a)') "Simulation box x-size:  ", xphystot, " cm / ", xphystot/PC, " pc"
    write(logu,'(1x,a,es12.5,a,f7.3,a)') "Simulation box y-size:  ", yphystot, " cm / ", yphystot/PC, " pc"
    write(logu,'(1x,a,es12.5,a,f7.3,a)') "Simulation box z-size:  ", zphystot, " cm / ", zphystot/PC, " pc"

    write(logu,*) ""
    if (mesh_method.eq.MESH_AUTO) then
      write(logu,'(1x,a)') "> Mesh initialization method: Automatic"
      write(logu,'(1x,a,i0)') "Max-level cells along x:  ", p_maxcells_x
      write(logu,'(1x,a,i0)') "Max-level cells along y:  ", p_maxcells_y
      write(logu,'(1x,a,i0)') "Max-level cells along z:  ", p_maxcells_z
    else if (mesh_method.eq.MESH_MANUAL) then
      write(logu,'(1x,a)') "> Mesh initialization method: Manual"
      write(logu,'(1x,a,i0)') " Number of root blocks along x:  ", p_nbrootx
      write(logu,'(1x,a,i0)') " Number of root blocks along y:  ", p_nbrooty
      write(logu,'(1x,a,i0)') " Number of root blocks along z:  ", p_nbrootz
      write(logu,'(1x,a,i0)') " Number of refinement levels:  ", p_maxlev
    end if

    write(logu,*) ""
    write(logu,'(1x,a,f7.1,a)') "Allowed RAM per process: ", RAM_per_proc, " MB"

    write(logu,*) ""
    write(logu,'(1x,a)') "> Mesh parameters"
    write(logu,'(1x,a,i0)') "Max blocks per processor:  ", nbMaxProc
    write(logu,'(1x,a,i0)') "Cells per block along x:   ", ncells_x
    write(logu,'(1x,a,i0)') "Cells per block along y:   ", ncells_y
    write(logu,'(1x,a,i0)') "Cells per block along z:   ", ncells_z
    write(logu,'(1x,a,i0)') "Ghost cell layer depth:    ", nghost
    write(logu,'(1x,a,f5.2)') "Refinement Threshold:  ", refineThres
    write(logu,'(1x,a,f5.2)') "Coarsening Threshold:  ", coarseThres

    write(logu,*) ""
    write(logu,'(1x,a)') "> Boundary Conditions"
    write(logu,'(1x,a,a)') "Left:   ", trim(bcname(bc_left))
    write(logu,'(1x,a,a)') "Right:  ", trim(bcname(bc_right))
    write(logu,'(1x,a,a)') "Front:  ", trim(bcname(bc_front))
    write(logu,'(1x,a,a)') "Back:   ", trim(bcname(bc_back))
    write(logu,'(1x,a,a)') "Bottom: ", trim(bcname(bc_bottom))
    write(logu,'(1x,a,a)') "Top:    ", trim(bcname(bc_top))

    write(logu,*) ""
    write(logu,'(1x,a)') "> Data Output"
    write(logu,'(1x,a,a)') "Data directory:     ", datadir
    write(logu,'(1x,a,a)') "Logging directory:  ", logdir
    write(logu,'(1x,a,a)') "Datafile template:  ", blockstpl
    write(logu,'(1x,a,a)') "Gridfile template:  ", gridtpl
    write(logu,'(1x,a,a)') "Statefile template: ", statetpl
    if (output_bin) then
      write(logu,'(1x,a)') "Output in internal format is ON"
    else
      write(logu,'(1x,a)') "Output in internal format is OFF"
    end if
    if (output_vtk) then
      write(logu,'(1x,a)') "Output in VisIt format is ON"
    else
      write(logu,'(1x,a)') "Output in VisIt format is OFF"
    end if
    if (units_type.eq.CODE_UNITS) write(logu,'(1x,a)') "Data dumped in Code Units"
    if (units_type.eq.PHYS_UNITS) write(logu,'(1x,a)') "Data dumped in Physical Units"

#ifdef BFIELD
    if ( mhd ) then
      write(logu,*) ""
      write(logu,'(1x,a)') "> Active magnetic field ENABLED (MHD)"
      write(logu,*) ""
      if (eight_wave)      &
      write(logu,'(1x,a)') "> div(B) constrained with 8 wave method'"
      if (enable_flux_cd) &
      write(logu,'(1x,a)') "> div(B) constrained with flux CD method'"
    else
      write(logu,*) ""
      write(logu,'(1x,a)') "> Passive magnetic field ENABLED"
    end if
#endif

    write(logu,*) ""
    write(logu,'(1x,a)') "> Hydro Solver"
    write(logu,'(1x,a,a)') "Type: ", trim(solvername(solver_type))
    if ( (solver_type.eq.SOLVER_HLL ).or.(solver_type.eq.SOLVER_HLLC).or. &
         (solver_type.eq.SOLVER_HLLE).or.(solver_type.eq.SOLVER_HLLD) ) then
      ! Check that two ghost cells are used for second-order solvers
      if (nghost.ne.2) then
        write(logu,*) "This solver requires TWO ghost cells!"
        write(logu,*) "Modify parameters.f90"
        write(logu,*) "***ABORTING***"
      end if
      write(logu,'(1x,a,a)') "Limiter: ", trim(limitername(limiter_type))
    end if
    write(logu,'(1x,a,i0)') "Hydro equations:  ", neqhydro
    write(logu,'(1x,a,i0)') "MHD equations:    ", neqmhd
    write(logu,'(1x,a,i0)') "Passive scalars:  ", npassive
    write(logu,'(1x,a,i0)') "Total equations:  ", neqtot
    write(logu,'(1x,a,f6.3)') "CFL parameter:  ", CFL
    write(logu,'(1x,a,f6.3)') "Artificial viscosity:  ", visc_eta

    write(logu,*) ""
    write(logu,'(1x,a)') "> Gas Parameters"
    write(logu,'(1x,a,f6.3)') "gamma = ", gamma
    write(logu,'(1x,a,f6.3)') "mu0   = ", mu0
    write(logu,'(1x,a,f6.3)') "mui   = ", mui
    write(logu,'(1x,a,f8.1)') "ion_thres = ", ion_thres

    ! Report unit scalings
    write(logu,*) ""
    write(logu,'(1x,a)') "> Unit scalings (1 code unit = ?)"
    write(logu,'(1x,a,es12.5,a)') "Length:   ", l_sc, " cm"
    write(logu,'(1x,a,es12.5,a)') "Density:  ", d_sc, " g cm^-3"
    write(logu,'(1x,a,es12.5,a)') "Velocity: ", v_sc, " cm s^-1"
    write(logu,'(1x,a,es12.5,a)') "Pressure: ", p_sc, " erg cm^-3"
    write(logu,'(1x,a,es12.5,a)') "Time:     ", t_sc, " s"

  end if
  ! Radiative cooling
  write(logu,*) ""
  if (cooling_type.eq.COOL_NONE) then
    if ( (verbosity > 0).and.(logged.or.(rank==master)) ) &
       write(logu,'(1x,a)') "> Radiative cooling is OFF"
  else
    if ( (verbosity > 0).and.(logged.or.(rank==master)) ) then
      if (verbosity > 0) write(logu,'(1x,a)') "> Radiative cooling is ON"
      if (verbosity > 0) write(logu,'(2x,a,a)') "Cooling Table: ", trim(cooling_file)
    end if
    if (cooling_type.eq.COOL_TABLE) then
      call loadcooldata ()
    else if (cooling_type.eq.COOL_TABLE_METAL) then
      if (npassive.lt.1) then
        write(logu,'(1x,a)') "At least one passive scalar is needed for &
        &metallicity dependent cooling!"
        write(logu,'(1x,a)') "Set npassive to at least 1 in parameters.f90"
        write(logu,'(1x,a)') "***ABORTING***"
        call clean_abort(ERROR_NOT_ENOUGH_PASSIVES)
        end if
      call loadcooldata_metal ()
    else if (cooling_type.eq.COOL_SCHURE) then
      write(logu,'(1x,a)') "Cooling w/Schure et al. XXXX table"
      call loadcooldata_schure ()
    end if
  end if

  ! Allocate memory and initialize big data arrays
  if ( (verbosity > 0).and.(logged.or.(rank==master)) ) then
    write(logu,*) ""
    write(logu,'(1x,a)') "> Allocating memory for big arrays ..."
    write(logu,*) ""
  end if

  ! Sanity check: abort if RAM will be insufficient to allocate big arrays
!  totMem = neqtot*(nxmax-nxmin+1)*(nymax-nymin+1)*(nymax-nymin+1)*(nbmaxProc*3+6)
!  if (mpi_real_kind.eq.MPI_DOUBLE_PRECISION) then
!    totMem = totMem*8
!  else if (mpi_real_kind.eq.MPI_REAL) then
!    totMem = totMem*4
!  end if
!  totMem = totMem/1024.0/1024.0
!  if (totMem.gt.maxMemProc) then
!    write(logu,'(a)') "Required memory would surpass maximum allowed memory!"
!    write(logu,'(a,f6.1,a)') "Required: ", totMem, " MB"
!    write(logu,'(a,f6.1,a)') "Allowed:  ", maxMemProc, " MB"
!    write(logu,'(a,f6.1,a)') "Either decrease nbMaxProc or increase maxMemProc in parameters.f90"
!    write(logu,*) "***ABORTING***"
!    call clean_abort (ERROR_INSUFFICIENT_RAM)
!  else
!    write(logu,'(1x,a,f6.1,a)') "Estimated required memory for big arrays: ", totMem, " MB"
!  end if

  ! Big data arrays

  allocate( U(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax), stat=istat)
  U(:,:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for U array!"
    write(logu,*) "***ABORTING***"
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  allocate( UP(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax), stat=istat)
  UP(:,:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for UP array!"
    write(logu,*) "***ABORTING***"
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  allocate( PRIM(nbMaxProc, neqtot, nxmin:nxmax, nymin:nymax, nzmin:nzmax), stat=istat)
  PRIM(:,:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for P array!"
    write(logu,*) "***ABORTING***"
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  ! Allocate memory and initialize auxiliary data arrays

  allocate( FC(neqtot, nxmin:nxmax, nxmin:nxmax, nxmin:nxmax), stat=istat)
  FC(:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for FC array!"
    write(logu,*) "***ABORTING***"
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  allocate( GC(neqtot, nxmin:nxmax, nxmin:nxmax, nxmin:nxmax), stat=istat)
  GC(:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for GC array!"
    write(logu,*) "***ABORTING***"
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  allocate( HC(neqtot, nxmin:nxmax, nxmin:nxmax, nxmin:nxmax), stat=istat)
  HC(:,:,:,:) = 0.0
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for HC array!"
    write(logu,*) "***ABORTING***"
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  ! Allocate and initialize global and local block registries

  allocate ( globalBlocks(nbMaxGlobal), stat=istat )
  globalBlocks(:) = -1
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for global block registry!"
    write(logu,*) "***ABORTING***"
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  allocate ( localBlocks(nbMaxProc), stat=istat )
  localBlocks(:) = -1
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for local block registry!"
    write(logu,*) "***ABORTING***"
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  ! Allocate and initialize global and local refinement flags lists

  allocate ( refineFlags(nbMaxGlobal), stat=istat )
  refineFlags(:) = -1
  if (istat.ne.0) then
    write(logu,'(a)') "Couldn't allocate memory for refinement flags!"
    write(logu,*) "***ABORTING***"
    call clean_abort (ERROR_NOALLOC_BIGARR)
  end if

  ! Calculate global block index range
  nbmin = rank*nbMaxProc + 1
  nbmax = rank*nbMaxProc + nbMaxProc

  if ( (verbosity > 0).and.(logged.or.(rank==master)) ) then
    ! Report allocation results
    write(logu,'(1x,a)') "Successfully allocated big arrays."
    write(logu,'(1x,a,f6.1,a)') "U:      ", sizeof(U)/(1024.*1024.), " MB"
    write(logu,'(1x,a,f6.1,a)') "UP:     ", sizeof(UP)/(1024.*1024.), " MB"
    write(logu,'(1x,a,f6.1,a)') "PRIM:   ", sizeof(PRIM)/(1024.*1024.), " MB"
    write(logu,'(1x,a,f6.1,a)') "Fluxes: ", (sizeof(FC)+sizeof(GC)+sizeof(HC))&
                                          /(1024.*1024.), " MB"
    totalsize = (sizeof(U)+sizeof(UP)+sizeof(PRIM)+sizeof(FC)+sizeof(GC)  &
    +sizeof(HC))/(1024.*1024.)
    write(logu,'(1x,a,f7.1,a,f7.1,a)') "Total: ", totalsize, " MB / ", &
    RAM_per_proc, "MB"
  endif
  if(.not.logged.and.rank==master) write(logu,'(a)') "> The above figures are per processor"
  ! =================================

  ! Initialize simulation state - (warm start does this later)
  time = 0.0
  it = 0
  nextout = 0

  ! initialize vaiables and modules defined by user
  !call initializeUserModule()

  ! =================================
  if (verbosity > 3) then
    write(logu,'(1x,a,a)') ""
    write(logu,'(1x,a,a)') "> Performed initializations and allocated big arrays in ", nicetoc(mark)
    write(logu,*) ""
  end if
  ! Barrier
  call mpi_barrier (mpi_comm_world, ierr)

end subroutine initmain

!==============================================================================
!> @brief Sets the initial flow conditions
!> @details This wrapper subroutine is called to set the initial flow
!! conditions on the mesh. It will first initialize the grid with a
!! uniform IC and then call userIC(), where user-defined custom
!! initial conditions can be specified.
subroutine initflow ()

  use parameters
  use globals
  use userconds
  implicit none
  if ( (verbosity > 0).and.(logged.or.(rank==master)) ) then
    write(logu,*) ""
    write(logu,*) "============================================"
    write(logu,'(1x,a)') " Setting Initial Conditions  ..."
    write(logu,*) "============================================"
    write(logu,*) ""
  end if
  ! IC (defined in user.f90)
  call setInitialCondition (U)

  if ( (verbosity > 0).and.(logged.or.(rank==master)) ) then
    write(logu,*) ""
    write(logu,'(1x,a)') "> Done setting ICs"
    write(logu,*) ""
    write(logu,'(1x,a)') "... and now for something completely different ..."
  end if

end subroutine initflow

!==============================================================================

!> @brief Performs a warm start. This requires reading the state file specified
!! in parameters.f90 and loading hydro data.
subroutine warmstart ()

  use parameters
  use globals
  use clean_quit, only : clean_abort
  use utils,      only : replace
  implicit none

  integer :: unitin, istat, noutput, l, nb, nblocks
  character(256) :: datadir_old
  character(256) :: blocksfile
  character(1) :: slash
  character(4) :: noutstr
  character(3) :: rankstr

  if ( (verbosity > 0).and.(logged.or.(rank==master)) ) then
    write(logu,*) ""
    write(logu,*) "============================================"
    write(logu,'(1x,a)') " Performing warm start ..."
    write(logu,*) "============================================"
    write(logu,*) ""
  end if

  ! Open state file
  if ( (verbosity > 0).and.(logged.or.(rank==master)) )  &
    write(logu,'(1x,a,a,a)') "Reading state file '", trim(warm_file), "' ..."
  unitin = 10 + rank
  open (unit=unitin, file=warm_file, status='old', iostat=istat)
  if (istat.ne.0) then
    write(logu,'(a,a,a)') "Could not open the file '", trim(warm_file), "' !"
    close(unitin)
    call clean_abort (ERROR_WARM_FILE_OPEN)
  end if

  ! Read simulation state variables and datadir
  read(unitin,'(es22.15,i8,i5)') time, it, noutput
  read(unitin,'(a)') datadir_old
  if ( (verbosity > 0).and.(logged.or.(rank==master)) ) then
    write(logu,*) ""
    write(logu,'(1x,a,i0,a)') "> Restarting simulation from output ", noutput, " ..."
    write(logu,'(1x,a,i0)') "Last iteration = ", it
    write(logu,'(1x,a,es22.15)') "Current time = ", time
    write(logu,*) ""
  end if

  time = time/t_sc
  nextout = noutput + 1

  close(unitin)

  ! Generate filename based on templates for data file
  l = len_trim(datadir_old)
  if (datadir(l:l).ne.'/') then
    slash = '/'
  else
    slash = ''
  end if
  write(rankstr,'(I3.3)') rank
  write(noutstr,'(I4.4)') noutput
  blocksfile = blockstpl
  call replace (blocksfile, 'XXX', rankstr)
  call replace (blocksfile, 'YYYY', noutstr)
  write(blocksfile,'(a)') trim(datadir_old) // trim(slash) // trim(blocksfile) // ".bin"

  ! Open data file
  if (verbosity > 2) write(logu,'(1x,a,a,a)') "Reading data file '", trim(blocksfile), "' ..."
  unitin = 10 + rank
  open (unit=unitin, file=blocksfile, status='old', access='stream', iostat=istat)
  if (istat.ne.0) then
    write(logu,'(a,a,a)') "Could not open the file '", trim(blocksfile), "' !"
    write(logu,'(a,a,a)') "Does the datadir '", trim(datadir_old), "' exist?"
    close(unitin)
    call clean_abort (ERROR_WARM_FILE_OPEN)
  end if

  ! Read data
  localBlocks(:) = -1
  nblocks = 0
  read(unitin) nbLocal
  do nb=1,nbLocal
    read(unitin) localBlocks(nb)
    if (verbosity > 2) write(logu,'(1x,a,i0,a)') "Loading block ", localBlocks(nb), " ..."
    read(unitin) U(nb,:,1:ncells_x,1:ncells_y,1:ncells_z)
    nblocks = nblocks + 1
  end do

  if (nbLocal.ne.nblocks) then
    write(logu,'(a,i0,a,i0,a)') "Read ", nblocks, " blocks from file, but expected ", nbLocal, " !"
    write(logu,'(a)') "***Aborting***"
    call clean_abort (ERROR_WARM_READ_BLOCKS)
  else
    write(logu,'(1x,a,i0,a,i0,a)') "Loaded ", nblocks, " blocks; expected ", nbLocal, "."
  end if

  close (unitin)

end subroutine warmstart

!===============================================================================
!===============================================================================

!> @brief Initializes root blocks and refines the mesh to prepare it for ICs
!> @details This is a main-level routine. Calculates the mesh geometry
!! (how many root blocks) and number of levels to meet the required resolution,
!! and refines the IC zone.
subroutine basegrid ()

  use parameters
  use globals
  use tictoc
  use clean_quit, only : clean_abort
  use amr, only : syncBlockLists
  implicit none

  ! Local variables
  integer :: maxlevx, maxlevy, maxlevz, smalldim
  integer :: nb, bID, ilev, x, y, z, mark
  real :: smallsize

  if (verbosity > 3) call tic(mark)

  if ( (verbosity > 0).and.(logged.or.(rank==master)) ) then
    write(logu,*) ""
    write(logu,*) "============================================"
    write(logu,'(1x,a)') " Initializing the base grid ..."
    write(logu,*) "============================================"
    write(logu,*) ""
  end if

  ! Every rank calculates the base grid geometry, based on the selected method
  if (mesh_method.eq.MESH_AUTO) then

    maxcells_x = p_maxcells_x
    maxcells_y = p_maxcells_y
    maxcells_z = p_maxcells_z

    ! Number of root blocks
    smalldim = min(maxcells_x,maxcells_y,maxcells_z)
    nbrootx = maxcells_x / smalldim
    nbrooty = maxcells_y / smalldim
    nbrootz = maxcells_z / smalldim

    ! Calculate number of levels
    maxlevx = int(log(float(maxcells_x)/nbrootx/ncells_x)/log(2.)+1)
    maxlevy = int(log(float(maxcells_y)/nbrooty/ncells_y)/log(2.)+1)
    maxlevz = int(log(float(maxcells_z)/nbrootz/ncells_z)/log(2.)+1)
    maxlev = max(maxlevx,maxlevy,maxlevz)


  else if (mesh_method.eq.MESH_MANUAL) then

    maxlev = p_maxlev
    nbrootx = p_nbrootx
    nbrooty = p_nbrooty
    nbrootz = p_nbrootz

    ! Number of max-resolution cells
    maxcells_x = ncells_x * 2**(maxlev-1)
    maxcells_y = ncells_y * 2**(maxlev-1)
    maxcells_z = ncells_z * 2**(maxlev-1)

  end if

  ! The maximum allowable number of levels is 10 at the moment, due to the
  ! fact that bIDs are handled as 4-byte (signed) integers. With 10 levels,
  ! one can have up to 14 root blocks.
  if (maxlev.gt.10) then
    write(logu,*) ""
    write(logu,'(1x,a)') "Desired finest-level resolution would require more than 10 levels!"
    write(logu,'(1x,a)') "***ABORTING***"
    call clean_abort (ERROR_TOO_MANY_LEVS)
  end if

  ! Check if physical sizes conform to desired max-resolution cell numbers
  smallsize = min(xphystot, yphystot, zphystot)
  if ((maxcells_x/smalldim.ne.xphystot/smallsize).or.&
      (maxcells_y/smalldim.ne.yphystot/smallsize).or.&
      (maxcells_z/smalldim.ne.zphystot/smallsize)) then
    write(logu,*)
    write(logu,'(a)') "ERROR: The physical size of the simulation box has"//&
      " a different aspect ratio than the requested number of cells!"
     write(logu,'(a,f4.1,a,f4.1,a,f4.1)') "Physical size aspect ratio: ",&
      xphystot/smallsize, " : ", yphystot/smallsize, " : ", zphystot/smallsize
    write(logu,'(a,i3,a,i3,a,i3)') "Cell number aspect ratio:   ",&
      maxcells_x/smalldim, " : ", maxcells_y/smalldim, " : ", maxcells_z/smalldim
    write(logu,'(a)') "> Modify parameters.f90"
    call clean_abort (ERROR_BASEGRID_BAD_ASPECT)
  else
    if ( (verbosity > 0).and.(logged.or.(rank==master)) )  then
      write(logu,'(1x,a,i2,a,i2,a,i2)') "Grid geometry (root blocks): ", &
      nbrootx, " x ", nbrooty, " x ", nbrootz
      write(logu,'(1x,a,i3)') "Number of root blocks: ", nbrootx*nbrooty*nbrootz
      write(logu,'(1x,a,i2)') "Number of refinement levels: ", maxlev
    end if
  end if
  ! Allocate and initialize grid spacing at each level (code units)
  allocate( dx(maxlev) )
  allocate( dy(maxlev) )
  allocate( dz(maxlev) )
  do ilev=1,maxlev
     dx(ilev) = xphystot / (ncells_x * nbrootx * 2.0**(ilev-1)) / l_sc
     dy(ilev) = yphystot / (ncells_y * nbrooty * 2.0**(ilev-1)) / l_sc
     dz(ilev) = zphystot / (ncells_z * nbrootz * 2.0**(ilev-1)) / l_sc
  end do

  ! Calculate total number of blocks per level
  allocate( nblockslev(maxlev) )
  do ilev=1,maxlev
    nblockslev(ilev) = nbrootx*nbrooty*nbrootz*8**(ilev-1)
  end do

  ! Calculate number of blocks along each dimension at each level
  allocate( nbx(maxlev) )
  allocate( nby(maxlev) )
  allocate( nbz(maxlev) )
  do ilev=1,maxlev
    nbx(ilev) = nbrootx*2**(ilev-1)
    nby(ilev) = nbrooty*2**(ilev-1)
    nbz(ilev) = nbrootz*2**(ilev-1)
  end do

  ! Calculate the bIDs of the first and last blocks at each level
  allocate( minID(maxlev) )
  allocate( maxID(maxlev) )
  minID(1) = 1
  maxID(1) = nblockslev(1)
  do ilev=2,maxlev
    minID(ilev) = maxID(ilev-1) + 1
    maxID(ilev) = maxID(ilev-1) + nblockslev(ilev)
  end do

  nbLocal = 0
  nbActive = nbrootx*nbrooty*nbrootz

  ! ==========================
  ! The following is only performed in cold starts

  if (.not.dowarm) then

    ! Activate root blocks and initialize block registry
    if (verbosity > 2) write(logu,'(1x,a)') "> Creating root blocks ..."

    ! Calculate bID of root blocks and register them in master's block list
    if (rank.eq.master) then
      nb = 1
      do z = 1,nbrootz
        do y = 1,nbrooty
          do x = 1,nbrootx

            ! bID is sensitive to numbering order
            bID = 1+(x-1)+(y-1)*nbrootx+(z-1)*nbrootx*nbrooty

            ! Register block
            localBlocks(nb) = bID
            nb = nb + 1
            nbLocal = nbLocal + 1

          end do
        end do
      end do
    end if

    ! Synchronize block lists
    call syncBlockLists ()

  else

    if (verbosity > 2) write(logu,*) ""
    if (verbosity > 2) write(logu,'(1x,a)') "> Skipping root block creation (warm start)"

  end if

  ! ==========================

  if (verbosity > 3) then
    write(logu,*) ""
    write(logu,'(1x,a,a)') "> Created base grid in ", nicetoc(mark)
    write(logu,*) ""
  end if

end subroutine basegrid

!===============================================================================

end module init
