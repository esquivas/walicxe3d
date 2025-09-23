module hrate

  implicit none

contains

  !=====================================================================
  !> @brief High level wrapper to apply cooling
  !> @details High level wrapper to apply cooling
  !! @n  parametrized cooling curve, uses the ionization state of
  !! hydrogen and ties the O I and II to it
  subroutine updateNeutralFraction()

    use parameters, only : t_sc, nbMaxProc, verbosity
    use globals   , only : dt, localBlocks, logu
    use tictoc
    implicit none
    real    :: dt_seconds
    integer :: mark, nb, bID

    dt_seconds = dt*t_sc

    if (verbosity > 3) call tic(mark)
    if (verbosity > 1) then
      write(logu,*) ""
      write(logu,'(1x,a)') "============================================"
      write(logu,'(1x,a)') " Updating Neutral Fraction ..."
      write(logu,'(1x,a)') "============================================"
      write(logu,*) ""
    end if

    do nb=1,nbMaxProc
      bID = localBlocks(nb)
      if (bID.ne.-1) then

        call solve_h_rate(nb, dt_seconds)

      end if
    end do

    if (verbosity > 3) write(logu,*) ""
    if (verbosity > 3) write(logu,'(1x,a,a)') "> Update completed in ",        &
                                               nicetoc(mark)


  end subroutine updateNeutralFraction

  !======================================================================
  !> @brief calculates the recombination rate (case B)
  !> @details calculates the recombination rate (case B)
  !> @param real8 [in] T : Temperature K
  function alpha(T)

    implicit none

    real (kind=8) :: alpha
    real (kind=8), intent(in) :: T

    alpha=2.55d-13*(1.d4/T)**0.79

  end function alpha

  !======================================================================
  !> @brief calculates the collisional ionization rate
  !> @details calculates the collisional ionization rate
  !> @param real8[in] T : Temperature K
  function colf(T)

    implicit none

    real (kind=8) :: colf
    real (kind=8), intent(in) :: T

    colf=5.83d-11*sqrt(T)*exp(-157828./T)

  end function colf

  !=======================================================================
  !> @brief Updates the ionization fraction using Hrate eqn.
  !> @param real [in] dt        : timestep (seconds)
  !> @param real [in] uu(neq)   : conserved variables in one cell
  !> @param real [in] prim(neq) : primitive variables in one cell
  subroutine solve_h_rate(bIndx,dt)

    use parameters,  only : ncells_x, ncells_y, ncells_z, inH0, mu0,       &
                           cooling_type, TauRT, l_sc, rad_transfer
    use globals,     only : U, PRIM, dz
    use constants
    use hydro_core
    use amr,         only : meshlevel
    use radTransfer, only : F0, a0
    implicit none

    integer, intent(in) :: bIndx
    real, intent(in)    :: dt
    integer             :: i, j, k, ilev
    real                :: T
    real                :: dtau , dl, flux !, dV
    real (kind=8)       :: dh, y0, g0, e, y1
    real (kind=8)       :: fpn
    real (kind=8)       :: col,rec,a,b,c,d
    !      xi - neutral carbon abundance (for non-zero electron density)
    real (kind=8), parameter ::  xi=1.e-4

    call meshlevel(bIndx, ilev )
    dl   = dz(ilev)*l_sc       !  [ cm ]
    !dV   = (dz(ilev)*l_sc)**3  !  [ cm^{-3} ]

    do k=1,ncells_z
      do j=1,ncells_y
        do i=1,ncells_x

          ! Calculate temperature of this cell
          call calcTemp (PRIM(bIndx,:,i,j,k), T)

          col = colf ( real(T,8) )       !# collisional ionization rate
          rec = alpha( real(T,8) )       !# rad. recombination rate

          dh  = real( PRIM(bIndx,  1  ,i,j,k)/ mu0, 8 ) !# H density
          y0  = real( PRIM(bIndx,inH0, i,j,k)/ dh , 8 ) !# neutral H fraction

           !# ionizing flux per nucleus
          if (rad_transfer) then
            dtau = a0 * dz(ilev) * l_sc * dh * y0
            flux = F0*exp(-PRIM(bIndx,TauRT,i,j,k-1)) ! x dA is cancelled below
            !fpn = phi / dh
            !fpn  = flux * (1.0 - exp(-dtau) ) / dV / dh**2 /y0
            fpn = flux * (1.0 - exp(-dtau) ) / dl / dh**2 /y0    !  dA/dV = 1/dl
          else
            fpn = 0.
          end if
          !fpn = 0.

          ! solve for the new neutral fraction using the analytical
          ! solution (see notes)
          a=rec+col
          b=-((2.0+xi)*rec+(1.0+xi)*col+fpn)
          c=(1.0+xi)*rec
          d=sqrt(b**2 - 4.0*a*c)
          g0=(2.0*a*y0+b+d)/(2.0*a*y0+b-d)
          e=exp( -d*dh*real(dt,8) )

          y1=(-b-d*(1.0+g0*e)/(1.0-g0*e))/(2.0*a) !# the new neutral fraction
          y1=min(y1,1.0-xi)
          y1=max(y1,0.0)

          !  update the density of neutrals
          PRIM(bIndx,inH0,i,j,k) = PRIM(bIndx,1,i,j,k)*real(y1)/mu0

          !  if Cool_H is enabled, U(y_H) is updated after cooling)
          if (cooling_type /= COOL_H) then
            U (bIndx,inH0,i,j,k) = PRIM(bIndx,inH0,i,j,k)
          end if

        end do
      end do
    end do

  end subroutine solve_h_rate

  !======================================================================

end module hrate
