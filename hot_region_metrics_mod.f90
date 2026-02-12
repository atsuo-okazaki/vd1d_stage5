module hot_region_metrics_mod
  use kind_params, only : dp, i4b
  implicit none
  private
  public :: init_hot_region_metrics, update_hot_region_metrics_after_commit

  real(dp), save :: Tcrit_hot = 1.0e4_dp
  integer(i4b), save :: iu_hot = -1
  logical, save :: is_init = .false.

contains

  subroutine init_hot_region_metrics(Tcrit_in, filename)
    use mod_global, only : t0, R0
    real(dp), intent(in), optional :: Tcrit_in
    character(len=*), intent(in), optional :: filename

    if (present(Tcrit_in)) Tcrit_hot = Tcrit_in

    if (present(filename)) then
      open(newunit=iu_hot, file=filename, status='replace', action='write')
      write(iu_hot,'("#  t0 =", 1pe14.6, ", R0 =", 1pe14.6, " [cm]")') t0, R0
      write(iu_hot,'(a)') '#     t/t0       r_hot_outer    r_island_max  n(island)  has_inner_hot'
    else
      iu_hot = -1
    end if

    is_init = .true.
  end subroutine init_hot_region_metrics


  pure real(dp) function crossing_r_logT(r0, T0, r1, T1, Tcrit) result(rc)
    ! Return r where T(r)=Tcrit assuming logT varies linearly with r between (r0,T0) and (r1,T1).
    ! If interpolation is ill-posed, fall back to r0 (conservative).
    real(dp), intent(in) :: r0, T0, r1, T1, Tcrit
    real(dp) :: x0, x1, xc, denom
    real(dp), parameter :: tinyT = 1.0e-300_dp

    rc = r0

    if (r1 <= r0) return
    if (Tcrit <= 0.0_dp) return

    ! Guard against non-positive temperatures
    x0 = log(max(T0, tinyT))
    x1 = log(max(T1, tinyT))
    xc = log(Tcrit)

    denom = (x1 - x0)
    if (abs(denom) <= 1.0e-300_dp) then
      ! Nearly flat in logT: cannot locate crossing reliably
      rc = r0
      return
    end if

    rc = r0 + (r1 - r0) * (xc - x0) / denom

    ! Clamp inside segment
    if (rc < r0) rc = r0
    if (rc > r1) rc = r1
  end function crossing_r_logT


  subroutine update_hot_region_metrics_after_commit(t_nd, r_nd, Tmid)
    use units_disk_mod, only : t_dim, r_dim
    real(dp), intent(in) :: t_nd
    real(dp), intent(in) :: r_nd(:)
    real(dp), intent(in) :: Tmid(:)

    integer(i4b) :: nr, i, i_end
    logical :: has_inner_hot
    integer(i4b) :: n_island
    real(dp) :: r_hot_outer, r_island_max
    integer(i4b) :: istart, iend

    if (.not. is_init) return

    nr = size(r_nd)
    if (size(Tmid) /= nr) return

    has_inner_hot = (Tmid(1) >= Tcrit_hot)

    ! ----------------------------
    ! (1) Inner-connected hot zone
    ! ----------------------------
    if (.not. has_inner_hot) then
      r_hot_outer = 0.0_dp
      i_end = 0
    else
      i_end = 1
      do i = 2, nr
        if (Tmid(i) >= Tcrit_hot) then
          i_end = i
        else
          exit
        end if
      end do

      if (i_end >= nr) then
        ! Hot all the way to outer boundary
        r_hot_outer = r_nd(nr)
      else
        ! Edge-interpolated crossing between i_end (hot) and i_end+1 (cold)
        r_hot_outer = crossing_r_logT( r_nd(i_end),   Tmid(i_end), &
                                       r_nd(i_end+1), Tmid(i_end+1), Tcrit_hot )
      end if
    end if

    ! ----------------------------
    ! (2) Detached hot islands
    ! ----------------------------
    n_island     = 0
    r_island_max = 0.0_dp

    i = max(i_end+1, 2)
    do while (i <= nr)

      if (Tmid(i) >= Tcrit_hot) then
        n_island = n_island + 1
        istart = i

        ! walk to the end of this island (last hot index = iend)
        iend = i
        do while (iend+1 <= nr)
           if (Tmid(iend+1) >= Tcrit_hot) then
              iend = iend + 1
           else
              exit
           end if
        end do

        ! compute edge-interpolated outer boundary of this island
        if (iend >= nr) then
          r_island_max = max(r_island_max, r_nd(nr))
        else
          r_island_max = max(r_island_max, crossing_r_logT( r_nd(iend),   Tmid(iend), &
                                                            r_nd(iend+1), Tmid(iend+1), Tcrit_hot ))
        end if

        i = iend + 1
      else
        i = i + 1
      end if

    end do

    write(iu_hot,'(3(1pe14.6, 1x), i6, 5x, i6)') &
         t_nd, r_hot_outer, r_island_max, n_island, merge(1,0,has_inner_hot)

        end subroutine update_hot_region_metrics_after_commit

end module hot_region_metrics_mod
