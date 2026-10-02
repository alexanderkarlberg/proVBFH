!-----------------------------------------------------------------------
! Per-point cache of matrix elements for the on-the-fly scale variations
! (cs_scales): within one phase-space point the matrix elements do not
! depend on mu_R, mu_F, so the extra scale points reuse the values of the
! first. Entries are keyed by integer labels (routine, flavours, line, ...)
! and the exact input momenta (and scale, where it enters), and are
! invalidated for the next point by a generation counter. Off (memo_on
! false) the callers evaluate directly, so a run without scale variations
! is unchanged.
!-----------------------------------------------------------------------
module cs_memo
  implicit none
  private
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter, public :: memo_nkey = 10, memo_nmom = 29, memo_nval = 17
  integer, parameter :: nslot = 32768
  logical, public, save :: memo_on = .false.
  integer(8), public, save :: memo_hits = 0, memo_misses = 0
  integer, save :: gen = 1
  integer, save :: sgen(nslot) = 0
  integer, save :: skey(memo_nkey, nslot)
  real(dp), save :: smom(memo_nmom, nslot), sval(memo_nval, nslot)
  public :: memo_new_point, memo_get, memo_put
contains

  subroutine memo_new_point()
    gen = gen + 1
    if (gen == huge(gen)) then
       sgen = 0; gen = 1
    endif
  end subroutine memo_new_point

  ! look up (key, mom(1:nm)); if found, val(1:nv) is returned; if not, slot
  ! is where memo_put stores it
  logical function memo_get(key, mom, nm, val, nv, slot) result(found)
    integer, intent(in) :: key(memo_nkey), nm, nv
    real(dp), intent(in) :: mom(nm)
    real(dp), intent(out) :: val(nv)
    integer, intent(out) :: slot
    integer(8) :: h
    integer :: i
    h = 0
    do i = 1, memo_nkey
       h = h*1000003_8 + key(i)
       h = iand(h, 2147483647_8)
    enddo
    do i = 1, nm, 3
       h = ieor(h*31_8, iand(transfer(mom(i), 1_8), 2147483647_8))
       h = iand(h, 2147483647_8)
    enddo
    do i = 0, nslot - 1
       slot = 1 + int(modulo(h + i, int(nslot, 8)))
       if (sgen(slot) /= gen) then
          found = .false.
          memo_misses = memo_misses + 1
          return
       endif
       if (all(skey(:,slot) == key)) then
          if (all(smom(1:nm,slot) == mom)) then
             val = sval(1:nv,slot)
             found = .true.
             memo_hits = memo_hits + 1
             return
          endif
       endif
    enddo
    stop 'cs_memo: table full'
  end function memo_get

  subroutine memo_put(slot, key, mom, nm, val, nv)
    integer, intent(in) :: slot, key(memo_nkey), nm, nv
    real(dp), intent(in) :: mom(nm), val(nv)
    sgen(slot) = gen
    skey(:,slot) = key
    smom(1:nm,slot) = mom
    sval(1:nv,slot) = val
  end subroutine memo_put
end module cs_memo
