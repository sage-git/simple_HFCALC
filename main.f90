PROGRAM simple_RHF
    implicit none
    integer :: IOP, N
    double precision :: R, ZETA1, ZETA2, ZA, ZB

    IOP = 2
    N = 3
    R = 1.4632d0
    ZETA1 = 2.0925d0
    ZETA2 = 1.24d0
    ZA = 2.0d0
    ZB = 1.0d0
    call HFCALC(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
    stop

end program simple_RHF

subroutine HFCALC(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
  implicit none
  integer,intent(in) :: IOP, N
  double precision,intent(in) :: R, ZETA1, ZETA2, ZA, ZB
  if(IOP /= 0)then
    write(*, '(A, I0, A, F5.2, A, F5.2)') "STO-", N, "G FOR ATOMC NUMBERS ", ZA, " AND ", ZB
  endif
  call INTGRL(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
  call COLECT(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
  call SCF(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
end subroutine HFCALC
