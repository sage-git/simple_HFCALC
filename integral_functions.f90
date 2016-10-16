
module integral_functions
  use constants
  implicit none
contains
  double precision function F0(ARG)
    implicit none
    double precision, intent(in) :: ARG
    if ( ARG < 1.0d-6)then
      F0 = 1.0d0 - ARG/3.0d0
    else
#if defined USE_MY_ERF
      F0 = sqrt(PI/ARG)*my_ERF(sqrt(ARG))/2.0d0
#else
      F0 = sqrt(PI/ARG)*erf(sqrt(ARG))/2.0d0
#endif
    endif
  end function F0

  double precision function my_ERF(ARG)
    implicit none
    double precision, intent(in) :: ARG
    double precision :: T, TN, POLY
    integer :: I
    double precision, parameter :: A(5) = (/0.254829592d0, -0.284496736d0, 1.421413741d0, -1.453152027d0, 1.061405429d0/)
    double precision, parameter :: P = 0.3275911d0
    T = 1.0d0/(1.0d0 + P*ARG)
    TN = T
    POLY = A(1)*TN
    do I = 2, 5
      TN = TN*T
      POLY = POLY + A(I)*TN
    enddo
    my_ERF = 1.d0 - POLY*exp( - ARG**2)
  end function my_ERF

  double precision function S(A, B, RAB2)
    implicit none
    double precision, intent(in) :: A, B, RAB2
    S = (PI/(A + B))**1.5d0*exp(-A*B*RAB2/(A + B))
  end function S

  double precision function T(A, B, RAB2)
    implicit none
    double precision, intent(in) :: A, B, RAB2
    double precision :: AB, C
    AB = A + B
    C = A*B/AB
    T = C * (3.0d0 - 2.0d0*C*RAB2) * (PI/AB)**1.5d0*exp( - C*RAB2)
  end function T

  double precision function V(A, B, RAB2, RCP2, ZC)
    implicit none
    double precision, intent(in) :: A, B, RAB2, RCP2, ZC
    double precision :: AB
    AB = A + B
    V = - 2.0d0*PI/AB*F0(AB*RCP2)*exp(-A*B*RAB2/AB)*ZC
  end function V

  double precision function TWOE(A, B, C, D, RAB2, RCD2, RPQ2)
    implicit none
    double precision, intent(in) :: A, B, C, D, RAB2, RCD2, RPQ2
    double precision :: AB, CD, ABCD
    AB = A + B
    CD = C + D
    ABCD = AB + CD
    TWOE = 2d0*(PI**2.5d0)/(AB*CD*sqrt(ABCD))*F0(AB*CD*RPQ2/ABCD)*exp(-A*B*RAB2/AB-C*D*RCD2/CD)
  end function TWOE

end module integral_functions
