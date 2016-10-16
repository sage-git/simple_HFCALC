subroutine INTGRL(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
  use integral_functions
  use constants
  use INT
  implicit none
  integer,intent(in) :: IOP, N
  double precision, intent(in) :: R, ZETA1, ZETA2, ZA, ZB
  double precision :: COEF(3, 3), EXPON(3, 3), D1(3), A1(3), D2(3), A2(3)
  double precision :: R2
  double precision :: RAP, RAP2, RBP, RBP2
  double precision :: RAQ, RAQ2, RBQ, RBQ2
  double precision :: RPQ, RPQ2
  integer :: I, J, K, L
  data COEF / 1.0d0, 0.0d0, 0.d0, &
            & 0.678914d0, 0.430129d0, 0.0d0, &
            & 0.444635d0, 0.535328d0, 0.154329d0 /
  data EXPON / 0.27095d0, 0.0d0, 0.d0, &
            & 0.151623d0, 0.851819d0, 0.0d0, &
            & 0.109818d0, 0.405771d0, 02.22766d0/

  R2 = R*R
  do I = 1, N
    A1(I) = EXPON(I, N)*(ZETA1**2)
    D1(I) = COEF(I, N)*((2.0d0*A1(I)/PI)**0.75d0)
    A2(I) = EXPON(I, N)*(ZETA2**2)
    D2(I) = COEF(I, N)*((2.0d0*A2(I)/PI)**0.75d0)
  enddo
  S12 = 0.0d0
  T11 = 0.0d0
  T12 = 0.0d0
  T22 = 0.0d0
  V11A = 0.0d0
  V12A = 0.0d0
  V22A = 0.0d0
  V11B = 0.0d0
  V12B = 0.0d0
  V22B = 0.0d0
  V1111 = 0.0d0
  V2111 = 0.0d0
  V2121 = 0.0d0
  V2211 = 0.0d0
  V2221 = 0.0d0
  V2222 = 0.0d0

   do I = 1, N
     do J = 1, N
       RAP = A2(J)*R/(A1(I) + A2(J))
       RAP2 = RAP**2
       RBP2 = (R - RAP)**2
       S12 = S12 + S(A1(I), A2(J), R2)*D1(I)*D2(J)
       T11 = T11 + T(A1(I), A1(J), 0.d0)*D1(I)*D1(J)
       T12 = T12 + T(A1(I), A2(J), R2)*D1(I)*D2(J)
       T22 = T22 + T(A2(I), A2(J), 0.d0)*D2(I)*D2(J)
       V11A = V11A + V(A1(I), A1(J), 0.0d0, 0.0d0, ZA)*D1(I)*D1(J)
       V12A = V12A + V(A1(I), A2(J), R2, RAP2, ZA)*D1(I)*D2(J)
       V22A = V22A + V(A2(I), A2(J), 0.0d0, R2, ZA)*D2(I)*D2(J)
       V11B = V11B + V(A1(I), A1(J), 0.0d0, R2, ZB)*D1(I)*D1(J)
       V12B = V12B + V(A1(I), A2(J), R2, RBP2, ZB)*D1(I)*D2(J)
       V22B = V22B + V(A2(I), A2(J), 0.0d0, 0.0d0, ZB)*D2(I)*D2(J)
     enddo
   enddo

   do I = 1, N
     do J = 1, N
       do K = 1, N
         do L = 1, N
           RAP = A2(I)*R/(A2(I) + A1(J))
           RBP = R - RAP
           RAQ = A2(K)*R/(A2(K) + A1(L))
           RBQ = R - RAQ
           RPQ = RAP - RAQ
           RAP2 = RAP * RAP
           RBP2 = RBP * RBP
           RAQ2 = RAQ * RAQ
           RBQ2 = RBQ * RBQ
           RPQ2 = RPQ * RPQ
           V1111 = V1111 + TWOE(A1(I), A1(J), A1(K), A1(L), 0.d0, 0.d0, 0.d0)*D1(I)*D1(J)*D1(K)*D1(L)
           V2111 = V2111 + TWOE(A2(I), A1(J), A1(K), A1(L), R2, 0.d0, RAP2)*D2(I)*D1(J)*D1(K)*D1(L)
           V2121 = V2121 + TWOE(A2(I), A1(J), A2(K), A1(L), R2, R2, RPQ2)*D2(I)*D1(J)*D2(K)*D1(L)
           V2211 = V2211 + TWOE(A2(I), A2(J), A1(K), A1(L), 0.d0, 0.d0, R2)*D2(I)*D2(J)*D1(K)*D1(L)
           V2221 = V2221 + TWOE(A2(I), A2(J), A2(K), A1(L), 0.d0, R2, RBQ2)*D2(I)*D2(J)*D2(K)*D1(L)
           V2222 = V2222 + TWOE(A2(I), A2(J), A2(K), A2(L), 0.d0, 0.d0, 0.d0)*D2(I)*D2(J)*D2(K)*D2(L)
         enddo
       enddo
     enddo
   enddo

   if (IOP /= 0)then
     write(*, *)
     write(*, '(A3, 5A11)')" ", "R          ", "ZETA1       ", "ZETA2      ", "S12        ", "T11        "
     write(*, '(5F11.6)') R, ZETA1, ZETA2, S12, T11
     write(*, '(A3, 5A11)')" ", "T12        ", "T22         ", "V11A       ", "V12A       ", "V22A       "
     write(*, '(5F11.6)') T12, T22, V11A, V12A, V22A
     write(*, '(A3, 5A11)')" ", "V11B       ", "V12B        ", "V22B       ", "V1111      ", "V2111      "
     write(*, '(5F11.6)') V11B, V12B, V22B, V1111, V2111
     write(*, '(A3, 4A11)')" ", "V2121      ", "V2211       ", "V2221      ", "V2222      "
     write(*, '(5F11.6)') V2121, V2211, V2221, V2222
     write(*, *)
   endif

 end subroutine INTGRL
