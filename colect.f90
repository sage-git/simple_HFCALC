SUBROUTINE COLECT(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
  USE MATRIX
  USE INT
  IMPLICIT NONE
  INTEGER, INTENT(in) :: IOP, N
  DOUBLE PRECISION, INTENT(in) :: R, ZETA1, ZETA2, ZA, ZB
  INTEGER :: I, J, K, L

  H(1, 1) = T11 + V11A + V11B
  H(1, 2) = T12 + V12A + V12B
  H(2, 1) = H(1, 2)
  H(2, 2) = T22 + V22A + V22B

  S(1, 1) = 1.0d0
  S(1, 2) = S12
  S(2, 1) = S12
  S(2, 2) = 1.0d0

  X(1, 1) = 1.0d0/SQRT(2.0d0*(1.0d0 + S12))
  X(2, 1) = X(1, 1)
  X(1, 2) = 1.0d0/SQRT(2.0d0*(1.0d0 - S12))
  X(2, 2) = - X(1, 2)

  XT(1, 1) = X(1, 1)
  XT(1, 2) = X(2, 1)
  XT(2, 1) = X(1, 2)
  XT(2, 2) = X(2, 2)

  TT(1, 1, 1, 1) = V1111
  TT(2, 1, 1, 1) = V2111
  TT(1, 2, 1, 1) = V2111
  TT(1, 1, 2, 1) = V2111
  TT(1, 1, 1, 2) = V2111
  TT(2, 1, 2, 1) = V2121
  TT(1, 2, 2, 1) = V2121
  TT(2, 1, 1, 2) = V2121
  TT(1, 2, 1, 2) = V2121
  TT(2, 2, 1, 1) = V2211
  TT(1, 1, 2, 2) = V2211
  TT(2, 2, 2, 1) = V2221
  TT(2, 2, 1, 2) = V2221
  TT(2, 1, 2, 2) = V2221
  TT(1, 2, 2, 2) = V2221
  TT(2, 2, 2, 2) = V2222

  IF(IOP /= 0)THEN
     CALL MATOUT(S, 2, 2, 2, 2, "S   ")
     CALL MATOUT(X, 2, 2, 2, 2, "X   ")
     CALL MATOUT(H, 2, 2, 2, 2, "H   ")
     WRITE(*, *) ""
     DO I = 1, 2
        DO J = 1, 2
           DO K = 1, 2
              DO L = 1, 2
                 WRITE(*, '(4(I2, A3), F16.10)')I, " ", J, " ", K, " ", L, " ", TT(I, J, K, L)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDIF

END SUBROUTINE COLECT
