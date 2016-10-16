subroutine FORMG
  use MATRIX
  implicit none
  integer :: I, J, K, L

  do I = 1, 2
    do J = 1, 2
      G(I, J) = 0.0d0
      do K = 1, 2
        do L = 1, 2
          G(I, J) = G(I, J) + P(K, L)*(TT(I, J, K, L) - 0.5d0*TT(I, L, K, J))
        enddo
      enddo
    enddo
  enddo
end subroutine FORMG

subroutine DIAG(F, C, E)
  use constants
  implicit none
  double precision,intent(in) :: F(2, 2)
  double precision, intent(out) :: C(2, 2), E(2, 2)
  double precision :: temp
  double precision :: sint, cost, sin2t
  double precision :: THETA
  THETA = PI/4.0d0
  if(abs(F(1, 1) - F(2, 2)) > 1.0d-20)then
    THETA = 0.5d0*atan(2.0d0*F(1, 2)/(F(1, 1) - F(2, 2)))
  endif
  sint = sin(THETA)
  cost = cos(THETA)
  sin2t = sin(2.d0*THETA)

  C(1, 1) = cost
  C(1, 2) = sint
  C(2, 1) = sint
  C(2, 2) = - cost

  E(1, 1) = F(1, 1)*cost**2 + F(2, 2)*sint**2 + F(1, 2)*sin2t
  E(2, 2) = F(2, 2)*cost**2 + F(1, 1)*sint**2 - F(1, 2)*sin2t
  E(1, 2) = 0.d0
  E(2, 1) = 0.d0

  if (E(2, 2) > E(1, 1)) then
    return
  endif

  temp = E(2, 2)
  E(2, 2) = E(1, 1)
  E(1, 1) = temp

  temp = C(1, 2)
  C(1, 2) = C(1, 1)
  C(1, 1) = temp

  temp = C(2, 2)
  C(2, 2) = C(2, 1)
  C(2, 1) = temp

end subroutine DIAG

subroutine MULT(A, B, C, IM, M)
  implicit none
  integer :: IM, M
  double precision :: A(IM, IM), B(IM, IM), C(IM, IM)
  integer :: I, J, K
  do I = 1, M
    do J = 1, M
      C(I, J) = 0.d0
      do K = 1, M
        C(I, J) = C(I, J) + A(I, K)*B(K, J)
      enddo
    enddo
  enddo

end subroutine MULT

subroutine MATOUT(A, IM, IN, M, N, LABEL)
    implicit none
    integer,intent(in) :: IM, IN, M, N
    character(len=*),intent(in) :: LABEL
    double precision,intent(in) :: A(IM, IN)
    integer :: IHIGH, LOW, I, J
    IHIGH = 0
    write(*, *)"THE ", LABEL, "ARRAY"
    do while (N - IHIGH > 0)
      LOW = IHIGH + 1
      IHIGH = IHIGH + 5
      IHIGH = min(IHIGH, N)
      write(*, '(A15)', advance='no')" "
      do I = LOW, IHIGH
        write(*, '(A10, I3, A6)', advance='no') " ", I, " "
      enddo
      write(*,*)
      do I = 1, N
        write(*, '(I10, A5, 5D19.10)')I, " ", (A(I, J), J=LOW, IHIGH)
      enddo
    enddo
end subroutine MATOUT
