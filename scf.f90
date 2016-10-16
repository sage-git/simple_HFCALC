subroutine SCF(IOP, N, R, ZETA1, ZETA2, ZA, ZB)
  use MATRIX
  use constants
  implicit none
  integer, intent(in) :: IOP, N
  double precision, intent(in) :: R, ZETA1, ZETA2, ZA, ZB

  integer :: ITER, I, J, K
  logical :: output
  double precision :: EN, DELTA, ENT

  double precision, parameter :: CRIT = 1.0d-4
  integer, parameter :: MAXIT = 25


  output = IOP >= 2
  ITER = 0
  do I = 1, 2
    do J = 1, 2
      P(I, J) = 0.0d0
    enddo
  enddo

  if(output) call MATOUT(P, 2, 2, 2, 2, "P    ")

  do while(.true.)
    ITER = ITER + 1
    if(output) then
      write(*, *)
      write(*, '(A, I2)')"START OF ITERATION NUMBER = ", ITER
      write(*, *)
    endif

    call FORMG
    if(output) call MATOUT(G, 2, 2, 2, 2, "G    ")

    do I = 1, 2
      do J = 1, 2
        F(I, J) = H(I, J) + G(I, J)
      enddo
    enddo

    EN = 0.d0
    do I = 1, 2
      do J = 1, 2
        EN = EN + 0.5d0*P(I, J)*(H(I, J) + F(I, J))
      enddo
    enddo
    if(output) call MATOUT(F, 2, 2, 2, 2, "F   ")
    if(output) then
      write(*, *)
      write(*, '(A, D20.12)')"ELECTRONIC ENERGY = ", EN
      write(*, *)
    endif

    call MULT(F, X, G, 2, 2)
    call MULT(XT, G, FPRIME, 2, 2)
    call DIAG(FPRIME, CPRIME, E)
    call MULT(X, CPRIME, C, 2, 2)

    do I = 1, 2
      do J = 1, 2
        OLDP(I, J) = P(I, J)
        P(I, J) = 0.d0
        do K = 1, 1
          P(I, J) = P(I, J) + 2.d0*C(I, K)*C(J, K)
        enddo
      enddo
    enddo

    if(output) then
      call MATOUT(FPRIME, 2, 2, 2, 2, "F'  ")
      call MATOUT(CPRIME, 2, 2, 2, 2, "C'  ")
      call MATOUT(E, 2, 2, 2, 2, "E    ")
      call MATOUT(C, 2, 2, 2, 2, "C    ")
      call MATOUT(P, 2, 2, 2, 2, "P    ")
    endif

    DELTA = 0.d0
    do I = 1, 2
      do J = 1, 2
        DELTA = DELTA + (P(I, J) - OLDP(I, J))**2
      enddo
    enddo
    DELTA = sqrt(DELTA/4.0d0)
    if(IOP /= 0)then
      write(*, *)
      write(*, '(A, F10.6)') "DELTA(CONVERGENCE OF DENSITY MATRIX) = ", DELTA
      write(*, *)
    endif

    if(DELTA < CRIT) exit
    if(ITER >= MAXIT)then
      write(*, *) "NO CONVERGENCE IN SCF"
      stop
    endif
  enddo ! do while(ITER < MAXIT)

  ENT = EN + ZA*ZB/R

  if (IOP /= 0) then
    write(*, *)
    write(*, *) "CALCULATION CONVERGED"
    write(*, *)
    write(*, '(A20, D20.12)') "ELECTRONIC ENERGY = ", EN
    write(*, '(A20, D20.12)') "TOTAL ENERGY = ", ENT
    write(*, *)
  endif

  if(IOP == 1)then
    call MATOUT(G, 2, 2, 2, 2, "G    ")
    call MATOUT(F, 2, 2, 2, 2, "F    ")
    call MATOUT(E, 2, 2, 2, 2, "E    ")
    call MATOUT(C, 2, 2, 2, 2, "C    ")
    call MATOUT(P, 2, 2, 2, 2, "P    ")
  endif

  call MULT(P, S, OLDP, 2, 2)

  if(IOP /= 0) call MATOUT(OLDP, 2, 2, 2, 2, "PS  ")

end subroutine SCF
