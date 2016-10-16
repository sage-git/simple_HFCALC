module int
  implicit none
  double precision :: S12, T11, T12, T22
  double precision :: V11A, V12A, V22A, V11B, V12B, V22B
  double precision :: V1111, V2111, V2121, V2211, V2221, V2222
end module int

module matrix
  implicit none
  double precision :: S(2, 2), X(2, 2), XT(2, 2)
  double precision :: H(2, 2), F(2, 2), G(2, 2), C(2, 2)
  double precision :: FPRIME(2, 2), CPRIME(2, 2)
  double precision :: P(2, 2), OLDP(2, 2)
  double precision :: TT(2, 2, 2, 2), E(2, 2)

end module matrix
