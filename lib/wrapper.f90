
module wrapper

use, intrinsic :: iso_c_binding

interface
    subroutine GetH2 (N,x,nn,r0,aPBC,mass,VAuFlag,Hp,dHp)
        INTEGER     :: N            ! number of Au atoms
        REAL(8)   :: x(:,:)       ! atom positions, 3 x N+2
        INTEGER   :: nn(:,:)      ! nearest neighbor atoms, N x 12
        REAL(8)     :: r0(3,12)     ! equilibrium positions of 12 nearest neighbors
        REAL(8) :: aPBC(3)      ! box dimensions
        REAL(8)  :: mass(:)      ! mass of nuclei
        INTEGER :: VAuFlag ! Flag for VAu--Au calculation
        REAL(8)    :: Hp(:)        ! neutral, ion and coupling energies of 3N+2 nuclei (3 X 1)
        REAL(8)    :: dHp(:,:)     ! neutral, ion and coupling forces of 3N+2 nuclei (3 X 3N+2)
    end subroutine
end interface

contains

subroutine get_nn(N, x, dnn, aPBC, nn) bind(c, name="get_nn")

    integer(kind=c_int), intent(in) :: N
    real(kind=c_double), intent(in) :: x(3,N+2)
    real(kind=c_double), intent(in) :: dnn
    real(kind=c_double), intent(in) :: aPBC(3)
    integer(kind=c_int), intent(out) :: nn(N,12)

    call GetNN3(N,x,dnn,aPBC,nn)

end subroutine

SUBROUTINE GetNN3(N,x,dnn,aPBC,nn)
! Calculates 1st nearest neighbor array for a set of atom positions in an fcc lattice
! Revised code based on coordinate transformation from 
! ([1 0 0],[0 1 0],[0 0 1]) --> ([1 0 -1]/sqrt(2)],[1 -4 1]/sqrt(6),[1 1 1]/sqrt(3))
! for use with [1 1 1] surface scattering simulations

IMPLICIT NONE

INTEGER :: N 			! number of atoms
REAL(8), INTENT(IN) :: x(:,:)	! atom positions, 3 x N
REAL(8), INTENT(IN) :: dnn 	! approximate equilibrium distance between 1st nearest neighbors
REAL(8), INTENT(IN) :: aPBC(3)	! supercell dimensions
REAL(8), PARAMETER :: RT2 = 1.41421356237310	! SQRT(2)
REAL(8), PARAMETER :: RT3 = 1.73205080756888	! SQRT(3)
REAL(8), PARAMETER :: RT6 = 2.44948974278318	! SQRT(6)
INTEGER :: nn(:,:) 		! nearest neighbors atoms, N x 12

REAL(8) :: r(3) 		! vector to nearest neighbor
REAL(8) :: rb(3)		! vector to nearest neighbor in Cartesian coordinates
REAL(8) :: TEMP(3)		
REAL(8)	:: U(3,3)	! Coordinate transormation matrix r(111) --> r(001)
INTEGER :: i,j,k,sint,jind
INTEGER :: s(3)

U = RESHAPE ( (/ -1d0/RT2,0d0,1d0/RT2, 1d0/RT6,-2d0/RT6,1d0/RT6, -1d0/RT3,-1d0/RT3,-1d0/RT3 /), (/3,3/) )
DO i = 1,N
  ! print*, x(:,i+2)
  DO j = 1,12
    nn(i,j) = 0
  END DO
  DO j = 1,N
    TEMP = (x(:,j+2) - x(:,i+2))/aPBC
    ! print*, TEMP, i, j
    TEMP = TEMP - floor(TEMP + 0.5)
    r = TEMP*aPBC
    IF (SUM(r*r) < dnn*dnn) THEN
      rb = MATMUL(U,r)
      s = NINT(rb/dnn)
      sint = s(1)*9+s(2)*3+s(3)
      SELECT CASE (sint)
        case (4)
          nn(i,1) = j
        case (2)
          nn(i,2) = j
        case (10)
          nn(i,3) = j
        case (-8)
          nn(i,4) = j
        case (12)
          nn(i,5) = j
        case (6)
          nn(i,6) = j
        case (-4)
          nn(i,7) = j
        case (-2)
          nn(i,8) = j
        case (-10)
          nn(i,9) = j
        case (8)
          nn(i,10) = j
        case (-12)
          nn(i,11) = j
        case (-6)
          nn(i,12) = j
      END SELECT
    END IF
  END DO
END DO

END SUBROUTINE GetNN3

subroutine get_energies_and_forces(N, x, nn, r0, aPBC, mass, VAuFlag, Hp, dHp) bind(c, name="get_energies_and_forces")

    integer(kind=c_int), intent(in)  :: N
    real(kind=c_double), intent(in)  :: x(3,N+2)
    integer(kind=c_int), intent(in)  :: nn(N,12)
    real(kind=c_double), intent(in)  :: r0(3,12)
    real(kind=c_double), intent(in)  :: aPBC(3)
    real(kind=c_double), intent(in)  :: mass(N+2)
    integer(kind=c_int), intent(in)  :: VAuFlag
    real(kind=c_double), intent(out) :: Hp(3)
    real(kind=c_double), intent(out) :: dHp(3,3*(N+2))

    call GetH2(N, x, nn, r0, aPBC, mass, VAuFlag, Hp, dHp)

end subroutine

end module
