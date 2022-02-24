SUBROUTINE GetH2(N,x,nn,r0,aPBC,mass,VAuFlag,Hp,dHp)

! Forces and Energies

IMPLICIT NONE

INTEGER, INTENT(IN)     :: N            ! number of Au atoms
REAL(8)     :: x(:,:)       ! atom positions, 3 x N+2
!REAL(8), INTENT(IN)     :: x(:,:)       ! atom positions, 3 x N+2
INTEGER, INTENT(IN)     :: nn(:,:)      ! nearest neighbor atoms, N x 12
REAL(8), INTENT(IN)     :: r0(3,12)     ! equilibrium positions of 12 nearest neighbors
REAL(8), INTENT(IN)     :: aPBC(3)      ! box dimensions
REAL(8), INTENT(IN)     :: mass(:)      ! mass of nuclei
INTEGER, INTENT(IN)     :: VAuFlag ! 1 to calculate VAu--Au, 0 to turn off calculation
REAL(8), INTENT(OUT)    :: Hp(:)        ! neutral, ion and coupling energies of 3(N+2) nuclei (3 X 1)
REAL(8), INTENT(OUT)    :: dHp(:,:)     ! neutral, ion and coupling forces of 3(N+2) nuclei (3 X 3(N+2))

REAL(8), PARAMETER :: Ang = 1.0d-10       ! 1 Angstrom in meters
REAL(8), PARAMETER :: eV = 1.60217646d-19    ! 1 eV in Joules
REAL(8), PARAMETER :: kJpermol = 1.660538863d-21 ! 1 kJ/mol in Joules
REAL(8), PARAMETER :: amu = 1.66053886d-27    ! 1 amu in kg

! DEFINE PARAMETERS FOR NEUTRAL SURFACE

REAL(8), PARAMETER    :: A_n=457000.052788d0*kJpermol             ! Au--O: exponential repulsion: A
REAL(8), PARAMETER     :: alpha_n=3.75257871821d0/Ang         ! Au--O: exponential repulsion: alpha
REAL(8), PARAMETER     :: Au_O_cutoff_n=10.0d0*Ang   ! Au--O: cutoff
REAL(8), PARAMETER     :: B_n=30788.8486039d0*kJpermol             ! Au--N: exponential repulsion: B
REAL(8), PARAMETER     :: beta_n=2.9728905911d0/Ang          ! Au--N: exponential repulsion: beta
REAL(8), PARAMETER     :: Au_N_cutoff_n=10.0d0*Ang   ! Au--N: cutoff
REAL(8), PARAMETER     :: F_n=638.5000000000d0*kJpermol             ! N--O: Morse: F
REAL(8), PARAMETER     :: delta_n=2.7430000000d0/Ang         ! N--O: Morse: delta
REAL(8), PARAMETER     :: rNO_e_n=1.1507700000d0*Ang         ! N--O: Morse: cutoff

! DEFINE PARAMETERS FOR ION SURFACE

REAL(8), PARAMETER     :: C_i=1.25581276843d0*Ang             ! Image potential: C
REAL(8), PARAMETER     :: D_i= 347.2225355d0*kJpermol*Ang             ! Image potential: D
REAL(8), PARAMETER     :: z0_i= 1.153606314d0*Ang            ! Image potential: z0
REAL(8), PARAMETER     :: A_i=457000.052788d0*kJpermol             ! Au--O: exponential repulsion: A
REAL(8), PARAMETER     :: alpha_i=3.75257871821d0/Ang         ! Au--O: exponential repulsion: alpha
REAL(8), PARAMETER     :: Au_O_cutoff_i=10.d0*Ang   ! Au--O: cutoff
REAL(8), PARAMETER     :: B_i=23.8597594272d0*kJpermol             ! Au--N: Morse: B
REAL(8), PARAMETER     :: beta_i=1.91014033785d0/Ang          ! Au--N: Morse: beta
REAL(8), PARAMETER     :: rN_sur_e_i=2.38958832878d0*Ang      ! Au--N: Morse: rN_sur_e
REAL(8), PARAMETER     :: Au_N_cutoff_i=10.d0*Ang   ! Au--N: cutoff
REAL(8), PARAMETER     :: F_i=495.9809807d0*kJpermol             ! N--O: Morse: F
REAL(8), PARAMETER     :: delta_i=2.47093477934d0/Ang         ! N--O: Morse: delta
REAL(8), PARAMETER     :: rNO_e_i=1.29289837288d0*Ang         ! N--O: Morse: cutoff
REAL(8), PARAMETER     :: KI=512.06425722d0*kJpermol         ! N--O: Morse: cutoff

! DEFINE PARAMETERS FOR COUPLING FUNCTION

REAL(8), PARAMETER     :: coup_a_N=-70.5259924491d0*kJpermol             ! Au--N: Exponential decay: a
REAL(8), PARAMETER     :: coup_b_N=0.00470023958504d0             ! Au--N: Exponential decay: b
REAL(8), PARAMETER     :: coup_beta_N=1.95982478112d0/Ang          ! Au--N: Exponential decay: beta
REAL(8), PARAMETER     :: Au_N_coupling_cutoff=10.d0*Ang ! Au--N: Exponential decay: cutoff
REAL(8), PARAMETER     :: coup_a_O=-16.7488672932d0*kJpermol             ! Au--O: Exponential decay: a
REAL(8), PARAMETER     :: coup_b_O=0.00617151653727d0            ! Au--O: Exponential decay: b
REAL(8), PARAMETER     :: coup_beta_O=1.35353579356d0/Ang          ! Au--O: Exponential decay: beta
REAL(8), PARAMETER     :: Au_O_coupling_cutoff=10.d0*Ang ! Au--O: Exponential decay: cutoff

! VARIABLES

INTEGER  :: i,j,k        ! counter
REAL(8)  :: rNO      ! bond length of NO
REAL(8)  :: zcom     ! z position of centre of mass of NO
REAL(8)  :: V_Au_Au,Vme  ! Energy of Au surface
REAL(8)  :: V_Au_N_n,V_Au_O_n,V_N_O_n   ! Components of "Neutral" energy
REAL(8)  :: V_Au_N_i,V_Au_O_i,V_N_O_i,V_image_pot   ! Components of "Ion" energy
REAL(8)  :: V_coupling  ! Coupling energy
REAL(8), ALLOCATABLE :: F_temp(:,:)
REAL(8), ALLOCATABLE :: F_NEUTRAL(:,:),F_ION(:,:),F_COUPLING(:,:) !"Neutral", "Ion" and "Coupling" forces (3 X N+2)
REAL(8) :: TIME

ALLOCATE(F_temp(3,N+2))
ALLOCATE(F_NEUTRAL(3,N+2),F_ION(3,N+2),F_COUPLING(3,N+2))

!--- BOND LENGTH AND CENTRE OF MASS OF NO

rNO = SUM((x(:,1)-x(:,2))**2)
rNO = DSQRT(rNO)

zcom = ((x(3,1)*mass(1))+(x(3,2)*mass(2)))/(mass(1)+mass(2))

!------------- CALCULATE ENERGIES ----------------------------

!--- ENERGY OF Au SURFACE

IF (VAuFlag == 1) THEN
CALL GetVNN2(V_Au_Au) 
END IF

!--- ENERGY OF THE "NEUTRAL" SURFACE

CALL ENERGY_Au_N_NEUTRAL(V_Au_N_n)
CALL ENERGY_Au_O_NEUTRAL(V_Au_O_n)
CALL ENERGY_N_O_NEUTRAL(V_N_O_n)

!WRITE(*,*) V_Au_N_n/eV, V_Au_O_n/eV, V_N_O_n/eV
Hp(1) = V_Au_N_n + V_Au_O_n + V_N_O_n + V_Au_Au

!--- ENERGY OF THE "ION" SURFACE


CALL ENERGY_Au_N_ION(V_Au_N_i)
CALL ENERGY_Au_O_ION(V_Au_O_i)
CALL ENERGY_N_O_ION(V_N_O_i)
CALL ENERGY_IMAGE_POTENTIAL(V_image_pot)

Hp(2) = V_Au_N_i + V_Au_O_i + V_N_O_i + V_image_pot + V_Au_Au + KI

!--- COUPLING FUNCTION

CALL ENERGY_COUPLING(V_coupling)

Hp(3) = V_coupling


!--------------- CALCULATE FORCES -----------------------------

!--- FORCES FROM "NEUTRAL" SURFACE:

CALL FORCE_Au_N_NEUTRAL(F_temp)
F_NEUTRAL = F_temp

CALL FORCE_Au_O_NEUTRAL(F_temp)
F_NEUTRAL = F_NEUTRAL+F_temp

CALL FORCE_N_O_NEUTRAL(F_temp)
F_NEUTRAL = F_NEUTRAL+F_temp


!--- FORCES FROM "ION" SURFACE:


CALL FORCE_Au_N_ION(F_temp) ! *****Expensive
F_ION = F_temp



CALL FORCE_Au_O_ION(F_temp) ! *****Expensive
F_ION = F_ION + F_temp




CALL FORCE_N_O_ION(F_temp)
F_ION = F_ION + F_temp




CALL FORCE_IMAGE_POTENTIAL(F_temp)
F_ION = F_ION + F_temp


!--- FORCES DUE TO COUPLING:

CALL FORCE_COUPLING(F_COUPLING) ! *****Expensive


!--- FORCES DUE TO Au--Au INTERACTION:

CALL GetFNN2(F_temp(:,3:N+2))


F_NEUTRAL(:,3:N+2) = F_NEUTRAL(:,3:N+2) - F_temp(:,3:N+2)
F_ION(:,3:N+2) = F_ION(:,3:N+2) - F_temp(:,3:N+2)

!----------------------------------------------------

dHp(1,:) = RESHAPE(F_NEUTRAL,(/3*(N+2)/))
dHp(2,:) = RESHAPE(F_ION,(/3*(N+2)/))
dHp(3,:) = RESHAPE(F_COUPLING,(/3*(N+2)/))

dHp = -dHp


DEALLOCATE(F_temp)
DEALLOCATE(F_NEUTRAL,F_ION,F_COUPLING)


CONTAINS

!---------------------------------------------------------

SUBROUTINE ENERGY_Au_N_NEUTRAL(V_Au_N_n)

IMPLICIT NONE

INTEGER               :: i
REAL(8)               :: delx_sur_N,dely_sur_N,delz_sur_N
REAL(8), INTENT(OUT)  :: V_Au_N_n
REAL(8)               :: dist_N_sur
REAL(8)               :: V,V_Au_N_cut

V_Au_N_n = 0.0d0 

do i = 1, N

  delx_sur_N = x(1,i+2) - x(1,1)
  dely_sur_N = x(2,i+2) - x(2,1)
  delz_sur_N = x(3,i+2) - x(3,1) 

  delx_sur_N = delx_sur_N - (aPBC(1)*DNINT(delx_sur_N/aPBC(1)))
  dely_sur_N = dely_sur_N - (aPBC(2)*DNINT(dely_sur_N/aPBC(2)))

  if (DABS(delx_sur_N).lt.Au_N_cutoff_n) then
    if (DABS(dely_sur_N).lt.Au_N_cutoff_n) then
      if (DABS(delz_sur_N).lt.Au_N_cutoff_n) then

        dist_N_sur = delx_sur_N*delx_sur_N + dely_sur_N*dely_sur_N + delz_sur_N*delz_sur_N
    
        if (dist_N_sur.lt.Au_N_cutoff_n*Au_N_cutoff_n) then
          dist_N_sur = DSQRT(dist_N_sur) 
          V = B_n*(dexp(-beta_n*dist_N_sur))
          V_Au_N_cut = B_n*(dexp(-beta_n*Au_N_cutoff_n))
          V_Au_N_n = V_Au_N_n + V - V_Au_N_cut 
        end if
 
      end if
    end if
  end if 

end do

END SUBROUTINE ENERGY_Au_N_NEUTRAL

!------------------------------------------------------------

SUBROUTINE ENERGY_Au_O_NEUTRAL(V_Au_O_n)

IMPLICIT NONE

INTEGER               :: i
REAL(8)               :: delx_sur_O,dely_sur_O,delz_sur_O
REAL(8), INTENT(OUT)  :: V_Au_O_n
REAL(8)               :: dist_O_sur
REAL(8)               :: V,V_Au_O_cut

V_Au_O_n = 0.0d0 

do i = 1, N 

  delx_sur_O = x(1,i+2) - x(1,2)
  dely_sur_O = x(2,i+2) - x(2,2)
  delz_sur_O = x(3,i+2) - x(3,2)

  delx_sur_O = delx_sur_O - (aPBC(1)*DNINT(delx_sur_O/aPBC(1)))
  dely_sur_O = dely_sur_O - (aPBC(2)*DNINT(dely_sur_O/aPBC(2)))

  if (DABS(delx_sur_O).lt.Au_O_cutoff_n) then
    if (DABS(dely_sur_O).lt.Au_O_cutoff_n) then
      if (DABS(delz_sur_O).lt.Au_O_cutoff_n) then

        dist_O_sur = delx_sur_O*delx_sur_O + dely_sur_O*dely_sur_O + delz_sur_O*delz_sur_O

        if (dist_O_sur.lt.Au_O_cutoff_n*Au_O_cutoff_n) then
          dist_O_sur = DSQRT(dist_O_sur)
          V = A_n*(dexp(-alpha_n*dist_O_sur))
          V_Au_O_cut = A_n*(dexp(-alpha_n*Au_O_cutoff_n))
          V_Au_O_n = V_Au_O_n + V - V_Au_O_cut 
        end if
 
      end if
    end if
  end if 

end do 

END SUBROUTINE ENERGY_Au_O_NEUTRAL

!---------------------------------------------------------

SUBROUTINE ENERGY_Au_N_ION(V_Au_N_i)

IMPLICIT NONE

INTEGER               :: i
REAL(8)               :: delx_sur_N,dely_sur_N,delz_sur_N
REAL(8), INTENT(OUT)  :: V_Au_N_i
REAL(8)               :: dist_N_sur
REAL(8)               :: cos_theta_sq,V,V_Au_N_cut

V_Au_N_i = 0.0d0 

do i = 1, N

  delx_sur_N = x(1,i+2)-x(1,1)
  dely_sur_N = x(2,i+2)-x(2,1)
  delz_sur_N = x(3,i+2)-x(3,1) 

  delx_sur_N = delx_sur_N - (aPBC(1)*DNINT(delx_sur_N/aPBC(1)))
  dely_sur_N = dely_sur_N - (aPBC(2)*DNINT(dely_sur_N/aPBC(2)))

  if (DABS(delx_sur_N).lt.Au_N_cutoff_i) then
    if (DABS(dely_sur_N).lt.Au_N_cutoff_i) then
      if (DABS(delz_sur_N).lt.Au_N_cutoff_i) then

        dist_N_sur = delx_sur_N*delx_sur_N + dely_sur_N*dely_sur_N + delz_sur_N*delz_sur_N

        if (dist_N_sur.lt.Au_N_cutoff_i*Au_N_cutoff_i) then

          dist_N_sur = DSQRT(dist_N_sur) 
          cos_theta_sq = (x(3,2) - x(3,1))*(x(3,2) - x(3,1))/(rNO*rNO)
          V = B_i*(dexp(-2.0d0*beta_i*(dist_N_sur-rN_sur_e_i))-&
          (2.0d0*cos_theta_sq*dexp(-beta_i*(dist_N_sur-rN_sur_e_i))))
          V_Au_N_cut = B_i*(dexp(-2.0d0*beta_i*(Au_N_cutoff_i-rN_sur_e_i))-&
          (2.0d0*cos_theta_sq*dexp(-beta_i*(Au_N_cutoff_i-rN_sur_e_i))))
          V_Au_N_i = V_Au_N_i + V - V_Au_N_cut 

        end if
      end if
    end if
 
  end if 

end do 

END SUBROUTINE ENERGY_Au_N_ION

!-------------------------------------------------------------

SUBROUTINE ENERGY_Au_O_ION(V_Au_O_i)

IMPLICIT NONE

INTEGER               :: i
REAL(8)               :: delx_sur_O,dely_sur_O,delz_sur_O
REAL(8), INTENT(OUT)  :: V_Au_O_i
REAL(8)               :: dist_O_sur
REAL(8)               :: V,V_Au_O_cut

V_Au_O_i = 0.0d0 

do i = 1, N 

  delx_sur_O = x(1,i+2)-x(1,2)
  dely_sur_O = x(2,i+2)-x(2,2)
  delz_sur_O = x(3,i+2)-x(3,2)

  delx_sur_O = delx_sur_O - (aPBC(1)*DNINT(delx_sur_O/aPBC(1)))
  dely_sur_O = dely_sur_O - (aPBC(2)*DNINT(dely_sur_O/aPBC(2)))

  if (DABS(delx_sur_O).lt.Au_O_cutoff_i) then
    if (DABS(dely_sur_O).lt.Au_O_cutoff_i) then
      if (DABS(delz_sur_O).lt.Au_O_cutoff_i) then

           dist_O_sur = delx_sur_O*delx_sur_O + dely_sur_O*dely_sur_O + delz_sur_O*delz_sur_O

        if (dist_O_sur.lt.Au_O_cutoff_i*Au_O_cutoff_i) then

          dist_O_sur = DSQRT(dist_O_sur)
          V = A_i*(dexp(-alpha_i*dist_O_sur))
          V_Au_O_cut = A_i*(dexp(-alpha_i*Au_O_cutoff_i))
          V_Au_O_i = V_Au_O_i + V - V_Au_O_cut 

        end if
      end if
    end if
 
  end if 

end do 

END SUBROUTINE ENERGY_Au_O_ION

!-----------------------------------------------------------------

SUBROUTINE ENERGY_IMAGE_POTENTIAL(V_image_pot)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: V_image_pot

V_image_pot = -D_i/DSQRT(C_i*C_i + (zcom - z0_i)*(zcom - z0_i))

END SUBROUTINE ENERGY_IMAGE_POTENTIAL

!-------------------------------------------------------------------

SUBROUTINE ENERGY_N_O_NEUTRAL(V_N_O_n)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: V_N_O_n

V_N_O_n = F_n*(1-dexp(-delta_n*(rNO - rNO_e_n)))*(1-dexp(-delta_n*(rNO - rNO_e_n)))

END SUBROUTINE ENERGY_N_O_NEUTRAL

!---------------------------------------------------------------------

SUBROUTINE ENERGY_N_O_ION(V_N_O_i)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: V_N_O_i

V_N_O_i = F_i*(1-dexp(-delta_i*(rNO - rNO_e_i)))*(1-dexp(-delta_i*(rNO - rNO_e_i)))

END SUBROUTINE ENERGY_N_O_ION

!-----------------------------------------------------------------------

SUBROUTINE ENERGY_COUPLING(V_coupling)

IMPLICIT NONE

INTEGER               :: i
REAL(8)               :: delx_sur_N,dely_sur_N,delz_sur_N
REAL(8)               :: delx_sur_O,dely_sur_O,delz_sur_O
REAL(8), INTENT(OUT)  :: V_coupling
REAL(8)               :: dist_N_sur,dist_O_sur
REAL(8)               :: V,V_cut

V_coupling = 0.0d0

do i = 1, N
  delx_sur_N = x(1,i+2) - x(1,1)
  dely_sur_N = x(2,i+2) - x(2,1)
  delz_sur_N = x(3,i+2) - x(3,1)

  delx_sur_O = x(1,i+2) - x(1,2)
  dely_sur_O = x(2,i+2) - x(2,2)
  delz_sur_O = x(3,i+2) - x(3,2)

  delx_sur_N = delx_sur_N - aPBC(1)*DNINT(delx_sur_N/aPBC(1))
  dely_sur_N = dely_sur_N - aPBC(2)*DNINT(dely_sur_N/aPBC(2))

  delx_sur_O = delx_sur_O - aPBC(1)*DNINT(delx_sur_O/aPBC(1))
  dely_sur_O = dely_sur_O - aPBC(2)*DNINT(dely_sur_O/aPBC(2))

  if (DABS(delx_sur_N).lt.Au_N_coupling_cutoff) then
    if (DABS(dely_sur_N).lt.Au_N_coupling_cutoff) then
      if (DABS(delz_sur_N).lt.Au_N_coupling_cutoff) then

        dist_N_sur = delx_sur_N*delx_sur_N + dely_sur_N*dely_sur_N + delz_sur_N*delz_sur_N
        dist_N_sur = DSQRT(dist_N_sur)
    
        if (dist_N_sur.lt.Au_N_coupling_cutoff) then
          V = coup_a_N/(1.0d0 + coup_b_N*dexp(coup_beta_N*dist_N_sur))
          V_cut = coup_a_N/(1.0d0 + coup_b_N*dexp(coup_beta_N*Au_N_coupling_cutoff))
          V_coupling = V_coupling + V - V_cut
        end if
 
      end if
    end if
  end if 

  if (DABS(delx_sur_O).lt.Au_O_coupling_cutoff) then
    if (DABS(dely_sur_O).lt.Au_O_coupling_cutoff) then
      if (DABS(delz_sur_O).lt.Au_O_coupling_cutoff) then

        dist_O_sur = delx_sur_O*delx_sur_O + dely_sur_O*dely_sur_O + delz_sur_O*delz_sur_O
        dist_O_sur = DSQRT(dist_O_sur)
    
        if (dist_O_sur.lt.Au_O_coupling_cutoff) then 

          V = coup_a_O/(1.0d0 + coup_b_O*dexp(coup_beta_O*dist_O_sur))
          V_cut = coup_a_O/(1.0d0 + coup_b_O*dexp(coup_beta_O*Au_O_coupling_cutoff))
          V_coupling = V_coupling + V - V_cut
       
        end if
 
      end if
    end if
  end if 

end do

END SUBROUTINE ENERGY_COUPLING

!--------------------------------------------------------------------------------

SUBROUTINE FORCE_Au_N_NEUTRAL(F_Au_N_n)

IMPLICIT NONE

INTEGER               :: i
REAL(8)               :: delx_sur_N,dely_sur_N,delz_sur_N
REAL(8), INTENT(OUT)  :: F_Au_N_n(:,:)
REAL(8)               :: dist_N_sur
REAL(8)               :: term,dV_xN,dV_yN,dV_zN,dV_xi,dV_yi,dV_zi

F_Au_N_n = 0.0d0
         
do i = 1, N 

  delx_sur_N = x(1,i+2) - x(1,1)
  dely_sur_N = x(2,i+2) - x(2,1)
  delz_sur_N = x(3,i+2) - x(3,1)

  delx_sur_N = delx_sur_N - (aPBC(1)*DNINT(delx_sur_N/aPBC(1)))
  dely_sur_N = dely_sur_N - (aPBC(2)*DNINT(dely_sur_N/aPBC(2)))

  if (DABS(delx_sur_N).lt.Au_N_cutoff_n) then
    if (DABS(dely_sur_N).lt.Au_N_cutoff_n) then
      if (DABS(delz_sur_N).lt.Au_N_cutoff_n) then

        dist_N_sur = delx_sur_N*delx_sur_N + dely_sur_N*dely_sur_N + delz_sur_N*delz_sur_N

        if (dist_N_sur.lt.Au_N_cutoff_n*Au_N_cutoff_n) then

          dist_N_sur = DSQRT(dist_N_sur) 
          term = -B_n*dexp(-beta_n*dist_N_sur)*beta_n/dist_N_sur

          dV_xi = term*delx_sur_N
          dV_yi = term*dely_sur_N
          dV_zi = term*delz_sur_N

          dV_xN = -dV_xi
          dV_yN = -dV_yi
          dV_zN = -dV_zi

          F_Au_N_n(1,1) = F_Au_N_n(1,1) - dV_xN
          F_Au_N_n(2,1) = F_Au_N_n(2,1) - dV_yN
          F_Au_N_n(3,1) = F_Au_N_n(3,1) - dV_zN

          F_Au_N_n(1,i+2) = F_Au_N_n(1,i+2) - dV_xi
          F_Au_N_n(2,i+2) = F_Au_N_n(2,i+2) - dV_yi
          F_Au_N_n(3,i+2) = F_Au_N_n(3,i+2) - dV_zi


          end if
   
        end if
      end if 
    end if  

end do

END SUBROUTINE FORCE_Au_N_NEUTRAL 

!--------------------------------------------------------------------------

SUBROUTINE FORCE_Au_O_NEUTRAL(F_Au_O_n)

IMPLICIT NONE

INTEGER               :: i
REAL(8)               :: delx_sur_O,dely_sur_O,delz_sur_O
REAL(8), INTENT(OUT)  :: F_Au_O_n(:,:)
REAL(8)               :: dist_O_sur
REAL(8)               :: term,dV_xO,dV_yO,dV_zO,dV_xi,dV_yi,dV_zi

F_Au_O_n = 0.0d0
         
do i = 1, N 

  delx_sur_O = x(1,i+2) - x(1,2)
  dely_sur_O = x(2,i+2) - x(2,2)
  delz_sur_O = x(3,i+2) - x(3,2)

  delx_sur_O = delx_sur_O - (aPBC(1)*DNINT(delx_sur_O/aPBC(1)))
  dely_sur_O = dely_sur_O - (aPBC(2)*DNINT(dely_sur_O/aPBC(2)))

  if (DABS(delx_sur_O).lt.Au_O_cutoff_n) then
    if (DABS(dely_sur_O).lt.Au_O_cutoff_n) then
      if (DABS(delz_sur_O).lt.Au_O_cutoff_n) then

        dist_O_sur = delx_sur_O*delx_sur_O + dely_sur_O*dely_sur_O + delz_sur_O*delz_sur_O

        if (dist_O_sur.lt.Au_O_cutoff_n*Au_O_cutoff_n) then

          dist_O_sur = DSQRT(dist_O_sur)

          term = -A_n*dexp(-alpha_n*dist_O_sur)*alpha_n/dist_O_sur 

          dV_xi = term*delx_sur_O
          dV_yi = term*dely_sur_O
          dV_zi = term*delz_sur_O

          dV_xO = -dV_xi
          dV_yO = -dV_yi
          dV_zO = -dV_zi

          F_Au_O_n(1,2) = F_Au_O_n(1,2) - dV_xO
          F_Au_O_n(2,2) = F_Au_O_n(2,2) - dV_yO
          F_Au_O_n(3,2) = F_Au_O_n(3,2) - dV_zO

          F_Au_O_n(1,i+2) = F_Au_O_n(1,i+2) - dV_xi
          F_Au_O_n(2,i+2) = F_Au_O_n(2,i+2) - dV_yi
          F_Au_O_n(3,i+2) = F_Au_O_n(3,i+2) - dV_zi

        end if
   
      end if
    end if 
  end if  

end do

END SUBROUTINE FORCE_Au_O_NEUTRAL 

!-----------------------------------------------------------------------------------

SUBROUTINE FORCE_Au_N_ION(F_Au_N_i)

IMPLICIT NONE

INTEGER               :: i
REAL(8)               :: delx_sur_N,dely_sur_N,delz_sur_N
REAL(8), INTENT(OUT)  :: F_Au_N_i(:,:)
REAL(8)               :: dist_N_sur
REAL(8)               :: r_diff,term1,term2,term3,term4 
REAL(8)               :: dV_xN_1,dV_xN_4,dV_xN_5,dV_xN,dV_yN_1,dV_yN_4,dV_yN_5,dV_yN,dV_zN_1,dV_zN_2,dV_zN_4,dV_zN_5,dV_zN
REAL(8)               :: dV_xO_1,dV_xO,dV_yO_1,dV_yO_2,dV_yO,dV_zO_1,dV_zO_2,dV_zO,dV_xi_3,dV_xi_4,dV_xi,dV_yi_3,dV_yi_4,dV_yi
REAL(8)               :: dV_zi_3,dV_zi_4,dV_zi
REAL(8)               :: dVc_xN_1,dVc_xN,dVc_yN_1,dVc_yN,dVc_zN_1,dVc_zN_2,dVc_zN
REAL(8)               :: dVc_xO_1,dVc_xO,dVc_yO_1,dVc_yO,dVc_zO_1,dVc_zO_2,dVc_zO,dVc_xi,dVc_yi,dVc_zi
         
F_Au_N_i = 0.0d0
 
do i = 1, N 

  delx_sur_N = x(1,i+2) - x(1,1)
  dely_sur_N = x(2,i+2) - x(2,1)
  delz_sur_N = x(3,i+2) - x(3,1)

  delx_sur_N = delx_sur_N - (aPBC(1)*DNINT(delx_sur_N/aPBC(1)))
  dely_sur_N = dely_sur_N - (aPBC(2)*DNINT(dely_sur_N/aPBC(2)))

  if (DABS(delx_sur_N).lt.Au_N_cutoff_i) then
    if (DABS(dely_sur_N).lt.Au_N_cutoff_i) then
      if (DABS(delz_sur_N).lt.Au_N_cutoff_i) then

        dist_N_sur = delx_sur_N*delx_sur_N + dely_sur_N*dely_sur_N + delz_sur_N*delz_sur_N

        if (dist_N_sur.lt.Au_N_cutoff_i*Au_N_cutoff_i) then 

          dist_N_sur = DSQRT(dist_N_sur)

          r_diff = dist_N_sur - rN_sur_e_i 

          term1 = x(3,2) - x(3,1)
          term2 = dexp(-beta_i*r_diff)
          term3 = dexp(-2.0d0*beta_i*r_diff)
          term4 = rNO**4.0d0

          dV_xN_1 = 4.0d0*term2*(x(1,1)-x(1,2))*term1*term1
          dV_xN_1 = dV_xN_1/term4
          dV_xN_4 = 2.0d0*term3*delx_sur_N*beta_i
          dV_xN_4 = dV_xN_4/dist_N_sur
          dV_xN_5 = 2.0d0*term2*delx_sur_N*term1*term1*beta_i
          dV_xN_5 = dV_xN_5/(dist_N_sur*rNO*rNO)
          dV_xN = B_i*(dV_xN_1+dV_xN_4-dV_xN_5) 

          dV_yN_1 = 4.0d0*term2*(x(2,1)-x(2,2))*term1*term1
          dV_yN_1 = dV_yN_1/term4
          dV_yN_4 = 2.0d0*term3*dely_sur_N*beta_i
          dV_yN_4 = dV_yN_4/dist_N_sur
          dV_yN_5 = 2.0d0*term2*dely_sur_N*term1*term1*beta_i
          dV_yN_5 = dV_yN_5/(dist_N_sur*rNO*rNO)
          dV_yN = B_i*(dV_yN_1+dV_yN_4-dV_yN_5)

          dV_zN_1 = 4.0d0*term2*(x(3,1)-x(3,2))*term1*term1
          dV_zN_1 = dV_zN_1/term4
          dV_zN_2 = 4.0d0*term2*term1
          dV_zN_2 = dV_zN_2/(rNO*rNO)
          dV_zN_4 = 2.0d0*term3*delz_sur_N*beta_i
          dV_zN_4 = dV_zN_4/dist_N_sur
          dV_zN_5 = 2.0d0*term2*delz_sur_N*term1*term1*beta_i
          dV_zN_5 = dV_zN_5/(dist_N_sur*rNO*rNO)
          dV_zN = B_i*(dV_zN_1+dV_zN_2+dV_zN_4-dV_zN_5)

          dV_xO_1 = 4.0d0*term2*(x(1,1)-x(1,2))*term1*term1
          dV_xO_1 = dV_xO_1/term4
          dV_xO = B_i*(-dV_xO_1) 

          dV_yO_1 = 4.0d0*term2*(x(2,1)-x(2,2))*term1*term1
          dV_yO_1 = dV_yO_1/term4
          dV_yO = B_i*(-dV_yO_1)

          dV_zO_1 = 4.0d0*term2*(x(3,1)-x(3,2))*term1*term1
          dV_zO_1 = dV_zO_1/term4
          dV_zO_2 = 4.0d0*term2*term1
          dV_zO_2 = dV_zO_2/(rNO*rNO)
          dV_zO = B_i*(-dV_zO_1-dV_zO_2) 

          dV_xi_3 = 2.0d0*term3*delx_sur_N*beta_i
          dV_xi_3 = dV_xi_3/dist_N_sur
          dV_xi_4 = 2.0d0*term2*delx_sur_N*term1*term1*beta_i
          dV_xi_4 = dV_xi_4/(dist_N_sur*rNO*rNO)
          dV_xi = B_i*(-dV_xi_3+dV_xi_4) 

          dV_yi_3 = 2.0d0*term3*dely_sur_N*beta_i
          dV_yi_3 = dV_yi_3/dist_N_sur
          dV_yi_4 = 2.0d0*term2*dely_sur_N*term1*term1*beta_i
          dV_yi_4 = dV_yi_4/(dist_N_sur*rNO*rNO)
          dV_yi = B_i*(-dV_yi_3+dV_yi_4)

          dV_zi_3 = 2.0d0*term3*delz_sur_N*beta_i
          dV_zi_3 = dV_zi_3/dist_N_sur
          dV_zi_4 = 2.0d0*term2*delz_sur_N*term1*term1*beta_i
          dV_zi_4 = dV_zi_4/(dist_N_sur*rNO*rNO)
          dV_zi = B_i*(-dV_zi_3+dV_zi_4)

!------   Calculating derivatives of cutoff function 

          r_diff = Au_N_cutoff_i - rN_sur_e_i

          term2 = dexp(-beta_i*r_diff)
          term3 = dexp(-2.0d0*beta_i*r_diff)

          dVc_xN_1 = 4.0d0*term2*(x(1,1)-x(1,2))*term1*term1
          dVc_xN_1 = dVc_xN_1/term4
          dVc_xN = B_i*(dVc_xN_1) 

          dVc_yN_1 = 4.0d0*term2*(x(2,1)-x(2,2))*term1*term1
          dVc_yN_1 = dVc_yN_1/term4
          dVc_yN = B_i*(dVc_yN_1)

          dVc_zN_1 = 4.0d0*term2*(x(3,1)-x(3,2))*term1*term1
          dVc_zN_1 = dVc_zN_1/term4
          dVc_zN_2 = 4.0d0*term2*term1
          dVc_zN_2 = dVc_zN_2/(rNO*rNO)
          dVc_zN = B_i*(dVc_zN_1+dVc_zN_2)

          dVc_xO_1 = 4.0d0*term2*(x(1,1)-x(1,2))*term1*term1
          dVc_xO_1 = dVc_xO_1/term4
          dVc_xO = B_i*(-dVc_xO_1) 

          dVc_yO_1 = 4.0d0*term2*(x(2,1)-x(2,2))*term1*term1
          dVc_yO_1 = dVc_yO_1/term4
          dVc_yO = B_i*(-dVc_yO_1)

          dVc_zO_1 = 4.0d0*term2*(x(3,1)-x(3,2))*term1*term1
          dVc_zO_1 = dVc_zO_1/term4
          dVc_zO_2 = 4.0d0*term2*term1
          dVc_zO_2 = dVc_zO_2/(rNO*rNO)
          dVc_zO = B_i*(-dVc_zO_1-dVc_zO_2)

          dVc_xi = 0.0d0
          dVc_yi = 0.0d0
          dVc_zi = 0.0d0

          F_Au_N_i(1,1) = F_Au_N_i(1,1) - dV_xN + dVc_xN
          F_Au_N_i(2,1) = F_Au_N_i(2,1) - dV_yN + dVc_yN
          F_Au_N_i(3,1) = F_Au_N_i(3,1) - dV_zN + dVc_zN

          F_Au_N_i(1,2) = F_Au_N_i(1,2) - dV_xO + dVc_xO
          F_Au_N_i(2,2) = F_Au_N_i(2,2) - dV_yO + dVc_yO
          F_Au_N_i(3,2) = F_Au_N_i(3,2) - dV_zO + dVc_zO

          F_Au_N_i(1,i+2) = F_Au_N_i(1,i+2) - dV_xi + dVc_xi
          F_Au_N_i(2,i+2) = F_Au_N_i(2,i+2) - dV_yi + dVc_yi
          F_Au_N_i(3,i+2) = F_Au_N_i(3,i+2) - dV_zi + dVc_zi

          end if
   
        end if
      end if 
    end if  

end do

END SUBROUTINE FORCE_Au_N_ION

!------------------------------------------------------------------

SUBROUTINE FORCE_Au_O_ION(F_Au_O_i)

IMPLICIT NONE

INTEGER               :: i
REAL(8)               :: delx_sur_O,dely_sur_O,delz_sur_O
REAL(8), INTENT(OUT)  :: F_Au_O_i(:,:)
REAL(8)               :: dist_O_sur
REAL(8)               :: term,dV_xO,dV_yO,dV_zO,dV_xi,dV_yi,dV_zi

F_Au_O_i = 0.0d0

do i = 1, N 

  delx_sur_O = x(1,i+2)-x(1,2)
  dely_sur_O = x(2,i+2)-x(2,2)
  delz_sur_O = x(3,i+2)-x(3,2)

  delx_sur_O = delx_sur_O - (aPBC(1)*DNINT(delx_sur_O/aPBC(1)))
  dely_sur_O = dely_sur_O - (aPBC(2)*DNINT(dely_sur_O/aPBC(2)))

  if (DABS(delx_sur_O).lt.Au_O_cutoff_i) then
    if (DABS(dely_sur_O).lt.Au_O_cutoff_i) then
      if (DABS(delz_sur_O).lt.Au_O_cutoff_i) then

        dist_O_sur = delx_sur_O*delx_sur_O + dely_sur_O*dely_sur_O + delz_sur_O*delz_sur_O

        if (dist_O_sur.lt.Au_O_cutoff_i*Au_O_cutoff_i) then

          dist_O_sur = DSQRT(dist_O_sur) 
          term = -A_i*dexp(-alpha_i*dist_O_sur)*alpha_i/dist_O_sur
          dV_xi = term*delx_sur_O
          dV_yi = term*dely_sur_O
          dV_zi = term*delz_sur_O

          dV_xO = -dV_xi
          dV_yO = -dV_yi
          dV_zO = -dV_zi

          F_Au_O_i(1,2) = F_Au_O_i(1,2) - dV_xO
          F_Au_O_i(2,2) = F_Au_O_i(2,2) - dV_yO
          F_Au_O_i(3,2) = F_Au_O_i(3,2) - dV_zO

          F_Au_O_i(1,i+2) = F_Au_O_i(1,i+2) - dV_xi
          F_Au_O_i(2,i+2) = F_Au_O_i(2,i+2) - dV_yi
          F_Au_O_i(3,i+2) = F_Au_O_i(3,i+2) - dV_zi

        end if
      end if
    end if 

  end if  

end do

END SUBROUTINE FORCE_Au_O_ION 

!-------------------------------------------------------------------
         
SUBROUTINE FORCE_IMAGE_POTENTIAL(F_image_pot) 

IMPLICIT NONE 

REAL(8), INTENT(OUT)  :: F_image_pot(:,:)
REAL(8)               :: term

F_image_pot = 0.0d0

term = (D_i*(zcom - z0_i))/(C_i*C_i+(zcom - z0_i)*(zcom - z0_i))**1.5d0
F_image_pot(3,1) = -term * mass(1)/(mass(1) + mass(2))
F_image_pot(3,2) = -term * mass(2)/(mass(1) + mass(2)) 

END SUBROUTINE FORCE_IMAGE_POTENTIAL
 
!-------------------------------------------------------------------

SUBROUTINE FORCE_N_O_NEUTRAL(F_N_O_n)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: F_N_O_n(:,:)
REAL(8)               :: r_diff,term

F_N_O_n = 0.0d0

r_diff = rNO - rNO_e_n

term = 2.0d0*F_n*delta_n*(dexp(-2.0d0*delta_n*(rNO - rNO_e_n))-dexp(-delta_n*(rNO - rNO_e_n)))
F_N_O_n(1,1) = -term*(x(1,2)-x(1,1))/rNO
F_N_O_n(1,2) = -F_N_O_n(1,1)
F_N_O_n(2,1) = -term*(x(2,2)-x(2,1))/rNO
F_N_O_n(2,2) = -F_N_O_n(2,1)
F_N_O_n(3,1) = -term*(x(3,2)-x(3,1))/rNO
F_N_O_n(3,2) = -F_N_O_n(3,1) 

END SUBROUTINE FORCE_N_O_NEUTRAL

!--------------------------------------------------------------------

SUBROUTINE FORCE_N_O_ION(F_N_O_i)

IMPLICIT NONE

REAL(8), INTENT(OUT)  :: F_N_O_i(:,:)
REAL(8)               :: r_diff,term

F_N_O_i = 0.0d0

r_diff = rNO - rNO_e_i

term = 2.0d0*F_i*delta_i*(dexp(-2.0d0*delta_i*(rNO - rNO_e_i))-dexp(-delta_i*(rNO - rNO_e_i)))



F_N_O_i(1,1) = -term*(x(1,2)-x(1,1))/rNO
F_N_O_i(1,2) = -F_N_O_i(1,1)
F_N_O_i(2,1) = -term*(x(2,2)-x(2,1))/rNO
F_N_O_i(2,2) = -F_N_O_i(2,1)
F_N_O_i(3,1) = -term*(x(3,2)-x(3,1))/rNO
F_N_O_i(3,2) = -F_N_O_i(3,1) 

END SUBROUTINE FORCE_N_O_ION

!---------------------------------------------------------------------

SUBROUTINE FORCE_COUPLING(F_COUPLING)

IMPLICIT NONE

INTEGER               :: i
REAL(8)               :: delx_sur_N,dely_sur_N,delz_sur_N
REAL(8)               :: delx_sur_O,dely_sur_O,delz_sur_O
REAL(8), INTENT(OUT)  :: F_COUPLING(:,:)
REAL(8)               :: dist_N_sur,dist_O_sur
REAL(8)               :: r_diff,term1,term2,term3
REAL(8)               :: dV_xN,dV_yN,dV_zN,dV_xO,dV_yO,dV_zO,dV_xi,dV_yi,dV_zi

F_COUPLING = 0.0d0

do i = 1, N 

  delx_sur_N = x(1,i+2) - x(1,1)
  dely_sur_N = x(2,i+2) - x(2,1)
  delz_sur_N = x(3,i+2) - x(3,1)

  delx_sur_N = delx_sur_N - aPBC(1)*floor(delx_sur_N/aPBC(1)+.5)
  dely_sur_N = dely_sur_N - aPBC(2)*floor(dely_sur_N/aPBC(2)+.5)

!  delx_sur_N = delx_sur_N - aPBC(1)*DNINT(delx_sur_N/aPBC(1))
!  dely_sur_N = dely_sur_N - aPBC(2)*DNINT(dely_sur_N/aPBC(2))

  if (DABS(delx_sur_N).lt.Au_N_coupling_cutoff) then
    if (DABS(dely_sur_N).lt.Au_N_coupling_cutoff) then
      if (DABS(delz_sur_N).lt.Au_N_coupling_cutoff) then

        dist_N_sur = delx_sur_N*delx_sur_N + dely_sur_N*dely_sur_N + delz_sur_N*delz_sur_N
    
        if (dist_N_sur.lt.Au_N_coupling_cutoff*Au_N_coupling_cutoff) then 

        dist_N_sur = DSQRT(dist_N_sur)

          term1 = dexp(coup_beta_N*dist_N_sur)
          term2 = (1.0d0 + coup_b_N*term1)*(1.0d0 + coup_b_N*term1)*dist_N_sur
          term3 = coup_a_N*coup_b_N*coup_beta_N/term2

          dV_xi = -term1*delx_sur_N*term3
          dV_yi = -term1*dely_sur_N*term3
          dV_zi = -term1*delz_sur_N*term3

          dV_xN = term1*delx_sur_N*term3
          dV_yN = term1*dely_sur_N*term3
          dV_zN = term1*delz_sur_N*term3

          dV_xO = 0.0d0
          dV_yO = 0.0d0
          dV_zO = 0.0d0

          F_COUPLING(1,1) = F_COUPLING(1,1) - dV_xN
          F_COUPLING(2,1) = F_COUPLING(2,1) - dV_yN
          F_COUPLING(3,1) = F_COUPLING(3,1) - dV_zN

          F_COUPLING(1,2) = F_COUPLING(1,2) - dV_xO
          F_COUPLING(2,2) = F_COUPLING(2,2) - dV_yO
          F_COUPLING(3,2) = F_COUPLING(3,2) - dV_zO

          F_COUPLING(1,i+2) = F_COUPLING(1,i+2) - dV_xi
          F_COUPLING(2,i+2) = F_COUPLING(2,i+2) - dV_yi
          F_COUPLING(3,i+2) = F_COUPLING(3,i+2) - dV_zi

        end if
       end if
      end if
 
    end if 

  delx_sur_O = x(1,i+2) - x(1,2)
  dely_sur_O = x(2,i+2) - x(2,2)
  delz_sur_O = x(3,i+2) - x(3,2)

  delx_sur_O = delx_sur_O - aPBC(1)*DNINT(delx_sur_O/aPBC(1))
  dely_sur_O = dely_sur_O - aPBC(2)*DNINT(dely_sur_O/aPBC(2))

  if (DABS(delx_sur_O).lt.Au_O_coupling_cutoff) then
    if (DABS(dely_sur_O).lt.Au_O_coupling_cutoff) then
      if (DABS(delz_sur_O).lt.Au_O_coupling_cutoff) then

        dist_O_sur = delx_sur_O*delx_sur_O + dely_sur_O*dely_sur_O + delz_sur_O*delz_sur_O
    
        if (dist_O_sur.lt.Au_O_coupling_cutoff*Au_O_coupling_cutoff) then 

          dist_O_sur = DSQRT(dist_O_sur)

          term1 = dexp(coup_beta_O*dist_O_sur)
          term2 = (1.0d0 + coup_b_O*term1)*(1.0d0 + coup_b_O*term1)*dist_O_sur
          term3 = coup_a_O*coup_b_O*coup_beta_O/term2

          dV_xi = -term1*delx_sur_O*term3
          dV_yi = -term1*dely_sur_O*term3
          dV_zi = -term1*delz_sur_O*term3

          dV_xN = 0.0d0
          dV_yN = 0.0d0
          dV_zN = 0.0d0

          dV_xO = term1*delx_sur_O*term3
          dV_yO = term1*dely_sur_O*term3
          dV_zO = term1*delz_sur_O*term3

          F_COUPLING(1,1) = F_COUPLING(1,1) - dV_xN
          F_COUPLING(2,1) = F_COUPLING(2,1) - dV_yN
          F_COUPLING(3,1) = F_COUPLING(3,1) - dV_zN

          F_COUPLING(1,2) = F_COUPLING(1,2) - dV_xO
          F_COUPLING(2,2) = F_COUPLING(2,2) - dV_yO
          F_COUPLING(3,2) = F_COUPLING(3,2) - dV_zO

          F_COUPLING(1,i+2) = F_COUPLING(1,i+2) - dV_xi
          F_COUPLING(2,i+2) = F_COUPLING(2,i+2) - dV_yi
          F_COUPLING(3,i+2) = F_COUPLING(3,i+2) - dV_zi

        end if
      end if
    end if
  end if 

end do

END SUBROUTINE FORCE_COUPLING

!------------------------------------------------------------------------------------
     
SUBROUTINE GetVNN2(Vme)

! Calculates potential energy for generalized force model
! Revised code based on coordinate transformation from 
! ([1 0 0],[0 1 0],[0 0 1]) --> ([1 0 -1]/sqrt(2)],[1 -4 1]/sqrt(6),[1 1 1]/sqrt(3))
! for use with [1 1 1] surface scattering simulations

IMPLICIT NONE

REAL(8), INTENT(OUT)	:: Vme		! Potential energy

REAL(8), PARAMETER :: RT2 = 1.41421356237310d0	! DSQRT(2)
REAL(8), PARAMETER :: RT3 = 1.73205080756888d0	! DSQRT(3)
REAL(8), PARAMETER :: RT6 = 2.44948974278318d0	! DSQRT(6)
REAL(8)			:: alpha(3) = (/-4.94d0, 17.15d0, 19.4d0/)     ! generalized force parameters
REAL(8)			:: r(3)		! Vector to nearest neighbors
REAL(8) 		:: TEMP(3)
REAL(8)			:: FmatOld(3,3)
REAL(8)			:: Fmat1(3,3),Fmat2(3,3),Fmat3(3,3),Fmat4(3,3),Fmat5(3,3),Fmat6(3,3)
					! Force matrices for each nearest neighbor position
REAL(8)			:: U(3,3),TU(3,3)	! Coordinate transormation matrix r(111) --> r(001)
INTEGER 		:: i,j,k
INTEGER			:: init = 0

SAVE init,Fmat1,Fmat2,Fmat3,Fmat4,Fmat5,Fmat6

IF (init .EQ. 0) THEN
U = RESHAPE ( (/ -1d0/RT2,0d0,1d0/RT2, 1d0/RT6,-2d0/RT6,1d0/RT6, -1d0/RT3,-1d0/RT3,-1d0/RT3 /), (/3,3/) )
TU = TRANSPOSE(U)
FmatOld = RESHAPE ( (/ alpha(1),0d0,0d0,0d0,alpha(2),alpha(3),0d0,alpha(3),alpha(2) /), (/3,3/) )
Fmat1 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(1),0d0,0d0,0d0,alpha(2),-alpha(3),0d0,-alpha(3),alpha(2) /), (/3,3/) )
Fmat2 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),0d0,alpha(3),0d0,alpha(1),0d0,alpha(3),0d0,alpha(2) /), (/3,3/) )
Fmat3 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),0d0,-alpha(3),0d0,alpha(1),0d0,-alpha(3),0d0,alpha(2) /), (/3,3/) )
Fmat4 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),alpha(3),0d0,alpha(3),alpha(2),0d0,0d0,0d0,alpha(1) /), (/3,3/) )
Fmat5 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),-alpha(3),0d0,-alpha(3),alpha(2),0d0,0d0,0d0,alpha(1) /), (/3,3/) )
Fmat6 = MATMUL(TU,MATMUL(FmatOld,U))
init = 1
END IF


Vme = 0
DO i = 1,N
   DO j = 1,12
     IF (nn(i,j) /= 0) THEN
      TEMP = (x(:,nn(i,j)+2)-x(:,i+2))/aPBC
      TEMP = TEMP - floor(TEMP + .5d0)
      r = TEMP*aPBC-r0(:,j)
      SELECT CASE (j)
         CASE (1)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat1,r))
         CASE (7)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat1,r))
         CASE (2)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat2,r))
	 CASE (8)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat2,r))
	 CASE (3)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat3,r))
	 CASE (9)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat3,r))
	 CASE (4)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat4,r))
	 CASE (10)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat4,r))
	 CASE (5)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat5,r))
	 CASE (11)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat5,r))
	 CASE (6)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat6,r))
	 CASE (12)
	 Vme = Vme + DOT_PRODUCT(r,MATMUL(Fmat6,r))
      END SELECT
      END IF
   END DO
END DO
Vme = Vme/4.0d0
END SUBROUTINE GetVNN2

!-------------------------------------------------------------------

SUBROUTINE GetFNN2(F)

! Calculates forces for generalized force model
! Revised code based on coordinate transformation from 
! ([1 0 0],[0 1 0],[0 0 1]) --> ([1 0 -1]/sqrt(2)],[1 -4 1]/sqrt(6),[1 1 1]/sqrt(3))
! for use with [1 1 1] surface scattering simulations
! Definitions for nearest neighbors, forces, etc. can be found in G.H. Begbie, Proc. Roy. Soc. Lon.,
! (1947) 180-208.

IMPLICIT NONE

REAL(8), INTENT(OUT)	:: F(:,:)	! forces on atoms

REAL(8), PARAMETER :: RT2 = 1.41421356237310d0	! DSQRT(2)
REAL(8), PARAMETER :: RT3 = 1.73205080756888d0	! DSQRT(3)
REAL(8), PARAMETER :: RT6 = 2.44948974278318d0	! DSQRT(6)
REAL(8)			:: alpha(3) = (/-4.94d0, 17.15d0, 19.4d0/)    ! generalized force parameters
REAL(8)			:: r(3)		! vector to nearest neighbors
REAL(8) 		:: TEMP(3)
REAL(8)			:: FmatOld(3,3)
REAL(8)			:: Fmat1(3,3),Fmat2(3,3),Fmat3(3,3),Fmat4(3,3),Fmat5(3,3),Fmat6(3,3)
					! Force matrices for each nearest neighbor position
REAL(8)			:: U(3,3),TU(3,3)	! Coordinate transormation matrix r(111) --> r(001)
INTEGER 		:: i,j,k
INTEGER			:: init = 0

SAVE init,Fmat1,Fmat2,Fmat3,Fmat4,Fmat5,Fmat6

IF (init .EQ. 0) THEN
U = RESHAPE ( (/ -1d0/RT2,0d0,1d0/RT2, 1d0/RT6,-2d0/RT6,1d0/RT6, -1d0/RT3,-1d0/RT3,-1d0/RT3 /), (/3,3/) )
TU = TRANSPOSE(U)
FmatOld = RESHAPE ( (/ alpha(1),0d0,0d0,0d0,alpha(2),alpha(3),0d0,alpha(3),alpha(2) /), (/3,3/) )
Fmat1 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(1),0d0,0d0,0d0,alpha(2),-alpha(3),0d0,-alpha(3),alpha(2) /), (/3,3/) )
Fmat2 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),0d0,alpha(3),0d0,alpha(1),0d0,alpha(3),0d0,alpha(2) /), (/3,3/) )
Fmat3 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),0d0,-alpha(3),0d0,alpha(1),0d0,-alpha(3),0d0,alpha(2) /), (/3,3/) )
Fmat4 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),alpha(3),0d0,alpha(3),alpha(2),0d0,0d0,0d0,alpha(1) /), (/3,3/) )
Fmat5 = MATMUL(TU,MATMUL(FmatOld,U))
FmatOld = RESHAPE ( (/ alpha(2),-alpha(3),0d0,-alpha(3),alpha(2),0d0,0d0,0d0,alpha(1) /), (/3,3/) )
Fmat6 = MATMUL(TU,MATMUL(FmatOld,U))
init = 1
END IF


DO i = 1,N
   F(:,i) = 0
   DO j = 1,12
     IF (nn(i,j) /= 0) THEN
      TEMP = (x(:,nn(i,j)+2)-x(:,i+2))/aPBC
      TEMP = TEMP - floor(TEMP + .5d0)
      r = TEMP*aPBC-r0(:,j)
      SELECT CASE (j)
         CASE (1)
	 F(:,i) = F(:,i) - MATMUL(Fmat1,r)
         CASE (7)
	 F(:,i) = F(:,i) - MATMUL(Fmat1,r)
         CASE (2)
	 F(:,i) = F(:,i) - MATMUL(Fmat2,r)
	 CASE (8)
	 F(:,i) = F(:,i) - MATMUL(Fmat2,r)
	 CASE (3)
	 F(:,i) = F(:,i) - MATMUL(Fmat3,r)
	 CASE (9)
	 F(:,i) = F(:,i) - MATMUL(Fmat3,r)
	 CASE (4)
	 F(:,i) = F(:,i) - MATMUL(Fmat4,r)
	 CASE (10)
	 F(:,i) = F(:,i) - MATMUL(Fmat4,r)
	 CASE (5)
	 F(:,i) = F(:,i) - MATMUL(Fmat5,r)
	 CASE (11)
	 F(:,i) = F(:,i) - MATMUL(Fmat5,r)
	 CASE (6)
	 F(:,i) = F(:,i) - MATMUL(Fmat6,r)
	 CASE (12)
	 F(:,i) = F(:,i) - MATMUL(Fmat6,r)
      END SELECT
      END IF
   END DO
END DO

END SUBROUTINE GetFNN2

!-----------------------------------------------------------------------

END SUBROUTINE GetH2