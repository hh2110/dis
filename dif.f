!-----------------------------------------------------------------------
! SUBROUTINE DEFINING TEMPORAL EVOLUTION OF THE SYSTEM (CONCENTRATION)
!-----------------------------------------------------------------------
      SUBROUTINE DIF(STEPID)
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES   
      INTEGER :: STEPID
      INTEGER :: IX,IY,IZ             !! POSITION INDICES
      INTEGER :: I,J                  !! DUMMY INDICES FOR IJK DIRECTIONS
      REAL :: CONC_0(NX,NY,NZ)        !! FOR THE IMPLICIT TRAPEZOIDAL...
                                      !! ...INTEGRATOR-predictor
      REAL :: CONC_1(NX,NY,NZ)        !! FOR THE IMPLICIT TRAPEZOIDAL
                                      !! ...INTEGRATOR-corrector
      REAL :: DIFFC                   !! CONC0-CONC1-REQUIRED 4 TOL TEST
      REAL :: CGC(NX,NY,NZ)           !! variable to store CALCGC result
      REAL :: GCONC(NX,NY,NZ)         !!THE EQUIV. OF DFDETA FOR ETA PROPAGATION
      INTEGER :: KK                   !!COUNTER FOR TOLERANCE TEST FOR
      INTEGER :: STATUS               !!FOR EXIT CALL during CONC EVO
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      CGC = CALCGC(CONC)
      CONC_0 = CONC + (CGC*DT)
      KK=0;
      DO 
        KK=KK+1
        CONC_1 = CONC + 
     &   ((DT/2.0)*(CGC+CALCGC(CONC_0)))
        DIFFC=MAXVAL(ABS(CONC_0-CONC_1))
!        PRINT *, 'diffc is ', DIFFC
        IF (DIFFC < DELTA) THEN
          EXIT
        ELSE IF(ISNAN(MAXVAL(ABS(CONC_0-CONC_1)))) THEN
          PRINT *, 'conc is NaN'
          STATUS = 0 
          CALL EXIT(STATUS)  
        END IF 
        CONC_0 = CONC_1
      END DO
      CONC = CONC_1
!
      CONTAINS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! FUNCTION TO CALCULATE DF/DC
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      FUNCTION CALCGC(CONCC)
! DEFINE VARIABLES NEEDED IN FUNCTION
      REAL :: CONCC(NX,NY,NZ)
      REAL :: CALCGC(NX,NY,NZ)
      REAL :: CONCT(NX+4,NY+4,NZ+4) !!TEMPORARY VECTOR TO ENFORCE PBCs
      REAL :: SUMMT(NX+4,NY+4,NZ+4)
      REAL :: SUMM(NX,NY,NZ)
      REAL :: GDIFFN(NX,NY,NZ)
      INTEGER :: IX, IY, IZ
      REAL :: TS0, TS1, TS2, TS3, TS4, TS5
      REAL :: TE0, TE1, TE2, TE3, TE4, TE5
! CALCILATING DIFFERENCE BETWEEN GS's (1-SISF CONC & 2-NORMAL CONC)
      CALL CPU_TIME(TS0)
      CALL CPU_TIME(TS1)
      GDIFF = CALCGD(ETA(1,1,:,:,:),ETA(1,2,:,:,:),ETA(1,3,:,:,:))
      CALL CPU_TIME(TE1)
! NORMAILISING CONSTANT TO for GDIFFN
      DELTAGSISF = 50.0/(2.0598*12.5*1.5*1000.0)
! NORMALISED CHANGE IN GS
      CALL CPU_TIME(TS2)
      GDIFFN = GDIFF/DELTAGSISF
      CALL CPU_TIME(TE2)
!      PRINT *, MAXVAL(GDIFF)
! EQUIL CONC - DEPENDS ON DEFORMATION IN X,Y,Z
      CALL CPU_TIME(TS3)
      C0=( (PHI)*((GDIFFN*CS)+((1-GDIFFN)*CGP)) ) + ( (1-PHI)*CG )
      CALL CPU_TIME(TE3)
! MATRIX (N+4)^3 - ENSURES PERIODIC BCs
      CONCT=0.0
      CONCT(3:NX+2, 3:NY+2, 3:NZ+2) = CONCC
      CONCT(1:2   , 3:NY+2, 3:NZ+2) = CONCC(NX-1:NX,    1:NY,    1:NZ) 
      CONCT(3:NX+2, 1:2   , 3:NZ+2) = CONCC(   1:NX, NY-1:NY,    1:NZ) 
      CONCT(3:NX+2, 3:NY+2,    1:2) = CONCC(   1:NX,    1:NY, NZ-1:NZ)
      CONCT(NX+3:NX+4, 3:NY+2   ,    3:NZ+2) = CONCC(1:2 , 1:NY, 1:NZ) 
      CONCT(3:NX+2   , NY+3:NY+4,    3:NZ+2) = CONCC(1:NX, 1:2 , 1:NZ) 
      CONCT(3:NX+2   , 3:NY+2   , NZ+3:NZ+4) = CONCC(1:NX, 1:NY,  1:2)
! CALCULATING ONE, TWO, THREE - TERMS IN CONC EVOLUTION EQN - SEE
! diffusion.pdf
      CALL CPU_TIME(TS4)
      ONE = (-1.0/(2.0*A2))*(1/COSH((A1-CONCC)/A2))*PHI*GDIFF
      TWO = 2.0*BETA1*(CONCC-C0)*EXP(-((CONCC-C0)**2))
      CALL CPU_TIME(TE4)
      CALL CPU_TIME(TS5)
      THREE = 2.0*EG*((1.0/(DDX**2))* 
     &      (CONCT(4:NX+3,3:NY+2,3:NZ+2)+CONCT(2:NX+1,3:NY+2,3:NZ+2)+
     &       CONCT(3:NX+2,4:NY+3,3:NZ+2)+CONCT(3:NX+2,2:NY+1,3:NZ+2)+
     &       CONCT(3:NX+2,3:NY+2,4:NZ+3)+CONCT(3:NX+2,3:NY+2,2:NZ+1)-
     &       (6.0*CONCC)))
      CALL CPU_TIME(TE5)
! SUMMING 1,2,3 UP
      SUMM=ONE+TWO-THREE
! APPLYING PBC TO SUMM FOR DERIVATIVE CALC 
      SUMMT=0.0
      SUMMT(3:NX+2, 3:NY+2, 3:NZ+2) = SUMM
      SUMMT(1:2   , 3:NY+2, 3:NZ+2) = SUMM(NX-1:NX,    1:NY,    1:NZ) 
      SUMMT(3:NX+2, 1:2   , 3:NZ+2) = SUMM(   1:NX, NY-1:NY,    1:NZ) 
      SUMMT(3:NX+2, 3:NY+2,    1:2) = SUMM(   1:NX,    1:NY, NZ-1:NZ)
      SUMMT(NX+3:NX+4, 3:NY+2   ,    3:NZ+2) = SUMM(1:2 , 1:NY, 1:NZ) 
      SUMMT(3:NX+2   , NY+3:NY+4,    3:NZ+2) = SUMM(1:NX, 1:2 , 1:NZ) 
      SUMMT(3:NX+2   , 3:NY+2   , NZ+3:NZ+4) = SUMM(1:NX, 1:NY,  1:2)
! APPLYING DERIVATIVE TO SUMMT - CALCULATING dc/dt
      CALCGC = MG*
     &      ((1.0/(DDX**2))* 
     &      (SUMMT(4:NX+3,3:NY+2,3:NZ+2)+SUMMT(2:NX+1,3:NY+2,3:NZ+2)+
     &       SUMMT(3:NX+2,4:NY+3,3:NZ+2)+SUMMT(3:NX+2,2:NY+1,3:NZ+2)+
     &       SUMMT(3:NX+2,3:NY+2,4:NZ+3)+SUMMT(3:NX+2,3:NY+2,2:NZ+1)-
     &       (6.0*SUMM)))
! CHECKING IF dc/dt IS REAL
      IF(ISNAN(MAXVAL(ABS(CALCGC)))) THEN
        PRINT *, 'conc_calcgc is NaN - please check dif.f file L105'
        STATUS = 0 
        CALL EXIT(1)  
      END IF
      CALL CPU_TIME(TE0)
! FINISHED
      WRITE(1011,*) TE0-TS0,TE1-TS1,TE2-TS2,TE3-TS3,TE4-TS4,TE5-TS5
      RETURN
      END FUNCTION 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! FUNCTION TO CALCULATE DIFFERENCE BETWEEN TWO GSs
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      FUNCTION CALCGD(EM1,EM2,EM3)
      REAL :: CALCGD(NX,NY,NZ)
      REAL :: EM1(NX,NY,NZ), EM2(NX,NY,NZ), EM3(NX,NY,NZ)
      REAL :: E1, E2, E3
      INTEGER :: IX,IY,IZ
!$OMP PARALLEL DO PRIVATE(IY,IZ,E1,E2,E3)
      DO IX=1,NX !!
      DO IY=1,NY
      DO IZ=1,NZ
        E1=EM1(IX,IY,IZ)
        E2=EM2(IX,IY,IZ)
        E3=EM3(IX,IY,IZ)
        CALCGD(IX,IY,IZ)= PHI(IX,IY,IZ)*(-(GC(3,1)+
     &        GC(3,2)*(COS(PI*(E1-E2))+COS(PI*(E2-E3))+
     &          COS(PI*(E3-E1))) +
     &        GC(3,3)*(COS(PI*(2*E1-E2-E3))+
     &          COS(PI*(2*E2-E3-E1))+
     &          COS(PI*(2*E3-E1-E2))) +
     &        GC(3,4)*(COS(TWOPI*(E1-E2))+COS(TWOPI*(E2-E3))+
     &          COS(TWOPI*(E3-E1))) +
     &        GC(3,5)*(COS(PI*(3*E1-E2-2*E3))+COS(PI*(3*E1-2*E2-E3))+
     &          COS(PI*(3*E2-E3-2*E1))+COS(PI*(3*E2-2*E3-E1))+
     &          COS(PI*(3*E3-E1-2*E2))+COS(PI*(3*E3-2*E1-E2)))+
     &        GC(3,6)*(COS(3*PI*(E1-E2))+COS(3*PI*(E2-E3))+
     &          COS(3*PI*(E3-E1))) +
     &        GC(3,7)*(COS(TWOPI*(2*E1-E2-E3))+
     &          COS(TWOPI*(2*E2-E3-E1))+
     &          COS(TWOPI*(2*E3-E1-E2))) +
     &        GC(3,8)*(COS(4*PI*(E1-E2))+COS(4*PI*(E2-E3))+
     &          COS(4*PI*(E3-E1))) +
     &        GC(3,10)*(SIN(PI*(E1-E2))+SIN(PI*(E2-E3))+
     &          SIN(PI*(E3-E1))) +
     &        GC(3,11)*(SIN(TWOPI*(E1-E2))+SIN(TWOPI*(E2-E3))+
     &          SIN(TWOPI*(E3-E1))) +
     &        GC(3,12)*(SIN(PI*(2*E1-3*E2+E3))+SIN(PI*(3*E1-2*E2-E3))+
     &          SIN(PI*(-2*E1-E2+3*E3))+SIN(PI*(E1+2*E2-3*E3))+
     &          SIN(PI*(-3*E1+E2+2*E3))+SIN(PI*(-E1+3*E2-2*E3)))+
     &        GC(3,13)*(SIN(3*PI*(E1-E2))+SIN(3*PI*(E2-E3))+
     &          SIN(3*PI*(E3-E1)))+
     &        GC(3,14)*(SIN(4*PI*(E1-E2))+SIN(4*PI*(E2-E3))+
     &          SIN(4*PI*(E3-E1))))
     &        +
     &        (GC(2,1)+
     &        GC(2,2)*(COS(PI*(E1-E2))+COS(PI*(E2-E3))+
     &          COS(PI*(E3-E1))) +
     &        GC(2,3)*(COS(PI*(2*E1-E2-E3))+
     &          COS(PI*(2*E2-E3-E1))+
     &          COS(PI*(2*E3-E1-E2))) +
     &        GC(2,4)*(COS(TWOPI*(E1-E2))+COS(TWOPI*(E2-E3))+
     &          COS(TWOPI*(E3-E1))) +
     &        GC(2,5)*(COS(PI*(3*E1-E2-2*E3))+COS(PI*(3*E1-2*E2-E3))+
     &          COS(PI*(3*E2-E3-2*E1))+COS(PI*(3*E2-2*E3-E1))+
     &          COS(PI*(3*E3-E1-2*E2))+COS(PI*(3*E3-2*E1-E2)))+
     &        GC(2,6)*(COS(3*PI*(E1-E2))+COS(3*PI*(E2-E3))+
     &          COS(3*PI*(E3-E1))) +
     &        GC(2,7)*(COS(TWOPI*(2*E1-E2-E3))+
     &          COS(TWOPI*(2*E2-E3-E1))+
     &          COS(TWOPI*(2*E3-E1-E2))) +
     &        GC(2,8)*(COS(4*PI*(E1-E2))+COS(4*PI*(E2-E3))+
     &          COS(4*PI*(E3-E1))) +
     &        GC(2,10)*(SIN(PI*(E1-E2))+SIN(PI*(E2-E3))+
     &          SIN(PI*(E3-E1))) +
     &        GC(2,11)*(SIN(TWOPI*(E1-E2))+SIN(TWOPI*(E2-E3))+
     &          SIN(TWOPI*(E3-E1))) +
     &        GC(2,12)*(SIN(PI*(2*E1-3*E2+E3))+SIN(PI*(3*E1-2*E2-E3))+
     &          SIN(PI*(-2*E1-E2+3*E3))+SIN(PI*(E1+2*E2-3*E3))+
     &          SIN(PI*(-3*E1+E2+2*E3))+SIN(PI*(-E1+3*E2-2*E3)))+
     &        GC(2,13)*(SIN(3*PI*(E1-E2))+SIN(3*PI*(E2-E3))+
     &          SIN(3*PI*(E3-E1)))+
     &        GC(2,14)*(SIN(4*PI*(E1-E2))+SIN(4*PI*(E2-E3))+
     &          SIN(4*PI*(E3-E1)))))
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO 
      RETURN
      END FUNCTION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! FUNCTION TO CALC STANDARD DEVIATION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      FUNCTION STDEV(X)
        REAL :: X(NX,NY,NZ)
        REAL :: STDEV
        REAL :: avg
! PERFORM CALCULATION
        avg = SUM(X)/(NXYZ)
        STDEV=SQRT((SUM((X-avg)**2))/(NXYZ))
        RETURN
      END FUNCTION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! FUNCTION TO CALCULATE THE FREE ENERGY RELATED TO CONCENTRATION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      FUNCTION CALCCONCFE(CON)
! DECLARE VARAIBLES
        REAL :: CALCCONCFE(NX,NY,NZ), CON(NX,NY,NZ)
        REAL :: CONCT(NX+4,NY+4,NZ+4)
! EXPRESSION FOR FREE ENERGY - SEE DIFFUSION DOC
! f = kappa + alpha*dG + fB + gradientEnergy
! ENSURING PBCs
        CONCT(3:NX+2, 3:NY+2, 3:NZ+2) = CON
        CONCT(1:2   , 3:NY+2, 3:NZ+2) = CON(NX-1:NX,    1:NY,    1:NZ) 
        CONCT(3:NX+2, 1:2   , 3:NZ+2) = CON(   1:NX, NY-1:NY,    1:NZ) 
        CONCT(3:NX+2, 3:NY+2,    1:2) = CON(   1:NX,    1:NY, NZ-1:NZ)
        CONCT(NX+3:NX+4, 3:NY+2   ,    3:NZ+2) = CON(1:2 , 1:NY, 1:NZ) 
        CONCT(3:NX+2   , NY+3:NY+4,    3:NZ+2) = CON(1:NX, 1:2 , 1:NZ) 
        CONCT(3:NX+2   , 3:NY+2   , NZ+3:NZ+4) = CON(1:NX, 1:NY,  1:2)
! INDIVIDUAL TERMS OF FREE ENERGY
        ALPDG=GDIFF*(0.5 + 0.5*TANH((A1-CON)/A2))
        FBENERGY=BETA1*(1-EXP(-((CON-C0)**2)) )
        CONCGE=EG*(((1/(2.0*DDX))*
     &    (CONCT(4:NX+3,3:NY+2,3:NZ+2)-CONCT(2:NX+1,3:NY+2,3:NZ+2)+
     &     CONCT(3:NX+2,4:NY+3,3:NZ+2)-CONCT(3:NX+2,2:NY+1,3:NZ+2)+
     &     CONCT(3:NX+2,3:NY+2,4:NZ+3)-CONCT(3:NX+2,3:NY+2,2:NZ+1)))**2)
! 
        CALCCONCFE=ALPDG + FBENERGY + CONCGE
        RETURN
      END FUNCTION
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      END SUBROUTINE DIF

