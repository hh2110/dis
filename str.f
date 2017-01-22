!-----------------------------------------------------------------------
! SUBROUTINE TO CALCULATE THE BPQ MATRIX
!-----------------------------------------------------------------------
      SUBROUTINE BMA
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES
      REAL :: OMEGA(3,3)          !! OMEGA
      REAL :: IOMEGA(3,3)         !! INVERSE OF OMEGA
      INTEGER :: I,J,K,L          !! TENSOR INDICES
      INTEGER :: IX,IY,IZ         !! POSITION INDICES
      INTEGER :: SA,SMA           !! INDICES OF SLIP PLANES AND DIRECTIONS
      INTEGER :: SAP,SMAP         !! SA', SMA'
      REAL :: C0IJKL(3,3,3,3)     !! COMPLIANCE MATRIX
      REAL :: TP1, TP2, TP3       !! TRANSFORMATION COSINES
      REAL :: LIIP,LJJP,LKKP,LLLP
      REAL :: N1(3,3), N2(3,3)    !! IDENTITY AND TRANSFORMATION MATRICES
      INTEGER :: IP,JP,KP,LP,KK   !! I',J',K',L' AND KK
! SET UP OPENMP STUFF     
      REAL(KIND=8) :: OMP_GET_WTIME      
! LOG THE TIME      
      BPQSTART=OMP_GET_WTIME()       
! ASSIGN INITIAL VALUES TO COMPLIANCE MATRIX
      DO I=1,3
      DO J=1,3
      DO K=1,3
      DO L=1,3
        C0IJKL(I,J,K,L)=0.0
      END DO
      END DO
      END DO
      END DO
!      
! ASSIGN ELASTIC CONSTANTS TO MATRIX 
      C0IJKL(1,1,1,1) =C11
      C0IJKL(2,2,2,2) =C11
      C0IJKL(3,3,3,3) =C11
      C0IJKL(1,1,2,2) =C12
      C0IJKL(2,2,3,3) =C12
      C0IJKL(3,3,1,1) =C12
      C0IJKL(2,2,1,1) =C12
      C0IJKL(3,3,2,2) =C12
      C0IJKL(1,1,3,3) =C12
      C0IJKL(1,2,1,2) =C44
      C0IJKL(2,3,2,3) =C44
      C0IJKL(3,1,3,1) =C44
      C0IJKL(2,1,2,1) =C44
      C0IJKL(1,2,2,1) =C44
      C0IJKL(2,1,1,2) =C44
      C0IJKL(3,2,3,2) =C44
      C0IJKL(3,2,2,3) =C44
      C0IJKL(2,3,3,2) =C44
      C0IJKL(1,3,1,3) =C44
      C0IJKL(1,3,3,1) =C44
      C0IJKL(3,1,1,3) =C44
!      
      IF (ROTFLAG == 0) THEN
        DO I=1,3
        DO J=1,3
        DO K=1,3
        DO L=1,3
          CIJKL(I,J,K,L)=C0IJKL(I,J,K,L)
        END DO
        END DO
        END DO
        END DO
      ELSE IF (ROTFLAG == 1) THEN  
! DEFINE IDENTITY AND TRANSFORMATION MATRICES
        CALL MAT(N1,1.,0.,0.,0.,1.,0.,0.,0.,1.)
!        
        DO I=1,3
        DO J=1,3
          N2(I,J)=ROT(I,J)
        END DO
        END DO
!        
        DO IP=1,3
        DO JP=1,3
        DO KP=1,3
        DO LP=1,3
          CIJKL(IP,JP,KP,LP) =0.0
          DO I=1,3
          DO J=1,3
          DO K=1,3
          DO L=1,3
! CALCULATE DOT PRODUCT N1(I)*N2(I') AND N1(J)*N2(J')
            LIIP =0.0
            LJJP =0.0
            LKKP =0.0
            LLLP =0.0
            DO KK=1,3
              LIIP=LIIP+N1(I,KK)*N2(IP,KK)
              LJJP=LJJP+N1(J,KK)*N2(JP,KK)
              LKKP=LKKP+N1(K,KK)*N2(KP,KK)
              LLLP=LLLP+N1(L,KK)*N2(LP,KK)
            END DO
            CIJKL(IP,JP,KP,LP)=CIJKL(IP,JP,KP,LP)+
     &        LIIP*LJJP*LKKP*LLLP*C0IJKL(I,J,K,L)
          END DO
          END DO
          END DO
          END DO
        END DO
        END DO
        END DO
        END DO
      END IF  
! DEFINE THE INTERPLANAR SPACING (FCC)
      SD =A0/SQRT3 
! DEFINITION TRANSFORMATION STRAIN E0
      DO I=1,3
      DO J=1,3
        DO SA=1,NP
        DO SMA=1,NQ
          E0(SA,SMA,I,J) =(BV(SA,SMA,I)*NV(SA,J)+BV(SA,SMA,J)*NV(SA,I))/
     &       SQRT( BV(SA,SMA,1)**2+BV(SA,SMA,2)**2+BV(SA,SMA,3)**2 )/2
        END DO
        END DO
      END DO
      END DO
!
      PRINT *,'>>NOTE - THE TRANSFORMATION STRAIN MATRICES'
      PRINT *,'>>NOTE - SLPLANE \\ SLPDIR \\ STRAIN'
      DO SA=1,NP
        PRINT '(I3)',SA
        DO SMA=1,NQ
          PRINT '(I3)',SMA
          PRINT'(F9.5,F9.5,F9.5)',
     &          E0(SA,SMA,1,1),E0(SA,SMA,1,2),E0(SA,SMA,1,3)
          PRINT'(F9.5,F9.5,F9.5)',
     &          E0(SA,SMA,2,1),E0(SA,SMA,2,2),E0(SA,SMA,2,3)
          PRINT'(F9.5,F9.5,F9.5)',
     &          E0(SA,SMA,3,1),E0(SA,SMA,3,2),E0(SA,SMA,3,3)
          PRINT *,' '
        END DO
        PRINT *,' '
      END DO
! DEFINITION TRANSFORMATION STRESS S0
      DO SA=1,NP
      DO SMA=1,NQ
        DO I=1,3
        DO J=1,3
          S0(SA,SMA,I,J) =0.0
          DO L=1,3
          DO K=1,3
            S0(SA,SMA,I,J)=S0(SA,SMA,I,J)+ 
     &        CIJKL(I,J,K,L)*E0(SA,SMA,K,L)
          END DO
          END DO
        END DO
        END DO
      END DO
      END DO
!      PRINT *, S0(SA,SMA,1,1),S0(SA,SMA,1,2),S0(SA,SMA,1,3)
!      PRINT *, S0(SA,SMA,2,1),S0(SA,SMA,2,2),S0(SA,SMA,2,3)
!      PRINT *, S0(SA,SMA,3,1),S0(SA,SMA,3,2),S0(SA,SMA,3,3)
! CALCULATION OF BPQ MATRIX
!
!$OMP PARALLEL DO PRIVATE(IY,IZ,OMEGA,IOMEGA)
      DO IX=1,NX
      DO IY=1,NY
      DO IZ=1,NZ
! CALCULATE OMEGA
        IF(IX  ==  1 .AND. IY  ==  1 .AND. IZ  ==  1) THEN
          DO I=1,3
          DO J=1,3
            IOMEGA(I,J)=0.0
          END DO
          END DO
        ELSE
          DO I=1,3
          DO J=1,3
            IOMEGA(I,J)=0.0
            DO K=1,3
            DO L=1,3
              IOMEGA(I,J)=IOMEGA(I,J)+
     &          CIJKL(I,K,L,J)*G(K,IX,IY,IZ)*G(L,IX,IY,IZ)/
     &          G2(IX,IY,IZ)
            END DO
            END DO
          END DO
          END DO
        END IF
!        
        IF(IX  ==  1 .AND. IY  ==  1 .AND. IZ  ==  1) THEN
          DO I=1,3
          DO J=1,3
            OMEGA(I,J) =0.0
          END DO
          END DO
        ELSE
          CALL INV(IOMEGA,OMEGA)
        END IF
!        
        IF(IX  ==  1 .AND. IY  ==  1 .AND. IZ  ==  1) THEN
          DO SA =1,NP
          DO SMA=1,NQ
          DO SAP=1,NP
          DO SMAP=1,NQ
            BPQ(SA,SMA,SAP,SMAP,IX,IY,IZ) =0.0
          END DO
          END DO
          END DO
          END DO        
        ELSE
          DO SA =1,NP
          DO SMA=1,NQ
          DO SAP=1,NP
          DO SMAP=1,NQ
            BPQ(SA,SMA,SAP,SMAP,IX,IY,IZ)=0.0            
            DO I=1,3
            DO J=1,3
            DO K=1,3
            DO L=1,3
              BPQ(SA,SMA,SAP,SMAP,IX,IY,IZ)=
     &          BPQ(SA,SMA,SAP,SMAP,IX,IY,IZ)+
     &          CIJKL(I,J,K,L)*E0(SA,SMA,I,J)*E0(SAP,SMAP,K,L)-
     &          G(I,IX,IY,IZ)*S0(SA,SMA,I,J)*OMEGA(J,K)*
     &          S0(SAP,SMAP,K,L)*G(L,IX,IY,IZ)/G2(IX,IY,IZ)
            END DO
            END DO
            END DO
            END DO
          END DO
          END DO
          END DO
          END DO
        END IF
!        
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO
!
! CALCULATE BPQ INITIALISATION TIME   
      BPQFINISH=OMP_GET_WTIME()
      BPQRUN=BPQFINISH-BPQSTART
      IF(BPQRUN < 60.0) THEN
      PRINT '(" >>NOTE - BPQ INITIALISATION TIME = ",F8.5," SECONDS.")',
     & BPQRUN
      ELSE IF(BPQRUN >= 60.0 .AND. BPQRUN<3600.0) THEN
      PRINT '(" >>NOTE - BPQ INITIALISATION TIME = ",F8.5," MINUTES.")',
     & BPQRUN/60.0
      ELSE IF(BPQRUN>=3600.0) THEN
      PRINT '(" >>NOTE - BPQ INITIALISATION TIME = ",ES14.5," HOURS.")',
     & BPQRUN/3600.0
      END IF
!      
      END SUBROUTINE BMA
!-----------------------------------------------------------------------
! SUBROUTINE TO DEFINE APPLIED STRESS MATRIX
!-----------------------------------------------------------------------
      SUBROUTINE STR
      USE VAR
      IMPLICIT NONE
! DEFINE VARIABLES
      REAL :: AS(3,3)          !! STRESS COMPONENT MATRIX 
      REAL :: INVROT(3,3)      !! INVERSE OF ROTATION MATRIX
      REAL :: S_CRYST(3,3)     !! APPLIED STRESS MATRIX WRT TRUE CRYSTAL AXES 
      REAL :: TP1,TP2,TP3      !! TRANSFORMATION COSINES
      INTEGER :: I,J           !! TENSOR INDICES
! ASSIGN COMPONENTS TO AS
      CALL MAT(AS,S11,S12,S13,S21,S22,S23,S31,S32,S33)
! ASSIGN STRESS MAGNITUDE
      DO I=1,3
      DO J=1,3
        AS(I,J) =AS(I,J)*S_APP_MAG
      END DO
      END DO
!      
      DO I=1,3
      DO J=1,3
        S_APP(I,J) =AS(I,J)
      END DO
      END DO      
! PRINT APPLIED STRESS MATRIX TO SCREEN
      PRINT *,''
      PRINT *,'>>NOTE - APPLIED STRESS MATRIX (REDUCED)'
      PRINT '(F8.5,F8.5,F8.5)',S_APP(1,1),S_APP(1,2),S_APP(1,3)
      PRINT '(F8.5,F8.5,F8.5)',S_APP(2,1),S_APP(2,2),S_APP(2,3)
      PRINT '(F8.5,F8.5,F8.5)',S_APP(3,1),S_APP(3,2),S_APP(3,3)
      PRINT *,''
      PRINT *,'>>NOTE - APPLIED STRESS MATRIX (MPa)'
      PRINT '(F8.0,F8.0,F8.0)',S_APP(1,1)*C44X*1000.0,
     & S_APP(1,2)*C44X*1000.0,S_APP(1,3)*C44X*1000.0
      PRINT '(F8.0,F8.0,F8.0)',S_APP(2,1)*C44X*1000.0,
     & S_APP(2,2)*C44X*1000.0,S_APP(2,3)*C44X*1000.0
      PRINT '(F8.0,F8.0,F8.0)',S_APP(3,1)*C44X*1000.0,
     & S_APP(3,2)*C44X*1000.0,S_APP(3,3)*C44X*1000.0
! CALCULATE APPLIED STRESS MATRIX FOR ORIGINAL CRYSTAL AXES (BEFORE ROTATION)     
      IF (ROTFLAG == 1) THEN
        CALL INV(ROT,INVROT)
        CALL MUL(S_CRYST,INVROT,S_APP)
      PRINT *,''
      PRINT *,'>>NOTE - APPLIED STRESS W.R.T. CRYSTAL AXES (REDUCED)'
      PRINT '(F8.5,F8.5,F8.5)',S_CRYST(1,1),S_CRYST(1,2),S_CRYST(1,3)
      PRINT '(F8.5,F8.5,F8.5)',S_CRYST(2,1),S_CRYST(2,2),S_CRYST(2,3)
      PRINT '(F8.5,F8.5,F8.5)',S_CRYST(3,1),S_CRYST(3,2),S_CRYST(3,3)
      PRINT *,''
      PRINT *,'>>NOTE - APPLIED STRESS W.R.T. CRYSTAL AXES (MPa)'
      PRINT '(F8.0,F8.0,F8.0)',S_CRYST(1,1)*C44X*1000.0,
     & S_CRYST(1,2)*C44X*1000.0,S_CRYST(1,3)*C44X*1000.0
      PRINT '(F8.0,F8.0,F8.0)',S_CRYST(2,1)*C44X*1000.0,
     & S_CRYST(2,2)*C44X*1000.0,S_CRYST(2,3)*C44X*1000.0
      PRINT '(F8.0,F8.0,F8.0)',S_CRYST(3,1)*C44X*1000.0,
     & S_CRYST(3,2)*C44X*1000.0,S_CRYST(3,3)*C44X*1000.0      
      END IF
      END SUBROUTINE STR    
