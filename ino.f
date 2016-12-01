!-----------------------------------------------------------------------
! SUBROUTINE TO DEFINE VARIABLE PARAMETERS
!-----------------------------------------------------------------------
      SUBROUTINE PAR
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES
!
! READ PARAMETER VALUES FROM EXTERNAL INPUT FILE (par.inp)
      OPEN(UNIT=11, FILE='par.inp', STATUS='UNKNOWN')
      READ(11,100) FFTP,NOMP,OUTGE,OUTFE,OUTETAS,INITBENCH,
     &         NX,NY,NZ,BORDER,T0,DT,OUT1,OUTDIM,SLICE_AXIS,SLICE_POS,
     &         AA,C11X,C12X,C44X,DXD,S_APP_MAG,
     &         S11,S12,S13,S21,S22,S23,S31,S32,S33,
     &         BETA,NSS,ROTFLAG,AXX(1),AXX(2),AXX(3),
     &         AXY(1),AXY(2),AXY(3),AXZ(1),AXZ(2),AXZ(3),
     &         NPHIS,GSURF(1),GSURF(2),NIN
! PARAMETER LIST
100   FORMAT(19X, I3/,
     &       19X, I3/,
     &       19X, I3/,     
     &       19X, I3/,
     &       19X, I3/,     
     &       19X, I3/,
     &       26X/,    
     &       19X, I4/,
     &       19X, I4/,
     &       19X, I4/,
     &       19X, I3/,
     &       19X, F6.2/,
     &       19X, F6.5/,
     &       19X, I5/,
     &       19X, I3/, 
     &       19X, I3/,
     &       19X, I4/,
     &       26X/,
     &       19X, F6.4/,
     &       19X, F6.2/,
     &       19X, F6.2/,
     &       19X, F6.2/,
     &       19X, F6.3/,
     &       19X, F6.5/,
     &       19X, F6.2/,
     &       19X, F6.2/,
     &       19X, F6.2/,
     &       19X, F6.2/,
     &       19X, F6.2/,
     &       19X, F6.2/,
     &       19X, F6.2/,
     &       19X, F6.2/,
     &       19X, F6.2/,     
     &       19X, F6.5/,
     &       26X/,     
     &       19X, I3/,
     &       19X, I3/,
     &       19X, F2.0, F2.0, F2.0/,  
     &       19X, F2.0, F2.0, F2.0/, 
     &       19X, F2.0, F2.0, F2.0/,
     &       26X/,       
     &       19X, I3/,
     &       19X, I3/,     
     &       19X, I3/,
     &       26X/,       
     &       19X, I3)

      CLOSE(11)
! PRINT INPUT VARIABLES TO SCREEN
      PRINT *,''
      PRINT *,'>>NOTE - CHECK SIMULATION PARAMETERS'
      PRINT '(" FFT PLAN TYPE     =",I8)', FFTP
      PRINT '(" OPENMP THREADS    =",I8)', NOMP
      PRINT '(" OUTPUT GRAD. E.   =",I8)', OUTGE
      PRINT '(" OUTPUT FAULT E.   =",I8)', OUTFE
      PRINT '(" OUTPUT ETAS       =",I8)', OUTETAS
      PRINT '(" INIT. BENCHMARK   =",I8)', INITBENCH
      PRINT '(" ---------------------------")'        
      PRINT '(" SYSTEM SIZE X     =",I8)', NX
      PRINT '(" SYSTEM SIZE Y     =",I8)', NY
      PRINT '(" SYSTEM SIZE Z     =",I8)', NZ
      PRINT '(" BORDER            =",I8)', BORDER
      PRINT '(" TOTAL TIME        =",F8.2)',T0
      PRINT '(" TIME STEP         =",F8.5)',DT
      PRINT '(" NUMBER OF OUTPUTS =",I8)',OUT1
      PRINT '(" 2D OR 3D OUTPUT?  =",I8)', OUTDIM
      PRINT '(" 2D SLICE AXIS     =",I8)', SLICE_AXIS
      PRINT '(" 2D SLICE POSITION =",I8)', SLICE_POS
      PRINT '(" ---------------------------")'     
      PRINT '(" LATTICE PARAMETER =",F8.5)',AA
      PRINT '(" C11               =",F8.2)',C11X
      PRINT '(" C12               =",F8.2)',C12X
      PRINT '(" C44               =",F8.2)',C44X
      PRINT '(" SCALING FACTOR    =",F8.3)',DXD      
      PRINT '(" APPLIED STRESS    =",F8.5)',S_APP_MAG
      PRINT '(" S11               =",F8.2)',S11
      PRINT '(" S12               =",F8.2)',S12
      PRINT '(" S13               =",F8.2)',S13
      PRINT '(" S21               =",F8.2)',S21
      PRINT '(" S22               =",F8.2)',S22
      PRINT '(" S23               =",F8.2)',S23
      PRINT '(" S31               =",F8.2)',S31
      PRINT '(" S32               =",F8.2)',S32
      PRINT '(" S33               =",F8.2)',S33       
      PRINT '(" BETA (GRAD COEF.) =",F8.5)',BETA
      PRINT '(" ---------------------------")'        
      PRINT '(" SLIP SYSTEM       =",I8)',NSS
      PRINT '(" ROTATE CRYSTAL?   =",I8)',ROTFLAG 
      PRINT '(" ROTATED X-AXIS    =",F3.0,F3.0,F3.0)'
     &,AXX(1),AXX(2),AXX(3)
      PRINT '(" ROTATED Y-AXIS    =",F3.0,F3.0,F3.0)'
     &,AXY(1),AXY(2),AXY(3) 
      PRINT '(" ROTATED Z-AXIS    =",F3.0,F3.0,F3.0)'
     &,AXZ(1),AXZ(2),AXZ(3)
      PRINT '(" ---------------------------")'
      PRINT '(" NO. OF PHASES     =",I8)',NPHIS      
      PRINT '(" GAMMA SURFACE G   =",I8)',GSURF(1)
      PRINT '(" GAMMA SURFACE GP  =",I8)',GSURF(2)
      PRINT '(" ---------------------------")'        
      PRINT '(" INITIAL CONDITION =",I8)',NIN
! VERIFY PARAMETERS
!      
      IF(INITBENCH /= 0 .AND. INITBENCH /=1) THEN
      STOP '>>ERROR - INVALID INITBENCH IN PAR.INP, ENTER 0-NO or 1-YES'
      END IF
! 
      IF(ROTFLAG /= 0 .AND. ROTFLAG /=1) THEN
      STOP '>>ERROR - INVALID ROTFLAG IN PAR.INP - ENTER 0-NO or 1-YES'
      END IF
!      
      IF(OUTDIM /= 2 .AND. OUTDIM /=3) THEN
      STOP '>>ERROR - INVALID OUTDIM IN PAR.INP - ENTER 2 OR 3 ONLY'
      END IF
!      
      IF(SLICE_AXIS /= 1 .AND. SLICE_AXIS /=2 
     & .AND. SLICE_AXIS /=3) THEN
      STOP '>>ERROR - INVALID SLIVE_AX IN PAR.INP, ENTER 1,2 OR 3 ONLY'
      END IF
!      
! DEFINE MATERIAL ELASTIC CONSTANTS
        C11=C11X/C44X
        C12=C12X/C44X
        C44=C44X/C44X
! PRINT ELASTIC COSNTANTS TO SCREEN
      PRINT *,''
      PRINT *,'>>NOTE - ELASTIC CONSTANTS, REDUCED and SCALED VALUES'
      PRINT '(" C11= ",F6.2," GPa,",F8.2,",",F8.2)',C11X,C11,C11/DXD
      PRINT '(" C12= ",F6.2," GPa,",F8.2,",",F8.2)',C12X,C12,C12/DXD
      PRINT '(" C44= ",F6.2," GPa,",F8.2,",",F8.2)',C44X,C44,C44/DXD
        C11=C11/DXD
        C12=C12/DXD
        C44=C44/DXD
!                
      END SUBROUTINE PAR
!-----------------------------------------------------------------------
! SUBROUTINE TO OUTPUT GRADIENT ENERGY
!-----------------------------------------------------------------------
      SUBROUTINE OGE(STEPID)
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES     
      INTEGER :: STEPID
      INTEGER :: IX,IY,IZ
      INTEGER :: SA,SMA,SAP,SMAP
      REAL :: B_ALPHA, B_ALPHAP, TP
! CALL GRAD-ETA
        CALL GRE
!
!$OMP PARALLEL DO PRIVATE(IY,IX)
      DO IZ=1,NZ
      DO IY=1,NY
      DO IX=1,NX
        TPR(IX,IY,IZ) =0.0
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO 
!
!
!$OMP PARALLEL DO PRIVATE(IY,IX,SA,SMA,SAP)
!$OMP+ PRIVATE(B_ALPHA,B_ALPHAP,TP) REDUCTION(+:TPR)
      DO IZ=1,NZ
      DO IY=1,NY
      DO IX=1,NX
        DO SA =1,NP
        DO SMA=1,NQ
        DO SAP =1,NP
        DO SMAP=1,NQ
          B_ALPHA=SQRT( BV(SA,SMA,1)**2+BV(SA,SMA,2)**2 
     &      +BV(SA,SMA,3)**2 )
!     
          B_ALPHAP=SQRT( BV(SAP,SMAP,1)**2+BV(SAP,SMAP,2)**2 
     &      +BV(SAP,SMAP,3)**2 )
!     
          TP=(BV(SA,SMA,1)*BV(SAP,SMAP,1)+
     &      BV(SA,SMA,2)*BV(SAP,SMAP,2)+
     &      BV(SA,SMA,3)*BV(SAP,SMAP,3))/B_ALPHA/B_ALPHAP
!     
          TPR(IX,IY,IZ)=TPR(IX,IY,IZ)+
     &      ( (NV(SA,2)*DER(SA,SMA,3,IX,IY,IZ)-
     &      NV(SA,3)*DER(SA,SMA,2,IX,IY,IZ))*
     &      (NV(SAP,2)*DER(SAP,SMAP,3,IX,IY,IZ)-
     &      NV(SAP,3)*DER(SAP,SMAP,2,IX,IY,IZ))+
     &      (NV(SA,3)*DER(SA,SMA,1,IX,IY,IZ)-
     &      NV(SA,1)*DER(SA,SMA,3,IX,IY,IZ))*
     &      (NV(SAP,3)*DER(SAP,SMAP,1,IX,IY,IZ)-
     &      NV(SAP,1)*DER(SAP,SMAP,3,IX,IY,IZ))+
     &      (NV(SA,1)*DER(SA,SMA,2,IX,IY,IZ)-
     &      NV(SA,2)*DER(SA,SMA,1,IX,IY,IZ))*
     &      (NV(SAP,1)*DER(SAP,SMAP,2,IX,IY,IZ)-
     &      NV(SAP,2)*DER(SAP,SMAP,1,IX,IY,IZ)))*TP
        END DO
        END DO
        END DO
        END DO
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO       
!
!$OMP PARALLEL DO PRIVATE(IY,IX)
      DO IZ=1,NZ
      DO IY=1,NY
      DO IX=1,NX
         TPR(IX,IY,IZ)=TPR(IX,IY,IZ)/2.0*BETA
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO 
!
! OUTPUT GE
      IF (OUTDIM == 2) THEN
        IF (SLICE_AXIS == 1) THEN      
          WRITE(20,*) 'T=',STEPID*DT
          WRITE(20,*) ((TPR(SLICE_POS,IY,IZ),IZ=1,NZ),IY=1,NY)
        ELSE IF (SLICE_AXIS == 2) THEN      
          WRITE(20,*) 'T=',STEPID*DT
          WRITE(20,*) ((TPR(IX,SLICE_POS,IZ),IZ=1,NZ),IX=1,NX)
        ELSE IF (SLICE_AXIS == 3) THEN      
          WRITE(20,*) 'T=',STEPID*DT
          WRITE(20,*) ((TPR(IX,IY,SLICE_POS),IY=1,NY),IX=1,NX)
        END IF  
      ELSE IF (OUTDIM == 3) THEN
        WRITE(20,*) 'T=',STEPID*DT
        WRITE(20,*) (((TPR(IX,IY,IZ),IZ=1,NZ),IY=1,NY),IX=1,NX)
      END IF  
!
      END SUBROUTINE OGE      
!-----------------------------------------------------------------------
! SUBROUTINE TO OUTPUT FAULT ('CRYSTALLINE') ENERGY
!-----------------------------------------------------------------------
      SUBROUTINE OFE(STEPID)
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES      
      INTEGER :: STEPID       !! SIMULATION ITERATION
      INTEGER :: IX,IY,IZ     !! POSITION INDICES
      INTEGER ::  SA,SMA      !! SLIP SYSTEM, SLIP DIRECTION
      REAL :: FE(NX,NY,NZ) !! GSF FAULT ENERGY
      REAL :: E1,E2,E3     !! TEMPORARY ETAS
! ASSIGN INITIAL VALUE TO FE     
      DO IX=1,NX
      DO IY=1,NY
      DO IZ=1,NZ
        FE(IX,IY,IZ) =0.0
      END DO
      END DO
      END DO
! 
      DO SA =1,NP
!$OMP PARALLEL DO PRIVATE(IY,IX,E1,E2,E3,PHID)      
      DO IZ=1,NZ
      DO IY=1,NY
      DO IX=1,NX
        E1=ETA(SA,1,IX,IY,IZ)
        E2=ETA(SA,2,IX,IY,IZ)
        E3=ETA(SA,3,IX,IY,IZ)
!-----------------------------------------------------------------------
! CRYSTALLINE ENERGY FOR GAMMA PRIME PHASE      
!-----------------------------------------------------------------------
        PHID=INT(PHI(IX,IY,IZ)+1)
!-----------------------------------------------------------------------
! PURE FCC SCHEME #1 WITH [1 LAYER]
!-----------------------------------------------------------------------        
        IF(GSF_EQ(PHID) == 10) THEN
          FE(IX,IY,IZ) =FE(IX,IY,IZ) + (
     &      GC(PHID,1)+
     &      GC(PHID,4)*(COS(TWOPI*(E1-E2))+COS(TWOPI*(E2-E3))+
     &        COS(TWOPI*(E3-E1))) +
     &      GC(PHID,7)*(COS(TWOPI*(2*E1-E2-E3))+
     &        COS(TWOPI*(2*E2-E3-E1))+
     &        COS(TWOPI*(2*E3-E1-E2))) +
     &      GC(PHID,8)*(COS(4*PI*(E1-E2))+COS(4*PI*(E2-E3))+
     &        COS(4*PI*(E3-E1))) +
     &      GC(PHID,9)*(COS(TWOPI*(3*E1-E2-2*E3))+
     &        COS(TWOPI*(3*E1-2*E2-E3))+
     &        COS(TWOPI*(3*E2-E3-2*E1))+COS(TWOPI*(3*E2-2*E3-E1))+
     &        COS(TWOPI*(3*E3-E1-2*E2))+COS(TWOPI*(3*E3-2*E1-E2)))+     
     &      GC(PHID,11)*(SIN(TWOPI*(E1-E2))+SIN(TWOPI*(E2-E3))+
     &        SIN(TWOPI*(E3-E1))) +
     &      GC(PHID,14)*(SIN(4*PI*(E1-E2))+SIN(4*PI*(E2-E3))+
     &        SIN(4*PI*(E3-E1)))+
     &      GC(PHID,15)*(SIN(TWOPI*(2*E1-3*E2+E3))+
     &        SIN(TWOPI*(3*E1-2*E2-E3))+
     &        SIN(TWOPI*(-2*E1-E2+3*E3))+SIN(TWOPI*(E1+2*E2-3*E3))+
     &        SIN(TWOPI*(-3*E1+E2+2*E3))+SIN(TWOPI*(-E1+3*E2-2*E3))))
!-----------------------------------------------------------------------
! PURE L1_2 ORDERED COMPOUND SCHEME #1A [1 LAYER]
!-----------------------------------------------------------------------      
        ELSE IF(GSF_EQ(PHID) == 11) THEN
          FE(IX,IY,IZ) =FE(IX,IY,IZ) + (
     &      GC(PHID,1)+
     &      GC(PHID,2)*(COS(PI*(E1-E2))+COS(PI*(E2-E3))+                                             
     &        COS(PI*(E3-E1))) +
     &      GC(PHID,3)*(COS(PI*(2*E1-E2-E3))+
     &        COS(PI*(2*E2-E3-E1))+
     &        COS(PI*(2*E3-E1-E2))) +
     &      GC(PHID,4)*(COS(TWOPI*(E1-E2))+COS(TWOPI*(E2-E3))+
     &        COS(TWOPI*(E3-E1))) +
     &      GC(PHID,5)*(COS(PI*(3*E1-E2-2*E3))+COS(PI*(3*E1-2*E2-E3))+
     &        COS(PI*(3*E2-E3-2*E1))+COS(PI*(3*E2-2*E3-E1))+
     &        COS(PI*(3*E3-E1-2*E2))+COS(PI*(3*E3-2*E1-E2)))+
     &      GC(PHID,6)*(COS(3*PI*(E1-E2))+COS(3*PI*(E2-E3))+
     &        COS(3*PI*(E3-E1))) +
     &      GC(PHID,8)*(COS(4*PI*(E1-E2))+COS(4*PI*(E2-E3))+
     &        COS(4*PI*(E3-E1))) +    
     &      GC(PHID,10)*(SIN(PI*(E1-E2))+SIN(PI*(E2-E3))+
     &        SIN(PI*(E3-E1))) +
     &      GC(PHID,11)*(SIN(TWOPI*(E1-E2))+SIN(TWOPI*(E2-E3))+
     &        SIN(TWOPI*(E3-E1))) +
     &      GC(PHID,12)*(SIN(PI*(2*E1-3*E2+E3))+SIN(PI*(3*E1-2*E2-E3))+
     &        SIN(PI*(-2*E1-E2+3*E3))+SIN(PI*(E1+2*E2-3*E3))+
     &        SIN(PI*(-3*E1+E2+2*E3))+SIN(PI*(-E1+3*E2-2*E3)))+
     &      GC(PHID,13)*(SIN(3*PI*(E1-E2))+SIN(3*PI*(E2-E3))+
     &        SIN(3*PI*(E3-E1)))+
     &      GC(PHID,14)*(SIN(4*PI*(E1-E2))+SIN(4*PI*(E2-E3))+
     &        SIN(4*PI*(E3-E1)))
     &      )
!-----------------------------------------------------------------------
! PURE L1_2 ORDERED COMPOUND SCHEME #1B [1 LAYER]
!-----------------------------------------------------------------------      
        ELSE IF(GSF_EQ(PHID) == 12) THEN
          FE(IX,IY,IZ) =FE(IX,IY,IZ) + (
     &      GC(PHID,1)+
     &      GC(PHID,2)*(COS(PI*(E1-E2))+COS(PI*(E2-E3))+                                             
     &        COS(PI*(E3-E1))) +
     &      GC(PHID,3)*(COS(PI*(2*E1-E2-E3))+
     &        COS(PI*(2*E2-E3-E1))+
     &        COS(PI*(2*E3-E1-E2))) +
     &      GC(PHID,4)*(COS(TWOPI*(E1-E2))+COS(TWOPI*(E2-E3))+
     &        COS(TWOPI*(E3-E1))) +
     &      GC(PHID,5)*(COS(PI*(3*E1-E2-2*E3))+COS(PI*(3*E1-2*E2-E3))+
     &        COS(PI*(3*E2-E3-2*E1))+COS(PI*(3*E2-2*E3-E1))+
     &        COS(PI*(3*E3-E1-2*E2))+COS(PI*(3*E3-2*E1-E2)))+
     &      GC(PHID,6)*(COS(3*PI*(E1-E2))+COS(3*PI*(E2-E3))+
     &        COS(3*PI*(E3-E1))) +
     &      GC(PHID,7)*(COS(TWOPI*(2*E1-E2-E3))+
     &        COS(TWOPI*(2*E2-E3-E1))+
     &        COS(TWOPI*(2*E3-E1-E2))) +
     &      GC(PHID,8)*(COS(4*PI*(E1-E2))+COS(4*PI*(E2-E3))+
     &        COS(4*PI*(E3-E1))) +   
     &      GC(PHID,10)*(SIN(PI*(E1-E2))+SIN(PI*(E2-E3))+
     &        SIN(PI*(E3-E1))) +
     &      GC(PHID,11)*(SIN(TWOPI*(E1-E2))+SIN(TWOPI*(E2-E3))+
     &        SIN(TWOPI*(E3-E1))) +
     &      GC(PHID,12)*(SIN(PI*(2*E1-3*E2+E3))+SIN(PI*(3*E1-2*E2-E3))+
     &        SIN(PI*(-2*E1-E2+3*E3))+SIN(PI*(E1+2*E2-3*E3))+
     &        SIN(PI*(-3*E1+E2+2*E3))+SIN(PI*(-E1+3*E2-2*E3)))+
     &      GC(PHID,13)*(SIN(3*PI*(E1-E2))+SIN(3*PI*(E2-E3))+
     &        SIN(3*PI*(E3-E1)))+
     &      GC(PHID,14)*(SIN(4*PI*(E1-E2))+SIN(4*PI*(E2-E3))+
     &        SIN(4*PI*(E3-E1)))
     &      )
!-----------------------------------------------------------------------
! PURE FCC SCHEME #2 WITH [2 LAYERS]
!-----------------------------------------------------------------------     
        ELSE IF(GSF_EQ(PHID) == 20) THEN  
          FE(IX,IY,IZ) =FE(IX,IY,IZ)+ (
     &      GC(PHID,1)+ 
     &      GC(PHID,4)*(COS(2*(E1-E2)*PI)+COS(2*(E2-E3)*PI)+ 
     &        COS(2*(-E1+E3)*PI))+ 
     &      GC(PHID,5)*(COS(2*(2*E1-E2-E3)*PI)+COS(2*(-E1+2*E2-E3)*PI)+ 
     &        COS(2*(-E1-E2+2*E3)*PI))+ 
     &      GC(PHID,6)*(COS(3*(E1-E2)*PI)+COS(3*(E2-E3)*PI)+ 
     &        COS(3*(-E1+E3)*PI))+ 
     &      GC(PHID,7)*(COS(3*(2*E1-E2-E3)*PI)+COS(3*(-E1+2*E2-E3)*PI)+ 
     &        COS(3*(-E1-E2+2*E3)*PI))+ 
     &      GC(PHID,8)*(COS(4*(E1-E2)*PI)+COS(4*(E2-E3)*PI)+ 
     &        COS(4*(-E1+E3)*PI))+ 
     &      GC(PHID,11)*(SIN(2*(E1-E2)*PI)+SIN(2*(E2-E3)*PI)+ 
     &        SIN(2*(-E1+E3)*PI))+ 
     &      GC(PHID,13)*(SIN(4*(E1-E2)*PI)+SIN(4*(E2-E3)*PI)+ 
     &        SIN(4*(-E1+E3)*PI)) 
     &      )
!-----------------------------------------------------------------------
! PURE L1_2 ORDERED COMPOUND SCHEME #2B [2 LAYERS] SOLVER 4
!-----------------------------------------------------------------------     
        ELSE IF(GSF_EQ(PHID) == 21) THEN  
          FE(IX,IY,IZ) =FE(IX,IY,IZ)+ (
     &      GC(PHID,1)+ 
     &      GC(PHID,2)*(COS((E1-E2)*PI)+COS((E2-E3)*PI)+ 
     &        COS((-E1+E3)*PI))+ 
     &      GC(PHID,3)*(COS((2*E1-E2-E3)*PI)+COS((-E1+2*E2-E3)*PI)+ 
     &        COS((-E1-E2+2*E3)*PI))+ 
     &      GC(PHID,4)*(COS(2*(E1-E2)*PI)+COS(2*(E2-E3)*PI)+ 
     &        COS(2*(-E1+E3)*PI))+ 
     &      GC(PHID,5)*(COS(2*(2*E1-E2-E3)*PI)+COS(2*(-E1+2*E2-E3)*PI)+ 
     &        COS(2*(-E1-E2+2*E3)*PI))+ 
     &      GC(PHID,6)*(COS(3*(E1-E2)*PI)+COS(3*(E2-E3)*PI)+ 
     &        COS(3*(-E1+E3)*PI))+ 
     &      GC(PHID,7)*(COS(3*(2*E1-E2-E3)*PI)+COS(3*(-E1+2*E2-E3)*PI)+ 
     &        COS(3*(-E1-E2+2*E3)*PI))+ 
     &      GC(PHID,8)*(COS(4*(E1-E2)*PI)+COS(4*(E2-E3)*PI)+ 
     &        COS(4*(-E1+E3)*PI))+ 
     &      GC(PHID,9)*(COS(5*(E1-E2)*PI)+COS(5*(E2-E3)*PI)+ 
     &        COS(5*(-E1+E3)*PI))+ 
     &      GC(PHID,10)*(SIN((E1-E2)*PI)+SIN((E2-E3)*PI)+ 
     &        SIN((-E1+E3)*PI))+ 
     &      GC(PHID,11)*(SIN(2*(E1-E2)*PI)+SIN(2*(E2-E3)*PI)+ 
     &        SIN(2*(-E1+E3)*PI))+ 
     &      GC(PHID,12)*(SIN(3*(E1-E2)*PI)+SIN(3*(E2-E3)*PI)+ 
     &        SIN(3*(-E1+E3)*PI))+ 
     &      GC(PHID,13)*(SIN(4*(E1-E2)*PI)+SIN(4*(E2-E3)*PI)+ 
     &        SIN(4*(-E1+E3)*PI))+ 
     &      GC(PHID,14)*(SIN(5*(E1-E2)*PI)+SIN(5*(E2-E3)*PI)+ 
     &        SIN(5*(-E1+E3)*PI)))     
        END IF
!
      END DO
      END DO
      END DO
!$OMP END PARALLEL DO      
      END DO
! OUTPUT FE
      IF (OUTDIM ==2) THEN
        IF (SLICE_AXIS == 1) THEN        
          WRITE(21,*) 'T=',STEPID*DT
          WRITE(21,*) ((FE(SLICE_POS,IY,IZ),IZ=1,NZ),IY=1,NY)
        ELSE IF (SLICE_AXIS == 2) THEN        
          WRITE(21,*) 'T=',STEPID*DT
          WRITE(21,*) ((FE(IX,SLICE_POS,IZ),IZ=1,NZ),IX=1,NX)
        ELSE IF (SLICE_AXIS == 3) THEN        
          WRITE(21,*) 'T=',STEPID*DT
          WRITE(21,*) ((FE(IX,IY,SLICE_POS),IY=1,NY),IX=1,NX)
        END IF  
      ELSE IF (OUTDIM == 3) THEN
        WRITE(21,*) 'T=',STEPID*DT
        WRITE(21,*) (((FE(IX,IY,IZ),IZ=1,NZ),IY=1,NY),IX=1,NX)      
      END IF
! 
      END SUBROUTINE OFE
!-----------------------------------------------------------------------
! SUBROUTINE TO OUTPUT FIELD VARIABLES
!-----------------------------------------------------------------------
      SUBROUTINE OET(STEPID)
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES
      INTEGER :: STEPID
      INTEGER :: IX,IY,IZ
      INTEGER :: SA,SMA
      INTEGER :: UNITCOUNT   !! OUTPUT FILE NUMBER (UNIT)        
!
      UNITCOUNT=24
!
      IF (OUTDIM == 2) THEN
        IF (SLICE_AXIS == 1) THEN        
          DO SA=1,NP
          DO SMA=1,NQ
            WRITE(UNITCOUNT,*) 'T=',STEPID*DT
            WRITE(UNITCOUNT,*) ((ETA(SA,SMA,SLICE_POS,IY,IZ)
     &      ,IZ=1,NZ),IY=1,NY)
            UNITCOUNT=UNITCOUNT+1
          END DO
          END DO
        ELSE IF (SLICE_AXIS == 2) THEN        
          DO SA=1,NP
          DO SMA=1,NQ
            WRITE(UNITCOUNT,*) 'T=',STEPID*DT
            WRITE(UNITCOUNT,*) ((ETA(SA,SMA,IX,SLICE_POS,IZ)
     &      ,IZ=1,NZ),IX=1,NX)
            UNITCOUNT=UNITCOUNT+1
          END DO
          END DO
        ELSE IF (SLICE_AXIS == 3) THEN        
          DO SA=1,NP
          DO SMA=1,NQ
            WRITE(UNITCOUNT,*) 'T=',STEPID*DT
            WRITE(UNITCOUNT,*) ((ETA(SA,SMA,IX,IY,SLICE_POS)
     &      ,IY=1,NY),IX=1,NX)
            UNITCOUNT=UNITCOUNT+1
          END DO
          END DO
        END IF  
!        
      ELSE IF (OUTDIM == 3) THEN
        DO SA=1,NP
        DO SMA=1,NQ
        WRITE(UNITCOUNT,*) 'T=',STEPID*DT
        WRITE(UNITCOUNT,*) (((ETA(SA,SMA,IX,IY,IZ),IZ=1,NZ)
     &  ,IY=1,NY),IX=1,NX)
        UNITCOUNT=UNITCOUNT+1
        END DO
        END DO     
      END IF
!
      END SUBROUTINE OET
!-----------------------------------------------------------------------
! SUBROUTINE TO OPEN ALL OUTPUT FILES
!-----------------------------------------------------------------------
      SUBROUTINE OPN
      USE VAR
      IMPLICIT NONE
      INTEGER :: SA,SMA      !! SLIP PLANE AND DIRECTION IDs
      INTEGER :: UNITCOUNT   !! OUTPUT FILE NUMBER (UNIT)
!
      IF (OUTGE==0) THEN
      ELSE IF (OUTGE==1) THEN
      OPEN(UNIT=20,FILE='ge.dis',STATUS='UNKNOWN')
      END IF
!
      IF (OUTFE==0) THEN
      ELSE IF (OUTFE==1) THEN      
      OPEN(UNIT=21,FILE='fe.dis',  STATUS='UNKNOWN')
      END IF
!
      UNITCOUNT=24
!
      IF (OUTETAS==0) THEN
      ELSE IF (OUTETAS==1) THEN 
!      
        DO SA=1,NP
        DO SMA=1,NQ
        WRITE (FILENAME,'(A3,I1,I1,A4)')'eta',SA,SMA,'.dis'
        OPEN (UNIT=UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')
        UNITCOUNT=UNITCOUNT+1
        END DO
        END DO
 !     
      END IF
      END SUBROUTINE OPN
!-----------------------------------------------------------------------
! SUBROUTINE TO CLOSE ALL OPEN OUTPUT FILES
!-----------------------------------------------------------------------
      SUBROUTINE CLS
      USE VAR
      IMPLICIT NONE
      INTEGER :: SA,SMA      !! SLIP PLANE AND DIRECTION IDs
      INTEGER :: UNITCOUNT   !! OUTPUT FILE NUMBER (UNIT)      
! CLOSE ALL FILES
!
      IF (OUTGE==0) THEN
        PRINT *,'>>NOTE - GRADIENT ENERGY WILL NOT BE SAVED'
      ELSE IF (OUTGE==1) THEN
      CLOSE(20)
      END IF
!
      IF (OUTFE==0) THEN
        PRINT *,'>>NOTE - CRYSTAL ENERGY WILL NOT BE SAVED'      
      ELSE IF (OUTFE==1) THEN      
      CLOSE(21)
      END IF
!
      UNITCOUNT=24
!
      IF (OUTETAS==0) THEN
        PRINT *,'>>NOTE - FIELD VARIABLES ETA WILL NOT BE SAVED'      
      ELSE IF (OUTETAS==1) THEN
!      
        DO SA=1,NP
        DO SMA=1,NQ
        CLOSE (UNITCOUNT)
        UNITCOUNT=UNITCOUNT+1
        END DO
        END DO
!        
      END IF
      END SUBROUTINE CLS
!-----------------------------------------------------------------------
! SUBROUTINE TO OUTPUT FINAL ETAS
!-----------------------------------------------------------------------
      SUBROUTINE FET
      USE VAR
      IMPLICIT NONE
      INTEGER :: IX,IY,IZ
      INTEGER :: SA,SMA
      INTEGER :: UNITCOUNT   !! OUTPUT FILE NUMBER (UNIT)            
!
      UNITCOUNT=50
!
      IF (OUTDIM == 2) THEN
        IF (SLICE_AXIS == 1) THEN
          DO SA=1,NP
          DO SMA=1,NQ
            WRITE (FILENAME,'(A3,I1,I1,A10)')'eta',SA,SMA,'_final.dis'
            OPEN (UNIT=UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')        
            WRITE(UNITCOUNT,*) ((ETA(SA,SMA,SLICE_POS,IY,IZ)
     &      ,IZ=1,NZ),IY=1,NY)
            CLOSE(UNITCOUNT)
            UNITCOUNT=UNITCOUNT+1
          END DO
          END DO
        ELSE IF (SLICE_AXIS == 2) THEN
          DO SA=1,NP
          DO SMA=1,NQ
            WRITE (FILENAME,'(A3,I1,I1,A10)')'eta',SA,SMA,'_final.dis'
            OPEN (UNIT=UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')        
            WRITE(UNITCOUNT,*) ((ETA(SA,SMA,IX,SLICE_POS,IZ)
     &      ,IZ=1,NZ),IX=1,NX)
            CLOSE(UNITCOUNT)
            UNITCOUNT=UNITCOUNT+1
          END DO
          END DO
        ELSE IF (SLICE_AXIS == 3) THEN
          DO SA=1,NP
          DO SMA=1,NQ
            WRITE (FILENAME,'(A3,I1,I1,A10)')'eta',SA,SMA,'_final.dis'
            OPEN (UNIT=UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')        
            WRITE(UNITCOUNT,*) ((ETA(SA,SMA,IX,IY,SLICE_POS)
     &      ,IY=1,NY),IX=1,NX)
            CLOSE(UNITCOUNT)
            UNITCOUNT=UNITCOUNT+1
          END DO
          END DO
        END IF  
!  
      ELSE IF (OUTDIM == 3) THEN
        DO SA=1,NP
        DO SMA=1,NQ
        WRITE (FILENAME,'(A3,I1,I1,A10)')'eta',SA,SMA,'_final.dis'
        OPEN (UNIT=UNITCOUNT, FILE=FILENAME, FORM = 'FORMATTED')        
        WRITE(UNITCOUNT,*) (((ETA(SA,SMA,IX,IY,IZ),IZ=1,NZ)
     &  ,IY=1,NY),IX=1,NX)
        CLOSE(UNITCOUNT)
        UNITCOUNT=UNITCOUNT+1
        END DO
        END DO
      END IF
!      
      END SUBROUTINE FET      