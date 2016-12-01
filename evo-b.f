!-----------------------------------------------------------------------
! SUBROUTINE DEFINING TEMPORAL EVOLUTION OF THE SYSTEM
!-----------------------------------------------------------------------
      SUBROUTINE EVO(STEPID)
      USE VAR 
      IMPLICIT NONE
! DEFINE VARIABLES   
      INTEGER :: STEPID
      REAL :: DFDETA(NX,NY,NZ)        !! DETA/DT
      COMPLEX :: DFKDETA(NX,NY,NZ)    !! DETA/DT IN K-SPACE
      REAL :: TP,B_ALPHA,B_ALPHAP     !! B MAGNITUDES
      REAL :: BB(NP,NQ,NP,NQ)         !! COSINE BETWEEN B
      INTEGER :: SA,SMA               !! SLIP PLANE, SLIP DIRECTION
      INTEGER :: SAP, SMAP            !! SA', SMA'
      INTEGER :: IX,IY,IZ             !! POSITION INDICES
      INTEGER :: I,J                  !! DUMMY INDICES FOR IJK DIRECTIONS
      REAL :: E1,E2,E3                !! TEMPORARY ETAS FOR GSF CALC.
! UPDATE ETA IN K-SPACE
      DO SA=1,NP
      DO SMA=1,NQ
!$OMP PARALLEL DO PRIVATE(IY,IZ)      
        DO IX=1,NX
        DO IY=1,NY
        DO IZ=1,NZ
          TPR(IX,IY,IZ)=ETA(SA,SMA,IX,IY,IZ)
        END DO
        END DO
        END DO
!$OMP END PARALLEL DO        
      CALL FFT_R2C(TPR,TPK)
!$OMP PARALLEL DO PRIVATE(IY,IZ)        
        DO IX=1,NX
        DO IY=1,NY
        DO IZ=1,NZ
          ETAK(SA,SMA,IX,IY,IZ)=TPK(IX,IY,IZ)                          
        END DO
        END DO
        END DO
!$OMP END PARALLEL DO        
      END DO
      END DO
! CALCULATE B MAGNITUDES AND BB COSINE
      DO SA=1,NP
      DO SMA=1,NQ
        B_ALPHA=SQRT(BV(SA,SMA,1)**2+BV(SA,SMA,2)**2
     &    +BV(SA,SMA,3)**2)
        DO SAP=1,NP
        DO SMAP=1,NQ
          B_ALPHAP=SQRT(BV(SAP,SMAP,1)**2+BV(SAP,SMAP,2)**2
     &      +BV(SAP,SMAP,3)**2)
          BB(SA,SMA,SAP,SMAP)=(
     &      BV(SA,SMA,1)*BV(SAP,SMAP,1)+
     &      BV(SA,SMA,2)*BV(SAP,SMAP,2)+
     &      BV(SA,SMA,3)*BV(SAP,SMAP,3))/
     &      (B_ALPHA*B_ALPHAP)
        END DO
        END DO
      END DO
      END DO
! BEGIN CALCULATING ETAS AND THEIR FUNCTIONS     
      DO SA=1,NP !! BEGIN MAIN ITERATIVE LOOP
      DO SMA=1,NQ
!$OMP PARALLEL DO PRIVATE(IY,IZ,TP)      
        DO IX=1,NX !! BEGIN IX,IY,IZ LOOP 1
        DO IY=1,NY
        DO IZ=1,NZ
          DFKDETA(IX,IY,IZ)=(0.,0.) !!STARTING VALUE
! INTRODUCE GRADIENT ENERGY COMPONENT
          DO SAP=1,NP
          DO SMAP=1,NQ
            TP=0.0
            DO I=1,3
            DO J=1,3
              TP=TP+NV(SA,I)*NV(SAP,J)*G(I,IX,IY,IZ)*G(J,IX,IY,IZ)
            END DO
            END DO
            TP=TP-G2(IX,IY,IZ)*
     &        (NV(SA,1)*NV(SAP,1)+NV(SA,2)*NV(SAP,2)+NV(SA,3)*NV(SAP,3))
            DFKDETA(IX,IY,IZ)=DFKDETA(IX,IY,IZ)+
     &        TP*ETAK(SAP,SMAP,IX,IY,IZ)*BB(SA,SMA,SAP,SMAP)
          END DO
          END DO
          DFKDETA(IX,IY,IZ)=-DFKDETA(IX,IY,IZ)*BETA
! INTRODUCE ELASTIC ENERGY COMPONENT
          DO SAP=1,NP
          DO SMAP=1,NQ
            DFKDETA(IX,IY,IZ)=DFKDETA(IX,IY,IZ)+
     &        BPQ(SA,SMA,SAP,SMAP,IX,IY,IZ)*
     &        ETAK(SAP,SMAP,IX,IY,IZ)
          END DO
          END DO
        END DO
        END DO
        END DO !! END IX,IY,IZ LOOP 1
!$OMP END PARALLEL DO 
!
! FFT ETAS BACK TO REAL SPACE
      CALL INVFFT_C2R(DFKDETA,DFDETA)
! INTRODUCE FAULT ENERGY COMPONENT
!$OMP PARALLEL DO PRIVATE(IY,IZ,E1,E2,E3,PHID)  
        DO IX=1,NX !! BEGIN IX,IY,IZ LOOP 2
        DO IY=1,NY
        DO IZ=1,NZ
          E1=ETA(SA,SMA,IX,IY,IZ)
          E2=ETA(SA,MOD(SMA,NQ)+1,IX,IY,IZ)
          E3=ETA(SA,MOD(SMA+1,NQ)+1,IX,IY,IZ)
!-----------------------------------------------------------------------
! CRYSTALLINE ENERGY GRADIENT FOR FCC 111 PLANES
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! PURE FCC SCHEME #1 WITH [1 LAYER]
!-----------------------------------------------------------------------
          PHID=INT(PHI(IX,IY,IZ)+1)
          IF(GSF_EQ(PHID) == 10) THEN 
            DFDETA(IX,IY,IZ)=DFDETA(IX,IY,IZ)+
     &        (PI*(
     &        -GC(PHID,4)*(2*SIN(TWOPI*(E1-E2)) 
     &          +2*SIN(TWOPI*(E1-E3)) )
     &        -GC(PHID,7)*(4*SIN(TWOPI*(2*E1-E2-E3))
     &          +2*SIN(TWOPI*(E1-2*E2+E3)) 
     &          +2*SIN(TWOPI*(E1+E2-2*E3)))
     &        -GC(PHID,8)*(4*SIN(4*PI*(E1-E2))+4*SIN(4*PI*(E1-E3)))
     &        -GC(PHID,9)*(4*SIN(TWOPI*(2*E1-3*E2+E3))
     &          +6*SIN(TWOPI*(3*E1-2*E2-E3))
     &          +4*SIN(TWOPI*(2*E1+E2-3*E3))
     &          +2*SIN(TWOPI*(E1+2*E2-3*E3))
     &          +6*SIN(TWOPI*(3*E1-E2-2*E3))
     &          +2*SIN(TWOPI*(E1-3*E2+2*E3)))     
     &        +GC(PHID,11)*(2*COS(TWOPI*(E1-E2))-2*COS(TWOPI*(E1-E3)))
     &        +GC(PHID,14)*(4*COS(4*PI*(E1-E2))-4*COS(4*PI*(E1-E3)))
     &        +GC(PHID,15)*(4*COS(TWOPI*(2*E1-3*E2+E3))
     &          +6*COS(TWOPI*(3*E1-2*E2-E3))
     &          -4*COS(TWOPI*(2*E1+E2-3*E3))
     &          +2*COS(TWOPI*(E1+2*E2-3*E3))
     &          -6*COS(TWOPI*(3*E1-E2-2*E3))
     &          -2*COS(TWOPI*(E1-3*E2+2*E3))))) 
!-----------------------------------------------------------------------
! PURE L1_2 ORDERED COMPOUND SCHEME #1A [1 LAYER]
!-----------------------------------------------------------------------     
          ELSE IF(GSF_EQ(PHID) == 11) THEN 
            DFDETA(IX,IY,IZ)=DFDETA(IX,IY,IZ)+
     &        (PI*(
     &        -GC(PHID,2)*(SIN(PI*(E1-E2))+SIN(PI*(E1-E3))) 
     &        -GC(PHID,3)*(2*SIN(PI*(2*E1-E2-E3))
     &          +SIN(PI*(E1-2*E2+E3)) 
     &          +SIN(PI*(E1+E2-2*E3)))
     &        -GC(PHID,4)*(2*SIN(TWOPI*(E1-E2)) 
     &          +2*SIN(TWOPI*(E1-E3)) )
     &        -GC(PHID,5)*(2*SIN(PI*(2*E1-3*E2+E3))
     &          +3*SIN(PI*(3*E1-2*E2-E3))
     &          +2*SIN(PI*(2*E1+E2-3*E3))
     &          +SIN(PI*(E1+2*E2-3*E3))
     &          +3*SIN(PI*(3*E1-E2-2*E3))
     &          +SIN(PI*(E1-3*E2+2*E3)))
     &        -GC(PHID,6)*(3*SIN(3*PI*(E1-E2))+3*SIN(3*PI*(E1-E3)))
     &        -GC(PHID,8)*(4*SIN(4*PI*(E1-E2))+4*SIN(4*PI*(E1-E3)))    
     &        +GC(PHID,10)*(COS(PI*(E1-E2))-COS(PI*(E1-E3)))
     &        +GC(PHID,11)*(2*COS(TWOPI*(E1-E2))-2*COS(TWOPI*(E1-E3)))
     &        +GC(PHID,12)*(2*COS(PI*(2*E1-3*E2+E3))
     &          +3*COS(PI*(3*E1-2*E2-E3))
     &          -2*COS(PI*(2*E1+E2-3*E3))
     &          +COS(PI*(E1+2*E2-3*E3))
     &          -3*COS(PI*(3*E1-E2-2*E3))
     &          -COS(PI*(E1-3*E2+2*E3)))
     &        +GC(PHID,13)*(3*COS(3*PI*(E1-E2))-3*COS(3*PI*(E1-E3)))
     &        +GC(PHID,14)*(4*COS(4*PI*(E1-E2))-4*COS(4*PI*(E1-E3)))
     &        ))
!-----------------------------------------------------------------------
! PURE L1_2 ORDERED COMPOUND SCHEME #1B [1 LAYER]
!-----------------------------------------------------------------------      
          ELSE IF(GSF_EQ(PHID) == 12) THEN 
            DFDETA(IX,IY,IZ)=DFDETA(IX,IY,IZ)+
     &        (PI*(
     &        -GC(PHID,2)*(SIN(PI*(E1-E2))+SIN(PI*(E1-E3))) 
     &        -GC(PHID,3)*(2*SIN(PI*(2*E1-E2-E3))
     &          +SIN(PI*(E1-2*E2+E3)) 
     &          +SIN(PI*(E1+E2-2*E3)))
     &        -GC(PHID,4)*(2*SIN(TWOPI*(E1-E2)) 
     &          +2*SIN(TWOPI*(E1-E3)) )
     &        -GC(PHID,5)*(2*SIN(PI*(2*E1-3*E2+E3))
     &          +3*SIN(PI*(3*E1-2*E2-E3))
     &          +2*SIN(PI*(2*E1+E2-3*E3))
     &          +SIN(PI*(E1+2*E2-3*E3))
     &          +3*SIN(PI*(3*E1-E2-2*E3))
     &          +SIN(PI*(E1-3*E2+2*E3)))
     &        -GC(PHID,6)*(3*SIN(3*PI*(E1-E2))+3*SIN(3*PI*(E1-E3)))
     &        -GC(PHID,7)*(4*SIN(TWOPI*(2*E1-E2-E3))
     &          +2*SIN(TWOPI*(E1-2*E2+E3)) 
     &          +2*SIN(TWOPI*(E1+E2-2*E3)))     
     &        -GC(PHID,8)*(4*SIN(4*PI*(E1-E2))+4*SIN(4*PI*(E1-E3)))    
     &        +GC(PHID,10)*(COS(PI*(E1-E2))-COS(PI*(E1-E3)))
     &        +GC(PHID,11)*(2*COS(TWOPI*(E1-E2))-2*COS(TWOPI*(E1-E3)))
     &        +GC(PHID,12)*(2*COS(PI*(2*E1-3*E2+E3))
     &          +3*COS(PI*(3*E1-2*E2-E3))
     &          -2*COS(PI*(2*E1+E2-3*E3))
     &          +COS(PI*(E1+2*E2-3*E3))
     &          -3*COS(PI*(3*E1-E2-2*E3))
     &          -COS(PI*(E1-3*E2+2*E3)))
     &        +GC(PHID,13)*(3*COS(3*PI*(E1-E2))-3*COS(3*PI*(E1-E3)))
     &        +GC(PHID,14)*(4*COS(4*PI*(E1-E2))-4*COS(4*PI*(E1-E3)))
     &        ))
!-----------------------------------------------------------------------
! PURE FCC SCHEME #2 WITH [2 LAYERS]
!-----------------------------------------------------------------------     
          ELSE IF(GSF_EQ(PHID) == 20) THEN  
            DFDETA(IX,IY,IZ)=DFDETA(IX,IY,IZ)+(
     &        GC(PHID,4)*(-2*PI*SIN(2*(E1-E2)*PI)+
     &          2*PI*SIN(2*(-E1+E3)*PI))+ 
     &        GC(PHID,5)*(-4*PI*SIN(2*(2*E1-E2-E3)*PI)+
     &          2*PI*SIN(2*(-E1+2*E2-E3)*PI)+ 
     &          2*PI*SIN(2*(-E1-E2+2*E3)*PI))+ 
     &        GC(PHID,6)*(-3*PI*SIN(3*(E1-E2)*PI)+
     &          3*PI*SIN(3*(-E1+E3)*PI))+ 
     &        GC(PHID,7)*(-6*PI*SIN(3*(2*E1-E2-E3)*PI)+
     &          3*PI*SIN(3*(-E1+2*E2-E3)*PI)+ 
     &          3*PI*SIN(3*(-E1-E2+2*E3)*PI))+
     &        GC(PHID,8)*(-4*PI*SIN(4*(E1-E2)*PI)+
     &          4*PI*SIN(4*(-E1+E3)*PI))+ 
     &        GC(PHID,11)*(2*PI*COS(2*(E1-E2)*PI)-
     &          2*PI*COS(2*(-E1+E3)*PI))+ 
     &        GC(PHID,13)*(4*PI*COS(4*(E1-E2)*PI)-
     &          4*PI*COS(4*(-E1+E3)*PI)))   
          ELSE IF(GSF_EQ(PHID) == 21) THEN
!-----------------------------------------------------------------------
! PURE L1_2 ORDERED COMPOUND SCHEME #2B [2 LAYERS] SOLVER 4
!-----------------------------------------------------------------------          
            DFDETA(IX,IY,IZ)=DFDETA(IX,IY,IZ)+(
     &        GC(PHID,2)*(-(PI*SIN((E1-E2)*PI))+PI*SIN((-E1+E3)*PI))+ 
     &        GC(PHID,3)*(-2*PI*SIN((2*E1-E2-E3)*PI)+
     &          PI*SIN((-E1+2*E2-E3)*PI)+ 
     &          PI*SIN((-E1-E2+2*E3)*PI))+ 
     &        GC(PHID,4)*(-2*PI*SIN(2*(E1-E2)*PI)+
     &          2*PI*SIN(2*(-E1+E3)*PI))+ 
     &        GC(PHID,5)*(-4*PI*SIN(2*(2*E1-E2-E3)*PI)+
     &          2*PI*SIN(2*(-E1+2*E2-E3)*PI)+ 
     &          2*PI*SIN(2*(-E1-E2+2*E3)*PI))+ 
     &        GC(PHID,6)*(-3*PI*SIN(3*(E1-E2)*PI)+
     &          3*PI*SIN(3*(-E1+E3)*PI))+ 
     &        GC(PHID,7)*(-6*PI*SIN(3*(2*E1-E2-E3)*PI)+
     &          3*PI*SIN(3*(-E1+2*E2-E3)*PI)+ 
     &          3*PI*SIN(3*(-E1-E2+2*E3)*PI))+
     &        GC(PHID,8)*(-4*PI*SIN(4*(E1-E2)*PI)+
     &          4*PI*SIN(4*(-E1+E3)*PI))+ 
     &        GC(PHID,9)*(-5*PI*SIN(5*(E1-E2)*PI)+
     &          5*PI*SIN(5*(-E1+E3)*PI))+ 
     &        GC(PHID,10)*(PI*COS((E1-E2)*PI)-PI*COS((-E1+E3)*PI))+ 
     &        GC(PHID,11)*(2*PI*COS(2*(E1-E2)*PI)-
     &          2*PI*COS(2*(-E1+E3)*PI))+ 
     &        GC(PHID,12)*(3*PI*COS(3*(E1-E2)*PI)-
     &          3*PI*COS(3*(-E1+E3)*PI))+ 
     &        GC(PHID,13)*(4*PI*COS(4*(E1-E2)*PI)-
     &          4*PI*COS(4*(-E1+E3)*PI))+ 
     &        GC(PHID,14)*(5*PI*COS(5*(E1-E2)*PI)-
     &          5*PI*COS(5*(-E1+E3)*PI)))     
          END IF
! INTRODUCE WORK COMPONENT DONE BY APPLIED STRESS
          DO J=1,3
          DO I=1,3
            DFDETA(IX,IY,IZ)=DFDETA(IX,IY,IZ)-
     &        S_APP(I,J)*E0(SA,SMA,I,J)
          END DO
          END DO
        END DO
        END DO
        END DO !! END IX,IY,IZ LOOP 2
!$OMP END PARALLEL DO 
!
! BEGIN SYSTEM EVOLUTION
!
!$OMP PARALLEL DO PRIVATE(IY,IX) 
        DO IZ=1+BORDER,NZ-BORDER
        DO IY=1+BORDER,NY-BORDER
        DO IX=1+BORDER,NX-BORDER
          ETA(SA,SMA,IX,IY,IZ)=ETA(SA,SMA,IX,IY,IZ) 
     &      -DT*DFDETA(IX,IY,IZ)
        END DO
        END DO
        END DO
!$OMP END PARALLEL DO         
      END DO
      END DO !! END MAIN ITERATIVE LOOP      
      END SUBROUTINE EVO